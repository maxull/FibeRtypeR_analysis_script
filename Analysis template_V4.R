
#!/usr/bin/env Rscript
# =============================================================================
# Template for exploring relationship between muscle fiber type (FibeRtypeR)
# and insulin resistance
#
# Author: Max Ullrich
# Last updated: 2026-01-16
#
# What collaborators edit: ONLY the CONFIG section below.
# Everything else runs automatically with one command:   run_all()
# =============================================================================

# -----------------------------------------------------------------------------
# 1) Deterministic behavior --- DO NOT CHANGE
# -----------------------------------------------------------------------------
Sys.setenv(TZ = "UTC")
Sys.setlocale(category = "LC_ALL", locale = "C")
suppressWarnings({
        RNGkind(kind = "Mersenne-Twister",
                normal.kind = "Inversion",
                sample.kind = "Rounding")
})
set.seed(20240101)
options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)

# -----------------------------------------------------------------------------
# 2) CONFIG (EDIT)
# -----------------------------------------------------------------------------

# --- Input paths
metadata_path  <- "/Users/maxullrich/Documents/GitHub/FiberTypeR-applied/data/testdata/sim_metadata.csv"
fibertype_path <- "/Users/maxullrich/Documents/GitHub/FiberTypeR-applied/data/testdata/sim_fibertype.csv"

# --- Study identifier (used in file names). If NULL, derived from metadata filename.
study_id <- "PMID40683252"

# --- Outcome to model (choose 1 present in data): "HOMA_IR", "M_value", or "Matsuda_index"
# If you have multiple outcomes, run the script separately for each.
outcome_var <- "HOMA_IR"

# --- Optional: filter to baseline/Pre values (set both to NULL to skip)
baseline_filter_col <- "Condition"
baseline_filter_val <- "Pre"

# --- Output root (folders created if missing)
output_dir <- paste0(study_id, "_outputs_4")  # creates {figures, results, logs}

# --- Column mapping: collaborator's original names -> standardized names below
# Fill only what exists; leave others NULL (handled gracefully).
column_map <- list(
        RNAseq_ID       = "GSE",                 # REQUIRED: must match fibertype table IDs
        BMI             = "d14.bmi",
        HOMA_IR         = "homair",
        M_value         = "m",
        Matsuda_index   = NULL,
        age             = "age",
        sex             = "sex",
        waist_hip_ratio = "d14.waist.hip",
        fat_percentage  = "total.whole.body.fat.1",
        # Manual mapping required for disease models (E):
        disease_state   = "disease_state"                   # e.g., "group", "status", "disease"
)

# --- FibeRtypeR table mapping (final standardized names)
fiber_id_col   <- "RNAseq_ID"
fiber_fast_col <- "fast"
fiber_slow_col <- "slow"

# -----------------------------------------------------------------------------
# 3) renv bootstrap (NO manual version pinning)
# -----------------------------------------------------------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
        install.packages("renv", repos = "https://cloud.r-project.org")
}
if (!file.exists("renv.lock")) {
        message("No renv.lock found — creating a new environment.")
        renv::init()
} else {
        message("renv.lock found — restoring environment for reproducibility.")
        renv::restore(prompt = FALSE)
}

# -----------------------------------------------------------------------------
# 4) Libraries
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
        library(tidyverse)
        library(readr)
        library(broom)
        library(sandwich)
        library(lmtest)
        library(ggplot2)
        library(purrr)
        library(tibble)
})

# -----------------------------------------------------------------------------
# 5) INTERNALS (DO NOT EDIT)
# -----------------------------------------------------------------------------

make_dirs <- function(root) {
        figs <- file.path(root, "figures")
        res  <- file.path(root, "results")
        logs <- file.path(root, "logs")
        dir.create(root, showWarnings = FALSE, recursive = TRUE)
        dir.create(figs, showWarnings = FALSE, recursive = TRUE)
        dir.create(res,  showWarnings = FALSE, recursive = TRUE)
        dir.create(logs, showWarnings = FALSE, recursive = TRUE)
        list(figures = figs, results = res, logs = logs)
}

safe_study_id <- function(meta_path, user_id = NULL) {
        if (!is.null(user_id) && nzchar(user_id)) return(user_id)
        base <- tools::file_path_sans_ext(basename(meta_path))
        gsub("[^A-Za-z0-9_\\-]+", "_", base)
}

rename_columns <- function(df, map) {
        for (std_name in names(map)) {
                src <- map[[std_name]]
                if (!is.null(src) && src %in% names(df)) {
                        df <- dplyr::rename(df, !!std_name := all_of(src))
                }
        }
        df
}

coerce_numeric <- function(df, cols) {
        cols <- intersect(cols, names(df))
        if (length(cols) == 0) return(df)
        df %>% mutate(across(all_of(cols), ~ suppressWarnings(as.numeric(.x))))
}

# keep original formula text; do NOT drop interactions
prune_formula <- function(fm, data) {
        trm <- terms(fm)
        resp <- as.character(attr(trm, "variables")[[2]])
        labs <- attr(trm, "term.labels")
        
        if (length(labs) == 0) return(fm)
        
        # Helper: check if a column has variation (at least 2 unique non-NA values)
        has_variation <- function(col_name) {
                if (!col_name %in% names(data)) return(FALSE)
                val <- data[[col_name]]
                val <- val[!is.na(val)]
                if (length(val) < 2) return(FALSE)
                
                # Check unique values. 
                # Note: factor with 2 levels but only 1 used counts as 1 unique value here,
                # which is what we want (drop it to avoid lm crash).
                length(unique(val)) >= 2
        }
        
        keep_terms <- character(0)
        for (lbl in labs) {
                # Split interactions (e.g. "fast:sex" -> c("fast", "sex"))
                parts <- unlist(strsplit(lbl, ":", fixed = TRUE))
                # Keep term ONLY if all parts exist AND have variation
                if (all(vapply(parts, has_variation, logical(1)))) {
                        keep_terms <- c(keep_terms, lbl)
                }
        }
        
        if (length(keep_terms) == 0) {
                return(as.formula(paste(resp, "~ 1")))
        }
        as.formula(paste(resp, "~", paste(keep_terms, collapse = " + ")))
}

label_formula <- function(fm) paste0(deparse(fm), collapse = "")

# Robust predictions with HC3
predict_with_robust <- function(mod, newdata, refs) {
        trm <- stats::terms(mod)
        preds <- attr(trm, "term.labels")
        needed <- setdiff(preds, names(newdata))
        if (length(needed) > 0) {
                newdata <- dplyr::bind_cols(newdata, refs %>% dplyr::select(dplyr::any_of(needed)))
        }
        X <- model.matrix(stats::delete.response(trm), newdata)
        V <- sandwich::vcovHC(mod, type = "HC3")
        beta <- stats::coef(mod)
        common <- intersect(colnames(X), names(beta))
        X <- X[, common, drop = FALSE]
        beta <- beta[common]
        fit <- as.vector(X %*% beta)
        se  <- sqrt(rowSums((X %*% V[common, common, drop = FALSE]) * X))
        tibble::tibble(fit = fit, se = se)
}

# Count complete cases for given variables
n_complete <- function(data, vars) {
        vars <- intersect(vars, names(data))
        if (length(vars) == 0) return(0L)
        sum(stats::complete.cases(data[, vars, drop = FALSE]))
}

# Fit LM + HC3 and return tidy/glance with annotations and CIs
run_model <- function(formula, data, model_id, group, adjustment = NA_character_,
                      stratum_type = NA_character_, stratum_label = NA_character_,
                      extras = list()) {
        mod <- lm(formula, data)
        V <- sandwich::vcovHC(mod, type = "HC3")
        ct <- lmtest::coeftest(mod, vcov. = V)
        
        tidy_df <- broom::tidy(ct) %>%
                mutate(
                        conf.low  = estimate - 1.96 * std.error,
                        conf.high = estimate + 1.96 * std.error,
                        model_id = model_id,
                        model_group = group,
                        adjustment = adjustment,
                        stratum_type = stratum_type,
                        stratum_label = stratum_label,
                        n = stats::nobs(mod),
                        model_formula = paste0(deparse(formula(mod)), collapse = "")
                )
        
        gl <- broom::glance(mod) %>%
                mutate(
                        model_id = model_id,
                        model_group = group,
                        adjustment = adjustment,
                        stratum_type = stratum_type,
                        stratum_label = stratum_label,
                        n = stats::nobs(mod),
                        model_formula = paste0(deparse(formula(mod)), collapse = "")
                )
        
        if (length(extras)) {
                for (nm in names(extras)) {
                        tidy_df[[nm]] <- extras[[nm]]
                        gl[[nm]]      <- extras[[nm]]
                }
        }
        list(model = mod, tidy = tidy_df, glance = gl)
}

add_partial_r2 <- function(tidy_df, glance_df) {
        df_map <- glance_df %>% select(model_id, df.residual)
        tidy_df %>%
                left_join(df_map, by = "model_id") %>%
                mutate(partial_r2 = ifelse(!is.na(statistic) & !is.na(df.residual),
                                           (statistic^2) / (statistic^2 + df.residual),
                                           NA_real_))
}

# -------- Dynamic base builder: degrades gracefully if age/sex missing --------
build_formulas_base <- function(outcome, data, include_z = TRUE) {
        has <- function(v) v %in% names(data) && any(!is.na(data[[v]]))
        has_age <- has("age")
        has_sex <- has("sex")
        has_bmi <- has("BMI")
        has_fat <- has("fat_percentage")
        has_whr <- has("waist_hip_ratio")
        
        models <- list()
        add <- function(name, rhs) {
                models[[name]] <<- as.formula(paste(outcome, "~", rhs))
                if (include_z) {
                        rhs_z <- gsub("\\bfast\\b", "fast_z", rhs)
                        models[[paste0(name, "_z")]] <<- as.formula(paste(outcome, "~", rhs_z))
                }
        }
        
        # A: unadjusted
        add("A_unadj", "fast")
        
        # B: minimally adjusted (prefer age+sex; degrade if one missing)
        if (has_age && has_sex) {
                add("B_minimal", "fast + age + sex")
        } else if (has_age) {
                add("B_minimal_age_only", "fast + age")
        } else if (has_sex) {
                add("B_minimal_sex_only", "fast + sex")
        }
        
        # C: adiposity separately (BMI / fat% / WHR), with graceful degradation
        if (has_bmi) {
                if (has_age && has_sex) add("C_BMI", "fast + age + sex + BMI")
                else if (has_age)       add("C_BMI_age_only", "fast + age + BMI")
                else if (has_sex)       add("C_BMI_sex_only", "fast + sex + BMI")
                else                    add("C_BMI_no_demo", "fast + BMI")
        }
        if (has_fat) {
                if (has_age && has_sex) add("C_fatpct", "fast + age + sex + fat_percentage")
                else if (has_age)       add("C_fatpct_age_only", "fast + age + fat_percentage")
                else if (has_sex)       add("C_fatpct_sex_only", "fast + sex + fat_percentage")
                else                    add("C_fatpct_no_demo", "fast + fat_percentage")
        }
        if (has_whr) {
                if (has_age && has_sex) add("C_WHR", "fast + age + sex + waist_hip_ratio")
                else if (has_age)       add("C_WHR_age_only", "fast + age + waist_hip_ratio")
                else if (has_sex)       add("C_WHR_sex_only", "fast + sex + waist_hip_ratio")
                else                    add("C_WHR_no_demo", "fast + waist_hip_ratio")
        }
        
        models
}

# Map model_id to friendly adjustment label
adjustment_label_from_name <- function(model_name) {
        base <- gsub("_z$", "", model_name)
        dplyr::case_when(
                base == "A_unadj"                ~ "none",
                base == "B_minimal"              ~ "age+sex",
                base == "B_minimal_age_only"     ~ "age (sex missing)",
                base == "B_minimal_sex_only"     ~ "sex (age missing)",
                base == "C_BMI"                  ~ "age+sex+BMI",
                base == "C_BMI_age_only"         ~ "age+BMI",
                base == "C_BMI_sex_only"         ~ "sex+BMI",
                base == "C_BMI_no_demo"          ~ "BMI only",
                base == "C_fatpct"               ~ "age+sex+fat_percentage",
                base == "C_fatpct_age_only"      ~ "age+fat_percentage",
                base == "C_fatpct_sex_only"      ~ "sex+fat_percentage",
                base == "C_fatpct_no_demo"       ~ "fat_percentage only",
                base == "C_WHR"                  ~ "age+sex+waist_hip_ratio",
                base == "C_WHR_age_only"         ~ "age+waist_hip_ratio",
                base == "C_WHR_sex_only"         ~ "sex+waist_hip_ratio",
                base == "C_WHR_no_demo"          ~ "waist_hip_ratio only",
                TRUE                             ~ NA_character_
        )
}

# Age grouping helpers (18–30, 30–50, 50+; fallback <50 vs >=50)
add_age_groups <- function(df) {
        df %>%
                mutate(
                        age_group = case_when(
                                !is.na(age) & age < 30 ~ "18_30",
                                !is.na(age) & age < 50 ~ "30_50",
                                !is.na(age) & age >= 50 ~ "50_plus",
                                TRUE ~ NA_character_
                        ),
                        age_group2 = case_when(
                                !is.na(age) & age < 50 ~ "<50",
                                !is.na(age) & age >= 50 ~ ">=50",
                                TRUE ~ NA_character_
                        )
                )
}

# -----------------------------------------------------------------------------
# 6) MAIN PIPELINE
# -----------------------------------------------------------------------------
run_all <- function() {
        message("=== Starting analysis ===")
        
        # Prepare directories
        dirs <- make_dirs(output_dir)
        sid  <- safe_study_id(metadata_path, study_id)
        
        # Single required log file
        sess_file <- file.path(dirs$logs, paste0(sid, "_sessionInfo.txt"))
        sess_txt <- utils::capture.output(sessionInfo())
        writeLines(sess_txt, sess_file)
        
        # ------------------------ Load data ------------------------
        stopifnot(file.exists(metadata_path))
        stopifnot(file.exists(fibertype_path))
        meta_raw <- readr::read_csv(metadata_path, show_col_types = FALSE)
        fib_raw  <- readr::read_csv(fibertype_path, show_col_types = FALSE)
        colnames(fib_raw) <- c("RNAseq_ID", "slow", "fast")
        
        # ------------------------ Normalize/rename ------------------------
        meta <- rename_columns(meta_raw, column_map)
        
        # Optional baseline filter
        if (!is.null(baseline_filter_col) && !is.null(baseline_filter_val) &&
            baseline_filter_col %in% names(meta)) {
                meta <- meta %>% filter(.data[[baseline_filter_col]] == !!baseline_filter_val)
        }
        
        # Keep standardized columns that exist
        std_cols <- c("RNAseq_ID","BMI","HOMA_IR","M_value","Matsuda_index",
                      "age","sex","waist_hip_ratio","fat_percentage","disease_state")
        meta <- meta %>% select(any_of(std_cols))
        
        # Coerce numeric variables
        numeric_candidates <- c("BMI","HOMA_IR","M_value","Matsuda_index","age","waist_hip_ratio","fat_percentage")
        meta <- coerce_numeric(meta, numeric_candidates)
        
        # Normalize sex/disease_state as factors if present
        if ("sex" %in% names(meta)) meta$sex <- factor(meta$sex)
        if ("disease_state" %in% names(meta)) meta$disease_state <- factor(meta$disease_state)
        
        # Fibertype normalization
        if (!("RNAseq_ID" %in% names(fib_raw))) {
                cand <- intersect(names(fib_raw), c("RNAseq_ID","sample","Sample","ID","id","GSE"))
                if (length(cand) == 0) stop("Fibertype file must have an ID column (e.g., RNAseq_ID).")
                fib_raw <- fib_raw %>% rename(RNAseq_ID = all_of(cand[1]))
        }
        if (!("fast" %in% names(fib_raw))) {
                cand_fast <- intersect(names(fib_raw), c("fast","Fast","type_II","typeII","fast_twitch"))
                if (length(cand_fast) == 0) stop("Fibertype file must have a 'fast' (fast-twitch) column.")
                fib_raw <- fib_raw %>% rename(fast = all_of(cand_fast[1]))
        }
        if (!("slow" %in% names(fib_raw))) {
                cand_slow <- intersect(names(fib_raw), c("slow","Slow","type_I","typeI","slow_twitch"))
                if (length(cand_slow) > 0) {
                        fib_raw <- fib_raw %>% rename(slow = all_of(cand_slow[1]))
                } else {
                        fib_raw$slow <- NA_real_
                }
        }
        fib <- fib_raw %>%
                select(RNAseq_ID, slow, fast) %>%
                coerce_numeric(c("slow","fast"))
        
        # ------------------------ Merge & prepare ------------------------
        dat <- meta %>%
                inner_join(fib, by = "RNAseq_ID") %>%
                mutate(
                        fast_z = as.numeric(scale(fast)),
                        age_centered = if ("age" %in% names(.)) age - mean(age, na.rm = TRUE) else NA_real_
                ) %>%
                add_age_groups()
        
        # Validate outcome
        if (!(outcome_var %in% names(dat))) {
                stop(sprintf("Configured outcome_var '%s' is not present after renaming.", outcome_var))
        }
        
        results <- list()
        
        # ------------------------ Base models ------------------------
        base_forms <- build_formulas_base(outcome_var, dat, include_z = TRUE)
        for (nm in names(base_forms)) {
                # Prune formula removes terms with 0 variance (avoids 'contrasts' error)
                fm <- prune_formula(base_forms[[nm]], dat)
                
                # Check if we have predictors left (intercept only? skip if desired)
                vars_needed <- all.vars(fm)[-1]
                if (length(vars_needed) == 0 && !grepl("unadj", nm)) next
                
                # Safety checks
                if (!all(vars_needed %in% names(dat))) next
                if (n_complete(dat, c(outcome_var, vars_needed)) < 5) next
                
                adj_label <- adjustment_label_from_name(nm)
                
                results[[nm]] <- run_model(
                        formula = fm, data = dat,
                        model_id = nm,
                        group = if (grepl("^A_", nm)) "A" else if (grepl("^B_", nm)) "B" else "C",
                        adjustment = adj_label
                )
        }
        
        # ------------------------ Model D: Age Interaction -------------------
        n_for_age <- n_complete(dat, c(outcome_var, "fast", "age"))
        if (n_for_age >= 20) {
                fm_D <- as.formula(paste(outcome_var, "~ fast * age_centered"))
                fm_D <- prune_formula(fm_D, dat)
                vars_needed <- all.vars(fm_D)[-1]
                
                if (length(vars_needed) > 0 && n_complete(dat, c(outcome_var, vars_needed)) >= 20) {
                        results[["D_age_interaction"]] <- run_model(
                                formula = fm_D, data = dat,
                                model_id = "D_age_interaction",
                                group = "D",
                                adjustment = "interaction fast*age_centered"
                        )
                }
                
                # Stratified models
                run_strata_models <- function(df, label_col, group_labels) {
                        B     <- as.formula(paste(outcome_var, "~ fast + age + sex"))
                        C_fat <- as.formula(paste(outcome_var, "~ fast + age + sex + fat_percentage"))
                        C_bmi <- as.formula(paste(outcome_var, "~ fast + age + sex + BMI"))
                        C_whr <- as.formula(paste(outcome_var, "~ fast + age + sex + waist_hip_ratio"))
                        
                        for (glab in group_labels) {
                                sub <- df %>% filter(.data[[label_col]] == glab)
                                if (n_complete(sub, c(outcome_var, "fast")) < 20) next
                                
                                for (tag in c("B_minimal","C_fatpct","C_BMI","C_WHR")) {
                                        form <- switch(tag, B_minimal=B, C_fatpct=C_fat, C_BMI=C_bmi, C_WHR=C_whr)
                                        fm <- prune_formula(form, sub)
                                        vars_needed <- all.vars(fm)[-1]
                                        if (n_complete(sub, c(outcome_var, vars_needed)) < 20) next
                                        
                                        rid <- paste0("D_strata_", glab, "_", tag)
                                        results[[rid]] <- run_model(
                                                formula = fm, data = sub,
                                                model_id = rid,
                                                group = "D",
                                                adjustment = tag,
                                                stratum_type = "age",
                                                stratum_label = glab
                                        )
                                }
                        }
                }
                
                # Check age strata sufficiency
                tab3 <- table(dat$age_group)
                if (all(c("18_30","30_50","50_plus") %in% names(tab3)) && min(tab3) >= 20) {
                        run_strata_models(dat, "age_group", c("18_30","30_50","50_plus"))
                } else {
                        tab2 <- table(dat$age_group2)
                        if (all(c("<50",">=50") %in% names(tab2)) && min(tab2) >= 20) {
                                run_strata_models(dat, "age_group2", c("<50",">=50"))
                        }
                }
        }
        
        # ------------------------ Model E: Disease Interaction -------------------
        if ("disease_state" %in% names(dat)) {
                ds <- dat %>% filter(!is.na(disease_state))
                if (nrow(ds) >= 60 && length(unique(ds$disease_state)) >= 2) {
                        fm_E <- as.formula(paste(outcome_var, "~ fast * disease_state + age + sex"))
                        fm_E <- prune_formula(fm_E, ds)
                        vars_needed <- all.vars(fm_E)[-1]
                        
                        if (n_complete(ds, c(outcome_var, vars_needed)) >= 60) {
                                results[["E_disease_interaction"]] <- run_model(
                                        formula = fm_E, data = ds,
                                        model_id = "E_disease_interaction",
                                        group = "E",
                                        adjustment = "interaction fast*disease_state"
                                )
                        }
                        # Stratified
                        for (lev in unique(ds$disease_state)) {
                                sub <- ds %>% filter(disease_state == lev)
                                if (nrow(sub) < 20) next
                                # Only run if sub has variation? prune_formula handles it inside run_strata logic if we copied it fully, 
                                # but here we just manually loop.
                                
                                # Define forms
                                B     <- as.formula(paste(outcome_var, "~ fast + age + sex"))
                                C_fat <- as.formula(paste(outcome_var, "~ fast + age + sex + fat_percentage"))
                                
                                for (tag in c("B_minimal", "C_fatpct")) {
                                        form <- switch(tag, B_minimal=B, C_fatpct=C_fat)
                                        fm <- prune_formula(form, sub)
                                        vars_needed <- all.vars(fm)[-1]
                                        if (n_complete(sub, c(outcome_var, vars_needed)) < 20) next
                                        
                                        rid <- paste0("E_strata_", lev, "_", tag)
                                        results[[rid]] <- run_model(
                                                formula = fm, data = sub,
                                                model_id = rid,
                                                group = "E",
                                                adjustment = tag,
                                                stratum_type = "disease_state",
                                                stratum_label = lev
                                        )
                                }
                        }
                }
        }
        
        # ------------------------ Model F: Sex Interaction -------------------
        if ("sex" %in% names(dat)) {
                sx <- dat %>% filter(!is.na(sex))
                if (nrow(sx) >= 60 && length(unique(sx$sex)) >= 2) {
                        fm_F <- as.formula(paste(outcome_var, "~ fast * sex + age"))
                        fm_F <- prune_formula(fm_F, sx)
                        vars_needed <- all.vars(fm_F)[-1]
                        
                        if (n_complete(sx, c(outcome_var, vars_needed)) >= 60) {
                                results[["F_sex_interaction"]] <- run_model(
                                        formula = fm_F, data = sx,
                                        model_id = "F_sex_interaction",
                                        group = "F",
                                        adjustment = "interaction fast*sex"
                                )
                        }
                        
                        for (lev in unique(sx$sex)) {
                                sub <- sx %>% filter(sex == lev)
                                if (nrow(sub) < 30) next
                                
                                B <- as.formula(paste(outcome_var, "~ fast + age + sex")) # 'sex' will be pruned by prune_formula
                                C_fat <- as.formula(paste(outcome_var, "~ fast + age + sex + fat_percentage"))
                                
                                for (tag in c("B_minimal", "C_fatpct")) {
                                        form <- switch(tag, B_minimal=B, C_fatpct=C_fat)
                                        fm <- prune_formula(form, sub)
                                        vars_needed <- all.vars(fm)[-1]
                                        if (n_complete(sub, c(outcome_var, vars_needed)) < 30) next
                                        
                                        rid <- paste0("F_strata_", lev, "_", tag)
                                        results[[rid]] <- run_model(
                                                formula = fm, data = sub,
                                                model_id = rid,
                                                group = "F",
                                                adjustment = tag,
                                                stratum_type = "sex",
                                                stratum_label = lev
                                        )
                                }
                        }
                }
        }
        
        # ------------------------ Collect Results ------------------------
        if (length(results) == 0) stop("No valid models could be run.")
        tidy_df   <- bind_rows(lapply(results, `[[`, "tidy"))
        glance_df <- bind_rows(lapply(results, `[[`, "glance"))
        tidy_df   <- add_partial_r2(tidy_df, glance_df)
        
        # ------------------------ Predictions (Robust) -------------------
        is_fast_model <- function(res) {
                trm  <- terms(res$model)
                labs <- attr(trm, "term.labels")
                ("fast" %in% labs) && !any(grepl(":", labs, fixed = TRUE))
        }
        base_for_plot <- Filter(is_fast_model, results)
        pred_plot <- NULL
        
        if (length(base_for_plot) > 0) {
                fast_range <- range(dat$fast, na.rm = TRUE)
                fast_grid <- tibble(fast = seq(fast_range[1], fast_range[2], length.out = 100))
                
                # --- ROBUST REFS ---
                refs <- dat %>% summarise(
                        age = if ("age" %in% names(.)) mean(age, na.rm = TRUE) else NA_real_,
                        age_centered = if ("age_centered" %in% names(.)) mean(age_centered, na.rm = TRUE) else NA_real_,
                        BMI = if ("BMI" %in% names(.)) mean(BMI, na.rm = TRUE) else NA_real_,
                        waist_hip_ratio = if ("waist_hip_ratio" %in% names(.)) mean(waist_hip_ratio, na.rm = TRUE) else NA_real_,
                        fat_percentage = if ("fat_percentage" %in% names(.)) mean(fat_percentage, na.rm = TRUE) else NA_real_,
                        sex = if ("sex" %in% names(.)) levels(sex)[1] else NA_character_,
                        disease_state = if ("disease_state" %in% names(.)) levels(disease_state)[1] else NA_character_
                )
                
                pred_df <- imap_dfr(base_for_plot, function(out, nm) {
                        nd <- fast_grid
                        preds <- predict_with_robust(out$model, nd, refs)
                        bind_cols(nd, preds) %>%
                                mutate(
                                        model = out$tidy$model_formula[1],
                                        lwr = fit - 1.96 * se,
                                        upr = fit + 1.96 * se
                                )
                })
                pred_plot <- ggplot(pred_df, aes(x = fast, y = fit, color = model, fill = model)) +
                        geom_point(data = dat, aes(x = fast, y = .data[[outcome_var]]),
                                   color = "gray40", alpha = 0.6, size = 2, inherit.aes = FALSE) +
                        geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA) +
                        geom_line(linewidth = 1.1) +
                        theme_bw() +
                        labs(title = paste(outcome_var, "vs Fast Fibers"), y = outcome_var)
        }
        
        # ------------------------ Outputs ------------------------
        readr::write_csv(tidy_df, file.path(dirs$results, paste0(sid, "_results_tidy.csv")))
        readr::write_csv(glance_df, file.path(dirs$results, paste0(sid, "_results_glance.csv")))
        
        qc_pdf <- file.path(dirs$figures, paste0(sid, "_QC_Plots.pdf"))
        grDevices::pdf(qc_pdf)
        for(r in results) { plot(r$model); mtext(deparse(formula(r$model)), outer=TRUE, line=-2) }
        if(!is.null(pred_plot)) print(pred_plot)
        dev.off()
        
        message("=== Done ===")
}

# Auto-run if executed directly (not sourced)
if (sys.nframe() == 0) run_all()
