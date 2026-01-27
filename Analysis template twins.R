
# =============================================================================
# Mixed-Effects Analysis: Muscle Fiber Type vs Insulin Resistance in MZ Twins
# =============================================================================

library(renv)
renv::restore(prompt = FALSE)

library(tidyverse)
library(lme4)        # mixed models
library(clubSandwich) # cluster-robust SE for mixed models
library(broom.mixed) # tidy mixed-model outputs
library(broom)        # glance merging for metadata
library(purrr)

# -----------------------------------------------------------------------------
# 1. INPUTS
# -----------------------------------------------------------------------------

study_ID       <- "PMID.............."
metadata_path  <- ".............."
fibertype_path <- ".............."
output_folder  <- ".............."


# ---- IMPORTANT FOR TWIN MODELS ----
col_twinpair  <- "twinpair_ID"      # MUST EXIST in metadata
col_id        <- "GSE"
col_HOMA_IR   <- "HOMA_IR"
col_bmi       <- "BMI"
col_age       <- "age"
col_sex       <- "sex"
col_fat       <- "fat_percentage"
col_whr       <- "waist_hip_ratio"
col_dis       <- "disease_state"

if(!dir.exists(output_folder)) dir.create(output_folder)

# -----------------------------------------------------------------------------
# 2. LOAD + CLEAN DATA
# -----------------------------------------------------------------------------

meta <- read_csv(metadata_path, show_col_types = FALSE)
fib  <- read_csv(fibertype_path, show_col_types = FALSE)

meta_clean <- meta %>%
        rename(
                RNAseq_ID = all_of(col_id),
                twinpair_ID = all_of(col_twinpair),
                HOMA_IR = all_of(col_HOMA_IR),
                BMI = all_of(col_bmi),
                age = all_of(col_age),
                sex = all_of(col_sex),
                fat_pct = all_of(col_fat),
                whr = all_of(col_whr),
                disease_state = all_of(col_dis)
        )

colnames(fib) <- c("RNAseq_ID","slow","fast")

df <- inner_join(meta_clean, fib, by = "RNAseq_ID") %>%
        mutate(
                fast = fast * 100,
                fast = as.numeric(fast),
                fast_z = scale(fast),
                HOMA_IR = as.numeric(HOMA_IR),
                BMI = as.numeric(BMI),
                age = as.numeric(age),
                age_centered = scale(age, center = TRUE, scale = FALSE),
                sex = factor(sex),
                disease_state = factor(disease_state),
                twinpair_ID = factor(twinpair_ID)
        )

message("Rows: ", nrow(df))
message("Twin pairs: ", n_distinct(df$twinpair_ID))

# -----------------------------------------------------------------------------
# 3. MODEL VARIABLE SET GENERATION
# -----------------------------------------------------------------------------

exclude_vars <- c("HOMA_IR","RNAseq_ID","slow","fast","fast_z","Condition",
                  "disease_state","twinpair_ID","age")

indep_vars <- setdiff(colnames(df), exclude_vars)

main_pred <- "fast"
obesity_metrics <- c("BMI","fat_pct","whr")
others <- setdiff(indep_vars, c(main_pred, obesity_metrics))

other_combos <- map(0:length(others), ~ combn(others, .x, simplify = FALSE)) %>% flatten()

final_var_sets <- list()

for (o_var in obesity_metrics) {
        for (others_set in other_combos) {
                final_var_sets[[length(final_var_sets)+1]] <- c(main_pred, o_var, others_set)
        }
}

for (others_set in other_combos) {
        final_var_sets[[length(final_var_sets)+1]] <- c(main_pred, others_set)
}

#######################################################################################
# --------------------- check if sex or disease state should be checked --------------

# --------------------- Robust Interaction Checks ---------------------

# 1. Check for Sex Interaction (requires column existence AND n >= 30)
if ("sex" %in% colnames(df)) {
        sex_counts <- df %>% count(sex)
        # Ensure there are at least 2 levels and each has n >= 30
        if (nrow(sex_counts) >= 2 && min(sex_counts$n, na.rm = TRUE) >= 30) {
                message("Adding Sex interaction model...")
                final_var_sets[[length(final_var_sets) + 1]] <- c("fast * sex", "age_centered")
        } else {
                message("Skipping Sex interaction: insufficient n or levels.")
        }
}

# 2. Check for Disease State Interaction (requires column existence AND n >= 20)
if ("disease_state" %in% colnames(df)) {
        disease_counts <- df %>% count(disease_state)
        # Ensure there are at least 2 levels and each has n >= 20
        if (nrow(disease_counts) >= 2 && min(disease_counts$n, na.rm = TRUE) >= 20) {
                message("Adding Disease State interaction model...")
                # Note: including 'sex' here only if sex exists
                base_vars <- if("sex" %in% colnames(df)) c("fast * disease_state", "sex", "age_centered") else c("fast * disease_state", "age_centered")
                final_var_sets[[length(final_var_sets) + 1]] <- base_vars
        } else {
                message("Skipping Disease State interaction: insufficient n or levels.")
        }
}







# -----------------------------------------------------------------------------
# 4. Mixed Model Helper Function
# -----------------------------------------------------------------------------



get_mixed_results <- function(model_obj, cluster_id) {
        robust_vcov <- clubSandwich::vcovCR(model_obj, type = "CR2", cluster = cluster_id)
        coef_names  <- names(lme4::fixef(model_obj))
        robust_vcov <- robust_vcov[coef_names, coef_names, drop = FALSE]
        
        std_err <- sqrt(diag(robust_vcov))
        
        tibble::tibble(
                term      = coef_names,
                estimate  = lme4::fixef(model_obj),
                std.error = std_err,
                statistic = estimate / std.error,
                p.value   = 2 * stats::pt(-abs(statistic), df = stats::df.residual(model_obj)),
                conf.low  = estimate - 1.96 * std.error,
                conf.high = estimate + 1.96 * std.error
        )
}




is_valid_model <- function(vars, df) {
        
        # Only keep predictor names that exist
        vars <- vars[vars %in% colnames(df)]
        
        # At least 1 predictor must be present
        if (length(vars) < 1) return(FALSE)
        
        for (v in vars) {
                
                x <- df[[v]]
                
                # 1. Factor with only 1 unique level → invalid model
                if (is.factor(x)) {
                        if (nlevels(x) <= 1) return(FALSE)
                        next
                }
                
                # 2. Character should not be included in models at all
                if (is.character(x)) return(FALSE)
                
                # 3. Numeric: check variance
                if (is.numeric(x)) {
                        if (var(x, na.rm = TRUE) == 0) return(FALSE)
                        next
                }
                
                # 4. If any other weird class appears → skip model
                return(FALSE)
        }
        
        return(TRUE)
}


# Identify which random-effect variance components are ~0 (boundary)
get_singular_terms <- function(model, tol = 1e-8) {
        # VarCorr -> data.frame with columns: grp, var1, var2, vcov, sdcor
        # Diagonal entries (var2 = NA) correspond to variances (sdcor = std dev)
        vc <- as.data.frame(VarCorr(model))
        if (nrow(vc) == 0) return(NA_character_)
        
        diag_terms <- vc %>%
                dplyr::filter(is.na(var2)) %>%             # variances only
                dplyr::mutate(
                        # term label: (var1 | grp) e.g., (Intercept | twinpair_ID) or (fast | twinpair_ID)
                        re_term = paste0("(", var1, " | ", grp, ")"),
                        near_zero = sdcor < tol
                )
        
        # Return a semicolon-separated list of near-zero variance components
        nz <- diag_terms %>% dplyr::filter(near_zero) %>% dplyr::pull(re_term)
        if (length(nz) == 0) return(NA_character_) else return(paste(nz, collapse = "; "))
}



# -----------------------------------------------------------------------------
# 5. Run All Mixed Models
# -----------------------------------------------------------------------------



run_mixed <- function(vars, df) {
        
        f <- as.formula(
                paste0("HOMA_IR ~ ", paste(vars, collapse = " + "), " + (1 | twinpair_ID)")
        )
        
        mod <- lmer(f, data = df, REML = FALSE)
        
        # Singular diagnostics
        singular_flag  <- isSingular(mod, tol = 1e-4)
        singular_terms <- get_singular_terms(mod, tol = 1e-8)
        
        # Tidy with cluster-robust SE (your safe version)
        tidy_res <- get_mixed_results(mod, df$twinpair_ID) %>%
                dplyr::mutate(
                        model.formula = deparse1(f),
                        singular_fit  = singular_flag,
                        singular_terms = ifelse(singular_flag, singular_terms, NA_character_)
                )
        
        glance_res <- broom.mixed::glance(mod) %>%
                dplyr::mutate(
                        model.formula = deparse1(f),
                        singular_fit  = singular_flag,
                        singular_terms = ifelse(singular_flag, singular_terms, NA_character_)
                )
        
        list(tidy = tidy_res, glance = glance_res)
}




valid_var_sets <- keep(final_var_sets, ~ is_valid_model(.x, df))
all_model_outputs <- map(valid_var_sets, ~ run_mixed(.x, df))
all_model_outputs_z <- map(
        map(final_var_sets, ~ gsub("fast", "fast_z", .x)),
        ~ run_mixed(.x, df)
)

final_tidy_df <- map_df(all_model_outputs, "tidy") %>%
        bind_rows(map_df(all_model_outputs_z, "tidy"))

final_glance_df <- bind_rows(
        map_df(all_model_outputs, "glance"),
        map_df(all_model_outputs_z, "glance")
) %>%
        mutate(study_ID = study_ID) %>%
        relocate(study_ID)




# -----------------------------------------------------------------------------
# 6. DESCRIPTIVE STATISTICS 
# -----------------------------------------------------------------------------

target_vars <- c(
        "fast", "BMI", "HOMA_IR", "M_value", "Matsuda_index",
        "age", "whr", "fat_percentage"
)

factor_cols <- names(df)[vapply(df, is.factor, logical(1))]

# Exclude IDs (twinpair_ID creates 2-row groups)
group_factors <- setdiff(factor_cols, c("twinpair_ID", "RNAseq_ID"))

# Only numeric vars that exist
existing_num <- intersect(target_vars, names(df))
existing_num <- existing_num[vapply(df[existing_num], is.numeric, logical(1))]


.desc_block <- function(data, group_vars, vars) {
        
        if (length(vars) == 0) {
                warning("No numeric variables available for descriptives.")
                return(tibble())
        }
        
        # GROUPING (allowed selecting context)
        if (length(group_vars) > 0) {
                data <- data %>% dplyr::group_by(dplyr::across(dplyr::all_of(group_vars)))
        }
        
        # SUMMARISE
        out <- data %>%
                dplyr::summarise(
                        dplyr::across(
                                dplyr::all_of(vars),
                                .fns = list(
                                        n = ~ sum(!is.na(.)),
                                        mean = ~ mean(., na.rm = TRUE),
                                        sd = ~ sd(., na.rm = TRUE),
                                        median = ~ median(., na.rm = TRUE)
                                ),
                                .names = "{.col}_{.fn}"
                        ),
                        .groups = "drop"
                )
        
        if (nrow(out) == 0) return(tibble())
        
        # Determine columns to pivot — plain vector OK, but inside pivot_longer we must wrap in all_of()
        cols_to_pivot <- setdiff(names(out), group_vars)
        
        out %>%
                tidyr::pivot_longer(
                        cols = dplyr::all_of(cols_to_pivot),
                        names_to = c("variable", "stat"),
                        names_pattern = "^(.*)_(n|mean|sd|median)$"
                ) %>%
                tidyr::pivot_wider(names_from = stat, values_from = value) %>%
                dplyr::mutate(study_ID = study_ID) %>%
                dplyr::relocate(study_ID)
}


# A) Grouped descriptives
study_descriptives <- .desc_block(df, group_factors, existing_num)

# B) Overall descriptives (no grouping)
overall_descriptives <- .desc_block(df, character(0), existing_num) %>%
        dplyr::mutate(grouping = "Overall") %>%
        dplyr::relocate(grouping, .after = study_ID)

# -----------------------------------------------------------------------------
# 7.   Add partial r^2
# -----------------------------------------------------------------------------

# Build a partial R^2 table per model formula (term-specific df)
.compute_partial_r2_for_formula <- function(form_str, data, cluster) {
        f <- as.formula(form_str)
        mod <- lme4::lmer(f, data = data, REML = FALSE)
        
        Vcr2 <- clubSandwich::vcovCR(mod, type = "CR2", cluster = cluster)
        ct   <- clubSandwich::coef_test(mod, vcov = Vcr2, test = "Satterthwaite")
        term_names <- rownames(ct)
        
        tibble::tibble(
                model.formula = deparse1(f),
                term          = term_names,
                t_stat        = ct[,"tstat"],
                df_Satt       = ct[,"df_Satt"]
        ) %>%
                dplyr::mutate(
                        partial_r2 = dplyr::if_else(
                                term == "(Intercept)",
                                as.numeric(NA),
                                (t_stat^2) / (t_stat^2 + df_Satt)
                        )
                ) %>%
                dplyr::select(model.formula, term, partial_r2)
}

# Compute for all unique models and join
unique_formulas <- unique(final_glance_df$model.formula)

partial_r2_table <- purrr::map_dfr(
        unique_formulas,
        ~ .compute_partial_r2_for_formula(.x, data = df, cluster = df$twinpair_ID)
)

final_tidy_df <- final_tidy_df %>%
        dplyr::left_join(partial_r2_table, by = c("model.formula", "term"))



# -----------------------------------------------------------------------------
# 7. Export
# -----------------------------------------------------------------------------

write.csv(final_tidy_df,   paste0(output_folder, "/", study_ID, "_tidy_df.csv"))
write.csv(final_glance_df, paste0(output_folder, "/", study_ID, "_glance_df.csv"))


# -----------------------------------------------------------------------------
# 8. QC PLOTS (PDF) FOR ALL MODELS
# -----------------------------------------------------------------------------

message("Generating QC plots...")

# Helper: safe QC plotting for an lmer model


.plot_qc_for_model <- function(mod, f, df, study_ID, singular_flag) {
        # Extract residuals and fitted values
        res <- resid(mod)     # conditional residuals
        fit <- fitted(mod)    # conditional fitted
        obs <- getME(mod, "y")
        
        # Build a neat, wrapped title
        header_line <- paste0("Study: ", study_ID, "   |   Singular fit: ", singular_flag)
        # Wrap the formula to avoid overlap; adjust width to taste (80–100)
        formula_line <- paste0("Model: ", deparse1(f))
        wrapped_formula <- paste(strwrap(formula_line, width = 95), collapse = "\n")
        
        # Save current par and set margins
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar), add = TRUE)
        
        # Larger outer margin on top (oma[3]) for multi-line header
        # Slightly roomier inner margins so axes labels don't collide
        par(mfrow = c(2, 2),
            mar = c(4.6, 4.8, 2.4, 1.2),
            oma = c(0.6, 0.6, 5.8, 0.6))  # top outer margin increased
        
        ## 1) Residuals vs Fitted
        plot(fit, res,
             pch = 19, col = rgb(0, 0, 0, 0.45),
             xlab = "Fitted values",
             ylab = "Residuals",
             main = "Residuals vs Fitted")
        abline(h = 0, col = "red", lty = 2, lwd = 1.3)
        
        ## 2) Normal Q-Q of residuals
        qqnorm(res, pch = 19, col = rgb(0, 0, 0, 0.45),
               main = "Normal Q-Q (residuals)")
        qqline(res, col = "red", lty = 2, lwd = 1.3)
        
        ## 3) Histogram of residuals
        hist(res, breaks = "FD", col = "grey80", border = "white",
             main = "Histogram of residuals", xlab = "Residuals")
        
        ## 4) Observed vs Fitted
        plot(fit, obs,
             pch = 19, col = rgb(0, 0, 0, 0.45),
             xlab = "Fitted values",
             ylab = "Observed HOMA_IR",
             main = "Observed vs Fitted")
        abline(a = 0, b = 1, col = "blue", lty = 2, lwd = 1.3)
        
        ## Outer titles: one for header, one for wrapped formula
        mtext(header_line, outer = TRUE, cex = 0.9, line = 2.6, font = 2)
        # The wrapped formula may occupy multiple lines; place it just below the header
        mtext(wrapped_formula, outer = TRUE, cex = 0.85, line = 0.9)
}

        

# Build the sets that were actually used in modeling
valid_var_sets_z <- purrr::map(valid_var_sets, ~ gsub("fast", "fast_z", .x))

# Prepare PDF
qc_pdf_path <- file.path(output_folder, paste0(study_ID, "_model_QC.pdf"))
pdf(qc_pdf_path, width = 10, height = 8)

# Loop over raw 'fast' models
for (vars in valid_var_sets) {
        f <- as.formula(
                paste0("HOMA_IR ~ ", paste(vars, collapse = ' + '), " + (1 | twinpair_ID)")
        )
        # Fit and plot with safety
        try({
                mod <- lme4::lmer(f, data = df, REML = FALSE)
                singular_flag <- lme4::isSingular(mod, tol = 1e-4)
                .plot_qc_for_model(mod, f, df, study_ID, singular_flag)
        }, silent = TRUE)
}

# Loop over z-scored 'fast_z' models
for (vars in valid_var_sets_z) {
        f <- as.formula(
                paste0("HOMA_IR ~ ", paste(vars, collapse = ' + '), " + (1 | twinpair_ID)")
        )
        # Fit and plot with safety
        try({
                mod <- lme4::lmer(f, data = df, REML = FALSE)
                singular_flag <- lme4::isSingular(mod, tol = 1e-4)
                .plot_qc_for_model(mod, f, df, study_ID, singular_flag)
        }, silent = TRUE)
}

dev.off()

message("QC PDF written to: ", qc_pdf_path)

