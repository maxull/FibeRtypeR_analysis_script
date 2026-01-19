# =============================================================================
# Manual Analysis: Muscle Fiber Type vs Insulin Resistance
# =============================================================================

# -----------------------------------------------------------------------------
# 1. SETUP & CONFIGURATION
# -----------------------------------------------------------------------------
renv::restore(prompt = FALSE)
library(tidyverse)
library(sandwich) # For robust SE
library(lmtest)   # For robust SE
library(broom)    # For tidy results

# --- INPUTS ---
study_ID       <- "PMID40683252"
metadata_path  <- "/Users/maxullrich/Documents/GitHub/FiberTypeR-applied/data/testdata/sim_metadata.csv"
fibertype_path <- "/Users/maxullrich/Documents/GitHub/FiberTypeR-applied/data/testdata/sim_fibertype.csv"
output_folder  <- "SIMDATA_output_manual"

# --- COLUMN MAPPING (Change the RIGHT side to match your CSV headers) ---
col_id      <- "GSE"                    # ID column in metadata
col_HOMA_IR <- "HOMA_IR"                # HOMA_IR variable (e.g., HOMA_IR)
col_bmi     <- "BMI"                    # BMI
col_age     <- "age"                    # Age
col_sex     <- "sex"                    # Sex
col_fat     <- "fat_percentage"         # Fat % (optional)
col_whr     <- "waist_hip_ratio"
col_dis     <- "disease_state"

# --- Create Output Directory ---
if(!dir.exists(output_folder)) dir.create(output_folder)

# -----------------------------------------------------------------------------
# 2. DATA PREPARATION
# -----------------------------------------------------------------------------
message("Loading data...")

# 1. Read files
meta <- read_csv(metadata_path, show_col_types = FALSE)
fib  <- read_csv(fibertype_path, show_col_types = FALSE)

# 2. Standardize Names (Manual Step)
# We rename columns here so the formulas below are easy to read
meta_clean <- meta %>%
        rename(
                RNAseq_ID = all_of(col_id),
                HOMA_IR = all_of(col_HOMA_IR),
                age = all_of(col_age),
                sex = all_of(col_sex),
                BMI = all_of(col_bmi),
                fat_pct = all_of(col_fat),
                whr = all_of(col_whr),
                disease_state = all_of(col_dis)
                # Add others here if needed: e.g., fat_pct = all_of(col_fat)
        )

# 3. Clean Fibertype data (Rename columns to standard 'id', 'fast', 'slow')
#    Adjust "RNAseq_ID", "fast", "slow" below if your raw file differs slightly
colnames(fib) <- c("RNAseq_ID", "slow", "fast") 

# 4. Merge
df <- inner_join(meta_clean, fib, by = "RNAseq_ID") %>%
        mutate(
                # Convert to numeric to be safe
                HOMA_IR = as.numeric(HOMA_IR),
                fast    = as.numeric(fast),
                age     = as.numeric(age),
                BMI     = as.numeric(BMI),
                # Ensure sex is a factor (categorical)
                sex     = factor(sex),
                disease_state = factor(disease_state)
        ) %>% 
        mutate(fast = fast*100,
               age_centered = scale(age, center = T, scale = F))

message("Data loaded. Rows: ", nrow(df))
message("Unique sexes in data: ", paste(unique(df$sex), collapse=", "))
message("Unique disease states in data: ", paste(unique(df$disease_state), collapse=", "))

str(df)
# Define columns to exclude
exclude_vars <- c("HOMA_IR", "RNAseq_ID", "age","slow","fast", "Condition", "disease_state") 

# Get all column names except the excluded ones
indep_vars <- setdiff(colnames(df), exclude_vars)

# 1. Define your groups
main_pred <- "fast"
obesity_metrics <- c("BMI", "fat_pct", "whr")
# Other variables that can always be included (age, sex)
others <- setdiff(indep_vars, c(main_pred, obesity_metrics))

# 2. Create combinations for the "other" variables (age, sex) 
# including an empty set (0) so we can have models with just 'fast' and 'BMI'
other_combos <- map(0:length(others), ~ combn(others, .x, simplify = FALSE)) %>% flatten()

# 3. Build the final list of variable sets
final_var_sets <- list()

for (o_var in obesity_metrics) {
        for (others_set in other_combos) {
                # Combine: fast + ONE obesity metric + any combination of others
                final_var_sets[[length(final_var_sets) + 1]] <- c(main_pred, o_var, others_set)
        }
}

# Add models that have 'fast' and 'others' but NO obesity metrics at all
for (others_set in other_combos) {
        final_var_sets[[length(final_var_sets) + 1]] <- c(main_pred, others_set)
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

#######################################################################################

# Helper function to get Robust HC3 Errors (Standard in academic econ/epi)
get_results <- function(model_obj) {
        # Calculate Robust Standard Errors
        robust_vcov <- vcovHC(model_obj, type = "HC3")
        
        # Get tidy table
        broom::tidy(coeftest(model_obj, vcov. = robust_vcov)) %>%
                mutate(
                        conf.low = estimate - 1.96 * std.error,
                        conf.high = estimate + 1.96 * std.error,
                        model.formula = deparse1(formula(model_obj))
                        
                )
}

# 4. Run the models
all_model_outputs <- map(final_var_sets, function(vars) {
        
        form <- as.formula(paste("HOMA_IR ~", paste(vars, collapse = " + ")))
        mod <- lm(form, data = df)
        
        # 1. Get coefficient results (using your get_results function)
        tidy_res <- get_results(mod)
        
        # 2. Get model-level stats (glance) and include formula
        glance_res <- broom::glance(mod) %>% 
                mutate(model.formula = deparse1(formula(mod)))
        
        # Return both as a list
        list(tidy = tidy_res, glance = glance_res)
})

# 5. Run the models fast_z

df <- df %>% mutate(fast_z = scale(fast))

final_var_sets_z <- map(final_var_sets, ~ gsub("fast", "fast_z", .x))

all_model_outputs_z <- map(final_var_sets_z, function(vars) {
        
        form <- as.formula(paste("HOMA_IR ~", paste(vars, collapse = " + ")))
        mod <- lm(form, data = df)
        
        # 1. Get coefficient results (using your get_results function)
        tidy_res <- get_results(mod)
        
        # 2. Get model-level stats (glance) and include formula
        glance_res <- broom::glance(mod) %>% 
                mutate(model.formula = deparse1(formula(mod)))
        
        # Return both as a list
        list(tidy = tidy_res, glance = glance_res)
})



# 5. Extract into two separate data frames
# Combine all coefficient tables
final_tidy_df <- map_df(all_model_outputs, "tidy")
final_tidy_df <- bind_rows(final_tidy_df, map_df(all_model_outputs_z, "tidy"))


# Combine all model-fit (glance) tables
final_glance_df <- map_df(all_model_outputs, "glance")
final_glance_df <- bind_rows(final_glance_df, map_df(all_model_outputs_z, "glance"))

# Friendlier glance names
final_glance_df <- final_glance_df %>%
        rename(
                Residual.Standard.Error = sigma,
                model.p.value = p.value,
                Sum.of.squared.residuals = deviance
        ) %>% 
        mutate(study_ID = study_ID) %>% 
        relocate(study_ID)

# add partial r^2

final_tidy_df <- left_join(final_tidy_df, final_glance_df %>% dplyr::select(df.residual, model.formula), by = "model.formula") %>% 
        mutate(partial_r2 = ifelse(term == "(Intercept)", "NA", round((statistic^2) / (statistic^2 + df.residual), 4))) %>% 
        mutate(study_ID = study_ID) %>% 
        relocate(study_ID)



# ------------------------ Descriptive statistics ----------------------------


target_vars <- c("fast", "BMI", "HOMA_IR", "M_value", "Matsuda_index",
                 "age", "sex", "whr", "fat_percentage", "disease_state")

factor_cols <- names(df)[sapply(df, is.factor)]

existing <- intersect(setdiff(target_vars, factor_cols), names(df))

study_descriptives <- df %>%
        group_by(across(all_of(factor_cols))) %>% 
        summarise(across(all_of(existing), list(
                n      = ~ sum(!is.na(.)),
                mean   = ~ if (is.numeric(.)) mean(., na.rm = TRUE),
                sd   = ~ if (is.numeric(.)) sd(., na.rm = TRUE),
                median   = ~ if (is.numeric(.)) median(., na.rm = TRUE)
        )), .groups = "drop") %>%
        pivot_longer(
                cols = -all_of(factor_cols),
                names_to = c("variable", "stat"),
                names_pattern = "(.*)_(.*)"
        ) %>%
        pivot_wider(names_from = stat, values_from = value) %>% 
        mutate(study_ID = study_ID) %>% 
        relocate(study_ID)




# ------------------------ Save data ----------------------------

write.csv(final_glance_df, file = paste0(output_folder,"/", study_ID, "_glance_df.csv"))
write.csv(final_tidy_df, file = paste0(output_folder,"/", study_ID, "_tidy_df.csv"))
write.csv(study_descriptives, file = paste0(output_folder,"/", study_ID, "_study_descriptives.csv"))


# -----------------------------------------------------------------------------
# 6. MODEL QC PLOTS (PDF EXPORT)
# -----------------------------------------------------------------------------
message("Generating QC plots...")

qc_pdf_path <- file.path(output_folder, paste0(study_ID, "_model_QC.pdf"))

pdf(qc_pdf_path, width = 10, height = 8)

# Combine both sets of variables to plot all models run
all_sets_to_plot <- c(final_var_sets, final_var_sets_z)

for (vars in all_sets_to_plot) {
        # Re-fit the model (or you could store 'mod' in the map function earlier)
        form <- as.formula(paste("HOMA_IR ~", paste(vars, collapse = " + ")))
        mod  <- lm(form, data = df)
        
        # Set up a 2x2 grid for the 4 standard plots
        par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
        
        # Generate plots
        plot(mod, sub.caption = "")
        
        # Add the formula as a title across the top
        mtext(paste("Study:", study_ID, "Model:", deparse1(formula(mod))), outer = TRUE, cex = 0.8)
}

dev.off() # Close the PDF device

