library(specr)
library(tidyverse)

meta <- read_csv("/Users/maxullrich/Documents/GitHub/FiberTypeR-applied/data/testdata/sim_metadata.csv", show_col_types = FALSE)
fib  <- read_csv("/Users/maxullrich/Documents/GitHub/FiberTypeR-applied/data/testdata/sim_fibertype.csv", show_col_types = FALSE)

merged_data <- merge(meta %>% rename(RNAseq_ID = "GSE"), fib, by = "RNAseq_ID") %>% 
        mutate(fast = fast * 100,
               slow = slow *100) %>% 
        mutate(slow_z = scale(slow))

colnames(merged_data)



# Requires to keep full model
tidy_full <- function(x) {
        fit <- broom::tidy(x, conf.int = TRUE)
        fit$res <- list(x)  # Store model object
        return(fit)
}

# 1. Define setup with mapping
specs_setup <- setup(
        data = merged_data,
        x = "slow",
        y = "HOMA_IR",
        model = "lm",
        controls = c("BMI", "age","sex", "waist_hip_ratio", "fat_percentage"),
        # Format: "ColumnInData" = values_to_subset
        # We use the actual column names as the keys here
        subsets = list(disease_state = unique(merged_data$disease_state),
                       sex = unique(merged_data$sex)),
        fun1 = tidy_full
)

# 2. Rename the columns in the specs table BEFORE running specr
# This prevents the name collision during bootstrapping later
specs_setup$specs <- specs_setup$specs %>%
        rename(subset_disease = disease_state, 
               subset_sex = sex) %>%
        # Update the filter to use the new names
        filter(!(str_detect(subset_sex, "Male") & str_detect(controls, "sex"))) %>% 
        filter(!(str_detect(subset_sex, "Female") & str_detect(controls, "sex")))

# 3. Run models
# Now specr will run because it used the original names to subset, 
# but the output table will have the 'subset_' prefix to avoid collision.
results <- specr(specs_setup, keep_results = TRUE)


# Generate the curve
plot(results, choices = c("x", "y", "controls", "subsets")) + 
        theme_classic() +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())+
        labs(title = "Specification Curve: Fiber Type vs Insulin Resistance")


as.data.frame(results)

str(results)

# Convert to tibble first
results_tibble <- as_tibble(results)

# fit_r.squared is the r squared for slow in the model

################################################################################

# run bootstrapping

boot_models <- boot_null(results, specs_setup, n_samples = 1000)

summary(boot_models)
plot(boot_models)

