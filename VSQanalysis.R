# ===============================
# SSQ: ART ANOVA + Pairwise Post-Hoc with Cohen's dz
# ===============================

# --- 1) Install and Load Packages ---
# Ensure all necessary packages are installed
pkgs <- c("ARTool", "dplyr", "stringr", "readr", "emmeans", "tidyr", "purrr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(ARTool)
library(dplyr)
library(stringr)
library(readr)
library(emmeans)
library(tidyr)
library(purrr)

# --- 2) Load and Prepare Data ---
# Adjust the file path as needed
# Make sure to use forward slashes "/" in your path
file_path <- "C:/Users/ITXPC/OneDrive - The Pennsylvania State University/Amir_abb6024/ICCMJournal/Codes/VRUsability/SimulatorSicknessQuestionnaire.csv"
df <- read_csv(file_path, show_col_types = FALSE)

# Basic cleaning: trim spaces in column names and key columns
names(df) <- names(df) |> str_replace_all("  +", " ") |> str_trim()
df <- df |>
  mutate(
    HeadSet  = str_trim(HeadSet),
    Gender   = str_trim(Gender),
    UniqueID = as.factor(UniqueID) # IMPORTANT for Error() term
  )

# Set factor levels (Oculus will be the reference group)
df <- df |>
  mutate(
    HeadSet = factor(HeadSet, levels = c("Oculus", "HTC", "HP")),
    Gender  = factor(Gender)
  )

# --- 3) Define SSQ Items & Compute Scales ---
nausea_items <- c("General discomfort (physical)", "Increased salivation", "Sweating", "Nausea", "Stomach awareness (a feeling of discomfort that is just short of nausea.)", "Burping")
oculomotor_items <- c("General discomfort (physical)", "Fatigue", "Headache", "Eyestrain", "Difficulty focusing", "Difficulty concentrating", "Blurred vision")
disorientation_items <- c("Difficulty focusing", "Nausea", "Fullness of head", "Blurred vision", "Dizzy (eyes open)", "Dizzy (eyes closed)", "Vertigo (A sense of spinning experienced even when someone is perfectly still.)")

# Clean item names to match column headers
clean_names <- function(v) v |> str_replace_all("  +", " ") |> str_trim()
nausea_items         <- clean_names(nausea_items)
oculomotor_items     <- clean_names(oculomotor_items)
disorientation_items <- clean_names(disorientation_items)

all_items <- unique(c(nausea_items, oculomotor_items, disorientation_items))
present_items <- intersect(all_items, names(df))

# Ensure item columns are numeric
df <- df |>
  mutate(across(all_of(present_items), as.numeric))

# Function to safely sum subscales
sum_if_present <- function(d, items) {
  items <- intersect(items, names(d))
  if (length(items) == 0) return(rep(NA_real_, nrow(d)))
  rowSums(d[, items, drop = FALSE], na.rm = TRUE)
}

# Calculate subscale scores and weighted scores
df <- df |>
  mutate(
    Nausea         = sum_if_present(cur_data_all(), nausea_items),
    Oculomotor     = sum_if_present(cur_data_all(), oculomotor_items),
    Disorientation = sum_if_present(cur_data_all(), disorientation_items),
    Total          = Nausea + Oculomotor + Disorientation,
    Nausea_Weighted         = Nausea * 9.54,
    Oculomotor_Weighted     = Oculomotor * 7.58,
    Disorientation_Weighted = Disorientation * 13.92,
    Total_Weighted          = (Nausea + Oculomotor + Disorientation) * 3.74
  )

# --- 4) Helper Functions for Post-Hoc Analysis ---

#' Calculate Cohen's dz from a t-statistic and number of pairs.
cohen_dz_from_t <- function(t_value, n) {
  if (n <= 1) return(NA)
  return(t_value / sqrt(n))
}


#' Run ART ANOVA and perform pairwise post-hoc tests for the HeadSet factor.
#' This function calculates p-values and Cohen's dz for all pairwise comparisons.
#'
#' @param dv The name of the dependent variable (string).
#' @param data The dataframe containing the data.
#' @return Invisibly returns the ART model object.
run_art_posthoc_dz <- function(dv, data) {
  cat("\n", paste0(rep("=", 40), collapse = ""), "\n")
  cat("  ANALYSIS FOR DEPENDENT VARIABLE:", dv, "\n")
  cat(paste0(rep("=", 40), collapse = ""), "\n")
  
  # Filter data to ensure no missing values for the DV or factors
  d_filtered <- data %>%
    select(UniqueID, HeadSet, Gender, all_of(dv)) %>%
    drop_na()
  
  if (nrow(d_filtered) < 2) {
    cat("Not enough complete observations to run analysis for", dv, "\n")
    return(invisible(NULL))
  }
  
  # --- Run ART ANOVA ---
  form <- as.formula(paste0(dv, " ~ HeadSet*Gender + Error(UniqueID/HeadSet)"))
  m <- art(form, data = d_filtered)
  
  cat("\n--- ART ANOVA Summary ---\n")
  print(anova(m))
  
  # --- Post-Hoc Analysis for HeadSet ---
  cat("\n--- Post-Hoc: Pairwise Comparisons for HeadSet (Holm-adjusted) ---\n")
  
  # artlm() extracts the linear model component for post-hoc tests
  m_headset <- artlm(m, "HeadSet")
  
  # emmeans() computes estimated marginal means
  emm_headset <- emmeans(m_headset, ~ HeadSet)
  
  # *** THIS IS THE MODIFIED PART ***
  # Perform all pairwise comparisons instead of comparing to a reference.
  # The 'pairs()' function is a convenient wrapper for this.
  pairwise_comps <- pairs(emm_headset, adjust = "holm")
  summary_pairs <- summary(pairwise_comps, infer = TRUE)
  
  if (nrow(summary_pairs) > 0) {
    # --- Calculate Cohen's dz ---
    # Find 'n' - the number of subjects with data for at least two HeadSet conditions.
    n_paired <- d_filtered %>%
      group_by(UniqueID) %>%
      summarize(n_levels = n_distinct(HeadSet), .groups = "drop") %>%
      filter(n_levels >= 2) %>%
      nrow()
    
    if (n_paired > 0) {
      cat("Number of subjects in paired comparisons (n):", n_paired, "\n\n")
      # Calculate dz for each contrast and add it to the results table
      summary_pairs$cohen_dz <- cohen_dz_from_t(summary_pairs$t.ratio, n_paired)
      
      # Print the final table with p-values and effect sizes
      print(summary_pairs, digits = 4)
      
    } else {
      cat("Warning: Could not determine the number of paired subjects.\n")
      print(summary_pairs, digits = 4)
    }
  } else {
    cat("No contrasts to compute for HeadSet.\n")
  }
  
  cat("\nInterpretation Guide for Cohen's dz:\n")
  cat("  |dz| < 0.2  : Negligible effect\n")
  cat("  |dz| ~ 0.2  : Small effect\n")
  cat("  |dz| ~ 0.5  : Medium effect\n")
  cat("  |dz| ~ 0.8+ : Large effect\n")
  
  invisible(m)
}


# --- 5) Run the Analysis for Each Weighted Score ---
# List of the dependent variables you want to analyze
weighted_scores <- c("Nausea_Weighted", "Oculomotor_Weighted", "Disorientation_Weighted", "Total_Weighted")

# Loop through the scores and run the analysis for each one
# Note: The `ref_headset` argument is no longer needed.
results <- map(weighted_scores, ~ run_art_posthoc_dz(.x, data = df))
names(results) <- weighted_scores

cat("\n\n--- Analysis Complete ---\n")
