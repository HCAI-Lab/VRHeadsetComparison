## ===============================
## ART ANOVA + Posthocs + Effect Sizes
## ===============================
## install.packages(c("ARTool","readr","dplyr","tidyr","emmeans","rstatix","stringr"))
library(ARTool)
library(readr)
library(dplyr)
library(tidyr)
library(emmeans)
library(rstatix)
library(stringr)

# ===== 1) Load & prepare data =====
df <- read_csv("C:/Users/ITXPC/OneDrive - The Pennsylvania State University/Amir_abb6024/ICCMJournal/Codes/VRUsability/IgroupPresenceQuestionnaire.csv")

# Standardize HeadSet & core factors
df$HeadSet <- trimws(df$HeadSet)
df$HeadSet <- factor(df$HeadSet, levels = c("Oculus","HTC","HP"))  # Oculus as reference (level 1)
df$Gender  <- factor(df$Gender)

# Rename IPQ items to match your Python mapping
aw_col <- names(df)[startsWith(names(df), "How aware")][1]
rename_map <- c(
  'In the computer generated world I had a sense of "being there"' = 'GP',
  'Somehow I felt that the virtual world surrounded me.' = 'SP1',
  'I felt like I was just perceiving pictures.' = 'SP2',
  '(-) I did not feel present in the virtual space.' = 'SP3',
  'I had a sense of acting in the virtual space, rather than operating something outside in real world.' = 'SP4',
  'I felt present in the virtual space.' = 'SP5',
  'I was not aware of my real environment.' = 'INV2',
  'I still paid attention to the real environment.' = 'INV3',
  'How real did the virtual world seem to you?' = 'REAL1',
  'How much did your experience in the virtual environment seem consistent with your real world experience ?' = 'REAL2',
  'The virtual world seemed more realistic than the real world.' = 'REAL3'
)
names(df)[match(c(names(rename_map), aw_col), names(df))] <- c(unname(rename_map), 'INV1')

# ===== Reverse code =====
# If your scale is 1..5, proper reverse is 6 - x. If it is 0..4, use 4 - x.
rev5 <- function(x) ifelse(is.na(x), NA_real_, 6 - as.numeric(x))
df <- df %>%
  mutate(
    SP2_rev  = rev5(SP2),
    SP3_rev  = rev5(SP3),
    SP4_rev  = 6 - SP4,
    INV1_rev = rev5(INV1),
    INV3_rev = rev5(INV3)
  )

# ===== Subscales =====
df <- df %>%
  mutate(
    GeneralPresence = GP,
    SpatialPresence = rowMeans(across(c(GP,SP1, SP2_rev, SP3_rev, SP4_rev, SP5)), na.rm = TRUE),
    Involvement    = rowMeans(across(c(INV1_rev, INV2, INV3_rev)), na.rm = TRUE),
    Realism        = rowMeans(across(c(REAL1, REAL2, REAL3)), na.rm = TRUE)
  )

# Ensure subject is a factor (fixes your previous error)
df$UniqueID <- factor(df$UniqueID)

# Keep only rows needed
df <- df %>% filter(!is.na(UniqueID), !is.na(HeadSet), !is.na(Gender))

# ===== Helper: Cohen's dz from t (paired) =====
# dz = t / sqrt(n), where n is # of paired observations used in the comparison
cohen_dz_from_t <- function(t_value, n) t_value / sqrt(n)

# ===== Helper: run ART, posthocs vs Oculus, and effect sizes =====
run_art_with_effects <- function(dv) {
  cat("\n========================================\n")
  cat("DV:", dv, "\n")
  
  # ---- ART model (within: HeadSet; between: Gender; subject: UniqueID) ----
  form <- as.formula(paste0(dv, " ~ HeadSet * Gender + Error(UniqueID/HeadSet)"))
  m <- art(form, data = df)
  cat("\nART ANOVA (Type III-like on aligned responses):\n")
  print(anova(m))
  
  ## ---- Posthocs vs reference (Oculus) on aligned model ----
  # Main effect of HeadSet averaged over Gender
  m_headset_main <- artlm(m, "HeadSet")
  emm_main <- emmeans(m_headset_main, ~ HeadSet)
  pairs_vs_ref_main <- contrast(emm_main, method = "trt.vs.ctrl", ref = 1, adjust = "holm")
  cat("\nPost hoc: HeadSet vs Oculus (averaged over Gender), Holm-adjusted:\n")
  print(pairs_vs_ref_main)
  
  # Compute Cohen's dz from emmeans contrasts (paired design)
  # Estimate 'n' as the number of subjects with non-missing DV across HeadSet levels
  n_subj <- df %>%
    filter(!is.na(.data[[dv]])) %>%
    group_by(UniqueID) %>%
    summarize(n_levels = n_distinct(HeadSet)) %>%
    filter(n_levels == 3 | n_levels == 2) %>% # include subjects who saw at least the two levels in a given contrast
    nrow()
  
  # Extract t-ratios and compute dz
  tvals <- summary(pairs_vs_ref_main)$t.ratio
  comps <- as.character(summary(pairs_vs_ref_main)$contrast)
  dz    <- sapply(tvals, cohen_dz_from_t, n = n_subj)
  dz_tbl <- data.frame(contrast = comps, t_ratio = tvals, n_paired = n_subj, cohen_dz = dz)
  cat("\nCohen's dz (paired) from aligned-model t-ratios (approximate):\n")
  print(dz_tbl, row.names = FALSE, digits = 4)
  
  ## ---- Optional: HeadSet differences within each Gender ----
  m_hg <- artlm(m, "HeadSet:Gender")
  emm_hg <- emmeans(m_hg, ~ HeadSet | Gender)
  pairs_vs_ref_by_gender <- contrast(emm_hg, method = "trt.vs.ctrl", ref = 1, adjust = "holm")
  cat("\nPost hoc: HeadSet vs Oculus within each Gender, Holm-adjusted:\n")
  print(pairs_vs_ref_by_gender)
  
  # Effect sizes dz by gender (paired)
  by_gender <- summary(pairs_vs_ref_by_gender)
  # Estimate n per gender as subjects with >= 2 headset levels within that gender
  n_by_gender <- df %>%
    filter(!is.na(.data[[dv]])) %>%
    group_by(Gender, UniqueID) %>%
    summarize(n_levels = n_distinct(HeadSet), .groups = "drop") %>%
    filter(n_levels >= 2) %>%
    group_by(Gender) %>% summarise(n_paired = n(), .groups = "drop")
  dz_by_gender <- by_gender %>%
    select(Gender, contrast, t.ratio) %>%
    left_join(n_by_gender, by = "Gender") %>%
    mutate(cohen_dz = t.ratio / sqrt(n_paired))
  cat("\nCohen's dz (paired) within each Gender (approximate):\n")
  print(dz_by_gender, row.names = FALSE, digits = 4)
  
  ## ---- Robust nonparametric effect sizes on raw DV: rank-biserial (paired Wilcoxon) ----
  # Prepare wide form by subject (and optionally by gender) to do paired Wilcoxon: Oculus vs HTC, Oculus vs HP
  wide_all <- df %>%
    select(UniqueID, Gender, HeadSet, !!sym(dv)) %>%
    pivot_wider(names_from = HeadSet, values_from = !!sym(dv)) %>%
    filter(!is.na(Oculus))  # need reference present
  
  # Overall (ignoring Gender)
  rb_overall <- list()
  if ("HTC" %in% names(wide_all)) {
    d <- wide_all %>% filter(!is.na(HTC)) %>% mutate(diff = HTC - Oculus)
    if (nrow(d) > 1) {
      w <- wilcox_test(data = d, diff ~ 1, mu = 0, alternative = "two.sided")  # sign test style
      eff <- wilcox_effsize(data = d, diff ~ 1, ci = TRUE, ci.type = "perc")
      rb_overall[["HTC vs Oculus"]] <- list(test = w, effect = eff)
    }
  }
  if ("HP" %in% names(wide_all)) {
    d <- wide_all %>% filter(!is.na(HP)) %>% mutate(diff = HP - Oculus)
    if (nrow(d) > 1) {
      w <- wilcox_test(data = d, diff ~ 1, mu = 0, alternative = "two.sided")
      eff <- wilcox_effsize(data = d, diff ~ 1, ci = TRUE, ci.type = "perc")
      rb_overall[["HP vs Oculus"]] <- list(test = w, effect = eff)
    }
  }
  cat("\nRank-biserial (paired Wilcoxon) on RAW DV (overall):\n")
  print(lapply(rb_overall, function(x) list(test = x$test, effect = x$effect)))
  
  # By Gender
  cat("\nRank-biserial by Gender (paired Wilcoxon) on RAW DV:\n")
  for (g in levels(df$Gender)) {
    wg <- wide_all %>% filter(Gender == g)
    res_g <- list()
    if ("HTC" %in% names(wg)) {
      d <- wg %>% filter(!is.na(HTC)) %>% mutate(diff = HTC - Oculus)
      if (nrow(d) > 1) {
        w <- wilcox_test(data = d, diff ~ 1, mu = 0)
        eff <- wilcox_effsize(data = d, diff ~ 1, ci = TRUE, ci.type = "perc")
        res_g[["HTC vs Oculus"]] <- list(test = w, effect = eff)
      }
    }
    if ("HP" %in% names(wg)) {
      d <- wg %>% filter(!is.na(HP)) %>% mutate(diff = HP - Oculus)
      if (nrow(d) > 1) {
        w <- wilcox_test(data = d, diff ~ 1, mu = 0)
        eff <- wilcox_effsize(data = d, diff ~ 1, ci = TRUE, ci.type = "perc")
        res_g[["HP vs Oculus"]] <- list(test = w, effect = eff)
      }
    }
    cat(paste0("\n  Gender = ", g, ":\n"))
    print(lapply(res_g, function(x) list(test = x$test, effect = x$effect)))
  }
  
  invisible(m)
}

# ===== 2) Run for each subscale =====
m_GP  <- run_art_with_effects("GeneralPresence")
m_SP  <- run_art_with_effects("SpatialPresence")
m_INV <- run_art_with_effects("Involvement")
m_REA <- run_art_with_effects("Realism")
