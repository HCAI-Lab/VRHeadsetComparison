## ===============================
## ART ANOVA on VRSUQ Dimensions + Tukey ART-C + Cronbach Alpha
## ===============================

# install.packages(c("ARTool","readr","dplyr","tidyr","emmeans","rstatix","stringr","psych"))
library(ARTool)
library(readr)
library(dplyr)
library(tidyr)
library(rstatix)
library(stringr)
library(psych)

# ===== 1) Load & clean data =====
df <- read_csv("C:/Users/goatp/OneDrive - The Pennsylvania State University/Amir_abb6024/ICCMJournal/Codes/VRUsability/VRSUQ.csv")

# Map Likert responses
likert_map <- c(
  "strongly disagree" = 1,
  "disagree" = 2,
  "neutral" = 3,
  "agree" = 4,
  "strongly agree" = 5
)

item_cols <- names(df)[4:(ncol(df)-1)]  # all items between HeadSet..Gender
for (col in item_cols) {
  df[[col]] <- tolower(trimws(df[[col]]))
  df[[col]] <- unname(likert_map[df[[col]]])
}

# Reverse-code negative items
neg_cols <- grep("^\\(-\\)", item_cols, value = TRUE)
for (col in neg_cols) {
  df[[col]] <- 6 - df[[col]]
}

# Factor coding
df$HeadSet  <- factor(trimws(df$HeadSet), levels = c("Oculus","HTC","HP"))
df$Gender   <- factor(trimws(df$Gender))
df$UniqueID <- factor(df$UniqueID)

# ===== 2) Define VRSUQ dimensions =====
# Adjust item grouping if needed to match your scale definition
efficiency_items  <- c(
  "The system responded well to my manipulations as expected with no delays.",
  "I think this virtual reality system provides clear feedback on my manipulations.",
  "(-) I kept making errors/mistakes while using the virtual reality system.",
  "I could clearly understand the information presented within the virtual environment.",
  "I think it is easy to correct errors made during virtual reality experiences."
)
satisfaction_items <- c(
  "I think this system is user-friendly, straightforward to learn, and designed in such a way that most people will find it easy to adapt to.",
  "I enjoyed the virtual reality experience."
)
effectiveness_items <- c(
  "(-) I felt dizzy, motion sickness, or a headache while experiencing virtual reality.",
  "(-) While experiencing virtual reality, I felt mental burdens such as tension, frustration, and time pressure."
)

df <- df %>%
  mutate(
    Efficiency   = rowMeans(across(all_of(efficiency_items)), na.rm = TRUE),
    Satisfaction = rowMeans(across(all_of(satisfaction_items)), na.rm = TRUE),
    Effectiveness = rowMeans(across(all_of(effectiveness_items)), na.rm = TRUE),
    Overall      = rowMeans(across(all_of(item_cols)), na.rm = TRUE)
  )

# ===== 3) Cronbach’s alpha =====
cat("\n--- Cronbach's alpha for subscales ---\n")
print(psych::alpha(df[efficiency_items])$total)
print(psych::alpha(df[satisfaction_items])$total)
print(psych::alpha(df[effectiveness_items])$total)
print(psych::alpha(df[item_cols])$total)  # overall

# ===== 4) Helper function: run ART + Tukey =====
run_art <- function(dv, data=df) {
  cat("\n========================================\n")
  cat("DV:", dv, "\n")
  
  form <- as.formula(paste0(dv, " ~ HeadSet * Gender + Error(UniqueID/HeadSet)"))
  m <- art(form, data=data)
  a <- anova(m)
  print(a)
  
  # Effect sizes (partial eta²)
  es <- eta_squared(a, partial = TRUE)
  print(es)
  
  # Posthoc with Tukey ART-C
  cat("\n--- Tukey ART-C posthoc for HeadSet ---\n")
  post <- art.con(m, "HeadSet")
  print(summary(post, adjust="tukey"))
  
  # Simple effects if interaction is significant
  if (a["HeadSet:Gender","Pr(>F)"] < 0.05) {
    cat("\n--- Simple effects: HeadSet within Gender ---\n")
    post_int <- art.con(m, "HeadSet:Gender")
    print(summary(post_int, adjust="tukey"))
  }
  
  invisible(m)
}

# ===== 5) Run analyses =====
m_eff  <- run_art("Efficiency")
m_sat  <- run_art("Satisfaction")
m_effv <- run_art("Effectiveness")
m_ovr  <- run_art("Overall")
