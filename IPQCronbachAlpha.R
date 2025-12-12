# ===== Packages =====

library(readr)
library(dplyr)
library(ARTool)
library(readr)
library(dplyr)
library(tidyr)
library(emmeans)
library(rstatix)
library(stringr)
# ===== Load & prepare =====
df <- read_csv("C:/Users/ITXPC/OneDrive - The Pennsylvania State University/Amir_abb6024/ICCMJournal/Codes/VRUsability/IgroupPresenceQuestionnaire.csv")

df$HeadSet <- trimws(as.character(df$HeadSet))
df$HeadSet <- factor(df$HeadSet, levels = c("Oculus","HTC","HP"))

aw_col <- names(df)[startsWith(names(df), "How aware")][1]
names(df)[match(c(
  'In the computer generated world I had a sense of "being there"',
  'Somehow I felt that the virtual world surrounded me.',
  'I felt like I was just perceiving pictures.',
  '(-) I did not feel present in the virtual space.',
  'I had a sense of acting in the virtual space, rather than operating something outside in real world.',
  'I felt present in the virtual space.',
  aw_col,
  'I was not aware of my real environment.',
  'I still paid attention to the real environment.',
  'How real did the virtual world seem to you?',
  'How much did your experience in the virtual environment seem consistent with your real world experience ?',
  'The virtual world seemed more realistic than the real world.'
), names(df))] <- c("GP","SP1","SP2","SP3","SP4","SP5","INV1","INV2","INV3","REAL1","REAL2","REAL3")

# Reverse-score (1–5 Likert → use 6 - x)
df <- df |>
  mutate(
    SP2_rev  = 6 - SP2,
    SP3_rev  = 6 - SP3,
    SP4_rev  = 6 - SP4,
    INV1_rev = 6 - INV1,
    INV3_rev = 6 - INV3
  )

sp_items  <- dplyr::select(df,GP, SP1, SP3_rev, SP4_rev, SP5)
inv_items <- dplyr::select(df, INV1_rev, INV2, INV3_rev)
rea_items <- dplyr::select(df, REAL1, REAL2, REAL3)

# ===== Overall Cronbach's alpha =====
a_sp  <- psych::alpha(na.omit(sp_items))
a_inv <- psych::alpha(na.omit(inv_items))
a_rea <- psych::alpha(na.omit(rea_items))

cat("\nOverall alpha:\n")
print(list(
  SpatialPresence = a_sp$total$raw_alpha,
  Involvement     = a_inv$total$raw_alpha,
  Realism         = a_rea$total$raw_alpha
))

# ===== Item diagnostics (corrected item–total r, alpha-if-deleted) =====
cat("\nSpatial Presence item diagnostics:\n")
print(a_sp$item.stats[, c("r.drop","alpha.if.deleted")])

# ===== By-headset alpha (sample sizes may be small) =====
by_headset_alpha <- df |>
  group_by(HeadSet) |>
  group_modify(~{
    tibble(
      SpatialPresence = psych::alpha(na.omit(dplyr::select(.x, SP1, SP2_rev, SP3_rev, SP4, SP5)))$total$raw_alpha,
      Involvement     = psych::alpha(na.omit(dplyr::select(.x, INV1_rev, INV2, INV3_rev)))$total$raw_alpha,
      Realism         = psych::alpha(na.omit(dplyr::select(.x, REAL1, REAL2, REAL3)))$total$raw_alpha
    )
  })
print(by_headset_alpha)

# ===== (Optional) Ordinal alpha via polychoric correlations =====
# This is often preferable for Likert data; requires psych::polychoric
sp_poly <- psych::polychoric(na.omit(sp_items))$rho
inv_poly <- psych::polychoric(na.omit(inv_items))$rho
rea_poly <- psych::polychoric(na.omit(rea_items))$rho

cat("\nOrdinal alpha (polychoric):\n")
print(list(
  SpatialPresence = psych::alpha(sp_poly, n.obs = nrow(na.omit(sp_items)))$total$raw_alpha,
  Involvement     = psych::alpha(inv_poly, n.obs = nrow(na.omit(inv_items)))$total$raw_alpha,
  Realism         = psych::alpha(rea_poly, n.obs = nrow(na.omit(rea_items)))$total$raw_alpha
))
