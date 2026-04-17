# =============================================================================
# Supplementary Figure 5 — Behavioral Performance: Detection and Recognition
#                          Response Latencies
# =============================================================================
#
# Paper:  ReFOFA (NCOMMS-24-11279)
# Sample: N = 42 (21 EXP, 21 CONT) before outlier exclusion
#         N = 40 (20 per group) after excluding 2 EXP participants
#
# Analysis:
#   Two mixed-effects ANOVAs on response latencies (image frames):
#     ANOVA 1: DV = detection_frms
#     ANOVA 2: DV = recognition_frms
#   Within factors: Category (Faces, Animals), Regulation (FFA, OFA)
#   Between factor: groupID (exp, cont)
#
# Outlier exclusion: rstatix::identify_outliers (Q1 − 1.5×IQR / Q3 + 1.5×IQR)
#   on detection_frms, grouped by groupID × Category × Regulation.
#   Subjects flagged in any cell are excluded entirely.
#
# Output:
#   Left  — Bar plot: Detection latencies by Category and Training Session
#   Right — Bar plot: Recognition latencies by Category and Training Session
#   Both faceted by Group (EXP, CONT), with individual data points overlaid
#
# Data source: sourceData.xlsx, tab "SFig. 5"
#
# Dependencies:
#   install.packages(c("readxl", "rstatix", "dplyr", "ggplot2"))
#
# Author: Lucas Peek
# Last updated: March 2026
# =============================================================================

rm(list = ls())

# --- 1. Libraries ------------------------------------------------------------

library(readxl)      # read_excel() for source data
library(rstatix)     # identify_outliers(), anova_test()
library(dplyr)       # data wrangling
library(ggplot2)     # plotting framework
library(colorspace)
library(extrafont)   # Times New Roman font support


# --- 2. Working directory & data ---------------------------------------------
# Set working directory to the script's location (works in RStudio)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {}
setwd("../..")  # move up to "code and source data" root

# Read from Source Data Excel (skip row 1 = description)
means <- read_excel("sourceData.xlsx", sheet = "SFig. 5", skip = 1)

cat("Data loaded:", nrow(means), "rows,", length(unique(means$subID)), "subjects\n")


# --- 3. Outlier exclusion ----------------------------------------------------
# Identify outliers on detection latencies (Q1 − 1.5×IQR / Q3 + 1.5×IQR)
# grouped by condition. Excludes entire subject if flagged in any cell.

outliers <- means %>%
  group_by(groupID, Category, Regulation) %>%
  identify_outliers(detection_frms)

subs2excl <- unique(outliers$subID)
means <- means[!(means$subID %in% subs2excl), ]

cat("Subjects excluded:", length(subs2excl), "\n")
cat("Remaining subjects:", length(unique(means$subID)), "\n")


# --- 4. ANOVA 1: Detection latencies ----------------------------------------

anova1 <- anova_test(
  data = means,
  dv = detection_frms,
  wid = subID,
  within = c(Category, Regulation),
  between = groupID,
  detailed = TRUE
)
cat("\n=== ANOVA 1: Detection (frms) ===\n")
print(get_anova_table(anova1))


# --- 5. ANOVA 2: Recognition latencies --------------------------------------

anova2 <- anova_test(
  data = means,
  dv = recognition_frms,
  wid = subID,
  within = c(Category, Regulation),
  between = groupID,
  detailed = TRUE
)
cat("\n=== ANOVA 2: Recognition (frms) ===\n")
print(get_anova_table(anova2))


# --- 6. Bar plots with individual data points --------------------------------

# Recode group labels for display: exp → EXP, sham → CONT
# (done after ANOVAs to avoid factor-level issues with anova_test)
means$groupID <- factor(means$groupID,
                        levels = c("exp", "sham"),
                        labels = c("EXP", "CONT"))

my_colors <- c("#2E994E", "#99201F")

# -- 6a. Detection bars --

se_det <- means %>%
  group_by(groupID, Category, Regulation) %>%
  summarize(se = sd(detection_frms) / sqrt(n()), .groups = "drop")

me_det <- means %>%
  group_by(groupID, Category, Regulation) %>%
  summarise(detection_frms = mean(detection_frms), .groups = "drop")

me_det <- me_det %>% left_join(se_det, by = c("groupID", "Category", "Regulation"))

p_det <- ggplot(me_det, aes(x = Category, y = detection_frms, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5, alpha = 0.9) +
  geom_jitter(data = means, aes(x = Category, y = detection_frms, color = Regulation),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5),
              size = 1.5, alpha = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = detection_frms - se, ymax = detection_frms + se),
                width = 0.2, position = position_dodge(0.5)) +
  facet_wrap(~ groupID, ncol = 2) +
  labs(x = "Category", y = "Detection (frms)", title = "Detection") +
  scale_fill_manual(values = my_colors, name = "Training Session") +
  scale_color_manual(values = darken(my_colors, 0.3)) +
  theme_classic(base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        legend.key.size = unit(1.5, "lines"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  coord_cartesian(ylim = c(60, NA))


# -- 6b. Recognition bars --

se_rec <- means %>%
  group_by(groupID, Category, Regulation) %>%
  summarize(se = sd(recognition_frms) / sqrt(n()), .groups = "drop")

me_rec <- means %>%
  group_by(groupID, Category, Regulation) %>%
  summarise(recognition_frms = mean(recognition_frms), .groups = "drop")

me_rec <- me_rec %>% left_join(se_rec, by = c("groupID", "Category", "Regulation"))

p_rec <- ggplot(me_rec, aes(x = Category, y = recognition_frms, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5, alpha = 0.9) +
  geom_jitter(data = means, aes(x = Category, y = recognition_frms, color = Regulation),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5),
              size = 1.5, alpha = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = recognition_frms - se, ymax = recognition_frms + se),
                width = 0.2, position = position_dodge(0.5)) +
  facet_wrap(~ groupID, ncol = 2) +
  labs(x = "Category", y = "Recognition (frms)", title = "Recognition") +
  scale_fill_manual(values = my_colors, name = "Training Session") +
  scale_color_manual(values = darken(my_colors, 0.3)) +
  theme_classic(base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        legend.key.size = unit(1.5, "lines"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  coord_cartesian(ylim = c(100, NA))

p_det
p_rec

# ggsave("detection_beh.png", plot = p_det, width = 871, height = 719, units = "px", dpi = 300)
# ggsave("recognition_beh.png", plot = p_rec, width = 871, height = 719, units = "px", dpi = 300)

# Descriptive stats for manuscript text
means %>%
  group_by(groupID) %>%
  summarise(M_det = mean(detection_frms), SD_det = sd(detection_frms),
            M_rec = mean(recognition_frms), SD_rec = sd(recognition_frms))

means %>%
  group_by(Category) %>%
  summarise(M_det = mean(detection_frms), SD_det = sd(detection_frms),
            M_rec = mean(recognition_frms), SD_rec = sd(recognition_frms))

# --- Post-hoc t-tests (Holm-corrected, family of 4) -------------------------
# Follow up significant main effects from both ANOVAs

# Detection: Group (collapsing across Category and Regulation)
t_det_group <- t.test(detection_frms ~ groupID, data = means, var.equal = FALSE)

# Detection: Category (collapsing across Group and Regulation)
t_det_cat <- t.test(detection_frms ~ Category, data = means, paired = TRUE)

# Recognition: Group
t_rec_group <- t.test(recognition_frms ~ groupID, data = means, var.equal = FALSE)

# Recognition: Category
t_rec_cat <- t.test(recognition_frms ~ Category, data = means, paired = TRUE)

# Collect p-values and apply Holm
p_raw <- c(det_group = t_det_group$p.value,
           det_cat   = t_det_cat$p.value,
           rec_group = t_rec_group$p.value,
           rec_cat   = t_rec_cat$p.value)
p_holm <- p.adjust(p_raw, method = "holm")

cat("\n=== Post-hoc t-tests (Holm-corrected, family of 4) ===\n")
for (i in seq_along(p_raw)) {
  cat(sprintf("%-12s  p.raw = %.4f  p.holm = %.4f\n", names(p_raw)[i], p_raw[i], p_holm[i]))
}

