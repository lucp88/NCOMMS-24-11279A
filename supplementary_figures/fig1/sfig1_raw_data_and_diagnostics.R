# =============================================================================
# Supplementary Figure 1 — Raw Data and Model Diagnostics for LMM 1
# =============================================================================
#
# Paper:  ReFOFA (NCOMMS-24-11279)
#
# This script generates all five panels of Supplementary Figure 1:
#   Panel A — Raw diffPSC data with group mean trajectories and SEM error bars
#   Panel B — Individual subject trajectories across runs
#   Panel C — Individual learning slopes (violin + boxplot + jitter)
#   Panel D — Residuals vs. fitted values (LMM 1 diagnostic)
#   Panel E — Observed vs. fitted values (LMM 1 diagnostic, r = 0.686)
#
# Data source: sourceData.xlsx, tab "Fig. 3" (same data as main Figure 3)
# Model:      LMM 1 (Run × Session × Group on diffPSC)
#
# Dependencies:
#   install.packages(c("readxl", "lme4", "lmerTest", "ggplot2", "dplyr",
#                       "cowplot"))
#
# Author: Lucas Peek
# Last updated: March 2026
# =============================================================================


# --- 1. Libraries ------------------------------------------------------------

library(readxl)
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(cowplot)


# --- 2. Load and prepare data -----------------------------------------------

# Set working directory to this script's location
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # Fallback: assumes the user is in the supplementary_figures/fig1/ directory
}
setwd("../..")  # Navigate up to 'code and source data' root

# Read source data (same as Fig. 3)
df <- read_excel("sourceData.xlsx", sheet = "Fig. 3", skip = 1)

# Re-scale DV to ensure scale attributes are available for lmer
df$diffROI_s <- scale(df$diffROI)

# Set factor levels and sum-to-zero contrasts
df$Regulation   <- factor(df$Regulation, levels = c("OFA", "FFA"))
df$groupID      <- factor(df$groupID, levels = c("Contr.", "Exp."))
df$Category     <- factor(df$Category, levels = c("Animals", "Faces"))
df$Training_Day <- factor(df$Training_Day, ordered = FALSE)

contrasts(df$Regulation)   <- c(-0.5, 0.5)
contrasts(df$groupID)      <- c(-0.5, 0.5)
contrasts(df$Category)     <- c(-0.5, 0.5)
contrasts(df$Training_Day) <- c(-0.5, 0.5)


# --- 3. Fit LMM 1 (same as Fig. 3) ------------------------------------------

md1 <- lmer(
  diffROI_s ~ runID + Regulation + Training_Day + groupID +
    Regulation * runID * groupID +
    (1 | subID) + (1 | stimID),
  data = df,
  REML = TRUE
)


# --- 4. Color scheme ---------------------------------------------------------

session_colors <- c("OFA" = "#0C4787", "FFA" = "#378AD3")


# =============================================================================
# PANEL A — Raw Data with Mean Trajectories
# =============================================================================

summary_data <- df %>%
  group_by(runID, Regulation, groupID) %>%
  summarise(
    mean_diff = mean(diffROI_s, na.rm = TRUE),
    se_diff   = sd(diffROI_s, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

pA <- ggplot() +
  # Individual trial data points
  geom_point(data = df,
             aes(x = runID, y = diffROI_s, color = Regulation),
             alpha = 0.1, size = 0.5) +
  # Group mean lines
  geom_line(data = summary_data,
            aes(x = runID, y = mean_diff, color = Regulation, group = Regulation),
            linewidth = 1.5) +
  # 95% CI error bars (±1.96 × SEM)
  geom_errorbar(data = summary_data,
                aes(x = runID,
                    ymin = mean_diff - 1.96 * se_diff,
                    ymax = mean_diff + 1.96 * se_diff,
                    color = Regulation),
                width = 0.2, alpha = 0.6) +
  facet_wrap(~ groupID, scales = "free_y") +
  scale_color_manual(values = session_colors,
                     labels = c("OFA Session", "FFA Session")) +
  scale_x_continuous(breaks = 1:7) +
  labs(x = "Run",
       y = "Differential PSC (FFA - OFA)",
       title = "Raw Data with Mean Trajectories",
       color = "Training Session") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 11, face = "bold"))


# =============================================================================
# PANEL B — Individual Subject Trajectories
# =============================================================================

subject_means <- df %>%
  group_by(subID, runID, Regulation, groupID) %>%
  summarise(mean_diff = mean(diffROI_s, na.rm = TRUE),
            .groups = "drop")

pB <- ggplot(subject_means,
             aes(x = runID, y = mean_diff,
                 color = Regulation,
                 group = interaction(subID, Regulation))) +
  # Thin individual lines
  geom_line(alpha = 0.3, linewidth = 0.5) +
  # Thick group mean
  stat_summary(aes(group = Regulation),
               fun = mean, geom = "line", linewidth = 2) +
  facet_wrap(~ groupID) +
  scale_color_manual(values = session_colors,
                     labels = c("OFA Session", "FFA Session")) +
  scale_x_continuous(breaks = 1:7) +
  labs(x = "Run",
       y = "Differential PSC (FFA - OFA)",
       title = "Individual Subject Trajectories",
       subtitle = "Thin lines = individuals, Thick lines = group means",
       color = "Training Session") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 11, face = "bold"))


# =============================================================================
# PANEL C — Individual Learning Slopes
# =============================================================================

# Compute per-subject learning slopes (beta of Run on diffPSC)
slopes <- df %>%
  group_by(subID, Regulation, groupID) %>%
  summarise(
    slope = coef(lm(diffROI_s ~ runID))[2],
    .groups = "drop"
  )

pC <- ggplot(slopes, aes(x = Regulation, y = slope, fill = Regulation)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
  facet_wrap(~ groupID) +
  scale_fill_manual(values = session_colors) +
  labs(x = "Training Session",
       y = "Learning Slope (beta per run)",
       title = "Individual Learning Slopes by Condition") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, face = "bold"))

# Print slope summary
slope_summary <- slopes %>%
  group_by(groupID, Regulation) %>%
  summarise(
    mean_slope = mean(slope),
    sd_slope   = sd(slope),
    n = n(),
    positive_slopes   = sum(slope > 0),
    percent_positive  = (sum(slope > 0) / n()) * 100,
    .groups = "drop"
  )
cat("\nLearning Slope Summary:\n")
print(slope_summary)


# =============================================================================
# PANELS D & E — Model Diagnostics (base R)
# =============================================================================

# Add residuals and fitted values to data
df$fitted    <- fitted(md1)
df$residuals <- residuals(md1)

# Correlation for panel E title
obs_fit_r <- round(cor(df$diffROI_s, df$fitted), 3)

# --- Save D & E as a single base R image ---
png("sfig1_panels_DE_diagnostics.png", width = 800, height = 1600, res = 150)
par(mfrow = c(2, 1), mar = c(5, 5, 4, 2))

# Panel D: Residuals vs Fitted
plot(df$fitted, df$residuals,
     main = "Residuals vs Fitted",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, cex = 0.5, col = rgb(0, 0, 0, 0.3))
abline(h = 0, col = "red", lty = 2)

# LOESS smooth with 95% CI band
lo <- loess(residuals ~ fitted, data = df, span = 0.75)
x_seq <- seq(min(df$fitted), max(df$fitted), length.out = 200)
pred <- predict(lo, newdata = data.frame(fitted = x_seq), se = TRUE)
polygon(c(x_seq, rev(x_seq)),
        c(pred$fit + 1.96 * pred$se.fit, rev(pred$fit - 1.96 * pred$se.fit)),
        col = adjustcolor("blue", alpha.f = 0.15), border = NA)
lines(x_seq, pred$fit, col = "blue", lwd = 2)
mtext("D", side = 3, line = 2, at = par("usr")[1], cex = 1.5, font = 2, adj = -0.2)

# Panel E: Observed vs Fitted
plot(df$diffROI_s, df$fitted,
     main = paste0("Observed vs Fitted (r = ", obs_fit_r, ")"),
     xlab = "Observed Values",
     ylab = "Fitted Values",
     pch = 19, cex = 0.5, col = rgb(0, 0, 1, 0.3))
abline(0, 1, col = "red", lty = 2, lwd = 2)
mtext("E", side = 3, line = 2, at = par("usr")[1], cex = 1.5, font = 2, adj = -0.2)

dev.off()
cat("\nPanels D & E saved: sfig1_panels_DE_diagnostics.png\n")


# =============================================================================
# SAVE PANELS A-C (ggplot-based)
# =============================================================================

# Save individual panels
ggsave("sfig1_panel_A_raw_data.pdf", pA, width = 10, height = 6)
ggsave("sfig1_panel_B_trajectories.pdf", pB, width = 10, height = 6)
ggsave("sfig1_panel_C_slopes.pdf", pC, width = 8, height = 6)
cat("Panels A, B, C saved as individual PDFs\n")

# Combined left column (A, B, C) using cowplot
left_column <- plot_grid(pA, pB, pC,
                         ncol = 1, nrow = 3,
                         rel_heights = c(1, 1, 1),
                         labels = c("A", "B", "C"),
                         label_size = 14)

ggsave("sfig1_panels_ABC_combined.png", left_column,
       width = 10, height = 18)
cat("Panels A-C combined saved: sfig1_panels_ABC_combined.png\n")

cat("\n=== Supplementary Figure 1 complete ===\n")
cat("Assemble final figure from:\n")
cat("  Left column:  sfig1_panels_ABC_combined.png\n")
cat("  Right column: sfig1_panels_DE_diagnostics.png\n")
