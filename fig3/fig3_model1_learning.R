# =============================================================================
# Figure 3 — Model 1: NFB Self-Regulation Learning
# =============================================================================
#
# Paper:  ReFOFA (NCOMMS-24-11279)
# DV:     diffPSC (FFA − OFA percent signal change during regulation), scaled
# Design: Linear mixed-effects model (LMM)
#         Fixed effects: Run (1–7) × Session (FFA/OFA) × Training Day × Group
#         Random effects: (1|Subject) + (1|Stimulus)
#
# This model tests whether participants in the experimental group learn to
# differentially regulate FFA vs OFA activity across successive NFB runs,
# compared to yoked controls receiving sham feedback.
#
# Output:
#   Left panel  — Dot-whisker plot of standardized fixed-effect coefficients
#   Right panel — Predicted marginal effects: Run × Session × Group
# =============================================================================

rm(list = ls())

# --- 1. Libraries ------------------------------------------------------------

library(readxl)      # read_excel() for source data
library(lme4)        # lmer() for mixed-effects models
library(lmerTest)    # p-values via Satterthwaite
library(sjPlot)      # plot_model() for predicted effects
library(ggplot2)     # plotting framework
library(dotwhisker)  # dwplot() for coefficient plots
library(dplyr)       # data wrangling
library(arm)         # standardize()
library(emmeans)     # emtrends() for simple slopes
library(egg)         # theme_article()
library(patchwork)   # combining plots


# --- 2. Load and prepare data -----------------------------------------------

# Set working directory to this script's location, then go up to 'code and source data' root
# This ensures the script works for anyone regardless of where the folder is on their machine
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # Fallback: assumes the user has opened R in the fig3/ directory
  # or sources the script from its own location
}
setwd("..")

# Read source data (already post-processed and labeled) from the Source Data file
# skip = 1 to skip the description row above the header
df <- read_excel("sourceData.xlsx", sheet = "Fig. 3", skip = 1)

# diffROI_s (scaled DV) is already in the source data, but we re-scale here
# to ensure the scale object attributes are available for lmer
df$diffROI_s <- scale(df$diffROI)

# Set factor levels and contrasts
# Reference levels: OFA session, Control group, Animals category
df$Regulation   <- factor(df$Regulation, levels = c("OFA", "FFA"))
df$groupID      <- factor(df$groupID, levels = c("Contr.", "Exp."))
df$Category     <- factor(df$Category, levels = c("Animals", "Faces"))
df$Training_Day <- factor(df$Training_Day, ordered = TRUE)

df$runID <- df$runID - mean(df$runID)

# Sum-to-zero contrasts (−0.5, +0.5) for interpretable main effects
# Training_Day uses polynomial contrasts (ordered factor) — consistent with fig5/fig6
contrasts(df$Regulation)   <- c(-0.5, 0.5)
contrasts(df$groupID)      <- c(-0.5, 0.5)
contrasts(df$Category)     <- c(-0.5, 0.5)


# --- 3. Fit the LMM ---------------------------------------------------------

# Model 1: NFB learning — does the Run × Session × Group interaction show
# differential learning between experimental and control participants?
md1 <- lmer(
  diffROI_s ~ runID + Regulation + Training_Day + groupID +
    Regulation * runID * groupID +
    (1 | subID) + (1 | stimID),
  data = df,
  REML = TRUE
)

# Print summary with Satterthwaite degrees of freedom
summary(md1)

# Construct standardized estimates and corresponding characteristics
md1_std <- standardize(md1)
summary(md1_std)

coefs <- as.data.frame(summary(md1_std)$coefficients)
ci <- confint(md1_std, method = "Wald")
ci <- ci[rownames(coefs), ]
cbind(coefs, ci)


# --- 4. Simple slopes analysis -----------------------------------------------
# Test the slope of Run (learning rate) for each Session × Group combination.
# This gives 4 slopes: Exp./FFA, Exp./OFA, Contr./FFA, Contr./OFA.
# Holm correction is applied across these 4 slopes.
# NOTE: we call emtrends() without adjust = "holm" so that emmeans uses
# Satterthwaite df (t.ratio), matching the LMM output. Holm correction
# is then applied manually via p.adjust().

slopes_by_group <- as.data.frame(summary(
  emtrends(md1, ~ Regulation * groupID, var = "runID", pbkrtest.limit = 4000),
  infer = c(TRUE, TRUE), level = 0.95
))
names(slopes_by_group)[grep("trend", names(slopes_by_group))] <- "slope"
slopes_by_group$p.holm <- p.adjust(slopes_by_group$p.value, method = "holm")
slopes_by_group$sig <- ifelse(slopes_by_group$p.holm < 0.001, "***",
                              ifelse(slopes_by_group$p.holm < 0.01, "**",
                                     ifelse(slopes_by_group$p.holm < 0.05, "*", "ns.")))
ratio_col <- ifelse("t.ratio" %in% names(slopes_by_group), "t.ratio", "z.ratio")
cat("\n=== Model 1: Learning — Simple slopes (Holm-corrected for 4 tests) ===\n")
print(slopes_by_group[, c("Regulation", "groupID", "slope", "SE", "df",
                          "lower.CL", "upper.CL", ratio_col,
                          "p.value", "p.holm", "sig")])


# Pairwise comparisons (for reference, not Holm-corrected)
emtrends(md1, pairwise ~ Regulation | groupID, var = "runID")
emtrends(md1, pairwise ~ groupID | Regulation, var = "runID")


# --- 5. Standardized coefficients -------------------------------------------
# md1_std is already created in Section 3 (line 87) via arm::standardize().
# Print summary here; the model object is reused for dwplot in Section 6.

print(summary(md1_std))


# --- 6. Plot — Left panel: Coefficient dot-whisker plot ----------------------
# Uses arm::standardize model object (md1_std) with by_2sd = FALSE so that
# plotted coefficients match the standardized values reported in the text.

color_coef <- c("#0F4889")

predictor_labels_std <- c(
  z.runID                          = "Run",
  c.Regulation                     = "Session",
  c.Training_Day                   = "Training Day",
  c.groupID                        = "Group",
  `z.runID:c.Regulation`           = "Run * Session",
  `c.Regulation:c.groupID`         = "Session * Group",
  `z.runID:c.groupID`              = "Run * Group",
  `z.runID:c.Regulation:c.groupID` = "Run * Session * Group"
)

plot_coef <- dwplot(
  list(md1_std),
  by_2sd       = FALSE,
  ci           = 0.95,
  dot_args     = list(size = 4, shape = 21, fill = "#fcd53a", color = "black", stroke = 1.2),
  whisker_args = list(size = 1, color = "black"),
  vline        = geom_vline(xintercept = 0, colour = "grey60", size = 1, linetype = 2),
  vars_order   = c(
    "z.runID", "c.Regulation", "c.Training_Day", "c.groupID",
    "z.runID:c.Regulation", "c.Regulation:c.groupID",
    "z.runID:c.groupID", "z.runID:c.Regulation:c.groupID"
  )
) %>%
  relabel_predictors(predictor_labels_std) +
  xlab("Standardized Coefficient (B)") +
  ylab("") +
  ggtitle("Model 1: Learning | DV: diffPSC (FFA\u2013OFA)") +
  theme_article() +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "none"
  )


# --- 7. Plot — Right panel: Predicted marginal effects -----------------------
# Uses raw (unstandardized) model for interpretable axis scales.

color_session <- c("#0C4787", "#378AD3")  # dark blue (OFA) | light blue (FFA)

plot_pred <- plot_model(
  md1, type = "pred", terms = c("runID", "Regulation", "groupID"),
  legend.title = "Session", show.legend = TRUE,
  colors = color_session
) +
  theme_sjplot2() +
  geom_line(size = 1) +
  xlab("Run") +
  ylab("diffPSC (FFA\u2013OFA, scaled)") +
  scale_x_continuous(breaks = seq(-3, 3, 1), labels = seq(1, 7, 1)) +
  ggtitle("Run \u00D7 Session \u00D7 Group") +
  theme(
    legend.justification = "top",
    plot.title   = element_text(size = 11),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 11)
  )


# --- 8. Combine panels -------------------------------------------------------

combined_plot <- wrap_elements(full = plot_coef) +
  wrap_elements(full = plot_pred) +
  plot_layout(ncol = 2, widths = c(2, 4))

print(combined_plot)
