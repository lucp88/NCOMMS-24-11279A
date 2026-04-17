# =============================================================================
# Supplementary Figure 8 — STS Control Analysis: Learning
# =============================================================================
#
# Paper:  ReFOFA (NCOMMS-24-11279)
# Sample: Experimental group only, all stimuli (faces + animals), N = 22,
#         2077 trial-level observations
#
# Model:
#   LMM 1 (STS): DV = regSTS_s (z-scored STS regulation activity)
#                 Predictors: Run × Session, Training Day
#
# This analysis mirrors Model 1 in the main text (Fig. 3) but replaces the
# differential FFA−OFA signal (diffROI_s) with regulation activity in the
# right STS — a co-activated but non-targeted face-processing region.
# The STS serves as a control to demonstrate that NFB effects are specific
# to the targeted ROIs and do not extend to STS.
#
# Note: Unlike the main Model 1, which includes both EXP and CONT groups
# (with a Group factor), STS data is available for the EXP group only.
# Therefore, the Group factor is not included in this model.
#
# Random effects: (1|subID) + (1|stimID)
#
# Output:
#   Left  — Dot-whisker coefficient plot
#   Right — Predicted Run × Session interaction
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
library(reghelper)   # simple_slopes()
library(egg)         # theme_article()
library(patchwork)   # combining plots


# --- 2. Working directory & data ---------------------------------------------
# Set working directory to this script's location, then go up to root
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # Fallback: assumes the user has opened R in the fig8/ directory
}
setwd("../..")  # move up to "code and source data" root

# Read source data (already post-processed and labeled) from the Source Data file
# skip = 1 to skip the description row above the header
df <- read_excel("sourceData.xlsx", sheet = "SFig. 8", skip = 1)

cat("Data loaded:", nrow(df), "trials,", length(unique(df$subID)), "subjects\n")

# Re-scale to ensure scale object attributes are available for lmer
df$regSTS_s <- scale(df$regSTS)


# --- 3. Factor coding --------------------------------------------------------
# Sum-to-zero (deviation) contrasts: −0.5 / +0.5
# Levels ordered so the contrast direction matches the manuscript

df$Regulation <- factor(df$Regulation, levels = c("OFA", "FFA"))
contrasts(df$Regulation) <- c(0.5, -0.5)

df$Training_Day <- factor(df$Training_Day, ordered = TRUE)


# --- 4. Model fitting --------------------------------------------------------
# STS learning model: does STS regulation activity change across runs
# and differ between training sessions?
# This mirrors Model 1 from the main text but with regSTS_s as DV
# and without the Group factor (EXP only).

cat("\n--- Fitting LMM 1 (STS): regSTS_s ~ Run × Session ---\n")
md1 <- lmer(regSTS_s ~ runID + Regulation + Training_Day +
              Regulation * runID +
              (1 | subID) + (1 | stimID),
            data = df, REML = TRUE)
print(summary(md1))


# --- 4a. Standardized coefficients -------------------------------------------
cat("\n--- Standardized coefficients (arm::standardize) ---\n")
md1_std <- standardize(md1)
print(summary(md1_std))


# --- 4b. Simple slopes via emtrends ------------------------------------------
cat("\n--- Simple slopes: Run slope by Session ---\n")
t_run <- as.data.frame(test(emtrends(md1, ~ Regulation, var = "runID")))
names(t_run)[grep("trend", names(t_run))] <- "slope"
t_run$predictor <- "runID"

# No Holm correction needed here: only 2 slopes (same as main Fig. 3
# where the learning slopes are the primary test)
ratio_col <- ifelse("t.ratio" %in% names(t_run), "t.ratio", "z.ratio")
cat("\nRun slope by Session (uncorrected):\n")
print(t_run[, c("predictor", "Regulation", "slope", "SE", ratio_col, "p.value")])


# --- 5. Coefficient plot (dot-whisker) ----------------------------------------
colorSpecs <- c('#0F4889')

mylabels <- c(
  runID              = "Run",
  Regulation1        = "Session",
  Training_Day.L     = "Training Day",
  `runID:Regulation1` = "Run * Session"
)

plot1 <- dwplot(
  list(md1),
  by_2sd       = TRUE,
  ci           = 0.95,
  dot_args     = list(size = 4, color = colorSpecs),
  whisker_args = list(size = 1, color = colorSpecs),
  vline        = geom_vline(xintercept = 0, colour = "grey60", size = 1, linetype = 2),
  vars_order   = c("runID", "Regulation1", "Training_Day.L", "runID:Regulation1")
) %>%
  relabel_predictors(mylabels) +
  xlab("Standardized Coefficient (B)") +
  ylab("") +
  ggtitle("LMM 1 (STS): Learning | DV: regSTS") +
  theme_blank() +
  theme_article()


# --- 6. Marginal effect plot: Run × Session -----------------------------------
colorSpecs_1 <- c('#0C4787', '#378AD3')  # dark blue | light blue

plot2 <- plot_model(
  md1, type = "pred", terms = c("runID", "Regulation"),
  legend.title = "Session", show.legend = TRUE,
  colors = colorSpecs_1
) +
  theme_sjplot2() +
  geom_line(size = 1) +
  xlab("Run") +
  ylab("STS Regulation Activity (z)") +
  scale_x_continuous(limits = c(1, 7), breaks = seq(1, 7, 1)) +
  ggtitle("Run by Session") +
  theme(
    legend.justification = "top",
    plot.title   = element_text(size = 11),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 11)
  )


# --- 7. Combine panels -------------------------------------------------------
plot1_wrapped <- wrap_elements(full = plot1)
plot2_wrapped <- wrap_elements(full = plot2)

final_figure <- plot1_wrapped + plot2_wrapped +
  plot_layout(ncol = 2, widths = c(2, 4))

print(final_figure)


# --- 8. Summary table ---------------------------------------------------------
tab_model(md1,
          dv.labels = "LMM 1 (STS): Regulation Activity",
          show.se = TRUE, show.std = "std2",
          string.est = "Estimate", string.se = "SE",
          title = "Linear Mixed Model: STS Regulation Activity (Learning)")
