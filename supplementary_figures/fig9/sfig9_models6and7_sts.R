# =============================================================================
# Supplementary Figure 9 — STS Control: Models 6 & 7 with STS
#                           Detection and Recognition Performance
# =============================================================================
#
# Paper:  ReFOFA (NCOMMS-24-11279)
# Sample: Experimental group only, face stimuli only (N = 22, 1494 trials)
#
# Models:
#   Model 6 (md6): DV = detection_frms (detection performance in frames)
#                  Predictors: FFA × OFA × Session + STS + STS × Session
#                              + Training Day
#
#   Model 7 (md7): DV = recognition_frms (recognition performance in frames)
#                  Predictors: FFA × OFA × Session + STS + STS × Session
#                              + Training Day + Detection_FRMS (covariate)
#
# This analysis mirrors Models 6 & 7 from the main text (Fig. 6) but adds
# STS task-evoked activity as an additional predictor alongside FFA and OFA,
# plus an STS × Session interaction. This tests whether the non-targeted STS
# region contributes to behavioral performance beyond the targeted FFA/OFA.
# STS is NOT included in the three-way FFA × OFA × Session interaction
# (as reviewed and accepted by reviewers).
#
# Key predictors: task_detect_FFA/OFA/STS_beta_s (Model 6),
#                 task_recogn_FFA/OFA/STS_beta_s (Model 7), all z-scored
# Random effects: (1|subID) + (1|stimID)
#
# Output:
#   Left panel  — Dot-whisker coefficients (Models 6 & 7 overlaid)
#   Right top   — Detection: predicted effects of FFA, OFA, STS, FFA×OFA
#   Right bottom— Recognition: predicted effects of FFA×Session, OFA×Session,
#                  FFA×OFA
# =============================================================================

rm(list = ls())

# --- 0. Settings --------------------------------------------------------------
# Set to TRUE to plot from standardized models; FALSE uses raw-scale predictors.
USE_STD <- FALSE

# --- 1. Libraries ------------------------------------------------------------

library(readxl)      # read_excel() for source data
library(lme4)        # lmer() for mixed-effects models
library(lmerTest)    # p-values via Satterthwaite
library(sjPlot)      # plot_model() for predicted effects
library(ggplot2)     # plotting framework
library(dotwhisker)  # dwplot() for coefficient plots
library(dplyr)       # data wrangling
library(tidyr)       # rename()
library(arm)         # standardize()
library(emmeans)     # emtrends() for simple slopes
library(reghelper)   # simple_slopes()
library(egg)         # theme_article()
library(patchwork)   # combining plots
library(cowplot)     # plot_grid()


# --- 2. Load and prepare data -----------------------------------------------
# Set working directory to this script's location, then go up to root
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # Fallback: assumes the user has opened R in the fig9/ directory
}
setwd("../..")  # move up to "code and source data" root

# Read source data from the Source Data file
# skip = 1 to skip the description row above the header
df <- read_excel("sourceData.xlsx", sheet = "SFig. 9", skip = 1)

cat("Data loaded:", nrow(df), "trials,", length(unique(df$subID)), "subjects\n")

# Re-scale to ensure scale object attributes are available for lmer
df$task_detect_FFA_beta_s <- scale(df$task_detect_FFA_beta)
df$task_detect_OFA_beta_s <- scale(df$task_detect_OFA_beta)
df$task_detect_STS_beta_s <- scale(df$task_detect_STS_beta)
df$task_recogn_FFA_beta_s <- scale(df$task_recogn_FFA_beta)
df$task_recogn_OFA_beta_s <- scale(df$task_recogn_OFA_beta)
df$task_recogn_STS_beta_s <- scale(df$task_recogn_STS_beta)


# --- 3. Prepare model-specific data frames -----------------------------------

# Model 6 (detection): rename task detection betas to generic FFA/OFA/STS
df_det <- df %>%
  dplyr::rename(
    FFA = task_detect_FFA_beta_s,
    OFA = task_detect_OFA_beta_s,
    STS = task_detect_STS_beta_s,
    Detection_FRMS = detection_frms
  )

# Model 7 (recognition): rename task recognition betas to generic FFA/OFA/STS
df_rec <- df %>%
  dplyr::rename(
    FFA = task_recogn_FFA_beta_s,
    OFA = task_recogn_OFA_beta_s,
    STS = task_recogn_STS_beta_s,
    Detection_FRMS = detection_frms,
    Recognition_FRMS = recognition_frms
  )

# Set factor levels and contrasts
# Regulation: OFA as first level, FFA as second
# Sum-to-zero contrast: OFA = -0.5, FFA = +0.5 (consistent with fig3/fig5/fig6)
df_det$Regulation   <- factor(df_det$Regulation, ordered = FALSE, levels = c("OFA", "FFA"))
df_det$Training_Day <- factor(df_det$Training_Day, ordered = TRUE)
contrasts(df_det$Regulation) <- c(-.5, .5)

df_rec$Regulation   <- factor(df_rec$Regulation, ordered = FALSE, levels = c("OFA", "FFA"))
df_rec$Training_Day <- factor(df_rec$Training_Day, ordered = TRUE)
contrasts(df_rec$Regulation) <- c(-.5, .5)


# --- 4. Fit the LMMs --------------------------------------------------------

# --- 4.0 STS-only models (reported in text, "when considered alone") ---------
# These test whether STS activity alone predicts behavioral performance,
# before adding FFA and OFA to the model.

cat("\n--- Fitting md0_det: STS-only → Detection ---\n")
md0_det <- lmer(
  Detection_FRMS ~ STS +
    Regulation + Training_Day + STS * Regulation +
    (1 | subID) + (1 | stimID),
  data = df_det,
  REML = TRUE
)
print(summary(md0_det))

cat("\n--- Fitting md0_rec: STS-only → Recognition ---\n")
md0_rec <- lmer(
  Recognition_FRMS ~ STS +
    Regulation + Training_Day + STS * Regulation +
    Detection_FRMS +
    (1 | subID) + (1 | stimID),
  data = df_rec,
  REML = TRUE
)
print(summary(md0_rec))

cat("\n--- STS-only standardized coefficients ---\n")
md0_det_std <- standardize(md0_det)
md0_rec_std <- standardize(md0_rec)
print(summary(md0_det_std))
print(summary(md0_rec_std))

# Full standardized stats with Wald 95% CIs for the supplementary text.
# std_with_ci() is defined in section 5; here we use an inline copy so
# this block is self-contained (section 4 runs before section 5).
.std_with_ci <- function(model, label) {
  std_obj <- standardize(model)
  coefs <- as.data.frame(summary(std_obj)$coefficients)
  ci <- confint(std_obj, method = "Wald")
  ci <- ci[rownames(coefs), ]
  result <- cbind(coefs, ci)
  cat(paste0("\n", label, ":\n"))
  print(result)
  return(result)
}

md0_det_std_ci <- .std_with_ci(md0_det, "STS-only Model: Detection performance")
md0_rec_std_ci <- .std_with_ci(md0_rec, "STS-only Model: Recognition performance")

tab_model(md0_det, md0_rec,
          dv.labels = c("Detection (STS only)", "Recognition (STS only)"),
          string.est = "Estimate", string.se = "SE",
          title = "STS-Only Models: Behavioral Performance")


# --- 4.1 Full models with FFA, OFA, and STS ---------------------------------

# Model 6: Detection performance (with STS)
# Same as main text Model 6, plus STS main effect and STS × Session
cat("\n--- Fitting Model 6 (with STS): Detection performance ---\n")
md6 <- lmer(
  Detection_FRMS ~ FFA + OFA + STS +
    Regulation + Training_Day +
    STS * Regulation +
    FFA * OFA * Regulation +
    (1 | subID) + (1 | stimID),
  data = df_det,
  REML = TRUE
)
summary(md6)

# Model 7: Recognition performance (with STS)
# Same as main text Model 7, plus STS main effect and STS × Session
cat("\n--- Fitting Model 7 (with STS): Recognition performance ---\n")
md7 <- lmer(
  Recognition_FRMS ~ FFA + OFA + STS +
    Regulation + Training_Day + Detection_FRMS +
    STS * Regulation +
    FFA * OFA * Regulation +
    (1 | subID) + (1 | stimID),
  data = df_rec,
  REML = TRUE
)
summary(md7)


# --- 5. Standardized coefficients & plot-ready models -------------------------

md6_std_obj <- standardize(md6)
md7_std_obj <- standardize(md7)
print(summary(md6_std_obj))
print(summary(md7_std_obj))

std_with_ci <- function(model, label) {
  std_obj <- standardize(model)
  coefs <- as.data.frame(summary(std_obj)$coefficients)
  ci <- confint(std_obj, method = "Wald")
  ci <- ci[rownames(coefs), ]
  result <- cbind(coefs, ci)
  cat(paste0("\n", label, ":\n"))
  print(result)
  return(result)
}

md6_std_ci <- std_with_ci(md6, "Model 6 (STS): Detection performance")
md7_std_ci <- std_with_ci(md7, "Model 7 (STS): Recognition performance")

# Select models for emtrends and plotting based on USE_STD flag
if (USE_STD) {
  md6_plot <- md6_std_obj
  md7_plot <- md7_std_obj
  var_FFA  <- "z.FFA"
  var_OFA  <- "z.OFA"
  var_STS  <- "z.STS"
  var_Reg  <- "c.Regulation"
  cat("\n>>> Using STANDARDIZED models for slopes and plots <<<\n")
} else {
  md6_plot <- md6
  md7_plot <- md7
  var_FFA  <- "FFA"
  var_OFA  <- "OFA"
  var_STS  <- "STS"
  var_Reg  <- "Regulation"
  cat("\n>>> Using RAW models for slopes and plots <<<\n")
}


# --- 6. Simple slopes analysis -----------------------------------------------
# 6a. STS main effect slopes (panels with STS)
cat("\n--- STS → Detection (main effect) ---\n")
print(test(emtrends(md6_plot, ~ 1, var = var_STS)))

cat("\n--- STS → Recognition (main effect) ---\n")
print(test(emtrends(md7_plot, ~ 1, var = var_STS)))

# 6b. STS × Session slopes
cat("\n--- STS × Session → Detection ---\n")
t_sts_det <- as.data.frame(test(emtrends(md6_plot, as.formula(paste("~", var_Reg)),
                                          var = var_STS)))
names(t_sts_det)[grep("trend", names(t_sts_det))] <- "slope"
print(t_sts_det)

cat("\n--- STS × Session → Recognition ---\n")
t_sts_rec <- as.data.frame(test(emtrends(md7_plot, as.formula(paste("~", var_Reg)),
                                          var = var_STS)))
names(t_sts_rec)[grep("trend", names(t_sts_rec))] <- "slope"
print(t_sts_rec)

# 6c. Model 6 — Detection: FFA × OFA slopes (panel c), same as main text
cat("\n\n=== Model 6 (Detection): FFA × OFA slopes ===\n")
ofa_at <- setNames(list(c(-3, 0, 3)), var_OFA)
t_c <- as.data.frame(test(emtrends(md6_plot, as.formula(paste("~", var_OFA)),
                                    var = var_FFA, at = ofa_at)))
names(t_c)[grep("trend", names(t_c))] <- "slope"
t_c$panel <- "c"
ratio_col <- ifelse("t.ratio" %in% names(t_c), "t.ratio", "z.ratio")
print(t_c[, c("panel", var_OFA, "slope", "SE", ratio_col, "p.value")])

# 6d. Model 7 — Recognition: FFA/OFA × Session slopes
cat("\n\n=== Model 7 (Recognition): FFA × Session ===\n")
t_d <- as.data.frame(test(emtrends(md7_plot, as.formula(paste("~", var_Reg)), var = var_FFA)))
names(t_d)[grep("trend", names(t_d))] <- "slope"
print(t_d)

cat("\n=== Model 7 (Recognition): OFA × Session ===\n")
t_e <- as.data.frame(test(emtrends(md7_plot, as.formula(paste("~", var_Reg)), var = var_OFA)))
names(t_e)[grep("trend", names(t_e))] <- "slope"
print(t_e)

# 6e. Model 7 — FFA × OFA interaction
cat("\n=== Model 7 (Recognition): FFA × OFA ===\n")
reg_formula <- as.formula(paste("~", var_OFA, "+", var_Reg))
t_f <- as.data.frame(test(emtrends(md7_plot, reg_formula,
                                    var = var_FFA, at = ofa_at)))
names(t_f)[grep("trend", names(t_f))] <- "slope"
print(t_f)


# --- 7. Plot — Left panel: Coefficient dot-whisker plot ----------------------

# Color scheme: green = detection, gold = recognition
color_det  <- "#597e59"
color_rec  <- "#dbbb5b"

predictor_labels_std <- c(
  z.FFA                        = "FFA",
  z.OFA                        = "OFA",
  z.STS                        = "STS",
  c.Regulation                 = "Session",
  c.Training_Day               = "Training Day",
  z.Detection_FRMS             = "Detection FRMS",
  `z.FFA:z.OFA`                = "FFA * OFA",
  `z.FFA:c.Regulation`         = "FFA * Session",
  `z.OFA:c.Regulation`         = "OFA * Session",
  `z.STS:c.Regulation`         = "STS * Session",
  `z.FFA:z.OFA:c.Regulation`   = "FFA * OFA * Session"
)


colorSpecs <- c(color_det, color_rec,
                color_det, color_rec,
                color_det, color_rec,
                color_det, color_rec,
                color_det, color_rec,
                color_det, color_rec,
                color_det, color_rec,
                color_det, color_rec,
                color_det, color_rec,
                color_det, color_rec)

plot_coef <- dwplot(
  list(md7_std_obj, md6_std_obj),
  by_2sd       = FALSE,
  ci           = 0.95,
  # dot_args     = list(size = 4, color = colorSpecs),
  dot_args     = list(size = 4, shape = 21, fill = colorSpecs, color = "black", stroke = 1.2),
  whisker_args = list(size = 1, color = colorSpecs),
  style        = "dotwhisker",
  vline        = geom_vline(xintercept = 0, colour = "grey60", size = 1, linetype = 2),
  vars_order   = c("z.FFA", "z.OFA", "z.STS", "c.Regulation", "c.Training_Day",
                    "z.FFA:z.OFA", "z.FFA:c.Regulation", "z.OFA:c.Regulation",
                    "z.STS:c.Regulation", "z.FFA:z.OFA:c.Regulation")
) %>%
  relabel_predictors(predictor_labels_std) +
  scale_color_manual(name = "Behavior",
                     labels = c("Detection", "Recognition"),
                     values = c(color_det, color_rec)) +
  xlab("Standardized Coefficient (B)") +
  ylab("") +
  ggtitle("Model 6 and 7 (with STS): Detection and Recognition") +
  theme_article()


# --- 8. Plot — Right panels: Predicted marginal effects ----------------------
# Uses md6_plot/md7_plot and var_FFA/var_OFA/var_STS set by USE_STD flag

# Color palettes
color_det_line <- c("#597e59", "#a2bca2")  # dark/light green
color_det_3    <- c("#a2bca2", "#7ba27b", "#597e59")
color_rec_line <- c("#dbbb5b", "#e4cd82")  # dark/light gold
color_rec_3    <- c("#e4cd82", "#d7b64a", "#dbbb5b")

ofa_3lvl <- paste0(var_OFA, " [-3, 0, 3]")

# --- Detection performance (top right) ---
# p1: STS → Detection FRMS
p1 <- plot_model(
  md6_plot, type = "pred", terms = c(var_STS),
  colors = color_det[1], show.legend = FALSE, grid = FALSE
) +
  theme_sjplot2() + geom_line(size = 1) +
  xlab(expression(beta)) + ylab("frms") +
  ggtitle("STS Detection Response") +
  theme(plot.title = element_text(size = 11),
        axis.title.y = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 11, face = "bold"))

# p2: OFA → Detection FRMS
p2 <- plot_model(
  md6_plot, type = "pred", terms = c(var_OFA),
  colors = color_det[1], show.legend = FALSE, grid = FALSE
) +
  theme_sjplot2() + geom_line(size = 1) +
  xlab(expression(beta)) + ylab("frms") +
  ggtitle("OFA Detection Response") +
  theme(plot.title = element_text(size = 11),
        axis.title.y = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 11, face = "bold"))

# p3: FFA × OFA → Detection FRMS
p3 <- plot_model(
  md6_plot, type = "pred", terms = c(var_FFA, ofa_3lvl),
  legend.title = "OFA Detection Resp.", colors = color_det_3,
  show.legend = FALSE, show.data = FALSE, grid = FALSE
) +
  theme_sjplot2() + geom_line(size = 1) +
  xlab(expression(beta)) + ylab("frms") +
  ggtitle("FFA Detection x OFA Detection") +
  theme(plot.title = element_text(size = 11),
        axis.title.y = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 11, face = "bold"))

# --- Recognition performance (bottom right) ---
# p4: FFA × Session → Recognition FRMS
p4 <- plot_model(
  md7_plot, type = "pred", terms = c(var_FFA, var_Reg),
  legend.title = "Session", colors = color_rec_line,
  show.legend = FALSE, grid = FALSE
) +
  theme_sjplot2() + geom_line(size = 1) +
  xlab(expression(beta)) + ylab("frms") +
  ggtitle("FFA Recognition x Session") +
  theme(plot.title = element_text(size = 11),
        axis.title.y = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 11, face = "bold"))

# p5: OFA × Session → Recognition FRMS
p5 <- plot_model(
  md7_plot, type = "pred", terms = c(var_OFA, var_Reg),
  legend.title = "OFA", colors = color_rec_line,
  show.legend = FALSE, show.data = FALSE, grid = FALSE
) +
  theme_sjplot2() + geom_line(size = 1) +
  xlab(expression(beta)) + ylab("frms") +
  ggtitle("OFA Recognition x Session") +
  theme(plot.title = element_text(size = 11),
        axis.title.y = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 11, face = "bold"))

# p6: FFA × OFA → Recognition FRMS
p6 <- plot_model(
  md7_plot, type = "pred", terms = c(var_FFA, ofa_3lvl),
  legend.title = "OFA Recognition", colors = color_rec_3,
  show.legend = FALSE, show.data = FALSE, grid = FALSE
) +
  theme_sjplot2() + geom_line(size = 1) +
  xlab(expression(beta)) + ylab("frms") +
  ggtitle("FFA Recognition by OFA Recognition") +
  theme(plot.title = element_text(size = 11),
        axis.title.y = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 11, face = "bold"))


# --- 9. Combine panels -------------------------------------------------------
left_col          <- plot_coef
top_right_left    <- cowplot::plot_grid(p1, p2, ncol = 1)
top_right         <- cowplot::plot_grid(top_right_left, p3, ncol = 2)
bottom_right_left <- cowplot::plot_grid(p4, p5, ncol = 1)
bottom_right      <- cowplot::plot_grid(bottom_right_left, p6, ncol = 2)

right_side <- cowplot::plot_grid(top_right, bottom_right, ncol = 1)
final_figure <- cowplot::plot_grid(left_col, right_side, ncol = 2, rel_widths = c(1, 2))
print(final_figure)


# --- 10. Summary tables (for supplementary materials) -----------------------
# Mirrors fig6 Section 10: two tables, one per model, each showing
# unstandardized coefficients alongside standardized coefficients (2-SD units).
#
# Standardized column is produced via an arm::rescale() workaround rather than
# passing arm::standardize() objects directly to tab_model() (sjPlot trips on
# the z./c. prefixed variable names). The workaround is mathematically
# equivalent to arm::standardize(standardize.y = FALSE).
#
# NOTE: only PREDICTORS are rescaled — never the DV. arm::standardize() leaves
# the DV on its original scale by default, and we match that here.
# Training_Day is an ordered factor and is also left alone (tab_model shows
# the linear contrast Training_Day.L, same as the unstandardized column).

df_det_std <- df_det
df_det_std$FFA            <- arm::rescale(df_det$FFA)
df_det_std$OFA            <- arm::rescale(df_det$OFA)
df_det_std$STS            <- arm::rescale(df_det$STS)
# Detection_FRMS is the DV in Model 6 — leave untouched
# Training_Day is an ordered factor — leave untouched
contrasts(df_det_std$Regulation) <- c(-.5, .5)

df_rec_std <- df_rec
df_rec_std$FFA              <- arm::rescale(df_rec$FFA)
df_rec_std$OFA              <- arm::rescale(df_rec$OFA)
df_rec_std$STS              <- arm::rescale(df_rec$STS)
df_rec_std$Detection_FRMS   <- arm::rescale(df_rec$Detection_FRMS)  # predictor (covariate) in Model 7
# Recognition_FRMS is the DV in Model 7 — leave untouched
# Training_Day is an ordered factor — leave untouched
contrasts(df_rec_std$Regulation) <- c(-.5, .5)

# Refit on rescaled data (same formula structure as section 4)
md6_std_tab <- lmer(
  Detection_FRMS ~ FFA + OFA + STS +
    Regulation + Training_Day +
    STS * Regulation +
    FFA * OFA * Regulation +
    (1 | subID) + (1 | stimID),
  data = df_det_std,
  REML = TRUE
)

md7_std_tab <- lmer(
  Recognition_FRMS ~ FFA + OFA + STS +
    Regulation + Training_Day + Detection_FRMS +
    STS * Regulation +
    FFA * OFA * Regulation +
    (1 | subID) + (1 | stimID),
  data = df_rec_std,
  REML = TRUE
)

# Pretty predictor labels
pred_labels_md6 <- c(
  `(Intercept)`         = "(Intercept)",
  FFA                   = "FFA Detection",
  OFA                   = "OFA Detection",
  STS                   = "STS Detection",
  Regulation1           = "Session",
  Training_Day.L        = "Training Day",
  `STS:Regulation1`     = "STS Detection × Session",
  `FFA:OFA`             = "FFA Detection × OFA Detection",
  `FFA:Regulation1`     = "FFA Detection × Session",
  `OFA:Regulation1`     = "OFA Detection × Session",
  `FFA:OFA:Regulation1` = "FFA Detection × OFA Detection × Session"
)

pred_labels_md7 <- c(
  `(Intercept)`         = "(Intercept)",
  FFA                   = "FFA Recognition",
  OFA                   = "OFA Recognition",
  STS                   = "STS Recognition",
  Regulation1           = "Session",
  Training_Day.L        = "Training Day",
  Detection_FRMS        = "Detection (frames)",
  `STS:Regulation1`     = "STS Recognition × Session",
  `FFA:OFA`             = "FFA Recognition × OFA Recognition",
  `FFA:Regulation1`     = "FFA Recognition × Session",
  `OFA:Regulation1`     = "OFA Recognition × Session",
  `FFA:OFA:Regulation1` = "FFA Recognition × OFA Recognition × Session"
)

# Table 1: Model 6 (STS) — unstandardized | standardized (refit)
tab_model(md6, md6_std_tab,
          dv.labels = c('Model 6 (STS): Detection (unstandardized)',
                        'Model 6 (STS): Detection (standardized)'),
          pred.labels = pred_labels_md6,
          show.se = TRUE,
          string.est = "Estimate", string.se = "SE",
          title = 'Model 6 (STS): Detection performance')

# Table 2: Model 7 (STS) — unstandardized | standardized (refit)
tab_model(md7, md7_std_tab,
          dv.labels = c('Model 7 (STS): Recognition (unstandardized)',
                        'Model 7 (STS): Recognition (standardized)'),
          pred.labels = pred_labels_md7,
          show.se = TRUE,
          string.est = "Estimate", string.se = "SE",
          title = 'Model 7 (STS): Recognition performance')


