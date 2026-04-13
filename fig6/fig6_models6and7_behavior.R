# =============================================================================
# Figure 6 — Models 6 & 7: Impact of NFB Self-Regulation on Behavioural
#                            Detection and Recognition Performance
# =============================================================================
#
# Paper:  ReFOFA (NCOMMS-24-11279)
# Sample: Experimental group only, face stimuli only (N = 22, 1494 trials)
#
# Models:
#   Model 6 (md6): DV = detection_frms (detection performance in frames)
#                  Predictors: FFA detection response × OFA detection response
#                              × Session, Training Day
#
#   Model 7 (md7): DV = recognition_frms (recognition performance in frames)
#                  Predictors: FFA recognition response × OFA recognition response
#                              × Session, Training Day, Detection_FRMS (covariate)
#
# Note: In the original raw script (model_6_and_7.R), the variable names are
#       swapped: md2 = Model 6 (detection), md1 = Model 7 (recognition).
#       This clean version uses md6/md7 to match paper numbering.
#
# Key predictors: task_detect_FFA/OFA_beta_s (Model 6),
#                 task_recogn_FFA/OFA_beta_s (Model 7), all z-scored
# Interactions:   FFA × OFA × Session (three-way)
# Random effects: (1|subID) + (1|stimID)
#
# Output:
#   Left panel  — Dot-whisker coefficients (Models 6 & 7 overlaid)
#   Right top   — Detection: predicted effects of FFA, OFA, FFA×OFA
#   Right bottom— Recognition: predicted effects of FFA×Session, OFA×Session,
#                  FFA×OFA
# =============================================================================

rm(list = ls())

# --- 0. Settings --------------------------------------------------------------
# Set to TRUE to plot and compute emtrends from standardized models (z-scored
# predictors. FALSE uses raw-scale predictors.
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
library(ggrastr)     # rasterise() for raw data overlay
library(ggExtra)     # ggMarginal() for density margins


# --- 2. Load and prepare data -----------------------------------------------

# Set working directory to this script's location, then go up to root
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # Fallback: assumes the user has opened R in the fig6/ directory
}
setwd("..")

# Read source data (already post-processed and labeled) from the Source Data file
# skip = 1 to skip the description row above the header
df <- read_excel("sourceData.xlsx", sheet = "Fig. 6", skip = 1)

# Re-scale to ensure scale object attributes are available for lmer
df$task_detect_FFA_beta_s <- scale(df$task_detect_FFA_beta)
df$task_detect_OFA_beta_s <- scale(df$task_detect_OFA_beta)
df$task_recogn_FFA_beta_s <- scale(df$task_recogn_FFA_beta)
df$task_recogn_OFA_beta_s <- scale(df$task_recogn_OFA_beta)


# --- 3. Prepare model-specific data frames -----------------------------------

# Model 6 (detection): rename task detection betas to generic FFA/OFA
df_det <- df %>%
  dplyr::rename(
    FFA = task_detect_FFA_beta_s,
    OFA = task_detect_OFA_beta_s,
    Detection_FRMS = detection_frms
  )

# Model 7 (recognition): rename task recognition betas to generic FFA/OFA
df_rec <- df %>%
  dplyr::rename(
    FFA = task_recogn_FFA_beta_s,
    OFA = task_recogn_OFA_beta_s,
    Detection_FRMS = detection_frms,
    Recognition_FRMS = recognition_frms
  )

# Set factor levels and contrasts
# Regulation: OFA as first level, FFA as second
# Sum-to-zero contrast: OFA = -0.5, FFA = +0.5 (consistent with fig3/fig5)
df_det$Regulation   <- factor(df_det$Regulation, ordered = FALSE, levels = c("OFA", "FFA"))
df_det$Training_Day <- factor(df_det$Training_Day, ordered = TRUE)
contrasts(df_det$Regulation) <- c(-.5, .5)

df_rec$Regulation   <- factor(df_rec$Regulation, ordered = FALSE, levels = c("OFA", "FFA"))
df_rec$Training_Day <- factor(df_rec$Training_Day, ordered = TRUE)
contrasts(df_rec$Regulation) <- c(-.5, .5)


# --- 4. Fit the LMMs --------------------------------------------------------

# Model 6: Detection performance
# Does task-evoked FFA/OFA activity during face detection predict
# detection speed (frames), modulated by training session?
md6 <- lmer(
  Detection_FRMS ~ FFA + OFA +
    Regulation + Training_Day +
    FFA * OFA * Regulation +
    (1 | subID) + (1 | stimID),
  data = df_det,
  REML = TRUE
)

summary(md6)

# Model 7: Recognition performance
# Does task-evoked FFA/OFA activity during face recognition predict
# recognition speed (frames), with detection speed as covariate?
md7 <- lmer(
  Recognition_FRMS ~ FFA + OFA +
    Regulation + Training_Day + Detection_FRMS +
    FFA * OFA +
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

md6_std_ci <- std_with_ci(md6, "Model 6: Detection performance")
md7_std_ci <- std_with_ci(md7, "Model 7: Recognition performance")

# Select models for emtrends and plotting based on USE_STD flag
if (USE_STD) {
  md6_plot <- md6_std_obj
  md7_plot <- md7_std_obj
  var_FFA  <- "z.FFA"
  var_OFA  <- "z.OFA"
  var_Reg  <- "c.Regulation"
  cat("\n>>> Using STANDARDIZED models for slopes and plots <<<\n")
} else {
  md6_plot <- md6
  md7_plot <- md7
  var_FFA  <- "FFA"
  var_OFA  <- "OFA"
  var_Reg  <- "Regulation"
  cat("\n>>> Using RAW models for slopes and plots <<<\n")
}


# --- 6. Simple slopes analysis (Holm-corrected) ------------------------------

# 6a. Main effect slopes — panels (a) and (b): no correction needed
cat("\n--- Panel (a): FFA → Detection (main effect, no correction) ---\n")
print(test(emtrends(md6_plot, ~ 1, var = var_FFA)))

cat("\n--- Panel (b): OFA → Detection (main effect, no correction) ---\n")
print(test(emtrends(md6_plot, ~ 1, var = var_OFA)))

# 6b. Model 6 — Detection: 3 post-hoc slopes (panel c), Holm-corrected
cat("\n\n=== Model 6 (Detection): Panel (c) FFA × OFA — Holm-corrected for 3 tests ===\n")
ofa_at <- setNames(list(c(-3, 0, 3)), var_OFA)
t_c <- as.data.frame(test(emtrends(md6_plot, as.formula(paste("~", var_OFA)),
                                    var = var_FFA, at = ofa_at)))
names(t_c)[grep("trend", names(t_c))] <- "slope"
t_c$panel <- "c"
t_c$p.holm <- p.adjust(t_c$p.value, method = "holm")
t_c$sig <- ifelse(t_c$p.holm < 0.001, "***",
            ifelse(t_c$p.holm < 0.01, "**",
            ifelse(t_c$p.holm < 0.05, "*", "ns.")))
ratio_col <- ifelse("t.ratio" %in% names(t_c), "t.ratio", "z.ratio")
print(t_c[, c("panel", var_OFA, "slope", "SE", ratio_col, "p.value", "p.holm", "sig")])

# 6c. Model 7 — Recognition: 10 post-hoc slopes (panels d, e, f), Holm-corrected
#     Panel f uses session-specific slopes (conditioned on Regulation) because
#     plot_model fixes Regulation at the first factor level (OFA session).
#     All session-specific slopes are corrected together as one family of 10.
cat("\n\n=== Model 7 (Recognition): Panels d, e, f — Holm-corrected for 10 tests ===\n")

# Panel (d): FFA slope by Session (2 slopes)
t_d <- as.data.frame(test(emtrends(md7_plot, as.formula(paste("~", var_Reg)), var = var_FFA)))
names(t_d)[grep("trend", names(t_d))] <- "slope"
t_d$panel <- "d"
t_d$predictor <- "FFA"

# Panel (e): OFA slope by Session (2 slopes)
t_e <- as.data.frame(test(emtrends(md7_plot, as.formula(paste("~", var_Reg)), var = var_OFA)))
names(t_e)[grep("trend", names(t_e))] <- "slope"
t_e$panel <- "e"
t_e$predictor <- "OFA"

# Panel (f): FFA slope at OFA = -3, 0, 3, by Session (6 slopes)
# plot_model shows the OFA-session view; we test both sessions for completeness
reg_formula <- as.formula(paste("~", var_OFA, "+", var_Reg))
t_f <- as.data.frame(test(emtrends(md7_plot, reg_formula,
                                    var = var_FFA, at = ofa_at)))
names(t_f)[grep("trend", names(t_f))] <- "slope"
t_f$panel <- "f"
t_f$predictor <- "FFA×OFA"

# Combine all 10 slopes and apply Holm correction across the family
combined_md7 <- bind_rows(t_d, t_e, t_f)
combined_md7$p.holm <- p.adjust(combined_md7$p.value, method = "holm")
combined_md7$sig <- ifelse(combined_md7$p.holm < 0.001, "***",
                    ifelse(combined_md7$p.holm < 0.01, "**",
                    ifelse(combined_md7$p.holm < 0.05, "*", "ns.")))

ratio_col <- ifelse("t.ratio" %in% names(combined_md7), "t.ratio", "z.ratio")

# Print with panel info so each slope is identifiable
print_cols <- c("panel", "predictor")
if (var_Reg %in% names(combined_md7)) print_cols <- c(print_cols, var_Reg)
if (var_OFA %in% names(combined_md7)) print_cols <- c(print_cols, var_OFA)
print_cols <- c(print_cols, "slope", "SE", ratio_col, "p.value", "p.holm", "sig")
print(combined_md7[, print_cols])

# 6d. Full simple_slopes decomposition (for reference, always on raw models)
simple_slopes(md6, levels = list(
  FFA        = c(-3, 0, 3, "sstest"),
  OFA        = c(-3, 0, 3, "sstest"),
  Regulation = c("FFA", "OFA", "sstest")
), confint = FALSE, ci.width = 0.95)

simple_slopes(md7, levels = list(
  FFA        = c(-3, 0, 3, "sstest"),
  OFA        = c(-3, 0, 3, "sstest"),
  Regulation = c("FFA", "OFA", "sstest")
), confint = FALSE, ci.width = 0.95)


# --- 7. Plot — Left panel: Coefficient dot-whisker plot ----------------------

# Predictor labels for the combined plot
predictor_labels_std <- c(
  z.FFA                      = "FFA",
  z.OFA                      = "OFA",
  c.Regulation               = "Session",
  c.Training_Day             = "Training Day",
  z.Detection_FRMS           = "Detection FRMS",
  `z.FFA:z.OFA`              = "FFA * OFA",
  `z.FFA:c.Regulation`       = "FFA * Session",
  `z.OFA:c.Regulation`       = "OFA * Session",
  `z.FFA:z.OFA:c.Regulation` = "FFA * OFA * Session"
)

# Color scheme: green = detection, gold = recognition
color_det  <- "#597e59"
color_rec  <- "#dbbb5b"

colorSpecs <- c(color_det, color_rec,
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
  vars_order   = c("z.FFA", "z.OFA", "c.Regulation", "c.Training_Day",
                   "z.FFA:z.OFA", "z.FFA:c.Regulation", "z.OFA:c.Regulation",
                   "z.FFA:z.OFA:c.Regulation")
) %>%
  relabel_predictors(predictor_labels_std) +
  scale_color_manual(name = "Behavior", labels = c("Detection", "Recognition")) +
  xlab("Standardized Coefficient (B)") +
  ylab("") +
  ggtitle("Model 6 and 7: Detection and Recognition Performance") +
  theme_article()


# --- 8. Plot — Right panels: Predicted marginal effects ----------------------
# Uses md6_plot/md7_plot and var_FFA/var_OFA set by USE_STD flag in section 5

# Color palettes
color_det_line <- c("#597e59", "#a2bca2")  # dark/light green
color_det_3    <- c("#a2bca2", "#7ba27b", "#597e59")  # 3-level OFA
color_rec_line <- c("#dbbb5b", "#e4cd82")  # dark/light gold
color_rec_3    <- c("#e4cd82", "#d7b64a", "#dbbb5b")  # 3-level OFA

# Build terms strings dynamically based on USE_STD flag
ofa_3lvl <- paste0(var_OFA, " [-3, 0, 3]")

# --- Detection performance (top right) ---

# p1: FFA → Detection FRMS
p1 <- plot_model(
  md6_plot, type = "pred", terms = c(var_FFA),
  legend.title = "Session", colors = color_det[1],
  show.legend = FALSE, grid = FALSE
) +
  theme_sjplot2() +
  geom_line(size = 1) +
  xlab(expression(beta)) +
  ylab("frms") +
  ggtitle("FFA Detection Response") +
  theme(
    plot.title   = element_text(size = 11),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold")
  )

# p2: OFA → Detection FRMS
p2 <- plot_model(
  md6_plot, type = "pred", terms = c(var_OFA),
  legend.title = "Session", colors = color_det[1],
  show.legend = FALSE, grid = FALSE
) +
  theme_sjplot2() +
  geom_line(size = 1) +
  xlab(expression(beta)) +
  ylab("frms") +
  ggtitle("OFA Detection Response") +
  theme(
    plot.title   = element_text(size = 11),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold")
  )

# p3: FFA × OFA interaction → Detection FRMS
p3 <- plot_model(
  md6_plot, type = "pred", terms = c(var_FFA, ofa_3lvl),
  legend.title = "OFA Detection Resp.", colors = color_det_3,
  show.legend = FALSE, show.data = FALSE, grid = FALSE
) +
  theme_sjplot2() +
  geom_line(size = 1) +
  xlab(expression(beta)) +
  ylab("frms") +
  ggtitle("FFA Detection x OFA Detection") +
  theme(
    plot.title   = element_text(size = 11),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold")
  )

# --- Recognition performance (bottom right) ---

# p4: FFA × Session → Recognition FRMS
p4 <- plot_model(
  md7_plot, type = "pred", terms = c(var_FFA, var_Reg),
  legend.title = "Session", colors = color_rec_line,
  show.legend = FALSE, grid = FALSE
) +
  theme_sjplot2() +
  geom_line(size = 1) +
  xlab(expression(beta)) +
  ylab("frms") +
  ggtitle("FFA Recognition x Session") +
  theme(
    plot.title   = element_text(size = 11),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold")
  )

# p5: OFA × Session → Recognition FRMS
p5 <- plot_model(
  md7_plot, type = "pred", terms = c(var_OFA, var_Reg),
  legend.title = "OFA", colors = color_rec_line,
  show.legend = FALSE, show.data = FALSE, grid = FALSE
) +
  theme_sjplot2() +
  geom_line(size = 1) +
  xlab(expression(beta)) +
  ylab("frms") +
  ggtitle("OFA Recognition x Session") +
  theme(
    plot.title   = element_text(size = 11),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold")
  )

# p6: FFA × OFA interaction → Recognition FRMS
p6 <- plot_model(
  md7_plot, type = "pred", terms = c(var_FFA, ofa_3lvl),
  legend.title = "OFA Recognition", colors = color_rec_3,
  show.legend = FALSE, show.data = FALSE, grid = FALSE
) +
  theme_sjplot2() +
  geom_line(size = 1) +
  xlab(expression(beta)) +
  ylab("frms") +
  ggtitle("FFA Recognition by OFA Recognition") +
  theme(
    plot.title   = element_text(size = 11),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold")
  )

# --- 9. Combine panels -------------------------------------------------------

# Sub-layouts matching original script structure
left_col          <- plot_coef
top_right_left    <- plot_grid(p1, p2, ncol = 1)
top_right         <- plot_grid(top_right_left, p3, ncol = 2)

bottom_right_left <- plot_grid(p4, p5, ncol = 1)
bottom_right      <- plot_grid(bottom_right_left, p6, ncol = 2)

right_side <- plot_grid(top_right, bottom_right, ncol = 1)
final_figure <- plot_grid(left_col, right_side, ncol = 2, rel_widths = c(1, 2))

print(final_figure)

# Alternative using patchwork (as in original script):
# finPlot <- (plot_coef | ((p1 / p2) | p3) / (p4 / p5 | p6))
# print(finPlot)


# --- 10. Summary tables (for supplementary materials) -----------------------
# Full model output tables for Models 6 and 7. Each table shows the
# unstandardized model alongside a formally standardized (refit) version so
# that standardized coefficients match those reported in the main text and
# dot-whisker plot (from arm::standardize).
#
# Note: tab_model's built-in show.std = "std2" applies POST-HOC rescaling
# (raw coefficient × 2·SD of predictor), which diverges from arm::standardize
# for interaction terms because arm::standardize actually refits the model on
# rescaled data. To keep text, figures, and table values consistent, we build
# "tab-friendly" standardized models by manually rescaling the data with
# arm::rescale() (same math as arm::standardize but preserves variable names,
# which tab_model requires). Deviation contrasts on Regulation match
# arm::standardize's c. coding so interaction signs align.

# --- Build rescaled data for tab_model ---
# NOTE: only rescale PREDICTORS — never the DV. arm::standardize() leaves the
# DV on its original scale by default (standardize.y = FALSE), so to match its
# output we must do the same here.
df_det_std <- df_det
df_det_std$FFA            <- arm::rescale(df_det$FFA)
df_det_std$OFA            <- arm::rescale(df_det$OFA)
# Detection_FRMS is the DV in Model 6 — leave untouched
contrasts(df_det_std$Regulation) <- c(-.5, .5)

df_rec_std <- df_rec
df_rec_std$FFA              <- arm::rescale(df_rec$FFA)
df_rec_std$OFA              <- arm::rescale(df_rec$OFA)
df_rec_std$Detection_FRMS   <- arm::rescale(df_rec$Detection_FRMS)  # predictor (covariate) in Model 7
# Recognition_FRMS is the DV in Model 7 — leave untouched
contrasts(df_rec_std$Regulation) <- c(-.5, .5)

# Refit on rescaled data (same formula structure as section 4)
md6_std_tab <- lmer(
  Detection_FRMS ~ FFA + OFA +
    Regulation + Training_Day +
    FFA * OFA * Regulation +
    (1 | subID) + (1 | stimID),
  data = df_det_std,
  REML = TRUE
)

md7_std_tab <- lmer(
  Recognition_FRMS ~ FFA + OFA +
    Regulation + Training_Day + Detection_FRMS +
    FFA * OFA +
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
  Regulation1           = "Session",
  Training_Day.L        = "Training Day",
  `FFA:OFA`             = "FFA Detection × OFA Detection",
  `FFA:Regulation1`     = "FFA Detection × Session",
  `OFA:Regulation1`     = "OFA Detection × Session",
  `FFA:OFA:Regulation1` = "FFA Detection × OFA Detection × Session"
)

pred_labels_md7 <- c(
  `(Intercept)`         = "(Intercept)",
  FFA                   = "FFA Recognition",
  OFA                   = "OFA Recognition",
  Regulation1           = "Session",
  Training_Day.L        = "Training Day",
  Detection_FRMS        = "Detection (frames)",
  `FFA:OFA`             = "FFA Recognition × OFA Recognition",
  `FFA:Regulation1`     = "FFA Recognition × Session",
  `OFA:Regulation1`     = "OFA Recognition × Session",
  `FFA:OFA:Regulation1` = "FFA Recognition × OFA Recognition × Session"
)

# Table 1: Model 6 — unstandardized | standardized (refit)
tab_model(md6, md6_std_tab,
          dv.labels = c('Model 6: Detection (unstandardized)',
                        'Model 6: Detection (standardized)'),
          pred.labels = pred_labels_md6,
          show.se = TRUE,
          string.est = "Estimate", string.se = "SE",
          title = 'Model 6: Detection performance')

# Table 2: Model 7 — unstandardized | standardized (refit)
tab_model(md7, md7_std_tab,
          dv.labels = c('Model 7: Recognition (unstandardized)',
                        'Model 7: Recognition (standardized)'),
          pred.labels = pred_labels_md7,
          show.se = TRUE,
          string.est = "Estimate", string.se = "SE",
          title = 'Model 7: Recognition performance')
