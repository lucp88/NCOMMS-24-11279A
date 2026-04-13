# =============================================================================
# Figure 5 — Models 2–5: Impact of NFB Self-Regulation on Face Detection
#                         and Recognition
# =============================================================================
#
# Paper:  ReFOFA (NCOMMS-24-11279)
# Sample: Experimental group only, face stimuli only (N = 22, 1494 trials)
#
# Models:
#   Model 2 (md1): DV = FFA detection response (task_detect_FFA_beta_s)
#   Model 3 (md3): DV = OFA detection response (task_detect_OFA_beta_s)
#   Model 4 (md2): DV = FFA recognition response (task_recogn_FFA_beta_s)
#   Model 5 (md4): DV = OFA recognition response (task_recogn_OFA_beta_s)
#
# Key predictors: regFFA_s, regOFA_s (NFB regulation PSC, z-scored)
# Interactions:   regFFA × regOFA × Session (three-way)
# Covariates:     runID, Training_Day, detection_frms,
#                 + detection beta for recognition models
# Random effects: (1|subID) + (1|stimID)
#
# Output:
#   Left panel  — Dot-whisker coefficients (Models 2 & 3: Detection)
#   Right panel — Dot-whisker coefficients (Models 4 & 5: Recognition)
#   Center      — Predicted marginal effects with marginal density plots
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
library(patchwork)   # plot composition
library(cowplot)     # plot_grid() for combining ggMarginal plots
library(ggExtra)     # ggMarginal() for marginal density plots
library(ggrastr)     # rasterise() for efficient raw data overlay


# --- 2. Working directory & data ---------------------------------------------
# Set working directory to the script's location (works in RStudio)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {}
setwd("..")  # move up to "code and source data" root

# Read from Source Data Excel (skip row 1 = description)
df <- read_excel("sourceData.xlsx", sheet = "Fig. 5", skip = 1)

cat("Data loaded:", nrow(df), "trials,", length(unique(df$subID)), "subjects\n")


# --- 3. Factor coding --------------------------------------------------------
# Sum-to-zero (deviation) contrasts: −0.5 / +0.5
# Levels ordered so the contrast direction matches the manuscript

df$Regulation <- factor(df$Regulation, levels = c("OFA", "FFA"))
contrasts(df$Regulation) <- c(-0.5, 0.5)

df$Training_Day <- factor(df$Training_Day, ordered = TRUE)

# Mean-center runID for consistent conditional effect interpretation
df$runID <- df$runID - mean(df$runID)

# Ensure scaled variables are numeric (read_excel may import as list-columns)
scale_cols <- c("regFFA_s", "regOFA_s", "diffROI_s",
                "task_detect_FFA_beta_s", "task_detect_OFA_beta_s",
                "task_recogn_FFA_beta_s", "task_recogn_OFA_beta_s",
                "task_FFA_beta_s", "task_OFA_beta_s")
for (col in scale_cols) {
  df[[col]] <- as.numeric(df[[col]])
}


# --- 4. Model fitting --------------------------------------------------------
# Note on model numbering:
#   Paper Model 2 → md1 (FFA detection DV)
#   Paper Model 3 → md3 (OFA detection DV)
#   Paper Model 4 → md2 (FFA recognition DV)
#   Paper Model 5 → md4 (OFA recognition DV)
# This follows the original analysis script naming convention.

cat("\n--- Fitting Model 2: FFA detection ---\n")
md1 <- lmer(task_detect_FFA_beta_s ~ regFFA_s + regOFA_s + Regulation +
              runID + Training_Day + detection_frms +
              regFFA_s * regOFA_s * Regulation +
              (1|subID) + (1|stimID),
            data = df, REML = TRUE)
print(summary(md1))

cat("\n--- Fitting Model 3: OFA detection ---\n")
md3 <- lmer(task_detect_OFA_beta_s ~ regFFA_s + regOFA_s + Regulation +
              runID + Training_Day + detection_frms +
              regFFA_s * regOFA_s * Regulation +
              (1|subID) + (1|stimID),
            data = df, REML = TRUE)
print(summary(md3))

cat("\n--- Fitting Model 4: FFA recognition ---\n")
md2 <- lmer(task_recogn_FFA_beta_s ~ regFFA_s + regOFA_s + Regulation +
              runID + Training_Day + detection_frms +
              task_detect_FFA_beta_s +
              regFFA_s * regOFA_s * Regulation +
              (1|subID) + (1|stimID),
            data = df, REML = TRUE)
print(summary(md2))

cat("\n--- Fitting Model 5: OFA recognition ---\n")
md4 <- lmer(task_recogn_OFA_beta_s ~ regFFA_s + regOFA_s + Regulation +
              runID + Training_Day + detection_frms +
              task_detect_OFA_beta_s +
              regFFA_s * regOFA_s * Regulation +
              (1|subID) + (1|stimID),
            data = df, REML = TRUE)
print(summary(md4))


# --- 4a. Standardized coefficients with CIs ----------------------------------
cat("\n--- Standardized coefficients (arm::standardize) with 95% Wald CIs ---\n")

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

md1_std_ci <- std_with_ci(md1, "Model 2: FFA detection")
md3_std_ci <- std_with_ci(md3, "Model 3: OFA detection")
md2_std_ci <- std_with_ci(md2, "Model 4: FFA recognition")
md4_std_ci <- std_with_ci(md4, "Model 5: OFA recognition")

# Keep arm::standardize model objects for dwplot (by_2sd = FALSE)
md1_std_obj <- standardize(md1)
md3_std_obj <- standardize(md3)
md2_std_obj <- standardize(md2)
md4_std_obj <- standardize(md4)


# --- 4b. Simple slopes via emtrends (Holm-corrected per model) ---------------
# For each model, we test 4 slopes: 2 predictors (regFFA, regOFA) × 2 sessions.
# Holm correction is applied across the 4 slopes within each model.
holm_correct_slopes <- function(model, model_label) {
  t_ffa <- as.data.frame(summary(emtrends(model, ~ Regulation, var = "regFFA_s"),
                                 infer = c(TRUE, TRUE), level = 0.95))
  t_ofa <- as.data.frame(summary(emtrends(model, ~ Regulation, var = "regOFA_s"),
                                 infer = c(TRUE, TRUE), level = 0.95))
  names(t_ffa)[grep("trend", names(t_ffa))] <- "slope"
  names(t_ofa)[grep("trend", names(t_ofa))] <- "slope"
  t_ffa$predictor <- "regFFA"
  t_ofa$predictor <- "regOFA"
  combined <- rbind(t_ffa, t_ofa)
  combined$p.holm <- p.adjust(combined$p.value, method = "holm")
  combined$sig <- ifelse(combined$p.holm < 0.001, "***",
                         ifelse(combined$p.holm < 0.01, "**",
                                ifelse(combined$p.holm < 0.05, "*", "ns.")))
  ratio_col <- ifelse("t.ratio" %in% names(combined), "t.ratio", "z.ratio")
  cat("\n===", model_label, "— Simple slopes (Holm-corrected for 4 tests) ===\n")
  print(combined[, c("predictor", "Regulation", "slope", "SE", "df",
                     "lower.CL", "upper.CL", ratio_col,
                     "p.value", "p.holm", "sig")])
  return(combined)
}

slopes_md1 <- holm_correct_slopes(md1, "Model 2: FFA detection")
slopes_md3 <- holm_correct_slopes(md3, "Model 3: OFA detection")
slopes_md2 <- holm_correct_slopes(md2, "Model 4: FFA recognition")
slopes_md4 <- holm_correct_slopes(md4, "Model 5: OFA recognition")


# --- 5. Coefficient plots (dot-whisker) --------------------------------------
# Uses arm::standardize model objects with by_2sd = FALSE so that plotted
# coefficients match the standardized values reported in the text.

# Color scheme: blue = FFA ROI, red = OFA ROI
colorSpecs <- c('#0F4889', '#D84545',
                '#0F4889', '#D84545',
                '#0F4889', '#D84545',
                '#0F4889', '#D84545',
                '#0F4889', '#D84545',
                '#0F4889', '#D84545',
                '#0F4889', '#D84545')

predictor_labels_std <- c(
  z.regFFA_s                              = "regFFA",
  z.regOFA_s                              = "regOFA",
  c.Regulation                            = "Session",
  `z.regFFA_s:z.regOFA_s`                 = "regFFA*regOFA",
  `z.regFFA_s:c.Regulation`               = "regFFA*Session",
  `z.regOFA_s:c.Regulation`               = "regOFA*Session",
  `z.regFFA_s:z.regOFA_s:c.Regulation`    = "regFFA*regOFA*Session"
)

variable_order_std <- c("z.regFFA_s", "z.regOFA_s", "c.Regulation",
                        "z.regFFA_s:z.regOFA_s",
                        "z.regFFA_s:c.Regulation", "z.regOFA_s:c.Regulation",
                        "z.regFFA_s:z.regOFA_s:c.Regulation")

# Detection: Models 2 & 3
p1 <- dwplot(list(md3_std_obj, md1_std_obj),
             by_2sd = FALSE, ci = 0.95,
             dot_args = list(size = 4, shape = 21, fill = colorSpecs, color = "black", stroke = 1.2),
             whisker_args = list(size = 1, color = colorSpecs),
             style = "dotwhisker",
             vline = geom_vline(xintercept = 0, colour = "grey60",
                                size = 1, linetype = 2),
             vars_order = variable_order_std) %>%
  relabel_predictors(predictor_labels_std) +
  scale_color_manual(name = "ROI", labels = c("FFA", "OFA")) +
  xlab("Standardized Coefficient (ß)") + ylab("") +
  ggtitle("Model 2 and 3: Detection") +
  theme_blank() +
  theme_article()

# Recognition: Models 4 & 5
p2 <- dwplot(list(md4_std_obj, md2_std_obj),
             by_2sd = FALSE, ci = 0.95,
             dot_args = list(size = 4, shape = 21, fill = colorSpecs, color = "black", stroke = 1.2),
             whisker_args = list(size = 1, color = colorSpecs),
             style = "dotwhisker",
             vline = geom_vline(xintercept = 0, colour = "grey60",
                                size = 1, linetype = 2),
             vars_order = variable_order_std) %>%
  relabel_predictors(predictor_labels_std) +
  scale_color_manual(name = "ROI", labels = c("FFA", "OFA")) +
  xlab("Standardized Coefficient (ß)") + ylab("") +
  ggtitle("Model 4 and 5: Recognition") +
  theme_blank() +
  theme_article()


# --- 6. Marginal effects with density margins --------------------------------
colorSpecs_ffa_2 <- c('#378AD3', '#0C4787')   # light/dark blue
colorSpecs_ofa_2 <- c('#FF8585', '#D74545')   # light/dark red
ylimits <- c(-1, 2)
densPlot_size <- 8

# -- 6a. DETECTION panels (pl1–pl4) --

# pl1: regFFA → FFA detection
pl1_base <- plot_model(md1, type = "pred", terms = c('regFFA_s', 'Regulation'),
                       legend.title = "Session", show.legend = FALSE,
                       colors = colorSpecs_ffa_2) +
  theme_sjplot2() + geom_line(size = 1) +
  rasterise(geom_point(data = df,
                       aes(x = regFFA_s, y = task_detect_FFA_beta_s, color = Regulation),
                       alpha = 0.15, size = 0.5, inherit.aes = FALSE)) +
  xlab('regFFA (PSC)') + ylab(expression(paste("Detection Resp ", (beta)))) +
  coord_cartesian(ylim = ylimits) + ggtitle('FFA')
pl1 <- ggMarginal(pl1_base, type = "density", groupColour = TRUE,
                  groupFill = TRUE, alpha = 0.3, size = densPlot_size)

# pl2: regFFA → OFA detection
pl2_base <- plot_model(md3, type = "pred", terms = c('regFFA_s', 'Regulation'),
                       legend.title = "Session", show.legend = FALSE,
                       colors = colorSpecs_ofa_2) +
  theme_sjplot2() + geom_line(size = 1) +
  rasterise(geom_point(data = df,
                       aes(x = regFFA_s, y = task_detect_OFA_beta_s, color = Regulation),
                       alpha = 0.15, size = 0.5, inherit.aes = FALSE)) +
  xlab('regFFA (PSC)') + ylab(expression(paste("Detection Resp ", (beta)))) +
  coord_cartesian(ylim = ylimits) + ggtitle('OFA')
pl2 <- ggMarginal(pl2_base, type = "density", groupColour = TRUE,
                  groupFill = TRUE, alpha = 0.3, size = densPlot_size)

# pl3: regOFA → FFA detection
pl3_base <- plot_model(md1, type = "pred", terms = c('regOFA_s', 'Regulation'),
                       legend.title = "Session", show.legend = FALSE,
                       colors = colorSpecs_ffa_2) +
  theme_sjplot2() + geom_line(size = 1) +
  rasterise(geom_point(data = df,
                       aes(x = regOFA_s, y = task_detect_FFA_beta_s, color = Regulation),
                       alpha = 0.15, size = 0.5, inherit.aes = FALSE)) +
  xlab('regOFA (PSC)') + ylab(expression(paste("Detection Resp ", (beta)))) +
  coord_cartesian(ylim = ylimits) + ggtitle('FFA')
pl3 <- ggMarginal(pl3_base, type = "density", groupColour = TRUE,
                  groupFill = TRUE, alpha = 0.3, size = densPlot_size)

# pl4: regOFA → OFA detection
pl4_base <- plot_model(md3, type = "pred", terms = c('regOFA_s', 'Regulation'),
                       legend.title = "Session", show.legend = FALSE,
                       colors = colorSpecs_ofa_2) +
  theme_sjplot2() + geom_line(size = 1) +
  rasterise(geom_point(data = df,
                       aes(x = regOFA_s, y = task_detect_OFA_beta_s, color = Regulation),
                       alpha = 0.15, size = 0.5, inherit.aes = FALSE)) +
  xlab('regOFA (PSC)') + ylab(expression(paste("Detection Resp ", (beta)))) +
  coord_cartesian(ylim = ylimits) + ggtitle('OFA')
pl4 <- ggMarginal(pl4_base, type = "density", groupColour = TRUE,
                  groupFill = TRUE, alpha = 0.3, size = densPlot_size)


# -- 6b. RECOGNITION panels (pl5–pl8) --

# pl5: regFFA → FFA recognition
pl5_base <- plot_model(md2, type = "pred", terms = c('regFFA_s', 'Regulation'),
                       legend.title = "Session", show.legend = FALSE,
                       colors = colorSpecs_ffa_2) +
  theme_sjplot2() + geom_line(size = 1) +
  rasterise(geom_point(data = df,
                       aes(x = regFFA_s, y = task_recogn_FFA_beta_s, color = Regulation),
                       alpha = 0.15, size = 0.5, inherit.aes = FALSE)) +
  xlab('regFFA (PSC)') + ylab(expression(paste("Recognition Resp ", (beta)))) +
  coord_cartesian(ylim = ylimits) + ggtitle('FFA')
pl5 <- ggMarginal(pl5_base, type = "density", groupColour = TRUE,
                  groupFill = TRUE, alpha = 0.3, size = densPlot_size)

# pl6: regFFA → OFA recognition
pl6_base <- plot_model(md4, type = "pred", terms = c('regFFA_s', 'Regulation'),
                       legend.title = "Session", show.legend = FALSE,
                       colors = colorSpecs_ofa_2) +
  theme_sjplot2() + geom_line(size = 1) +
  rasterise(geom_point(data = df,
                       aes(x = regFFA_s, y = task_recogn_OFA_beta_s, color = Regulation),
                       alpha = 0.15, size = 0.5, inherit.aes = FALSE)) +
  xlab('regFFA (PSC)') + ylab(expression(paste("Recognition Resp ", (beta)))) +
  coord_cartesian(ylim = ylimits) + ggtitle('OFA')
pl6 <- ggMarginal(pl6_base, type = "density", groupColour = TRUE,
                  groupFill = TRUE, alpha = 0.3, size = densPlot_size)

# pl7: regOFA → FFA recognition
pl7_base <- plot_model(md2, type = "pred", terms = c('regOFA_s', 'Regulation'),
                       legend.title = "Session", show.legend = FALSE,
                       colors = colorSpecs_ffa_2) +
  theme_sjplot2() + geom_line(size = 1) +
  rasterise(geom_point(data = df,
                       aes(x = regOFA_s, y = task_recogn_FFA_beta_s, color = Regulation),
                       alpha = 0.15, size = 0.5, inherit.aes = FALSE)) +
  xlab('regOFA (PSC)') + ylab(expression(paste("Recognition Resp ", (beta)))) +
  coord_cartesian(ylim = ylimits) + ggtitle('FFA')
pl7 <- ggMarginal(pl7_base, type = "density", groupColour = TRUE,
                  groupFill = TRUE, alpha = 0.3, size = densPlot_size)

# pl8: regOFA → OFA recognition
pl8_base <- plot_model(md4, type = "pred", terms = c('regOFA_s', 'Regulation'),
                       legend.title = "Session", show.legend = FALSE,
                       colors = colorSpecs_ofa_2) +
  theme_sjplot2() + geom_line(size = 1) +
  rasterise(geom_point(data = df,
                       aes(x = regOFA_s, y = task_recogn_OFA_beta_s, color = Regulation),
                       alpha = 0.15, size = 0.5, inherit.aes = FALSE)) +
  xlab('regOFA (PSC)') + ylab(expression(paste("Recognition Resp ", (beta)))) +
  coord_cartesian(ylim = ylimits) + ggtitle('OFA')
pl8 <- ggMarginal(pl8_base, type = "density", groupColour = TRUE,
                  groupFill = TRUE, alpha = 0.3, size = densPlot_size)


# --- 7. Assemble figure ------------------------------------------------------

# Detection half: regFFA panels | regOFA panels | coefficient plot
figure1 <- cowplot::plot_grid(
  cowplot::plot_grid(pl1, pl2, ncol = 1),   # regFFA → FFA/OFA detection
  cowplot::plot_grid(pl3, pl4, ncol = 1),   # regOFA → FFA/OFA detection
  NULL,                             # spacing
  p1,                               # dot-whisker
  ncol = 4,
  rel_widths = c(1, 1, 0.3, 1.2)
)

# Recognition half: coefficient plot | regFFA panels | regOFA panels
figure2 <- cowplot::plot_grid(
  p2,                               # dot-whisker
  NULL,                             # spacing
  cowplot::plot_grid(pl5, pl6, ncol = 1),   # regFFA → FFA/OFA recognition
  cowplot::plot_grid(pl7, pl8, ncol = 1),   # regOFA → FFA/OFA recognition
  ncol = 4,
  rel_widths = c(1.2, 0.3, 1, 1)
)


# --- 8. Summary table (for SFig. 2) ------------------------------------------
# Full unstandardized model output tables for the supplementary materials
# (SFig. 2). Raw coefficients are reported in the original units of the model
# (z-scored task beta responses) and include all covariates, random effects,
# and fit statistics. Standardized coefficients are reported in the main text
# and dot-whisker plots (from arm::standardize) for visual comparison across
# predictors on a common scale.

myDV <- c('Model 2: FFA detection', 'Model 3: OFA Detection',
          'Model 4: FFA Recognition', 'Model 5: OFA Recognition')

tab_model(md1, md3, md2, md4,
          dv.labels = myDV,
          show.se = TRUE,
          string.est = "Estimate", string.se = "SE",
          title = 'Linear Mixed Models: Self-Regulation → Task Responses')
