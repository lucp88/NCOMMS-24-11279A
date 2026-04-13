# =============================================================================
# Figure 6 — Output Log Generator
# =============================================================================
# Run this AFTER executing fig6_models6and7_behavior.R (sections 1–6).
# It captures all model summaries, emtrends, and Holm-corrected slopes
# to a single text file in the fig6/ directory.
#
# Requires: md6, md7, md6_std_obj, md7_std_obj, md6_plot, md7_plot,
#           var_FFA, var_OFA, var_Reg, USE_STD  (all set by the main script)
# =============================================================================

output_file <- "fig6_output_log.txt"
sink(output_file)

cat("=============================================================================\n")
cat("  Figure 6 — Models 6 & 7: Full Statistical Output\n")
cat("  ReFOFA (NCOMMS-24-11279)\n")
cat("  Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("  Mode:", ifelse(USE_STD, "STANDARDIZED (z-scored predictors)",
                                "RAW (original scale)"), "\n")
cat("=============================================================================\n")


# --- Model 6: LMM summary (raw) ---------------------------------------------
cat("\n\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")
cat("  MODEL 6 — Detection Performance (raw scale)\n")
cat("  DV: Detection_FRMS ~ FFA * OFA * Regulation + Training_Day\n")
cat("  Random: (1|subID) + (1|stimID)\n")
cat("─────────────────────────────────────────────────────────────────────────────\n\n")
print(summary(md6))


# --- Model 6: Standardized coefficients -------------------------------------
cat("\n\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")
cat("  MODEL 6 — Standardized Coefficients (arm::standardize)\n")
cat("─────────────────────────────────────────────────────────────────────────────\n\n")
print(summary(md6_std_obj))


# --- Model 7: LMM summary (raw) ---------------------------------------------
cat("\n\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")
cat("  MODEL 7 — Recognition Performance (raw scale)\n")
cat("  DV: Recognition_FRMS ~ FFA * OFA * Regulation + Training_Day\n")
cat("       + Detection_FRMS (covariate)\n")
cat("  Random: (1|subID) + (1|stimID)\n")
cat("─────────────────────────────────────────────────────────────────────────────\n\n")
print(summary(md7))


# --- Model 7: Standardized coefficients -------------------------------------
cat("\n\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")
cat("  MODEL 7 — Standardized Coefficients (arm::standardize)\n")
cat("─────────────────────────────────────────────────────────────────────────────\n\n")
print(summary(md7_std_obj))


# --- Panel (a): FFA main effect on detection ---------------------------------
cat("\n\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")
cat("  PANEL (a): FFA → Detection (main effect, no correction needed)\n")
cat("─────────────────────────────────────────────────────────────────────────────\n\n")
print(test(emtrends(md6_plot, ~ 1, var = var_FFA)))


# --- Panel (b): OFA main effect on detection ---------------------------------
cat("\n\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")
cat("  PANEL (b): OFA → Detection (main effect, no correction needed)\n")
cat("─────────────────────────────────────────────────────────────────────────────\n\n")
print(test(emtrends(md6_plot, ~ 1, var = var_OFA)))


# --- Panel (c): FFA × OFA interaction — Detection (Holm, 3 tests) -----------
cat("\n\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")
cat("  PANEL (c): FFA slope at OFA = -3, 0, 3 — Detection\n")
cat("  Holm-corrected for 3 tests (Model 6 family)\n")
cat("─────────────────────────────────────────────────────────────────────────────\n\n")

ofa_at <- setNames(list(c(-3, 0, 3)), var_OFA)
t_c <- as.data.frame(test(emtrends(md6_plot, as.formula(paste("~", var_OFA)),
                                    var = var_FFA, at = ofa_at)))
names(t_c)[grep("trend", names(t_c))] <- "slope"
t_c$p.holm <- p.adjust(t_c$p.value, method = "holm")
t_c$sig <- ifelse(t_c$p.holm < 0.001, "***",
            ifelse(t_c$p.holm < 0.01, "**",
            ifelse(t_c$p.holm < 0.05, "*", "ns.")))
ratio_col <- ifelse("t.ratio" %in% names(t_c), "t.ratio", "z.ratio")
print(t_c[, c(var_OFA, "slope", "SE", ratio_col, "p.value", "p.holm", "sig")])


# --- Panels (d, e, f): Recognition slopes (Holm, 7 tests) -------------------
cat("\n\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")
cat("  PANELS (d, e, f): Recognition post-hoc slopes\n")
cat("  Holm-corrected for 7 tests (Model 7 family)\n")
cat("    d = FFA by Session (2 slopes)\n")
cat("    e = OFA by Session (2 slopes)\n")
cat("    f = FFA at OFA = -3, 0, 3 (3 slopes)\n")
cat("─────────────────────────────────────────────────────────────────────────────\n\n")

t_d <- as.data.frame(test(emtrends(md7_plot, as.formula(paste("~", var_Reg)), var = var_FFA)))
names(t_d)[grep("trend", names(t_d))] <- "slope"
t_d$panel <- "d"
t_d$predictor <- "FFA"

t_e <- as.data.frame(test(emtrends(md7_plot, as.formula(paste("~", var_Reg)), var = var_OFA)))
names(t_e)[grep("trend", names(t_e))] <- "slope"
t_e$panel <- "e"
t_e$predictor <- "OFA"

t_f <- as.data.frame(test(emtrends(md7_plot, as.formula(paste("~", var_OFA)),
                                    var = var_FFA, at = ofa_at)))
names(t_f)[grep("trend", names(t_f))] <- "slope"
t_f$panel <- "f"
t_f$predictor <- "FFA×OFA"

combined_md7 <- dplyr::bind_rows(t_d, t_e, t_f)
combined_md7$p.holm <- p.adjust(combined_md7$p.value, method = "holm")
combined_md7$sig <- ifelse(combined_md7$p.holm < 0.001, "***",
                    ifelse(combined_md7$p.holm < 0.01, "**",
                    ifelse(combined_md7$p.holm < 0.05, "*", "ns.")))
ratio_col <- ifelse("t.ratio" %in% names(combined_md7), "t.ratio", "z.ratio")

print_cols <- c("panel", "predictor")
if (var_Reg %in% names(combined_md7)) print_cols <- c(print_cols, var_Reg)
if (var_OFA %in% names(combined_md7)) print_cols <- c(print_cols, var_OFA)
print_cols <- c(print_cols, "slope", "SE", ratio_col, "p.value", "p.holm", "sig")
print(combined_md7[, print_cols])


# --- Full simple_slopes decomposition (raw models, for reference) ------------
cat("\n\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")
cat("  SIMPLE SLOPES — Model 6 (Detection), raw scale\n")
cat("  Full decomposition at FFA/OFA = -3, 0, 3 × Regulation\n")
cat("  (For reference only — Holm-corrected emtrends above are authoritative)\n")
cat("─────────────────────────────────────────────────────────────────────────────\n\n")
print(simple_slopes(md6, levels = list(
  FFA        = c(-3, 0, 3, "sstest"),
  OFA        = c(-3, 0, 3, "sstest"),
  Regulation = c("FFA", "OFA", "sstest")
), confint = FALSE, ci.width = 0.95))

cat("\n\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")
cat("  SIMPLE SLOPES — Model 7 (Recognition), raw scale\n")
cat("  Full decomposition at FFA/OFA = -3, 0, 3 × Regulation\n")
cat("  (For reference only — Holm-corrected emtrends above are authoritative)\n")
cat("─────────────────────────────────────────────────────────────────────────────\n\n")
print(simple_slopes(md7, levels = list(
  FFA        = c(-3, 0, 3, "sstest"),
  OFA        = c(-3, 0, 3, "sstest"),
  Regulation = c("FFA", "OFA", "sstest")
), confint = FALSE, ci.width = 0.95))


cat("\n\n=== END OF LOG ===\n")
sink()

cat("Output saved to:", normalizePath(output_file), "\n")
