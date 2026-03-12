# ────────────────────────────────────────────
# Rating Curve Analysis — R10.004 Cabula
# Cagayan de Oro River
# ────────────────────────────────────────────
#
# This script fits a power-law rating curve to observed
# stage-discharge data using nonlinear least squares (NLS).
# It produces six diagnostic plots and a statistical summary:
#
#   1. Observed Rating Curve
#   2. Log-Log Transformation
#   3. Reliability Assessment
#   4. Residual Plot
#   5. Normal Q-Q Plot
#   6. Statistical Summary Table
#
# Output files are saved as PDFs to ~/Downloads/THESIS PLOTS/
#
# Requirements:
#   minpack.lm, ggplot2, qqplotr, gridExtra, grid
#
# Author: [Your Name]
# Date:   [Date]
# ────────────────────────────────────────────


# ────────────────────────
# 1. LOAD LIBRARIES
# ────────────────────────

graphics.off()  # Clear any existing plots before starting

# Install packages automatically if not yet available
if (!require(minpack.lm)) install.packages("minpack.lm")
if (!require(ggplot2))    install.packages("ggplot2")
if (!require(qqplotr))    install.packages("qqplotr")
if (!require(gridExtra))  install.packages("gridExtra")
if (!require(grid))       install.packages("grid")

library(minpack.lm)
library(ggplot2)
library(qqplotr)
library(gridExtra)
library(grid)


# ────────────────────────
# 2. OUTPUT FOLDER SETUP
# ────────────────────────

# All PDFs will be saved here
output_dir <- file.path(Sys.getenv("HOME"), "Downloads", "THESIS PLOTS")
if (!dir.exists(output_dir)) dir.create(output_dir)


# ────────────────────────
# 3. LOAD DATA
# ────────────────────────

cat("Please select your CSV file for CABULA (CAGAYAN DE ORO)...\n")

raw_data   <- read.csv(file.choose(), header = FALSE, stringsAsFactors = FALSE)
data_clean <- na.omit(raw_data[2:nrow(raw_data), ])  # Remove header row and NAs

# Convert columns to numeric and assign names
data_clean[, 1:2] <- sapply(data_clean[, 1:2], as.numeric)
colnames(data_clean) <- c("Gage_Height", "Discharge")

# Assign sequential dates starting from Jan 1, 2015
data_clean$Date <- seq(
  as.Date("2015-01-01"),
  by = "day",
  length.out = nrow(data_clean)
)


# ────────────────────────
# 4. OPTIMIZATION ENGINE
# ────────────────────────
# This function searches for the longest window of consecutive
# observations where a power-law curve fits within 5% error
# for at least 95% of data points.
#
# Power-law form:  Q = a * (h - h0)^b
#   Q  = Discharge (m³/s)
#   h  = Gage height (m)
#   h0 = Effective zero-flow stage (m)
#   a  = Coefficient
#   b  = Exponent

find_best_window <- function(df) {

  best_window <- NULL
  best_params <- NULL
  max_len     <- 0
  start_h0    <- min(df$Gage_Height) - 0.5  # Initial guess for h0

  for (i in 1:(nrow(df) - 60)) {
    for (j in (i + 60):nrow(df)) {

      tmp <- df[i:j, ]

      # Attempt NLS fit; skip if it fails to converge
      fit <- tryCatch({
        nlsLM(
          Discharge ~ a * (Gage_Height - h0)^b,
          data    = tmp,
          start   = list(a = 20, b = 2, h0 = start_h0),
          control = nls.lm.control(maxiter = 50)
        )
      }, error = function(e) NULL)

      if (!is.null(fit)) {

        # Compute percent deviation for each observation
        pred            <- predict(fit)
        pred[pred == 0] <- NA
        rel <- sum(
          abs(((tmp$Discharge - pred) / pred) * 100) <= 5,
          na.rm = TRUE
        ) / nrow(tmp) * 100

        # Keep this window if reliability >= 95% and longer than previous best
        if (rel >= 95 && (j - i) > max_len) {
          max_len     <- j - i
          best_window <- tmp
          best_params <- coef(fit)
        }
      }
    }
  }

  if (is.null(best_window)) stop("No window satisfied the 95% reliability criteria.")
  return(list(data = best_window, params = best_params))
}

cat("Optimizing Rating Curve for Cabula (Wait lang, Inday...)\n")
res       <- find_best_window(data_clean)
best_data <- res$data
co        <- res$params


# ────────────────────────
# 5. FINAL FIT & VARIABLES
# ────────────────────────

# Re-fit using the best window parameters
fit_final <- nlsLM(
  Discharge ~ a * (Gage_Height - h0)^b,
  data  = best_data,
  start = as.list(co)
)

# Generate a smooth curve for plotting (200 points across stage range)
stage_seq    <- seq(min(best_data$Gage_Height), max(best_data$Gage_Height), length.out = 200)
smooth_curve <- data.frame(
  Gage_Height = stage_seq,
  Q_Pred      = predict(fit_final, newdata = data.frame(Gage_Height = stage_seq))
)

# Add predicted values, residuals, and percent deviation to the dataset
best_data$Q_Pred             <- predict(fit_final)
best_data$Residual           <- best_data$Discharge - best_data$Q_Pred
best_data$Pct_Dev            <- ((best_data$Discharge - best_data$Q_Pred) / best_data$Q_Pred) * 100
best_data$Reliability_Status <- ifelse(abs(best_data$Pct_Dev) <= 5, "Within Range", "Outside Range")


# ────────────────────────
# 6. STATISTICS
# ────────────────────────

r_sq            <- cor(best_data$Discharge, best_data$Q_Pred)^2
n_obs           <- nrow(best_data)
pass_rate_final <- round(
  sum(best_data$Reliability_Status == "Within Range", na.rm = TRUE) / n_obs * 100,
  2
)


# ────────────────────────
# 7. EQUATION LABELS
# ────────────────────────

# Round fitted parameters for display
a_val  <- round(as.numeric(unname(co[1])), 3)
b_val  <- round(as.numeric(unname(co[2])), 2)   # 2 decimal places
h0_val <- round(as.numeric(unname(co[3])), 3)

# Build (h ± h0) term — sign depends on whether h0 is negative or positive
h0_abs     <- round(abs(h0_val), 3)
stage_term <- if (h0_val < 0) {
  bquote(italic(h) + .(h0_abs))
} else {
  bquote(italic(h) - .(h0_val))
}

# Equation labels for annotations
eq_label     <- bquote(italic(Q) == .(a_val) * (.(stage_term))^.(b_val))
log_eq_label <- bquote(italic(log)(italic(Q)) == .(b_val) %*% italic(log)(.(stage_term)) + italic(log)(.(a_val)))
r2_label     <- bquote(italic(R)^2 == .(round(r_sq, 4)))


# ────────────────────────
# 8. SHARED PLOT THEME
# ────────────────────────

common_theme <- theme_bw(base_size = 9, base_family = "serif") +
  theme(
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle     = element_text(hjust = 0.5, size = 7.5),
    legend.background = element_rect(color = "black", linewidth = 0.3),
    legend.margin     = margin(2, 4, 2, 4),
    legend.key.size   = unit(0.35, "cm"),
    legend.text       = element_text(size = 7),
    legend.title      = element_text(size = 7, face = "bold")
  )


# ────────────────────────
# 9. PLOTS
# ────────────────────────

# Shared annotation positions
x_anno <- min(best_data$Gage_Height, na.rm = TRUE)
y_anno <- max(c(best_data$Discharge, smooth_curve$Q_Pred), na.rm = TRUE)

# --- Plot 1: Observed Rating Curve ---
p1 <- ggplot() +
  geom_point(
    data  = best_data,
    aes(x = Gage_Height, y = Discharge, color = "Observations"),
    shape = 21, fill = "palegreen3", size = 1.2, alpha = 0.7
  ) +
  geom_line(
    data = smooth_curve,
    aes(x = Gage_Height, y = Q_Pred, color = "Rating Curve"),
    linewidth = 0.8
  ) +
  scale_color_manual(
    name   = "Legend",
    values = c("Observations" = "black", "Rating Curve" = "blue"),
    breaks = c("Observations", "Rating Curve")
  ) +
  guides(color = guide_legend(override.aes = list(
    shape = c(21, NA), fill = c("palegreen3", NA),
    linetype = c(0, 1), linewidth = c(0, 0.8)
  ))) +
  labs(
    title    = "Observed Rating Curve",
    subtitle = paste0(
      "Station: R10.004 (Cagayan de Oro River - Cabula)\n",
      "N = ", n_obs, " | Period: ", min(best_data$Date), " to ", max(best_data$Date)
    ),
    x = expression(Stage*","~italic(h)~(m)),
    y = expression(Discharge*","~italic(Q)~(m^3/s))
  ) +
  common_theme +
  theme(legend.position = c(0.8, 0.2)) +
  annotate("text", x = x_anno, y = y_anno,
           label = as.character(as.expression(eq_label)),
           parse = TRUE, hjust = 0, vjust = 1.2, size = 3.5, color = "blue", family = "serif") +
  annotate("text", x = x_anno, y = y_anno,
           label = as.character(as.expression(r2_label)),
           parse = TRUE, hjust = 0, vjust = 2.6, size = 3.5, family = "serif")


# --- Plot 2: Log-Log Transformation ---
p2 <- ggplot(best_data, aes(x = Gage_Height - h0_val, y = Discharge)) +
  geom_point(aes(color = "Observations"), shape = 21, fill = "palegreen3", size = 1.2, alpha = 0.7) +
  geom_smooth(
    aes(color = "Linear Fit"), method = "lm", linewidth = 0.8,
    se = TRUE, fill = "lightsteelblue", alpha = 0.3
  ) +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks(sides = "bl", size = 0.3) +
  scale_color_manual(
    name   = "Legend",
    values = c("Observations" = "black", "Linear Fit" = "red")
  ) +
  guides(color = guide_legend(order = 1, override.aes = list(
    shape = c(NA, 21), fill = c(NA, "palegreen3"),
    linetype = c(1, 0), linewidth = c(0.8, 0)
  ))) +
  labs(
    title    = "Log-Log Transformation",
    subtitle = paste0(
      "Station: R10.004 (Cagayan de Oro River - Cabula)\n",
      "N = ", n_obs, " | Period: ", min(best_data$Date), " to ", max(best_data$Date)
    ),
    x = expression(italic(Log)*"("*italic(h) - italic(h)[0]*")"),
    y = expression(italic(Log)*"("*italic(Q)*")")
  ) +
  common_theme +
  theme(legend.position = c(0.8, 0.2)) +
  annotate("text",
           x = min(best_data$Gage_Height - h0_val, na.rm = TRUE),
           y = max(best_data$Discharge, na.rm = TRUE),
           label = as.character(as.expression(log_eq_label)),
           parse = TRUE, hjust = 0, vjust = 1.2, size = 3.2, color = "red", family = "serif") +
  annotate("text",
           x = min(best_data$Gage_Height - h0_val, na.rm = TRUE),
           y = max(best_data$Discharge, na.rm = TRUE),
           label = as.character(as.expression(r2_label)),
           parse = TRUE, hjust = 0, vjust = 2.8, size = 3.5, family = "serif")


# --- Plot 3: Reliability Assessment ---
p3 <- ggplot(best_data, aes(x = Gage_Height, y = Pct_Dev)) +
  geom_hline(yintercept = c(5, -5), linetype = "dashed", color = "red",   linewidth = 0.6) +
  geom_hline(yintercept = 0,        linetype = "solid",  color = "black", linewidth = 0.4) +
  geom_point(aes(fill = Reliability_Status), shape = 21, color = "black", size = 1.2, alpha = 0.7) +
  scale_fill_manual(
    name   = "Legend",
    values = c("Within Range" = "palegreen3", "Outside Range" = "red"),
    breaks = c("Within Range", "Outside Range")
  ) +
  labs(
    title    = "Observed Rating Curve Reliability Assessment",
    subtitle = paste0(
      "Station: R10.004 (Cagayan de Oro River - Cabula)\n",
      "N = ", n_obs, " | Period: ", min(best_data$Date), " to ", max(best_data$Date),
      "\nReliability: ", pass_rate_final, "% | Threshold: \u00B15%"
    ),
    x = expression(Stage*","~italic(h)~(m)),
    y = "Percent Error (%)"
  ) +
  common_theme +
  theme(legend.position = c(0.85, 0.15))


# --- Plot 4: Residual Plot ---
p4 <- ggplot(best_data, aes(x = Q_Pred, y = Residual)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_smooth(aes(fill = "95% CI Band",     color = "95% CI Band"),     method = "loess", se = TRUE,  alpha = 0.25) +
  geom_smooth(aes(color = "LOESS Smoother", fill = "LOESS Smoother"),   method = "loess", se = FALSE, linewidth = 0.8) +
  geom_point(aes(fill = "Residuals", color = "Residuals"), shape = 21, size = 1.2, alpha = 0.7) +
  scale_fill_manual(
    name   = "Legend",
    values = c("Residuals" = "palegreen3", "LOESS Smoother" = "transparent", "95% CI Band" = "gray70"),
    breaks = c("Residuals", "LOESS Smoother", "95% CI Band")
  ) +
  scale_color_manual(
    name   = "Legend",
    values = c("Residuals" = "black", "LOESS Smoother" = "blue", "95% CI Band" = "transparent"),
    breaks = c("Residuals", "LOESS Smoother", "95% CI Band")
  ) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = c(21, NA, 22), fill = c("palegreen3", NA, "gray70"),
      color = c("black", "blue", NA), linetype = c(0, 1, 0),
      linewidth = c(0, 0.8, 0), alpha = c(0.7, 1, 0.4)
    )),
    color = "none"
  ) +
  labs(
    title    = "Residual Plot",
    subtitle = paste0(
      "Station: R10.004 (Cagayan de Oro River - Cabula)\n",
      "N = ", n_obs, " | Period: ", min(best_data$Date), " to ", max(best_data$Date)
    ),
    x = expression(Predicted~italic(Q)~(m^3/s)),
    y = expression(Residuals~(m^3/s))
  ) +
  common_theme +
  theme(legend.position = c(0.82, 0.18))


# --- Plot 5: Normal Q-Q Plot ---
p5 <- ggplot(data = best_data, mapping = aes(sample = Residual)) +
  stat_qq_band(
    aes(fill = "95% CI Band", color = "95% CI Band"),
    alpha = 0.3, band_type = "pointwise"
  ) +
  stat_qq_line(aes(color = "Normal Line", fill = "Normal Line"), linewidth = 0.8) +
  stat_qq_point(aes(fill = "Residuals", color = "Residuals"), shape = 21, size = 1.2, alpha = 0.7) +
  scale_fill_manual(
    name   = "Legend",
    values = c("Residuals" = "palegreen3", "Normal Line" = "transparent", "95% CI Band" = "gray70"),
    breaks = c("Residuals", "Normal Line", "95% CI Band")
  ) +
  scale_color_manual(
    name   = "Legend",
    values = c("Residuals" = "black", "Normal Line" = "red", "95% CI Band" = "transparent"),
    breaks = c("Residuals", "Normal Line", "95% CI Band")
  ) +
  guides(
    fill  = "none",
    color = guide_legend(override.aes = list(
      shape = c(21, NA, 22), fill = c("palegreen3", NA, "gray70"),
      color = c("black", "red", NA), linetype = c(0, 1, 0),
      linewidth = c(0, 0.8, 0), alpha = c(0.7, 1, 0.4)
    ))
  ) +
  labs(
    title    = "Normal Q-Q Plot",
    subtitle = paste0(
      "Station: R10.004 (Cagayan de Oro River - Cabula)\n",
      "N = ", n_obs, " | Period: ", min(best_data$Date), " to ", max(best_data$Date)
    ),
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  common_theme +
  theme(legend.position = c(0.82, 0.18))


# ────────────────────────
# 10. STATISTICAL SUMMARY TABLE
# ────────────────────────

# Compute summary statistics
min_stage <- min(best_data$Gage_Height)
max_stage <- max(best_data$Gage_Height)
df_val    <- df.residual(fit_final)
adj_r_sq  <- 1 - (1 - r_sq) * ((n_obs - 1) / df_val)
mse_val   <- sum(residuals(fit_final)^2) / df_val
se_val    <- summary(fit_final)$sigma

# Build summary data frame
stats_df <- data.frame(
  Parameter = c(
    "Min. Stage (m)", "Max. Stage (m)",
    "R-Squared (R²)", "Adjusted R²",
    "Degrees of Freedom (DF)",
    "Mean Squared Error (MSE)",
    "Standard Error (Se)"
  ),
  Value = c(
    sprintf("%.3f", min_stage),
    sprintf("%.3f", max_stage),
    sprintf("%.4f", r_sq),
    sprintf("%.4f", adj_r_sq),
    sprintf("%d",   df_val),
    sprintf("%.3f", mse_val),
    sprintf("%.3f", se_val)
  )
)

# Render as a formatted table grob
t1 <- tableGrob(
  stats_df, rows = NULL,
  theme = ttheme_default(
    core    = list(
      bg_params = list(fill = c("grey95", "white")),
      fg_params = list(fontfamily = "serif", fontsize = 12)
    ),
    colhead = list(
      fg_params = list(fontfamily = "serif", fontface = "bold", fontsize = 13)
    )
  )
)

title_grob <- textGrob(
  paste0("Statistical Summary: Station R10.004\n(Cagayan de Oro River - Cabula)"),
  gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "serif")
)

p6_table <- arrangeGrob(title_grob, t1, heights = c(0.2, 0.8))


# ────────────────────────
# 11. DISPLAY & EXPORT
# ────────────────────────

timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
cat("\nSaving plots to:", output_dir, "\n")
cat("NOTE: Each plot will appear in RStudio AND be saved to your Downloads folder.\n\n")

# Display all plots on screen
print(p1); Sys.sleep(1)
print(p2); Sys.sleep(1)
print(p3); Sys.sleep(1)
print(p4); Sys.sleep(1)
print(p5); Sys.sleep(1)
grid.newpage(); grid.draw(p6_table)

# Save each plot as a PDF
save_plot <- function(plot_obj, prefix) {
  fname <- paste0(prefix, "_Cabula_", timestamp, ".pdf")
  pdf(file.path(output_dir, fname), width = 5.8, height = 5.8)
  if (inherits(plot_obj, "ggplot")) {
    print(plot_obj)
  } else {
    grid.draw(plot_obj)
  }
  dev.off()
  cat("Saved:", fname, "\n")
}

save_plot(p1,       "1_RatingCurve")
save_plot(p2,       "2_LogLog")
save_plot(p3,       "3_Reliability")
save_plot(p4,       "4_Residual")
save_plot(p5,       "5_NormalQQ")
save_plot(p6_table, "6_StatSummary")

cat("\nSUCCESS! All 6 plots for Cabula saved with timestamp:", timestamp, "\n")
