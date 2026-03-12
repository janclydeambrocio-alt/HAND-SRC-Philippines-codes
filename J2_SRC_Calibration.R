# ────────────────────────────────────────────────────────────────
# SRC Calibration Plot with Best-Fit Overlay — Cabula
# Station: R10.004 (Cagayan de Oro River)
# ────────────────────────────────────────────────────────────────
#
# This script identifies the best-fitting Synthetic Rating Curve
# (SRC) from a sensitivity analysis output matrix and produces
# a calibration plot comparing the SRC against the Observed
# Rating Curve (ORC):
#
#     Q_ORC = 17.779 * (h - 1.025)^3.43
#
# The best SRC is selected by minimizing RMSE against the ORC.
# A power-law is then re-fitted to the best SRC data points,
# and goodness-of-fit statistics are computed and displayed.
#
# Output: Calibration plot with metrics panel (on screen)
#         Console summary of all statistics
#
# Author: [Your Name]
# Date:   [Date]
# ────────────────────────────────────────────────────────────────


# ────────────────────────
# 1. LOAD DATA
# ────────────────────────

cat("Select CABULA CSV file (Columns A to AE)...\n")
filename <- file.choose()
raw_df   <- read.csv(filename, header = FALSE, stringsAsFactors = FALSE)


# ────────────────────────
# 2. STATION PARAMETERS
# ────────────────────────
# Stage range and ORC equation for Station R10.004 (Cabula).
# Centerline depth H = Gage Height + offset_val

offset_val <- -1.025   # h0: zero-flow offset (m)
gage_min   <-  2.450   # Minimum gage height in valid range (m)
gage_max   <-  3.870   # Maximum gage height in valid range (m)

# ORC power-law: Q = coeff_A * H^power_B
coeff_A <- 17.779
power_B <-  3.43

depth_min <- gage_min + offset_val   # Minimum centerline depth
depth_max <- gage_max + offset_val   # Maximum centerline depth


# ────────────────────────
# 3. EXTRACT AND CLEAN DATA
# ────────────────────────
# CSV structure:
#   Row 1         → model depth values (columns 2 to 31)
#   Column 1      → Manning's n values (rows 3 onward)
#   Columns 2–31  → simulated Q per depth (rows 3 onward)

cat("Processing depth values from Row 1...\n")
h_model_all <- suppressWarnings(as.numeric(raw_df[1, 2:31]))

# Retain only depths within the valid stage range
valid_idx <- which(
  !is.na(h_model_all) &
  h_model_all >= depth_min &
  h_model_all <= depth_max
)

if (length(valid_idx) <= 1) stop("ERROR: No valid depths found within the defined range.")

h_subset_depth <- h_model_all[valid_idx]

# Extract Manning's n and Q matrix for valid depth columns
start_row   <- 3
end_row     <- min(74594, nrow(raw_df))
col_indices <- valid_idx + 1   # Offset by 1 (column 1 = Manning's n)

mannings_list <- suppressWarnings(as.numeric(raw_df[start_row:end_row, 1]))

q_raw    <- raw_df[start_row:end_row, col_indices]
q_matrix <- as.matrix(q_raw)
mode(q_matrix) <- "numeric"


# ────────────────────────
# 4. OBSERVED RATING CURVE (ORC)
# ────────────────────────

target_q <- coeff_A * (h_subset_depth)^power_B


# ────────────────────────
# 5. FIND BEST-FITTING SRC
# ────────────────────────
# Compute RMSE for each simulated curve against the ORC.
# Rows with more than half their values missing are penalized.

cat("Scanning", nrow(q_matrix), "simulated curves for best RMSE fit...\n")
all_rmse <- numeric(nrow(q_matrix))

for (i in 1:nrow(q_matrix)) {
  sim_q <- q_matrix[i, ]

  if (sum(!is.na(sim_q)) > (length(sim_q) / 2)) {
    diff        <- sim_q - target_q
    all_rmse[i] <- sqrt(mean(diff^2, na.rm = TRUE))
  } else {
    all_rmse[i] <- 999999   # Penalize incomplete rows
  }
}

best_idx         <- which.min(all_rmse)
best_n           <- mannings_list[best_idx]
best_curve_q     <- as.vector(q_matrix[best_idx, ])
actual_excel_row <- best_idx + (start_row - 1)   # For traceability to original CSV


# ────────────────────────
# 6. GOODNESS-OF-FIT STATISTICS
# ────────────────────────

rmse_val  <- all_rmse[best_idx]
mae_val   <- mean(abs(target_q - best_curve_q), na.rm = TRUE)
mean_obs  <- mean(target_q, na.rm = TRUE)

nrmse_val <- rmse_val / mean_obs     # Normalized RMSE
nmae_val  <- mae_val  / mean_obs     # Normalized MAE

# Percent Bias: positive = underestimate, negative = overestimate
pbias_val <- 100 * sum(target_q - best_curve_q, na.rm = TRUE) / sum(target_q, na.rm = TRUE)

# Nash-Sutcliffe Efficiency (NSE)
numerator_nse   <- sum((target_q - best_curve_q)^2, na.rm = TRUE)
denominator_nse <- sum((target_q - mean_obs)^2,     na.rm = TRUE)
nse_val <- 1 - (numerator_nse / denominator_nse)

# RMSE-Observations Standard Deviation Ratio (RSR)
rsr_val <- rmse_val / sd(target_q, na.rm = TRUE)

# Pass/Fail against standard performance thresholds
check_nrmse <- ifelse(nrmse_val      < 0.60, "PASS", "FAIL")
check_nmae  <- ifelse(nmae_val       < 0.60, "PASS", "FAIL")
check_pbias <- ifelse(abs(pbias_val) < 20,   "PASS", "FAIL")
check_nse   <- ifelse(nse_val        > 0.50, "PASS", "FAIL")
check_rsr   <- ifelse(rsr_val       <= 0.70, "PASS", "FAIL")


# ────────────────────────
# 7. FIT POWER-LAW TO BEST SRC
# ────────────────────────
# Re-fit Q = a * H^b to the best SRC using nls().
# Falls back to log-log linear regression if nls() fails.

valid_log_idx <- which(best_curve_q > 0 & h_subset_depth > 0)
h_fit <- h_subset_depth[valid_log_idx]
q_fit <- best_curve_q[valid_log_idx]

lm_start <- lm(log(q_fit) ~ log(h_fit))
a_start  <- exp(coef(lm_start)[1])
b_start  <- coef(lm_start)[2]

fit_model <- tryCatch({
  nls(
    q_fit ~ a * h_fit^b,
    start   = list(a = a_start, b = b_start),
    control = nls.control(maxiter = 200)
  )
}, error = function(e) {
  message("nls() failed — falling back to log-log lm: ", e$message)
  lm_start
})

if (inherits(fit_model, "nls")) {
  a_fit <- coef(fit_model)[["a"]]
  b_fit <- coef(fit_model)[["b"]]
} else {
  a_fit <- exp(coef(fit_model)[1])
  b_fit <- coef(fit_model)[2]
}

a_txt <- sprintf("%.3f", a_fit)
b_txt <- sprintf("%.3f", b_fit)

# Smooth curve for the fitted SRC
h_seq     <- seq(min(h_subset_depth), max(h_subset_depth), length.out = 200)
src_reg_q <- a_fit * h_seq^b_fit


# ────────────────────────
# 8. CALIBRATION PLOT
# ────────────────────────

par(
  family = "serif",
  mar    = c(3.5, 4.5, 3.5, 1.0),
  mgp    = c(2.5, 0.7, 0)
)

max_y <- max(c(target_q, best_curve_q, src_reg_q), na.rm = TRUE) * 1.15

# Initialize empty plot canvas
plot(
  h_subset_depth, target_q,
  type = "n",
  xlab = "", ylab = "",
  ylim = c(0, max_y),
  xlim = c(min(h_subset_depth), max(h_subset_depth)),
  axes = FALSE
)

# Grid, axes, and border
grid(col = "lightgray", lty = "dotted", lwd = 1.5)
axis(1, col = "black", tck = -0.02, cex.axis = 1.0)
axis(2, col = "black", tck = -0.02, las = 1, cex.axis = 1.0)
box()

# Axis labels
title(xlab = expression(bold("Centerline Water Depth, H (m)")), line = 2.0, cex.lab = 1.0)
title(ylab = expression(bold("Discharge, Q (m"^3*"/s)")),       line = 2.3, cex.lab = 1.0)

# Main title
mtext(
  "Station: R10.004 (Cagayan de Oro River - Cabula) Calibration",
  side = 3, line = 2.0, cex = 1.15, font = 2, adj = 0
)

# Subtitle: ORC equation | SRC fitted equation
ex_orc  <- expression(Q[ORC] == 17.779*(h - 1.025)^3.43)
ex_pipe <- "|"
ex_src  <- bquote(Q[SRC] == .(a_txt)*(H)^.(b_txt))

x_start <- par("usr")[1]
w_orc   <- strwidth(ex_orc,  units = "user", cex = 0.95)
w_gap   <- (par("usr")[2] - par("usr")[1]) * 0.03

mtext(ex_orc,  side = 3, line = 0.5, cex = 0.9, col = "blue",    adj = 0, at = x_start)
mtext(ex_pipe, side = 3, line = 0.5, cex = 0.9, col = "black",   adj = 0, at = x_start + w_orc + w_gap)
mtext(ex_src,  side = 3, line = 0.5, cex = 0.9, col = "darkred", adj = 0,
      at = x_start + w_orc + w_gap + strwidth(ex_pipe, units = "user", cex = 0.9) + w_gap)

# Plot curves and data points
lines(h_subset_depth, target_q,      col = "blue",    lwd = 2.0, lty = 4)   # ORC
points(h_subset_depth, best_curve_q, col = "salmon",  pch = 3,   cex = 0.9) # Best SRC points
lines(h_seq, src_reg_q,              col = "darkred", lwd = 2.0, lty = 4)   # Fitted SRC

# Main legend
legend(
  "topleft",
  legend  = c("ORC", "Best Simulated Data Points", "Fitted SRC"),
  col     = c("blue", "salmon", "darkred"),
  lty     = c(4, NA, 4),
  pch     = c(NA, 3, NA),
  lwd     = c(2.0, NA, 2.0),
  bty     = "o",
  bg      = "white",
  box.col = "black",
  cex     = 0.95,
  inset   = c(0.01, 0.01)
)

# Metrics panel (bottom-right, monospace font for alignment)
metrics_text <- c(
  sprintf("%-9s : %8.3f  (%s)", "NSE",   nse_val,                       check_nse),
  sprintf("%-9s : %8.3f  (%s)", "RSR",   rsr_val,                       check_rsr),
  sprintf("%-9s : %8s  (%s)",   "PBIAS", sprintf("%.2f%%", pbias_val),  check_pbias),
  sprintf("%-9s : %8.3f  (%s)", "NRMSE", nrmse_val,                     check_nrmse),
  sprintf("%-9s : %8.3f  (%s)", "nMAE",  nmae_val,                      check_nmae)
)

op <- par(family = "mono")
legend(
  "bottomright",
  legend    = metrics_text,
  bty       = "o",
  bg        = "white",
  box.col   = "black",
  cex       = 0.70,
  text.font = 1,
  inset     = c(0.01, 0.01),
  y.intersp = 1.00,
  x.intersp = 0.5,
  seg.len   = 0,
  pch       = NA,
  lty       = 0
)
par(op)


# ────────────────────────
# 9. CONSOLE SUMMARY
# ────────────────────────

actual_min_h <- min(h_subset_depth, na.rm = TRUE)
actual_max_h <- max(h_subset_depth, na.rm = TRUE)

cat("\n--------------------------------------------\n")
cat(" CABULA CALIBRATION RESULTS (Row Limit: 74594)\n")
cat("--------------------------------------------\n")
cat("Best Manning's n    :", sprintf("%.6f", best_n),        "\n")
cat("Best Index (R)      :", best_idx,                       "\n")
cat("Excel Row Number    :", actual_excel_row,                "\n")
cat("--------------------------------------------\n")
cat("Depth Range (H)     :", sprintf("%.3f", actual_min_h),
    "-", sprintf("%.3f", actual_max_h), "m\n")
cat("--------------------------------------------\n")
cat("NSE                 :", nse_val,   " (", check_nse,    ")\n")
cat("RSR                 :", rsr_val,   " (", check_rsr,    ")\n")
cat("PBIAS               :", pbias_val, "% (", check_pbias, ")\n")
cat("NRMSE               :", nrmse_val, " (", check_nrmse,  ")\n")
cat("nMAE                :", nmae_val,  " (", check_nmae,   ")\n")
cat("--------------------------------------------\n")
cat("Raw RMSE            :", sprintf("%.4f", rmse_val), "m3/s\n")
cat("Raw MAE             :", sprintf("%.4f", mae_val),  "m3/s\n")
cat("--------------------------------------------\n")
cat("SRC Fit Method      :", ifelse(inherits(fit_model, "nls"), "nls (direct)", "lm log-log (fallback)"), "\n")
cat("SRC Equation        : Q =", a_txt, "* H ^", b_txt, "\n")
cat("--------------------------------------------\n")
