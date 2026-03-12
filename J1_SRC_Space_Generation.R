# ────────────────────────────────────────────────────────────────
# Synthetic Rating Curve (SRC) Calibration Plot — Cabula
# Station: R10.004 (Cagayan de Oro River)
# ────────────────────────────────────────────────────────────────
#
# This script loads a matrix of simulated discharge curves (SRCs)
# and compares them against an Observed Rating Curve (ORC) defined
# by a known power-law equation:
#
#     Q = 17.779 * (h - 1.025)^3.43
#
# The script identifies the best-fitting SRC by minimizing RMSE,
# computes goodness-of-fit statistics, and produces a calibration
# plot showing the SRC space overlaid with the ORC.
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
# These values define the valid stage range and the ORC equation
# for Station R10.004 (Cabula).

offset_val <- -1.025   # h0: effective zero-flow offset (m)
gage_min   <-  2.450   # Minimum observed gage height (m)
gage_max   <-  3.870   # Maximum observed gage height (m)

# ORC power-law coefficients: Q = coeff_A * (h + offset_val)^power_B
coeff_A <- 17.779
power_B <-  3.43

# Convert gage height range to centerline depth (h = gage + offset)
depth_min <- gage_min + offset_val
depth_max <- gage_max + offset_val


# ────────────────────────
# 3. EXTRACT AND CLEAN DATA
# ────────────────────────
# Row 1  → model depth values (columns 2 to 31)
# Row 3+ → Manning's n (column 1) and simulated Q per depth (columns 2+)

cat("Processing depth values from Row 1...\n")
h_model_all <- suppressWarnings(as.numeric(raw_df[1, 2:31]))

# Keep only depth values that fall within the valid stage range
valid_idx <- which(
  !is.na(h_model_all) &
  h_model_all >= depth_min &
  h_model_all <= depth_max
)

if (length(valid_idx) <= 1) stop("ERROR: No valid depths found within the defined range.")

h_subset_depth <- h_model_all[valid_idx]  # Depths used for comparison

# Extract Manning's n values and the corresponding Q matrix
start_row   <- 3
end_row     <- min(74594, nrow(raw_df))   # Row limit for Cabula dataset
col_indices <- valid_idx + 1              # Offset by 1 since column 1 = Manning's n

mannings_list <- suppressWarnings(as.numeric(raw_df[start_row:end_row, 1]))

q_raw    <- raw_df[start_row:end_row, col_indices]
q_matrix <- as.matrix(q_raw)
mode(q_matrix) <- "numeric"


# ────────────────────────
# 4. OBSERVED RATING CURVE (ORC)
# ────────────────────────
# Compute target discharge values using the ORC equation
# at each valid depth in the subset.

target_q <- coeff_A * (h_subset_depth)^power_B


# ────────────────────────
# 5. FIND BEST-FITTING SRC
# ────────────────────────
# For each row (i.e., each Manning's n combination), compute the RMSE
# between the simulated Q and the ORC target. The row with the lowest
# RMSE is selected as the best-fitting Synthetic Rating Curve.

cat("Scanning", nrow(q_matrix), "simulated curves for best RMSE fit...\n")
all_rmse <- numeric(nrow(q_matrix))

for (i in 1:nrow(q_matrix)) {
  sim_q <- q_matrix[i, ]

  # Only compute RMSE if at least half the values are valid
  if (sum(!is.na(sim_q)) > (length(sim_q) / 2)) {
    diff        <- sim_q - target_q
    all_rmse[i] <- sqrt(mean(diff^2, na.rm = TRUE))
  } else {
    all_rmse[i] <- 999999  # Penalize rows with too many missing values
  }
}

best_idx         <- which.min(all_rmse)
best_n           <- mannings_list[best_idx]
best_curve_q     <- as.vector(q_matrix[best_idx, ])
actual_excel_row <- best_idx + (start_row - 1)  # For traceability back to the CSV


# ────────────────────────
# 6. GOODNESS-OF-FIT STATISTICS
# ────────────────────────

rmse_val  <- all_rmse[best_idx]
mae_val   <- mean(abs(target_q - best_curve_q), na.rm = TRUE)

mean_obs  <- mean(target_q, na.rm = TRUE)
nrmse_val <- rmse_val / mean_obs   # Normalized RMSE
nmae_val  <- mae_val  / mean_obs   # Normalized MAE

# Percent Bias (PBIAS): positive = model underestimates, negative = overestimates
pbias_val <- 100 * sum(target_q - best_curve_q, na.rm = TRUE) / sum(target_q, na.rm = TRUE)

# Nash-Sutcliffe Efficiency (NSE)
numerator_nse   <- sum((target_q - best_curve_q)^2, na.rm = TRUE)
denominator_nse <- sum((target_q - mean_obs)^2,     na.rm = TRUE)
nse_val <- 1 - (numerator_nse / denominator_nse)

# RMSE-Observations Standard Deviation Ratio (RSR)
rsr_val <- rmse_val / sd(target_q, na.rm = TRUE)

# Pass/Fail evaluation against standard thresholds
check_nrmse <- ifelse(nrmse_val      < 0.60, "PASS", "FAIL")
check_nmae  <- ifelse(nmae_val       < 0.60, "PASS", "FAIL")
check_pbias <- ifelse(abs(pbias_val) < 20,   "PASS", "FAIL")
check_nse   <- ifelse(nse_val        > 0.50, "PASS", "FAIL")
check_rsr   <- ifelse(rsr_val       <= 0.70, "PASS", "FAIL")


# ────────────────────────
# 7. FIT POWER-LAW TO BEST SRC
# ────────────────────────
# Re-fit a power-law curve to the best SRC using nls().
# Falls back to a log-log linear model if nls() fails to converge.

valid_log_idx <- which(best_curve_q > 0 & h_subset_depth > 0)
h_fit <- h_subset_depth[valid_log_idx]
q_fit <- best_curve_q[valid_log_idx]

# Use log-log regression as starting values for nls
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

# Extract fitted coefficients
if (inherits(fit_model, "nls")) {
  a_fit <- coef(fit_model)[["a"]]
  b_fit <- coef(fit_model)[["b"]]
} else {
  a_fit <- exp(coef(fit_model)[1])
  b_fit <- coef(fit_model)[2]
}

a_txt <- sprintf("%.3f", a_fit)
b_txt <- sprintf("%.3f", b_fit)

# Generate smooth regression line for the best SRC
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

# Set y-axis upper limit with 15% headroom
max_y <- max(c(target_q, best_curve_q, src_reg_q, as.vector(q_matrix)), na.rm = TRUE) * 1.15

# Initialize empty plot area
plot(
  h_subset_depth, target_q,
  type = "n",
  xlab = "", ylab = "",
  ylim = c(0, max_y),
  xlim = c(min(h_subset_depth), max(h_subset_depth)),
  axes = FALSE
)

# Add grid, axes, and box
grid(col = "lightgray", lty = "dotted", lwd = 1.5)
axis(1, col = "black", tck = -0.02, cex.axis = 1.0)
axis(2, col = "black", tck = -0.02, las = 1, cex.axis = 1.0)
box()

# Axis labels
title(xlab = expression(bold("Centerline Water Depth, H (m)")), line = 2.0, cex.lab = 1.0)
title(ylab = expression(bold("Discharge, Q (m"^3*"/s)")),       line = 2.3, cex.lab = 1.0)

# Plot title and ORC equation subtitle
mtext(
  "Station: R10.004 (Cagayan de Oro River - Cabula) SRC Space",
  side = 3, line = 2.0, cex = 1.15, font = 2, adj = 0
)

ex_orc <- expression(Q[ORC] == 17.779*(h - 1.025)^3.43)
mtext(ex_orc, side = 3, line = 0.5, cex = 0.9, col = "blue", adj = 0, at = par("usr")[1])

# Draw all simulated SRC curves (light red, semi-transparent)
n_curves <- nrow(q_matrix)

for (i in 1:n_curves) {
  sim_q <- as.vector(q_matrix[i, ])
  if (sum(!is.na(sim_q)) > (length(sim_q) / 2)) {
    lines(h_subset_depth, sim_q, col = adjustcolor("red", alpha.f = 0.08), lwd = 0.5)
  }
}

# Overlay the ORC on top of the SRC space
lines(h_subset_depth, target_q, col = "blue", lwd = 2.0, lty = 4)

# Legend
n_label <- paste0("SRCs (N = ", format(n_curves, big.mark = ","), " curves)")

legend(
  "topleft",
  legend  = c("ORC", n_label),
  col     = c("blue", adjustcolor("red", alpha.f = 0.5)),
  lty     = c(4, 1),
  lwd     = c(2.0, 1.5),
  bty     = "o",
  bg      = "white",
  box.col = "black",
  cex     = 0.95,
  inset   = c(0.01, 0.01)
)
