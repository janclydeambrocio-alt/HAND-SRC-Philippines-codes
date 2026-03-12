# ───────────────────────────────────────────────────────────────────
#  CABULA — ALL-PASS SRC ENVELOPE PLOT
#  Appendix J.4: All-Pass SRC Generation Code
#  Station: R10.004 (Cagayan de Oro River)
#
#  This script scans all simulated Synthetic Rating Curves (SRCs)
#  from a sensitivity analysis output and retains only those that
#  pass ALL five goodness-of-fit criteria simultaneously:
#
#    - NRMSE  < 0.60
#    - nMAE   < 0.60
#    - |PBIAS|< 20%
#    - NSE    > 0.50
#    - RSR   <= 0.70
#
#  The resulting "All-PASS" envelope is plotted against the ORC:
#      Q_ORC = 17.779 * (h - 1.025)^3.43
#
#  Upper and lower bounding curves are identified and their
#  fitted equations and metrics are displayed on the plot.
#
#  Author: [Your Name]
#  Date:   [Date]
# ───────────────────────────────────────────────────────────────────


# ───────────────────────────────────────────────────────────────────
#  1. LOAD DATA
# ───────────────────────────────────────────────────────────────────

cat("Select CABULA CSV file...\n")
filename <- file.choose()
raw_df   <- read.csv(filename, header = FALSE, stringsAsFactors = FALSE)


# ───────────────────────────────────────────────────────────────────
#  2. PARAMETERS
# ───────────────────────────────────────────────────────────────────

offset_val <- -1.025   # h0: effective zero-flow offset (m)
gage_min   <-  2.450   # Minimum gage height in valid range (m)
gage_max   <-  3.870   # Maximum gage height in valid range (m)

# ORC power-law coefficients: Q = coeff_A * H^power_B
coeff_A <- 17.779
power_B <-  3.43

# Convert gage height range to centerline depth
depth_min <- gage_min + offset_val
depth_max <- gage_max + offset_val


# ───────────────────────────────────────────────────────────────────
#  3. EXTRACT DATA
# ───────────────────────────────────────────────────────────────────

cat("Processing depths...\n")
h_model_all <- suppressWarnings(as.numeric(raw_df[1, 2:31]))

# Retain only depths within the valid stage range
valid_idx <- which(
  !is.na(h_model_all) &
  h_model_all >= depth_min &
  h_model_all <= depth_max
)

if (length(valid_idx) <= 1) stop("ERROR: No valid depths found in range.")

h_subset_depth <- h_model_all[valid_idx]

start_row <- 3
end_row   <- min(74594, nrow(raw_df))

mannings_list  <- suppressWarnings(as.numeric(raw_df[start_row:end_row, 1]))
col_indices    <- valid_idx + 1   # Offset by 1 (column 1 = Manning's n)
q_raw          <- raw_df[start_row:end_row, col_indices]
q_matrix       <- as.matrix(q_raw)
mode(q_matrix) <- "numeric"

# Compute ORC target values at each valid depth
target_q <- coeff_A * (h_subset_depth)^power_B
mean_obs  <- mean(target_q, na.rm = TRUE)


# ───────────────────────────────────────────────────────────────────
#  4. SCAN ALL ROWS — COLLECT ALL-PASS CURVES + THEIR METRICS
# ───────────────────────────────────────────────────────────────────

cat("Scanning", nrow(q_matrix), "rows for all-PASS curves...\n")

all_pass_curves  <- list()
all_pass_n       <- c()
all_pass_nse     <- c()
all_pass_metrics <- list()

h_seq <- seq(min(h_subset_depth), max(h_subset_depth), length.out = 200)

for (i in 1:nrow(q_matrix)) {
  sim_q <- q_matrix[i, ]

  # Skip rows with too many missing values
  if (sum(!is.na(sim_q)) <= length(sim_q) / 2) next

  # Compute goodness-of-fit metrics
  rmse_v  <- sqrt(mean((sim_q - target_q)^2, na.rm = TRUE))
  mae_v   <- mean(abs(target_q - sim_q), na.rm = TRUE)
  nrmse_v <- rmse_v / mean_obs
  nmae_v  <- mae_v  / mean_obs
  pbias_v <- 100 * sum(target_q - sim_q, na.rm = TRUE) / sum(target_q, na.rm = TRUE)
  nse_v   <- 1 - sum((target_q - sim_q)^2, na.rm = TRUE) /
                     sum((target_q - mean_obs)^2, na.rm = TRUE)
  rsr_v   <- rmse_v / sd(target_q, na.rm = TRUE)

  # All five criteria must PASS
  all_pass <- (nrmse_v  < 0.60) &
              (nmae_v   < 0.60) &
              (abs(pbias_v) < 20) &
              (nse_v    > 0.50) &
              (rsr_v   <= 0.70)

  if (!all_pass) next

  # Fit power-law to this passing curve
  vl <- which(sim_q > 0 & h_subset_depth > 0)
  if (length(vl) < 3) next

  lm_s <- lm(log(sim_q[vl]) ~ log(h_subset_depth[vl]))
  a_s  <- exp(coef(lm_s)[1])
  b_s  <- coef(lm_s)[2]

  fit_s <- tryCatch({
    nls(
      sim_q[vl] ~ a * h_subset_depth[vl]^b,
      start   = list(a = a_s, b = b_s),
      control = nls.control(maxiter = 200)
    )
  }, error = function(e) lm_s)

  if (inherits(fit_s, "nls")) {
    a_f <- coef(fit_s)[["a"]]
    b_f <- coef(fit_s)[["b"]]
  } else {
    a_f <- exp(coef(fit_s)[1])
    b_f <- coef(fit_s)[2]
  }

  src_q <- a_f * h_seq^b_f

  all_pass_curves[[length(all_pass_curves) + 1]] <- src_q
  all_pass_n   <- c(all_pass_n,   mannings_list[i])
  all_pass_nse <- c(all_pass_nse, nse_v)
  all_pass_metrics[[length(all_pass_metrics) + 1]] <- list(
    nse   = nse_v,
    rsr   = rsr_v,
    pbias = pbias_v,
    nrmse = nrmse_v,
    nmae  = nmae_v,
    a     = a_f,
    b     = b_f
  )
}

cat(sprintf("\nTotal all-PASS curves found: %d\n", length(all_pass_curves)))
if (length(all_pass_curves) == 0) stop("No all-PASS curves found.")


# ───────────────────────────────────────────────────────────────────
#  5. IDENTIFY UPPER AND LOWER BOUND BY MEAN Q
# ───────────────────────────────────────────────────────────────────

# Upper = highest Q at the maximum depth; Lower = lowest
q_at_max_h <- sapply(all_pass_curves, function(cv) cv[length(cv)])
upper_idx  <- which.max(q_at_max_h)
lower_idx  <- which.min(q_at_max_h)

m_upper <- all_pass_metrics[[upper_idx]]
m_lower <- all_pass_metrics[[lower_idx]]


# ───────────────────────────────────────────────────────────────────
#  6. PLOT
# ───────────────────────────────────────────────────────────────────

par(
  family = "serif",
  mar    = c(3.5, 4.5, 3.5, 1.0),
  mgp    = c(2.5, 0.7, 0)
)

all_q_vals <- c(target_q, unlist(all_pass_curves), coeff_A * h_seq^power_B)
max_y      <- max(all_q_vals, na.rm = TRUE) * 1.15

# Initialize empty plot canvas
plot(
  h_seq, h_seq * 0,
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

# Main title and subtitle
mtext(
  "Station: R10.004 (Cagayan de Oro River - Cabula)",
  side = 3, line = 2.0, cex = 1.15, font = 2, adj = 0
)

ex_orc <- expression(Q[ORC] == 17.779*(h - 1.025)^3.43)
mtext(ex_orc, side = 3, line = 0.5, cex = 0.9, col = "blue", adj = 0)

# Draw all-PASS SRC curves (semi-transparent steel blue)
for (k in seq_along(all_pass_curves)) {
  lines(
    h_seq, all_pass_curves[[k]],
    col = adjustcolor("steelblue", alpha.f = 0.25),
    lwd = 1.2, lty = 1
  )
}

# Overlay ORC on top
orc_q <- coeff_A * h_seq^power_B
lines(h_seq, orc_q, col = "blue", lwd = 2.2, lty = 4)

# Main legend
legend(
  "topleft",
  legend = c("ORC", sprintf("All-PASS SRCs (N = %d)", length(all_pass_curves))),
  col    = c("blue", adjustcolor("steelblue", alpha.f = 0.5)),
  lty    = c(4, 1),
  lwd    = c(2.2, 1.5),
  bty    = "o", bg = "white", box.col = "black",
  cex    = 0.90, inset = c(0.01, 0.01),
  x.intersp  = 0.8,
  text.width = strwidth("All-PASS SRCs (N = 99999)")
)

# Metrics panel — Upper (U) and Lower (L) bounds
metrics_text <- c(
  sprintf("%-9s   %-12s %-12s %s", "",      "Upper (U)",                  "Lower (L)",                  "Status"),
  sprintf("%-9s : %-12s %-12s %s", "NSE",   sprintf("%.3f",   m_upper$nse),   sprintf("%.3f",   m_lower$nse),   "PASS"),
  sprintf("%-9s : %-12s %-12s %s", "RSR",   sprintf("%.3f",   m_upper$rsr),   sprintf("%.3f",   m_lower$rsr),   "PASS"),
  sprintf("%-9s : %-12s %-12s %s", "PBIAS", sprintf("%.2f%%", m_upper$pbias), sprintf("%.2f%%", m_lower$pbias), "PASS"),
  sprintf("%-9s : %-12s %-12s %s", "NRMSE", sprintf("%.3f",   m_upper$nrmse), sprintf("%.3f",   m_lower$nrmse), "PASS"),
  sprintf("%-9s : %-12s %-12s %s", "nMAE",  sprintf("%.3f",   m_upper$nmae),  sprintf("%.3f",   m_lower$nmae),  "PASS"),
  "",
  sprintf("%-9s : Q = %.3f * H ^ %.3f", "SRC (U)", m_upper$a, m_upper$b),
  sprintf("%-9s : Q = %.3f * H ^ %.3f", "SRC (L)", m_lower$a, m_lower$b)
)

op <- par(family = "mono")
legend(
  "bottomright",
  legend    = metrics_text,
  bty       = "o", bg = "white", box.col = "black",
  cex       = 0.65, text.font = 1,
  inset     = c(0.01, 0.01), y.intersp = 1.00,
  x.intersp = 0.5, seg.len = 0, pch = NA, lty = 0
)
par(op)


# ───────────────────────────────────────────────────────────────────
#  7. CONSOLE OUTPUT
# ───────────────────────────────────────────────────────────────────

cat("────────────────────────────────────────────\n")
cat(sprintf(" Manning's n range (all-PASS): %.6f – %.6f\n",
            min(all_pass_n), max(all_pass_n)))
cat(sprintf(" NSE range        (all-PASS): %.4f – %.4f\n",
            min(all_pass_nse), max(all_pass_nse)))
cat(sprintf(" Upper bound n    : %.6f\n", all_pass_n[upper_idx]))
cat(sprintf(" Lower bound n    : %.6f\n", all_pass_n[lower_idx]))
cat("────────────────────────────────────────────\n")
