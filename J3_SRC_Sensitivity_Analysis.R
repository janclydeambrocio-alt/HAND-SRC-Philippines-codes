# ────────────────────────────────────────────────────────────────
# Spearman Sensitivity Analysis — Manning's n vs Discharge
# Station: R10.004 (Cagayan de Oro River - Cabula)
# ────────────────────────────────────────────────────────────────
#
# This script computes Spearman rank correlations (ρ) between
# each land cover's Manning's n value and the resulting discharge
# (Q) across a range of water depths.
#
# For each depth-specific sensitivity Excel file within the
# defined range (1.5 m to 2.8 m), the script:
#   1. Reads the Manning's n columns and Discharge_Q column
#   2. Computes Spearman ρ for each land cover class
#   3. Ranks classes by mean |ρ| (overall influence strength)
#   4. Produces a ρ vs. Depth plot with a summary panel
#
# Output: PNG and PDF saved to the same folder as input files
#
# Requirements: readxl, dplyr, stringr
#
# Author: [Your Name]
# Date:   [Date]
# ────────────────────────────────────────────────────────────────


# ────────────────────────
# 1. LOAD LIBRARIES
# ────────────────────────

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
})


# ────────────────────────
# 2. CONFIGURATION
# ────────────────────────

# Folder containing all Sensitivity_*.xlsx files
# Note: File paths shown are representative.
# The full production script includes complete path strings
# specific to the study area directory structure.
folder <- "/Users/janclydeambrocio/Desktop/THESISNEWESTERA/ERRORCHECKFEB24/CABULA/Excel Files (54)"

# Only files matching this naming pattern will be processed
file_pattern <- "^Sensitivity_[0-9]+(\\.[0-9]+)?\\.xlsx$"

# Depth range to include in the analysis (meters)
DEPTH_MIN <- 1.5
DEPTH_MAX <- 2.8

# Column name for discharge in the Excel files
Q_COL    <- "Discharge_Q"
AREA_COL <- "Area"   # Loaded but not used in correlation; kept for reference

# Mapping from Excel column names to readable land cover labels
id_map <- c(
  "n_ID_1"  = "Water",
  "n_ID_11" = "Rangelands",
  "n_ID_7"  = "Built Areas",
  "n_ID_5"  = "Crops",
  "n_ID_2"  = "Trees"
)

# Plot labels
station_heading <- "Station: R10.004 (Cagayan de Oro River - Cabula)"
sub_heading     <- "Sensitivity Analysis of the LAD Composite Manning's n Formulation"

# Line and point colors per land cover class
class_cols <- c(
  "Trees"       = "darkolivegreen4",
  "Water"       = "blue",
  "Crops"       = "darkgreen",
  "Built Areas" = "gray30",
  "Rangelands"  = "sienna3"
)

# Point symbols per land cover class
class_pch <- c(
  "Trees"       = 18,
  "Water"       = 16,
  "Crops"       = 3,
  "Built Areas" = 15,
  "Rangelands"  = 17
)

# Output file paths (saved alongside input files)
png_path <- file.path(folder, "Sensitivity_Spearman_rho_vs_Depth_Passed_1p5_to_2p8.png")
pdf_path <- file.path(folder, "Sensitivity_Spearman_rho_vs_Depth_Passed_1p5_to_2p8.pdf")


# ────────────────────────
# 3. HELPER FUNCTIONS
# ────────────────────────

# Extract the numeric depth value from a filename like "Sensitivity_1.5.xlsx"
extract_depth <- function(filename) {
  m <- str_match(basename(filename), "Sensitivity_([0-9]+(?:\\.[0-9]+)?)\\.xlsx")
  as.numeric(m[, 2])
}

# Compute Spearman correlation safely; returns NA if fewer than 3 valid pairs
safe_spearman <- function(x, y) {
  ok <- complete.cases(x, y)
  x  <- x[ok]
  y  <- y[ok]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

# Format a correlation value: 3 decimal places, scientific notation if near zero
fmt_r <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(
      abs(x) < 0.001,
      formatC(x, format = "e", digits = 2),
      formatC(x, format = "f", digits = 3)
    )
  )
}

# Read an Excel file and rename Manning's n columns to land cover labels
read_and_rename <- function(fp) {
  df <- read_excel(fp)
  for (old in names(id_map)) {
    if (old %in% names(df)) names(df)[names(df) == old] <- id_map[[old]]
  }
  df
}


# ────────────────────────
# 4. LOAD AND FILTER FILES
# ────────────────────────

files <- list.files(folder, pattern = file_pattern, full.names = TRUE)
if (length(files) == 0) stop("No Sensitivity_*.xlsx files found in the specified folder.")

# Build a table of files with their extracted depth values
file_tbl <- data.frame(
  file  = files,
  Depth = sapply(files, extract_depth),
  stringsAsFactors = FALSE
)
file_tbl <- file_tbl[order(file_tbl$Depth), ]

# Keep only files within the defined depth range
passed_tbl <- subset(file_tbl, !is.na(Depth) & Depth >= DEPTH_MIN & Depth <= DEPTH_MAX)
if (nrow(passed_tbl) == 0) stop("No files found within the specified depth range (1.5–2.8 m).")

depths  <- sort(unique(passed_tbl$Depth))
classes <- unname(id_map)


# ────────────────────────
# 5. COMPUTE SPEARMAN ρ PER DEPTH
# ────────────────────────
# For each file (one per depth), compute the Spearman correlation
# between each land cover's Manning's n and the discharge column.

cor_list <- vector("list", nrow(passed_tbl))

for (i in seq_len(nrow(passed_tbl))) {

  fp <- passed_tbl$file[i]
  d  <- passed_tbl$Depth[i]

  df <- read_and_rename(fp)

  # Check that all required columns are present
  missing_cols <- setdiff(c(classes, Q_COL), names(df))
  if (length(missing_cols) > 0) {
    warning(paste0("Missing columns in Depth = ", d, ": ", paste(missing_cols, collapse = ", ")))
    cor_list[[i]] <- NULL
    next
  }

  # Ensure numeric types for all relevant columns
  for (cc in intersect(c(classes, Q_COL, AREA_COL), names(df))) {
    df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
  }

  # Compute Spearman ρ for each land cover class at this depth
  tmp <- data.frame(
    Depth      = rep(d, length(classes)),
    Class      = classes,
    Spearman_r = NA_real_,
    stringsAsFactors = FALSE
  )

  for (j in seq_along(classes)) {
    cls               <- classes[j]
    tmp$Spearman_r[j] <- safe_spearman(df[[cls]], df[[Q_COL]])
  }

  cor_list[[i]] <- tmp
}

# Combine all depth results into one data frame
cor_by_depth <- dplyr::bind_rows(cor_list)
if (nrow(cor_by_depth) == 0) stop("Spearman table is empty. Check column names, numeric conversion, and data pairs.")

classes_present <- intersect(names(class_cols), unique(cor_by_depth$Class))


# ────────────────────────
# 6. RANK BY MEAN |ρ|
# ────────────────────────
# Average absolute Spearman ρ across all depths for each class.
# Higher mean |ρ| = stronger overall influence on discharge.

rank_tbl <- cor_by_depth %>%
  group_by(Class) %>%
  summarise(
    mean_abs_rho = mean(abs(Spearman_r), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Class %in% classes_present) %>%
  arrange(desc(mean_abs_rho)) %>%
  mutate(Rank = row_number())


# ────────────────────────
# 7. PLOT FUNCTION
# ────────────────────────

plot_rho_vs_depth <- function() {

  par(
    family = "serif",
    mar    = c(3.5, 4.5, 3.8, 1.0),
    mgp    = c(2.5, 0.7, 0)
  )

  # Initialize empty plot canvas
  plot(
    depths, rep(0, length(depths)),
    type = "n",
    xlab = "", ylab = "",
    ylim = c(-1, 1),
    xlim = c(min(depths), max(depths)),
    axes = FALSE
  )

  grid(col = "lightgray", lty = "dotted", lwd = 1.5)
  axis(1, at = depths, col = "black", tck = -0.02, cex.axis = 1.0)
  axis(2, col = "black", tck = -0.02, las = 1, cex.axis = 1.0)
  box()

  title(xlab = expression(bold("Water Depth, H (m)")),           line = 2.0, cex.lab = 1.0)
  title(ylab = expression(bold("Spearman correlation, " * rho)), line = 2.3, cex.lab = 1.0)

  mtext(station_heading, side = 3, line = 2.2, cex = 1.15, font = 2, adj = 0)
  mtext(sub_heading,     side = 3, line = 1.2, cex = 1.10,           adj = 0)

  # Reference line at ρ = 0
  abline(h = 0, col = "black", lwd = 1)

  # Draw lines, points, and ρ value labels for each land cover class
  for (cls in classes_present) {

    sub <- cor_by_depth[cor_by_depth$Class == cls, c("Depth", "Spearman_r")]
    sub <- sub[order(sub$Depth), ]
    x   <- sub$Depth
    y   <- sub$Spearman_r

    lines(x, y,  col = class_cols[cls], lwd = 2.0)
    points(x, y, col = class_cols[cls], pch = class_pch[cls], cex = 0.9)

    # Annotate each point with its ρ value
    text(x, y, labels = fmt_r(y), pos = 3, cex = 0.75,
         col = class_cols[cls], xpd = NA)
  }

  # ── SUMMARY PANEL (top-right) ──────────────────────────────────
  # Shows each class, its mean |ρ|, and its influence rank

  usr <- par("usr")

  px_right <- usr[2] - (usr[2] - usr[1]) * 0.02
  py_top   <- usr[4] - (usr[4] - usr[3]) * 0.03

  box_w  <- (usr[2] - usr[1]) * 0.30
  line_h <- (usr[4] - usr[3]) * 0.040
  box_h  <- line_h * (nrow(rank_tbl) + 2.0)

  px_left <- px_right - box_w
  py_bot  <- py_top - box_h

  rect(px_left, py_bot, px_right, py_top,
       col = "white", border = "black", lwd = 1.2)

  # Column x positions within the panel
  pad   <- box_w * 0.05
  x_cls <- px_left  + pad
  x_rho <- px_left  + box_w * 0.62
  x_rnk <- px_right - pad * 1.5

  # Column headers
  y_hdr <- py_top - line_h * 0.8
  text(x_cls, y_hdr, "Class",                            adj = c(0,   0.5), cex = 0.82, font = 2)
  text(x_rho, y_hdr, expression(bold(Mean~"|"*rho*"|")), adj = c(0.5, 0.5), cex = 0.80)
  text(x_rnk, y_hdr, "Rank",                             adj = c(0.8, 0.4), cex = 0.82, font = 2)

  # Separator line below headers
  seg_y <- y_hdr - line_h * 0.55
  segments(px_left + pad * 0.6, seg_y, px_right - pad * 0.6, seg_y,
           col = "gray60", lwd = 0.8)

  # Data rows — one per land cover class, colored to match plot lines
  for (i in seq_len(nrow(rank_tbl))) {
    cls <- rank_tbl$Class[i]
    yy  <- seg_y - line_h * (i - 0.2)

    text(x_cls, yy, cls,
         adj = c(0, 0.5), cex = 0.80, col = class_cols[cls])
    text(x_rho, yy, formatC(rank_tbl$mean_abs_rho[i], format = "f", digits = 3),
         adj = c(0.5, 0.5), cex = 0.80, col = class_cols[cls])
    text(x_rnk, yy, rank_tbl$Rank[i],
         adj = c(1.8, 0.5), cex = 0.80, col = class_cols[cls])
  }
}


# ────────────────────────
# 8. DISPLAY AND SAVE
# ────────────────────────

# Show on screen first
plot_rho_vs_depth()

# Save as high-resolution PNG
png(filename = png_path, width = 3300, height = 1800, res = 300)
plot_rho_vs_depth()
dev.off()

# Save as PDF
pdf(file = pdf_path, width = 11, height = 6.5, family = "serif")
plot_rho_vs_depth()
dev.off()

cat("\nSpearman plot saved:\n")
cat(" PNG:", png_path, "\n")
cat(" PDF:", pdf_path, "\n")
