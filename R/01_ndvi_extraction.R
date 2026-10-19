#Manuscript: Section 2.2.2, paragraphs 1-2 (bimonthly composites and handling of missing values) | run 1st
#Extraction of the bimonthly NDVI series (2000-2022) - Areas A and B ------------------------------
#Input: data/raw/ndvi_area_A and data/raw/ndvi_area_B (138 images each = 23 years x 6 bimonths),
#exported from Google Earth Engine. No separation between dry and rainy seasons: one series per area.
#All scripts are run from the project root; paths are relative to it.

#Handling of missing values (Landsat 7 SLC-off gaps, clouds and the 2012 gap):
#1. Images with more than 'max_na_img' % of NA within the polygon are discarded entirely (set to NA),
#   so that observed and interpolated pixels are not mixed in the same bimonth
#2. Pixels with more than 'max_na_px' % of NA in the series (after step 1) are removed
#3. Remaining NAs are filled by temporal interpolation for each pixel (bfast does not accept NA)
#4. A mask records which values were interpolated
#All pixels are exported; the mean series of each area (scripts 02-04) is computed from all of them.

#Control options --------------------------------------------------------------------------------
options(OutDec = ".")
rm(list = ls(all = TRUE))


# Packages -------------------------------------------------------------------------------------
for (p in c("terra", "zoo", "writexl")) {
  if (!require(p, character.only = TRUE)) { install.packages(p, dependencies = TRUE); library(p, character.only = TRUE) }
}

max_na_img <- 50 #maximum % of NA to keep an image
max_na_px <- 20 #maximum % of NA to keep a pixel


# Extraction function --------------------------------------------------------------------------
extract_ndvi <- function(folder, max_na_img, max_na_px) {

  #Images in chronological order (file names end in YYYY_MM_MM.tif)
  files <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
  period <- sub(".*_harmOLI_(\\d{4})_(\\d{2})_(\\d{2})\\.tif$", "\\1_\\2_\\3", basename(files))
  files <- files[order(period)]
  period <- sort(period)

  #Check the series: 23 years x 6 bimonths
  cat("\n=====", folder, ":", length(files), "images |", period[1], "to", period[length(period)], "\n")

  #Load the images and align them to the grid of the first one (extent, resolution and CRS)
  images <- lapply(files, rast)
  ref <- images[[1]]
  images <- lapply(images, function(r) {
    if (!compareGeom(r, ref, stopOnError = FALSE)) {
      if (crs(r) != crs(ref)) r <- project(r, ref)
      r <- resample(r, ref, method = "near")
    }
    r
  })

  #Stack and build the data frame (rows = pixels, columns = bimonths)
  stack <- rast(images)
  names(stack) <- period
  df <- as.data.frame(stack, xy = TRUE, na.rm = FALSE)

  #Remove pixels without any value (outside the area polygon)
  values <- as.matrix(df[, period])
  keep <- rowSums(!is.na(values)) > 0
  df <- df[keep, ]
  values <- values[keep, , drop = FALSE]
  na_original <- colMeans(is.na(values)) * 100
  cat("Pixels in the polygon:", nrow(values), "| original NA:", round(mean(is.na(values)) * 100, 1), "%\n")

  #1. Discard images with many NAs
  discarded <- na_original > max_na_img
  values[, discarded] <- NA
  cat("Discarded images (>", max_na_img, "% NA):", sum(discarded), "->",
      paste(period[discarded], collapse = ", "), "\n")

  #2. Remove pixels with many NAs
  na_px <- rowMeans(is.na(values)) * 100
  keep <- na_px <= max_na_px
  cat("Removed pixels (>", max_na_px, "% NA):", sum(!keep), "\n")
  df <- df[keep, ]
  values <- values[keep, , drop = FALSE]

  #3. Temporal interpolation for each pixel (series ends: nearest value)
  interpolated <- is.na(values) #4. mask of interpolated values
  gap <- apply(interpolated, 1, function(x) { r <- rle(x); if (any(r$values)) max(r$lengths[r$values]) else 0 })
  cat("Pixels kept:", nrow(values), "| interpolated values:", round(mean(interpolated) * 100, 1),
      "% | longest continuous gap:", max(gap), "bimonths\n")
  values <- t(apply(values, 1, zoo::na.approx, na.rm = FALSE, rule = 2))
  colnames(values) <- period
  df[, period] <- values
  mask <- data.frame(df[, c("x", "y")], interpolated, check.names = FALSE)

  #Summary per image (reported in the manuscript)
  summary_img <- data.frame(period = period,
                            na_original_pct = round(na_original, 1),
                            discarded = discarded,
                            interpolated_pct = round(colMeans(interpolated) * 100, 1))

  list(ndvi = df, interpolated = mask, summary = summary_img)
}


# Area A ---------------------------------------------------------------------------------------
ndvi_A <- extract_ndvi("data/raw/ndvi_area_A", max_na_img = max_na_img, max_na_px = max_na_px)
summary(as.vector(as.matrix(ndvi_A$ndvi[, -(1:2)])))

# Area B ---------------------------------------------------------------------------------------
ndvi_B <- extract_ndvi("data/raw/ndvi_area_B", max_na_img = max_na_img, max_na_px = max_na_px)
summary(as.vector(as.matrix(ndvi_B$ndvi[, -(1:2)])))


# Export NDVI data -----------------------------------------------------------------------------
#ndvi_bimonthly_area_X.xlsx -> sheets:
#  ndvi         : x, y + 138 bimonths (interpolated series, bfast input), all pixels
#  interpolated : x, y + 138 bimonths (TRUE = interpolated value)
#  summary      : original % NA per image, discarded images and % of interpolated values
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
writexl::write_xlsx(ndvi_A, "data/processed/ndvi_bimonthly_area_A.xlsx")
writexl::write_xlsx(ndvi_B, "data/processed/ndvi_bimonthly_area_B.xlsx")
