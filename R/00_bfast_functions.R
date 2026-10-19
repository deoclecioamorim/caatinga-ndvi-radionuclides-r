#Manuscript: helper functions (no section of their own) | used by 05_bfast_pixel_metrics.R
#Shared functions: extraction of the bimonthly NDVI and BFAST metrics per pixel ---------------------
#Same logic as 01_ndvi_extraction.R and 02_bfast_mean_series.R


# Extraction of all pixels of an area ----------------------------------------------------------
#Images with > max_na_img % of NA are discarded; pixels with > max_na_px % of NA are removed;
#remaining NAs are filled by temporal interpolation for each pixel
extract_all_pixels <- function(folder, max_na_img = 50, max_na_px = 20) {
  files <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
  period <- sub(".*_harmOLI_(\\d{4})_(\\d{2})_(\\d{2})\\.tif$", "\\1_\\2_\\3", basename(files))
  files <- files[order(period)]; period <- sort(period)
  images <- lapply(files, rast)
  ref <- images[[1]]
  images <- lapply(images, function(r) {
    if (!compareGeom(r, ref, stopOnError = FALSE)) {
      if (crs(r) != crs(ref)) r <- project(r, ref)
      r <- resample(r, ref, method = "near")
    }
    r
  })
  stack <- rast(images); names(stack) <- period
  df <- as.data.frame(stack, xy = TRUE, na.rm = FALSE)
  values <- as.matrix(df[, period])
  keep <- rowSums(!is.na(values)) > 0
  df <- df[keep, ]; values <- values[keep, , drop = FALSE]
  values[, colMeans(is.na(values)) * 100 > max_na_img] <- NA
  keep <- rowMeans(is.na(values)) * 100 <= max_na_px
  df <- df[keep, ]; values <- values[keep, , drop = FALSE]
  interpolated <- is.na(values)
  values <- t(apply(values, 1, zoo::na.approx, na.rm = FALSE, rule = 2))
  list(xy = df[, c("x", "y")], ndvi = values, interp = interpolated, period = period, crs = crs(ref))
}


# BFAST metrics of a series --------------------------------------------------------------------
#Default bfast(): harmonic seasonal model (order 3, fixed in the package)
pixel_metrics <- function(y, h = 0.15, max.iter = 1, freq = 6, interp = rep(FALSE, length(y)), breaks = FALSE) {
  n <- length(y)
  yt <- ts(y, start = c(2000, 1), frequency = freq)
  fit <- tryCatch(bfast(yt, h = h, season = "harmonic", max.iter = max.iter), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  o <- fit$output[[length(fit$output)]]
  Tt <- as.numeric(o$Tt); Nt <- as.numeric(o$Nt)
  time_yr <- as.numeric(time(yt))
  bp <- if (fit$nobp$Vt) integer(0) else o$Vt.bp
  start_seg <- c(1, bp + 1); end_seg <- c(bp, n)
  slope <- (Tt[end_seg] - Tt[start_seg]) / (time_yr[end_seg] - time_yr[start_seg])
  n_seg <- end_seg - start_seg + 1
  mag <- Tt[bp + 1] - Tt[bp]
  decline <- mag < 0
  St <- as.numeric(o$St)
  res <- data.frame(
    n_trend_breaks = length(bp),
    n_seasonal_breaks = if (fit$nobp$Wt) 0 else length(o$Wt.bp),
    slope_first_segment = slope[1], #NDVI/year
    end_first_segment = time_yr[end_seg[1]],
    slope_2000_2006 = unname(coef(lm(Tt[time_yr < 2006.25] ~ time_yr[time_yr < 2006.25]))[2]), #trend in 2000-2006 (NDVI/year)
    pct_time_negative_trend = sum(n_seg[slope < 0]) / n * 100,
    seasonal_amplitude = mean(tapply(St, floor(time_yr), function(s) max(s) - min(s))), #max - min of St in each year
    largest_decline = if (any(decline)) min(mag) else 0,
    year_largest_decline = if (any(decline)) time_yr[bp[which.min(mag)]] else NA,
    decline_2002_2006 = any(decline & time_yr[bp] >= 2002 & time_yr[bp] < 2007),
    decline_2011_2012 = any(decline & time_yr[bp] >= 2011 & time_yr[bp] < 2013),
    net_trend_change = Tt[n] - Tt[1],
    remainder_sd = sd(Nt),
    mean_ndvi = mean(y),
    pct_interpolated = mean(interp) * 100,
    n_interpolated_breaks = sum(interp[bp] | interp[bp + 1]) #trend breakpoints at interpolated composites
  )
  #breaks = TRUE: also returns each trend breakpoint (year and magnitude)
  if (!breaks) return(res)
  list(summary = res,
       breaks = if (length(bp) > 0) data.frame(year = time_yr[bp], magnitude = mag) else NULL)
}
