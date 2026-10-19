#Manuscript: Section 2.2.2 items (i)-(ii); Results 3.2 (mean series and difference A - B); Fig. 5 | run after 01
#BFAST applied to the mean bimonthly NDVI series (2000-2022) - Areas A and B ----------------------
#Input: data/processed/ndvi_bimonthly_area_A.xlsx and ndvi_bimonthly_area_B.xlsx
#       (created by 01_ndvi_extraction.R; sheet 'ndvi', all pixels)
#
#Steps:
#1. BFAST on the mean series of each area (mean of all pixels; Fig. 5)
#2. BFAST on the difference NDVI(A) - NDVI(B): it removes the regional climatic signal shared by both
#   areas, so breakpoints in this series indicate local divergence
#3. Linear trends for 2000-2006 fitted to the deseasonalised (STL) series
#Pixel-level metrics are computed in 05_bfast_pixel_metrics.R and compared in 06_bfast_spatial_resampling.R

#Control options --------------------------------------------------------------------------------
options(OutDec = ".")
rm(list = ls(all = TRUE))


# Packages -------------------------------------------------------------------------------------
for (p in c("bfast", "strucchange", "readxl", "writexl")) {
  if (!require(p, character.only = TRUE)) { install.packages(p, dependencies = TRUE); library(p, character.only = TRUE) }
}
source("R/00_publication_theme.R") #figure standard (font and print size)


# BFAST parameters (final) ---------------------------------------------------------------------
#Sensitivity tested in 03_bfast_sensitivity_h.R and 04_bfast_sensitivity_k.R.
#The default bfast() is used, which fixes the harmonic order at K = 3; with 6 obs/year the sine term
#of the 3rd harmonic is zero. With max.iter = 1, trend breakpoints are searched in Yt - St(STL).
freq <- 6 #6 bimonths per year
h <- 0.15 #minimum segment size (fraction of the series): 0.15 x 138 = ~21 bimonths (~3.5 years)
season_model <- "harmonic" #harmonic seasonal model (K = 3, fixed in bfast)
max_iter <- 1


# Data -----------------------------------------------------------------------------------------
read_area <- function(area) {
  ndvi <- read_excel(paste0("data/processed/ndvi_bimonthly_area_", area, ".xlsx"), sheet = "ndvi")
  as.matrix(ndvi[, -(1:2)])
}
ndvi <- list(A = read_area("A"), B = read_area("B"))
period <- colnames(ndvi$A)
n <- length(period)
cat("Observations per series:", n, "| start:", period[1], "| end:", period[n],
    "| pixels - A:", nrow(ndvi$A), "B:", nrow(ndvi$B), "\n")
mean_series <- lapply(ndvi, colMeans) #mean series of each area (all pixels)


# Function: BFAST and breakpoints of a series --------------------------------------------------
fit_bfast <- function(y) {
  yt <- ts(y, start = c(2000, 1), frequency = freq)
  fit <- bfast(yt, h = h, season = season_model, max.iter = max_iter)
  o <- fit$output[[length(fit$output)]]
  Tt <- as.numeric(o$Tt)
  time_yr <- as.numeric(time(yt))

  #Trend breakpoints (Vt); none if the OLS-MOSUM test is not significant
  bp <- if (fit$nobp$Vt) integer(0) else o$Vt.bp

  #Trend segments and their slopes (NDVI/year)
  start_seg <- c(1, bp + 1); end_seg <- c(bp, n)
  slope <- (Tt[end_seg] - Tt[start_seg]) / (time_yr[end_seg] - time_yr[start_seg])

  #Magnitude of each break: trend jump between the end of a segment and the start of the next one
  mag <- Tt[bp + 1] - Tt[bp]

  list(fit = fit,
       segments = data.frame(start = time_yr[start_seg], end = time_yr[end_seg], slope_per_year = slope),
       breaks = if (length(bp) > 0) data.frame(obs = bp, year = time_yr[bp], magnitude = mag,
                                               slope_before = slope[-length(slope)], slope_after = slope[-1]) else NULL,
       n_seasonal_breaks = if (fit$nobp$Wt) 0 else length(o$Wt.bp))
}


# 1. BFAST on the mean series of each area -----------------------------------------------------
fit_mean <- lapply(mean_series, fit_bfast)
fit_mean$A$segments; fit_mean$A$breaks
fit_mean$B$segments; fit_mean$B$breaks


# 2. BFAST on the difference A - B (control for climate) ---------------------------------------
fit_diff <- fit_bfast(mean_series$A - mean_series$B)
fit_diff$segments; fit_diff$breaks

#Structural change test (OLS-MOSUM) on the trend of the difference, as in bfast (Yt - St from STL)
diff_ts <- ts(mean_series$A - mean_series$B, start = c(2000, 1), frequency = freq)
diff_v <- diff_ts - stl(diff_ts, "periodic")$time.series[, "seasonal"]
time_diff <- time(diff_ts)
p_mosum_diff <- sctest(efp(diff_v ~ time_diff, h = h, type = "OLS-MOSUM"))$p.value
cat("OLS-MOSUM on the difference A - B: p =", round(p_mosum_diff, 3), "\n")


# 3. Linear trend 2000-2006 on the deseasonalised (STL) series ---------------------------------
#Complements BFAST: a long first segment can mask an early decline
trend_2000_2006 <- do.call(rbind, lapply(c("A", "B", "diff"), function(s) {
  y <- if (s == "diff") mean_series$A - mean_series$B else mean_series[[s]]
  yt <- ts(y, start = c(2000, 1), frequency = freq)
  v <- as.numeric(yt - stl(yt, "periodic")$time.series[, "seasonal"])
  tt <- as.numeric(time(yt)); sel <- tt < 2006.25
  cf <- summary(lm(v[sel] ~ tt[sel]))$coefficients
  data.frame(series = s, slope_per_year = cf[2, 1], std_error = cf[2, 2], p = cf[2, 4])
}))
trend_2000_2006


# Figures (Fig. 5A: Area A; Fig. 5B: Area B) ---------------------------------------------------
for (a in c("A", "B")) {
  png_pub(paste0("Fig5", a), width_double, 150, pointsize = 13) #plot.bfast has 4 panels (text x 0.66)
  plot(fit_mean[[a]]$fit, type = "components", main = paste("Area", a), ANOVA = TRUE)
  close_pub(paste0("Fig5", a))
}


# Export results -------------------------------------------------------------------------------
#bfast_mean_series_results.xlsx -> sheets:
#  segments        : trend segments and slopes of the mean series of A, B and the difference
#  breaks          : trend breakpoints of the mean series of A, B and the difference
#  trend_2000_2006 : linear trend 2000-2006 on the deseasonalised series
#  mosum_diff      : OLS-MOSUM test on the difference A - B
fits <- list(A = fit_mean$A, B = fit_mean$B, diff = fit_diff)
segments <- do.call(rbind, lapply(names(fits), function(s) data.frame(series = s, fits[[s]]$segments)))
breaks <- do.call(rbind, lapply(names(fits), function(s)
  if (!is.null(fits[[s]]$breaks)) data.frame(series = s, fits[[s]]$breaks)))
write_xlsx(list(segments = segments,
                breaks = if (is.null(breaks)) data.frame() else breaks,
                trend_2000_2006 = trend_2000_2006,
                mosum_diff = data.frame(test = "OLS-MOSUM (trend of A - B)", p = p_mosum_diff)),
           "data/processed/bfast_mean_series_results.xlsx")
