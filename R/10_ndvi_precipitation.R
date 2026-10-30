#Manuscript: NDVI x precipitation (reviewer request); Figs. A.6 and A.7 | run after 02
#Association between bimonthly precipitation (CHIRPS) and NDVI (2000-2022) - Areas A and B -------
#Inputs: data/raw/ndvi_area_A and data/raw/ndvi_area_B (138 GEE images each, observed values only)
#        data/raw/Cabrobo_CHIRPS_2000_2022 (138 bimonthly CHIRPS totals, 0.05 deg, region enclosing both areas)
#        data/processed/bfast_mean_series_results.xlsx (breakpoints, created by 02_bfast_mean_series.R)
#
#Steps:
#1. Bimonthly NDVI of each area = median of the valid pixels of each image (no interpolation),
#   with the % of the polygon covered by valid pixels; precipitation = mean of the CHIRPS cells
#   (one regional series applied to both areas)
#2. Bimonths with NDVI coverage below 'coverage_main' are set to NA; anomalies = value minus the
#   2000-2022 mean of the same bimonth (removes the mean seasonal cycle)
#3. Pearson correlation and linear regression NDVI anomaly ~ precipitation anomaly at lags 0, 1 and 2
#   bimonths, plus Spearman correlation (robust to extreme bimonths); p-values with Newey-West (HAC)
#   standard errors (on the ranks for Spearman); Benjamini-Hochberg over the 6 tests (2 areas x 3 lags)
#   of each scenario
#4. Uncertainty: circular moving block bootstrap on the complete calendar (blocks of 6 bimonths = 1 year),
#   95% percentile intervals of r, slope and R2; the observed values are the estimates
#5. Sensitivity: coverage threshold (50, 70, 90%), detrended anomalies (piecewise trend with the BFAST
#   break, so that the 2011-2012 break does not drive the correlation) and block length (3, 6, 12)
#6. Precipitation in the bimonth of each BFAST trend break and in the two previous bimonths
#Association only: no causal inference.

#Control options --------------------------------------------------------------------------------
options(OutDec = ".")
rm(list = ls(all = TRUE))


# Packages -------------------------------------------------------------------------------------
for (p in c("terra", "readxl", "writexl", "sandwich", "lmtest", "ggplot2")) {
  if (!require(p, character.only = TRUE)) { install.packages(p, dependencies = TRUE); library(p, character.only = TRUE) }
}
source("R/00_publication_theme.R") #figure standard (font and print size)


# Parameters -----------------------------------------------------------------------------------
coverage_main <- 70 #minimum % of valid pixels to keep a bimonth (main analysis)
coverage_sens <- c(50, 70, 90) #sensitivity of the coverage threshold (50 = criterion of script 01)
lags <- 0:2 #precipitation lags (bimonths)
n_boot <- 1000 #bootstrap replicates
block_main <- 6 #block length (bimonths) = 1 year
block_sens <- c(3, 6, 12)
min_pairs <- 12 #minimum number of pairs to compute the statistics
set.seed(20261020)

years <- 2000:2022
calendar <- data.frame(year = rep(years, each = 6), bimonth = rep(1:6, times = length(years)))
calendar$t <- seq_len(nrow(calendar))
calendar$time <- calendar$year + (calendar$bimonth - 1) / 6
bimonth_names <- c("Jan-Feb", "Mar-Apr", "May-Jun", "Jul-Aug", "Sep-Oct", "Nov-Dec")


# 1. Data --------------------------------------------------------------------------------------
#NDVI: median of the valid pixels of each image and % coverage of the polygon
ndvi_series <- function(folder) {
  files <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
  period <- sub(".*_harmOLI_(\\d{4})_(\\d{2})_(\\d{2})\\.tif$", "\\1_\\2_\\3", basename(files))
  files <- files[order(period)]
  period <- sort(period)
  stopifnot(length(files) == nrow(calendar))

  #Load the images and align them to the grid of the first one (as in 01_ndvi_extraction.R)
  images <- lapply(files, rast)
  ref <- images[[1]]
  images <- lapply(images, function(r) {
    if (!compareGeom(r, ref, stopOnError = FALSE)) {
      if (crs(r) != crs(ref)) r <- project(r, ref)
      r <- resample(r, ref, method = "near")
    }
    r
  })
  values <- terra::values(rast(images), mat = TRUE)
  values[values < -1 | values > 1] <- NA #declared NoData (-9999)
  values <- values[rowSums(!is.na(values)) > 0, , drop = FALSE] #pixels of the polygon

  data.frame(year = as.integer(substr(period, 1, 4)),
             bimonth = (as.integer(substr(period, 6, 7)) + 1) / 2,
             ndvi = apply(values, 2, function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)),
             coverage = colMeans(!is.na(values)) * 100)
}
ndvi <- list(A = ndvi_series("data/raw/ndvi_area_A"), B = ndvi_series("data/raw/ndvi_area_B"))

#CHIRPS: regional mean of the cells (file names end in YYYY_Mon_Mon.tif)
files <- list.files("data/raw/Cabrobo_CHIRPS_2000_2022", pattern = "\\.tif$", full.names = TRUE)
rain <- data.frame(year = as.integer(sub(".*_(\\d{4})_\\w{3}_\\w{3}\\.tif$", "\\1", basename(files))),
                   bimonth = match(sub(".*_\\d{4}_(\\w{3})_(\\w{3})\\.tif$", "\\1-\\2", basename(files)),
                                   bimonth_names),
                   rain_mm = sapply(files, function(f) global(rast(f), "mean", na.rm = TRUE)[1, 1]),
                   row.names = NULL)
stopifnot(nrow(rain) == nrow(calendar), !anyNA(rain), !anyDuplicated(rain[, c("year", "bimonth")]),
          all(rain$rain_mm >= 0))

#BFAST trend breaks of the mean series (obs = last bimonth before the trend jump)
breaks <- as.data.frame(read_excel("data/processed/bfast_mean_series_results.xlsx", sheet = "breaks"))


# 2. Series and anomalies ----------------------------------------------------------------------
anomaly <- function(x, bimonth) x - ave(x, bimonth, FUN = function(v) mean(v, na.rm = TRUE))
lag_vec <- function(x, k) if (k == 0) x else c(rep(NA_real_, k), head(x, -k))

series <- lapply(names(ndvi), function(a) {
  z <- merge(calendar, ndvi[[a]], by = c("year", "bimonth"), all.x = TRUE)
  z <- merge(z, rain, by = c("year", "bimonth"), all.x = TRUE)
  z <- z[order(z$t), ]
  z$area <- a
  z$rain_anom <- anomaly(z$rain_mm, z$bimonth)
  for (cv in coverage_sens) {
    y <- ifelse(is.finite(z$coverage) & z$coverage >= cv, z$ndvi, NA_real_)
    z[[paste0("ndvi_anom_cov", cv)]] <- anomaly(y, z$bimonth)
  }

  #Detrended anomalies: piecewise linear trend with the BFAST break of the area (level and slope
  #may change at the break) removed from NDVI; linear trend removed from precipitation
  bp <- breaks$obs[breaks$series == a]
  segment <- factor(z$t > bp)
  y <- z[[paste0("ndvi_anom_cov", coverage_main)]]
  z$ndvi_anom_detrended <- residuals(lm(y ~ z$t * segment, na.action = na.exclude))
  z$rain_anom_detrended <- residuals(lm(z$rain_anom ~ z$t, na.action = na.exclude))
  z
})
names(series) <- names(ndvi)

#Bimonths kept per coverage threshold
coverage_summary <- do.call(rbind, lapply(names(series), function(a) {
  z <- series[[a]]
  data.frame(area = a, threshold_pct = coverage_sens,
             bimonths_kept = sapply(coverage_sens, function(cv) sum(is.finite(z$ndvi) & z$coverage >= cv)),
             bimonths_total = nrow(z))
}))
coverage_summary


# 3-4. Statistics, HAC test and block bootstrap ------------------------------------------------
fit_stats <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < min_pairs || sd(x) == 0 || sd(y) == 0)
    return(c(r = NA_real_, rho = NA_real_, slope = NA_real_, r2 = NA_real_, n = length(x)))
  r <- cor(x, y)
  c(r = r, rho = cor(x, y, method = "spearman"), slope = unname(coef(lm(y ~ x))[2]), r2 = r^2, n = length(x))
}

#Circular moving block bootstrap on the complete calendar: blocks keep contiguous bimonths (and their
#NAs); lags are computed before resampling, so each pair keeps its precipitation-NDVI timing
block_index <- function(n, L) {
  starts <- sample.int(n, size = ceiling(n / L), replace = TRUE)
  idx <- unlist(lapply(starts, function(s) ((s - 1 + 0:(L - 1)) %% n) + 1))
  idx[seq_len(n)]
}

#Newey-West standard errors with automatic bandwidth (Newey and West, 1994); the Spearman test is the
#same HAC test applied to the ranks of the pairs (robust to extreme bimonths, e.g. Jan-Feb 2004)
hac_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  hac <- function(fit) lmtest::coeftest(fit, vcov. = sandwich::NeweyWest(fit, prewhite = FALSE, adjust = TRUE))
  ct <- hac(lm(y[ok] ~ x[ok]))
  ct_rank <- hac(lm(rank(y[ok]) ~ rank(x[ok])))
  c(se_hac = ct[2, 2], p_hac = ct[2, 4], p_hac_rho = ct_rank[2, 4])
}

analyse <- function(x, y, L) {
  obs <- fit_stats(x, y)
  boot <- t(replicate(n_boot, { ix <- block_index(length(x), L); fit_stats(x[ix], y[ix]) }))
  ci <- apply(boot[, c("r", "rho", "slope", "r2")], 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  data.frame(n_pairs = obs[["n"]], r = obs[["r"]], r_ci_low = ci[1, "r"], r_ci_high = ci[2, "r"],
             rho = obs[["rho"]], rho_ci_low = ci[1, "rho"], rho_ci_high = ci[2, "rho"],
             slope = obs[["slope"]], slope_ci_low = ci[1, "slope"], slope_ci_high = ci[2, "slope"],
             r2 = obs[["r2"]], r2_ci_low = ci[1, "r2"], r2_ci_high = ci[2, "r2"],
             t(hac_test(x, y)), block = L, valid_replicates = sum(is.finite(boot[, "r"])))
}

#Scenarios: main (coverage 70%), coverage 50 and 90%, detrended anomalies (coverage 70%)
scenarios <- list(main = c(y = paste0("ndvi_anom_cov", coverage_main), x = "rain_anom"),
                  coverage_50 = c(y = "ndvi_anom_cov50", x = "rain_anom"),
                  coverage_90 = c(y = "ndvi_anom_cov90", x = "rain_anom"),
                  detrended = c(y = "ndvi_anom_detrended", x = "rain_anom_detrended"))

results <- list()
for (s in names(scenarios)) for (a in names(series)) for (k in lags) {
  z <- series[[a]]
  res <- analyse(lag_vec(z[[scenarios[[s]]["x"]]], k), z[[scenarios[[s]]["y"]]], block_main)
  results[[length(results) + 1]] <- data.frame(scenario = s, area = a, lag_bimonths = k, res)
}
results <- do.call(rbind, results)
results$p_bh <- ave(results$p_hac, results$scenario, FUN = function(p) p.adjust(p, method = "BH"))
results$p_bh_rho <- ave(results$p_hac_rho, results$scenario, FUN = function(p) p.adjust(p, method = "BH"))
rownames(results) <- NULL
results[, c("scenario", "area", "lag_bimonths", "n_pairs", "r", "r_ci_low", "r_ci_high", "p_bh",
            "rho", "rho_ci_low", "rho_ci_high", "p_bh_rho")]

#Sensitivity to the block length (main scenario)
block_results <- list()
for (L in block_sens) for (a in names(series)) for (k in lags) {
  z <- series[[a]]
  res <- analyse(lag_vec(z$rain_anom, k), z[[scenarios$main["y"]]], L)
  block_results[[length(block_results) + 1]] <- data.frame(area = a, lag_bimonths = k,
                                                            res[, c("block", "r", "r_ci_low", "r_ci_high",
                                                                    "rho", "rho_ci_low", "rho_ci_high",
                                                                    "slope_ci_low", "slope_ci_high")])
}
block_results <- do.call(rbind, block_results)
block_results


# 6. Precipitation at the BFAST trend breaks ---------------------------------------------------
#For the bimonth of the break and the two previous ones: total, 2000-2022 mean of the same bimonth,
#anomaly, rank among the 23 years (1 = driest) and precipitation of the 12 months ending there
rain_12m <- stats::filter(rain$rain_mm[order(rain$year, rain$bimonth)], rep(1, 6), sides = 1)
rain_12m_mean <- mean(rain_12m, na.rm = TRUE)
break_rain <- list()
for (i in seq_len(nrow(breaks))) {
  a <- breaks$series[i]
  z <- series[[a]]
  for (k in lags) {
    j <- breaks$obs[i] - k
    same <- z$rain_mm[z$bimonth == z$bimonth[j]]
    clim <- mean(same)
    break_rain[[length(break_rain) + 1]] <- data.frame(
      area = a, break_year = z$year[breaks$obs[i]], break_bimonth = bimonth_names[z$bimonth[breaks$obs[i]]],
      break_magnitude = breaks$magnitude[i], lag_bimonths = k,
      year = z$year[j], bimonth = bimonth_names[z$bimonth[j]],
      rain_mm = z$rain_mm[j], climatology_mm = clim,
      anomaly_mm = z$rain_mm[j] - clim, anomaly_pct = 100 * (z$rain_mm[j] - clim) / clim,
      rank_driest = rank(same)[match(z$year[j], years)],
      rain_12m_mm = rain_12m[j], rain_12m_anomaly_pct = 100 * (rain_12m[j] - rain_12m_mean) / rain_12m_mean,
      ndvi_coverage_pct = z$coverage[j])
  }
}
break_rain <- do.call(rbind, break_rain)
break_rain


# Export results -------------------------------------------------------------------------------
#ndvi_precipitation_results.xlsx -> sheets:
#  series           : NDVI (median), coverage, precipitation and anomalies per area and bimonth
#  coverage_summary : bimonths kept per coverage threshold
#  correlation      : r, Spearman rho, slope, R2 (95% block bootstrap intervals), HAC p and BH p
#                     per scenario/area/lag
#  block_length     : sensitivity of the intervals to the block length (main scenario)
#  breaks_rain      : precipitation at the BFAST trend breaks and in the two previous bimonths
writexl::write_xlsx(list(series = do.call(rbind, series), coverage_summary = coverage_summary,
                         correlation = results, block_length = block_results, breaks_rain = break_rain),
                    "data/processed/ndvi_precipitation_results.xlsx")


# Fig. A6: NDVI anomaly x precipitation anomaly (main scenario) --------------------------------
lag_labels <- c("Lag 0 (same bimonth)", "Lag 1 (previous bimonth)", "Lag 2 (two bimonths before)")
fig_A6 <- do.call(rbind, lapply(names(series), function(a) do.call(rbind, lapply(lags, function(k) {
  z <- series[[a]]
  data.frame(area = paste("Area", a), lag = lag_labels[k + 1],
             x = lag_vec(z$rain_anom, k), y = z[[scenarios$main["y"]]])
}))))
fig_A6 <- fig_A6[is.finite(fig_A6$x) & is.finite(fig_A6$y), ]
fig_A6$lag <- factor(fig_A6$lag, levels = lag_labels)

labels_A6 <- subset(results, scenario == "main")
labels_A6 <- data.frame(area = paste("Area", labels_A6$area),
                        lag = factor(lag_labels[labels_A6$lag_bimonths + 1], levels = lag_labels),
                        label = sprintf("r = %.2f [%.2f, %.2f], p(BH) = %.3f\nρ = %.2f [%.2f, %.2f], p(BH) = %.3f",
                                        labels_A6$r, labels_A6$r_ci_low, labels_A6$r_ci_high, labels_A6$p_bh,
                                        labels_A6$rho, labels_A6$rho_ci_low, labels_A6$rho_ci_high,
                                        labels_A6$p_bh_rho))

g <- ggplot(fig_A6, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(size = 0.8, alpha = 0.6, colour = "grey30") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black", linewidth = 0.6) +
  geom_text(data = labels_A6, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.2,
            size = text_pub(7), family = "serif", lineheight = 0.9) +
  facet_grid(area ~ lag) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) + #room for the labels above the points
  labs(x = "Bimonthly precipitation anomaly (mm)", y = "NDVI anomaly") +
  theme_pub()
save_pub("FigA6", g, width = width_double, height = 120)


# Fig. A7: precipitation and NDVI anomalies with the BFAST trend breaks ------------------------
fig_rain <- series$A[, c("time", "rain_anom")]
fig_rain$panel <- "Precipitation anomaly (mm)"
fig_rain$sign <- ifelse(fig_rain$rain_anom >= 0, "positive", "negative")
fig_ndvi <- do.call(rbind, lapply(names(series), function(a)
  data.frame(time = series[[a]]$time, y = series[[a]][[scenarios$main["y"]]], area = paste("Area", a))))
fig_ndvi$panel <- "NDVI anomaly"
panels <- c("Precipitation anomaly (mm)", "NDVI anomaly")
fig_rain$panel <- factor(fig_rain$panel, levels = panels)
fig_ndvi$panel <- factor(fig_ndvi$panel, levels = panels)
fig_breaks <- data.frame(time = breaks$year, area = paste("Area", breaks$series))

g <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_col(data = fig_rain, aes(time, rain_anom, fill = sign), width = 1 / 6) +
  geom_line(data = fig_ndvi, aes(time, y, linetype = area), linewidth = 0.4, na.rm = TRUE) +
  geom_vline(data = fig_breaks, aes(xintercept = time, linetype = area), colour = "red3", linewidth = 0.5) +
  facet_wrap(~panel, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_fill_manual(values = c(positive = "steelblue", negative = "tan3"), guide = "none") +
  scale_linetype_manual(values = c("Area A" = "solid", "Area B" = "22"), name = NULL) +
  scale_x_continuous(breaks = seq(2000, 2022, 2), expand = c(0.01, 0)) +
  labs(x = "Year", y = NULL) +
  theme_pub() +
  theme(strip.placement = "outside", strip.background = element_blank(), legend.position = "top")
save_pub("FigA7", g, width = width_double, height = 120)
