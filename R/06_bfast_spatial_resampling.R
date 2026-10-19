#Manuscript: Section 2.2.2 (comparison with spatially thinned samples); Results 3.2; Figs. A.3 and A.4 | run after 05
#Comparison of pixel-level BFAST metrics between areas using repeated, spatially thinned samples ----
#Input: data/processed/bfast_pixel_metrics.xlsx (sheet 'pixels'; BFAST metrics for all pixels,
#       created by 05_bfast_pixel_metrics.R)
#
#Rationale: neighbouring pixels (30 m) are not independent (spatial autocorrelation inflates p-values),
#and a single random sample may depend on the seed.
#
#Steps:
#1. Correlogram (Moran's I by distance class) to assess the range of spatial autocorrelation (Fig. A.3)
#2. In each iteration, spatially thinned sampling (minimum distance d_min between pixels) in each area
#   (several values of d_min are accepted; only 90 m is used)
#3. Wilcoxon test (with Holm correction) and Cliff's delta between areas in each iteration
#4. Summary of the iterations: median and 95% interval of delta, % of iterations with adjusted p < 0.05 (Fig. A.4)

#Control options --------------------------------------------------------------------------------
options(OutDec = ".")
rm(list = ls(all = TRUE))


# Packages -------------------------------------------------------------------------------------
for (p in c("readxl", "writexl", "ggplot2")) {
  if (!require(p, character.only = TRUE)) { install.packages(p, dependencies = TRUE); library(p, character.only = TRUE) }
}
source("R/00_publication_theme.R") #figure standard (font and print size)


# Parameters -----------------------------------------------------------------------------------
d_min <- 90 #minimum distance between sampled pixels (m) = 3 pixels (reduces, but does not eliminate, autocorrelation)
n_iter <- 1000 #number of samples
seed <- 2026
dist_classes <- seq(0, 300, by = 30) #distance classes of the correlogram (m)

metrics <- c(n_trend_breaks = "Trend breaks (n)",
             largest_decline = "Largest abrupt decline",
             pct_time_negative_trend = "Time under negative trend (%)",
             slope_2000_2006 = "NDVI trend 2000-2006",
             net_trend_change = "Net trend change 2000-2022",
             seasonal_amplitude = "Seasonal amplitude",
             remainder_sd = "Remainder SD",
             mean_ndvi = "Mean NDVI",
             decline_2002_2006 = "Abrupt decline 2002-2006",
             decline_2011_2012 = "Abrupt decline 2011-2012")


# Metrics of all pixels ------------------------------------------------------------------------
pix <- as.data.frame(read_excel("data/processed/bfast_pixel_metrics.xlsx", sheet = "pixels"))
pix$decline_2002_2006 <- as.numeric(as.logical(pix$decline_2002_2006))
pix$decline_2011_2012 <- as.numeric(as.logical(pix$decline_2011_2012))
table(pix$area)


# 1. Correlogram: Moran's I by distance class --------------------------------------------------
moran_class <- function(z, D, lower, upper) {
  W <- (D > lower & D <= upper) * 1
  if (sum(W) == 0) return(NA)
  z <- z - mean(z)
  (length(z) / sum(W)) * sum(W * outer(z, z)) / sum(z^2)
}
correlogram <- do.call(rbind, lapply(c("A", "B"), function(a) {
  p <- pix[pix$area == a, ]
  D <- as.matrix(dist(p[, c("x", "y")]))
  do.call(rbind, lapply(c("remainder_sd", "slope_2000_2006", "largest_decline", "mean_ndvi"), function(v) {
    data.frame(area = a, metric = v,
               dist_upper = dist_classes[-1],
               moran = sapply(seq_along(dist_classes)[-1], function(k)
                 moran_class(p[[v]], D, dist_classes[k - 1], dist_classes[k])))
  }))
}))
cat("\n===== Moran's I by distance class (upper limit, m) =====\n")
print(reshape(correlogram, idvar = c("area", "metric"), timevar = "dist_upper", direction = "wide"), digits = 2)


# 2. Spatially thinned sampling ----------------------------------------------------------------
#Random order; a pixel is accepted if it is >= d_min from all pixels already accepted
thinned_sample <- function(xy, d_min) {
  ord <- sample(nrow(xy))
  accepted <- ord[1]
  for (i in ord[-1]) {
    if (all((xy[accepted, 1] - xy[i, 1])^2 + (xy[accepted, 2] - xy[i, 2])^2 >= d_min^2)) accepted <- c(accepted, i)
  }
  accepted
}

cliff <- function(a, b) mean(sign(outer(a, b, "-"))) #Cliff's delta: P(A > B) - P(A < B)
pA <- pix[pix$area == "A", ]; pB <- pix[pix$area == "B", ]
xyA <- as.matrix(pA[, c("x", "y")]); xyB <- as.matrix(pB[, c("x", "y")])

#|delta| < 0.147 negligible, < 0.33 small, < 0.474 medium, >= 0.474 large
cliff_class <- function(d) as.character(cut(abs(d), c(-Inf, 0.147, 0.33, 0.474, Inf),
                                            labels = c("negligible", "small", "medium", "large"), right = FALSE))

#Resampling for a given minimum distance
resample_metrics <- function(d) {
  set.seed(seed)
  t0 <- Sys.time()
  it_d <- do.call(rbind, lapply(seq_len(n_iter), function(it) {
    a <- pA[thinned_sample(xyA, d), ]
    b <- pB[thinned_sample(xyB, d), ]
    p <- sapply(names(metrics), function(v) wilcox.test(a[[v]], b[[v]], exact = FALSE)$p.value)
    p[is.na(p)] <- 1
    data.frame(d_min = d, iteration = it, metric = names(metrics), n_A = nrow(a), n_B = nrow(b),
               delta = sapply(names(metrics), function(v) cliff(a[[v]], b[[v]])),
               p = p, p_adjusted = p.adjust(p, method = "holm"),
               mean_A = sapply(names(metrics), function(v) mean(a[[v]])),
               mean_B = sapply(names(metrics), function(v) mean(b[[v]])),
               median_A = sapply(names(metrics), function(v) median(a[[v]])),
               median_B = sapply(names(metrics), function(v) median(b[[v]])))
  }))
  cat("\nd_min =", d, "m | time:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min",
      "| pixels per sample - A:", paste(range(it_d$n_A), collapse = "-"), "| B:", paste(range(it_d$n_B), collapse = "-"), "\n")
  it_d
}
iterations <- do.call(rbind, lapply(d_min, resample_metrics))


# 3. Summary of the iterations -----------------------------------------------------------------
resampling_summary <- do.call(rbind, lapply(d_min, function(dd) do.call(rbind, lapply(names(metrics), function(v) {
  d <- iterations[iterations$d_min == dd & iterations$metric == v, ]
  data.frame(d_min = dd, metric = v,
             n_A_median = median(d$n_A), n_B_median = median(d$n_B),
             mean_A = median(d$mean_A), mean_B = median(d$mean_B),
             median_A = median(d$median_A), median_B = median(d$median_B),
             delta_median = median(d$delta),
             delta_ci95_lower = unname(quantile(d$delta, 0.025)),
             delta_ci95_upper = unname(quantile(d$delta, 0.975)),
             effect = cliff_class(median(d$delta)),
             p_adjusted_median = median(d$p_adjusted),
             pct_iter_significant = mean(d$p_adjusted < 0.05) * 100,
             delta_all_pixels = cliff(pA[[v]], pB[[v]]))
}))))
cat("\n===== Summary of", n_iter, "spatially thinned samples per minimum distance =====\n")
print(resampling_summary[, c("d_min", "metric", "n_A_median", "n_B_median", "delta_median", "delta_ci95_lower",
                             "delta_ci95_upper", "effect", "pct_iter_significant", "delta_all_pixels")], digits = 2)


# Figures --------------------------------------------------------------------------------------
#Fig. A.4: distribution of Cliff's delta across samples, by minimum distance
iterations$label <- factor(metrics[iterations$metric], levels = rev(metrics))
iterations$d_label <- factor(paste0(iterations$d_min, " m"), levels = paste0(d_min, " m"))
g1 <- ggplot(iterations, aes(delta, label, fill = d_label)) +
  annotate("rect", xmin = -0.147, xmax = 0.147, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = 0.3, width = 0.7, position = position_dodge(0.8)) +
  scale_fill_manual(values = if (length(d_min) == 1) "#9ECAE1" else c("#DEEBF7", "#9ECAE1", "#3182BD")) +
  labs(x = "Cliff's delta (Area A vs. Area B)", y = NULL, fill = "Minimum distance",
       caption = if (length(d_min) == 1)
         paste0(format(n_iter, big.mark = ","), " spatially thinned samples (minimum distance ", d_min, " m); grey band: negligible effect")
       else paste0(format(n_iter, big.mark = ","), " spatially thinned samples per minimum distance; grey band: negligible effect")) +
  theme_pub() + theme(legend.position = if (length(d_min) == 1) "none" else "bottom")
save_pub("FigA4", g1, width_mid, 130)

#Fig. A.3: correlogram
correlogram$label <- metrics[correlogram$metric]
correlogram$area <- paste("Area", correlogram$area)
g2 <- ggplot(correlogram, aes(dist_upper, moran, colour = label)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_vline(xintercept = d_min, linetype = "dashed", colour = "grey50") +
  geom_line() + geom_point(size = 1.2) +
  facet_wrap(~ area) +
  labs(x = "Distance class (upper limit, m)", y = "Moran's I", colour = NULL) +
  theme_pub() + theme(legend.position = "bottom")
save_pub("FigA3", g2, width_double, 100)


# Export results -------------------------------------------------------------------------------
write_xlsx(list(summary = resampling_summary, correlogram = correlogram[, c("area", "metric", "dist_upper", "moran")],
                iterations = iterations[, setdiff(names(iterations), c("label", "d_label"))]),
           "data/processed/bfast_spatial_resampling.xlsx")
