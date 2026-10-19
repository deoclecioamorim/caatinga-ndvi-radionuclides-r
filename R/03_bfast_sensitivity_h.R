#Manuscript: Section 2.2.2 (sensitivity to h); Results 3.2 (robustness); Fig. A.1 | run after 01
#Sensitivity of BFAST to the parameter h - mean bimonthly NDVI series (2000-2022) ------------------
#h = minimum segment size, as a fraction of the series length (138 bimonths)
#Same settings as 02_bfast_mean_series.R, repeated for each value of h

#Control options --------------------------------------------------------------------------------
options(OutDec = ".")
rm(list = ls(all = TRUE))


# Packages -------------------------------------------------------------------------------------
for (p in c("bfast", "readxl", "writexl", "ggplot2")) {
  if (!require(p, character.only = TRUE)) { install.packages(p, dependencies = TRUE); library(p, character.only = TRUE) }
}
source("R/00_publication_theme.R") #figure standard (font and print size)


# Parameters -----------------------------------------------------------------------------------
freq <- 6 #6 bimonths per year
h_values <- c(0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45) #~7 to 62 bimonths
season_model <- "harmonic"
max_iter <- 1
year_agriculture <- 2003 #agricultural land use detected by the LULC analysis in Area A


# Data: mean series of each area (all pixels) and difference A - B ---------------------------
read_area <- function(area) {
  ndvi <- read_excel(paste0("data/processed/ndvi_bimonthly_area_", area, ".xlsx"), sheet = "ndvi")
  colMeans(as.matrix(ndvi[, -(1:2)]))
}
mean_series <- list(A = read_area("A"), B = read_area("B"))
series <- list(A = mean_series$A, B = mean_series$B, diff = mean_series$A - mean_series$B)
n <- length(mean_series$A)


# Function: BFAST of a series ------------------------------------------------------------------
fit_bfast <- function(y, h) {
  yt <- ts(y, start = c(2000, 1), frequency = freq)
  fit <- tryCatch(bfast(yt, h = h, season = season_model, max.iter = max_iter), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  o <- fit$output[[length(fit$output)]]
  Tt <- as.numeric(o$Tt); Nt <- as.numeric(o$Nt)
  time_yr <- as.numeric(time(yt))
  bp <- if (fit$nobp$Vt) integer(0) else o$Vt.bp
  start_seg <- c(1, bp + 1); end_seg <- c(bp, n)
  slope <- (Tt[end_seg] - Tt[start_seg]) / (time_yr[end_seg] - time_yr[start_seg])
  mag <- Tt[bp + 1] - Tt[bp]
  list(
    Tt = Tt, time_yr = time_yr,
    summary = data.frame(n_trend_breaks = length(bp),
                         n_seasonal_breaks = if (fit$nobp$Wt) 0 else length(o$Wt.bp),
                         slope_first_segment = slope[1], end_first_segment = time_yr[end_seg[1]],
                         remainder_sd = sd(Nt)),
    breaks = if (length(bp) > 0) data.frame(year = time_yr[bp], magnitude = mag,
                                            slope_before = slope[-length(slope)], slope_after = slope[-1]) else NULL
  )
}


# Loop over the values of h --------------------------------------------------------------------
res_series <- list(); res_trend <- list(); res_breaks <- list()
for (h in h_values) {
  cat("\n===== h =", h, "(minimum segment ~", round(h * n), "bimonths =", round(h * n / freq, 1), "years) =====\n")
  for (s in names(series)) {
    r <- fit_bfast(series[[s]], h)
    if (is.null(r)) { cat("Failed:", s, "\n"); next }
    q <- r$breaks
    cat(sprintf("%-4s breakpoints: %s\n", s,
                if (is.null(q)) "none" else paste(sprintf("%.2f (mag %.3f; %.3f -> %.3f/yr)",
                                                          q$year, q$magnitude, q$slope_before, q$slope_after), collapse = " | ")))
    res_series[[length(res_series) + 1]] <- data.frame(h = h, series = s, r$summary)
    if (!is.null(q)) res_breaks[[length(res_breaks) + 1]] <- data.frame(h = h, series = s, q)
    res_trend[[length(res_trend) + 1]] <- data.frame(h = h, series = s, time_yr = r$time_yr, Tt = r$Tt)
  }
}
res_series <- do.call(rbind, res_series)
res_breaks <- do.call(rbind, res_breaks)
res_trend <- do.call(rbind, res_trend)
cat("\n===== Summary of the mean series =====\n"); print(res_series, digits = 3)


# Figure A.1: trend (Tt) of the mean series and of the difference, with h grouped by result ----
#Values of h producing the same trend are drawn as a single curve

#Compact label: consecutive values of h_values become a range (e.g. 0.10–0.25)
label_h <- function(hs) {
  i <- sort(match(hs, h_values))
  blocks <- split(i, cumsum(c(1, diff(i) != 1)))
  parts <- sapply(blocks, function(b) {
    if (length(b) == 1) format(h_values[b], nsmall = 2) else paste0(format(h_values[b[1]], nsmall = 2), "–", format(h_values[b[length(b)]], nsmall = 2))
  })
  paste("h =", paste(parts, collapse = ", "))
}

groups <- do.call(rbind, lapply(split(res_trend, res_trend$series), function(d) {
  key <- tapply(d$Tt, d$h, function(x) paste(round(x, 6), collapse = ";"))
  ids <- split(as.numeric(names(key)), key)
  ids <- ids[order(-lengths(ids))] #group with most values of h first
  do.call(rbind, lapply(seq_along(ids), function(k) {
    dk <- d[d$h == ids[[k]][1], ]
    data.frame(series = d$series[1], time_yr = dk$time_yr, Tt = dk$Tt,
               group = k, label = label_h(ids[[k]]))
  }))
}))
groups$series <- factor(groups$series, levels = c("A", "B", "diff"),
                        labels = c("Area A", "Area B", "Difference (A - B)"))
#Legend inside each panel: labels stacked above the curves, to the right of the 2003 line
labels <- do.call(rbind, lapply(split(groups, groups$series), function(d) {
  if (nrow(d) == 0) return(NULL)
  r <- diff(range(d$Tt)); g_ <- sort(unique(d$group))
  data.frame(series = d$series[1], group = g_, label = d$label[match(g_, d$group)],
             time_yr = year_agriculture + 0.4, Tt = max(d$Tt) + r * 0.14 * rev(seq_along(g_)))
}))

g <- ggplot(groups, aes(time_yr, Tt, colour = factor(group))) +
  geom_vline(xintercept = year_agriculture, linetype = "dashed", colour = "grey40") +
  geom_line(data = groups[groups$group != 1, ], linewidth = 0.8) +
  #main group dashed and drawn on top, so it remains visible where it overlaps another curve
  geom_line(data = groups[groups$group == 1, ], linewidth = 0.8, linetype = "22") +
  geom_text(data = labels, aes(label = label), hjust = 0, size = text_pub(), family = "serif",
            fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ series, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "4" = "#E7298A")) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  labs(x = "Year", y = "NDVI trend component (Tt)") +
  theme_pub() + theme(legend.position = "none")
save_pub("FigA1", g, width_mid, 150)


# Export results -------------------------------------------------------------------------------
write_xlsx(list(mean_series = res_series, breaks = if (is.null(res_breaks)) data.frame() else res_breaks),
           "data/processed/bfast_sensitivity_h.xlsx")
