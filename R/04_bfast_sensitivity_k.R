#Manuscript: Section 2.2.2 (sensitivity to K); Results 3.2 (robustness); Fig. A.2 | run after 01
#Sensitivity of BFAST to the harmonic order (K) - mean bimonthly NDVI series (2000-2022) -----------
#
#The bfast() function fixes K = 3 (co, si, co2, si2, co3, si3). Here the algorithm is replicated
#(bfast 1.7.2, season = "harmonic", decomp = "stl") with an adjustable K and no other changes.
#
#Notes:
#- With 6 obs/year, the maximum K is 3 (Nyquist). The term sin(2*pi*3*t/6) = sin(pi*t) = 0, so it is
#  removed when k = f/2 (the original bfast() keeps it, but it has no effect)
#- With max.iter = 1, trend breakpoints are searched in Yt - St(STL), so K only affects the
#  seasonal and remainder components. With max.iter > 1, the harmonic St of order K feeds back
#  into the trend. Therefore max.iter = 1 and 10 are tested
#
#Steps:
#1. Validation: bfast_k(K = 3) must reproduce the breakpoints of the original bfast()
#2. Choice of K by BIC (linear trend + harmonic model of order K, without breakpoints)
#3. Mean series and difference A - B: K = 1, 2, 3 x max.iter = 1, 10 x h = 0.10, 0.15, 0.20, 0.25

#Control options --------------------------------------------------------------------------------
options(OutDec = ".")
rm(list = ls(all = TRUE))


# Packages -------------------------------------------------------------------------------------
for (p in c("bfast", "strucchange", "readxl", "writexl", "ggplot2")) {
  if (!require(p, character.only = TRUE)) { install.packages(p, dependencies = TRUE); library(p, character.only = TRUE) }
}
source("R/00_publication_theme.R") #figure standard (font and print size)


# Parameters -----------------------------------------------------------------------------------
freq <- 6 #6 bimonths per year
k_values <- 1:3
iter_values <- c(1, 10)
h_values <- c(0.10, 0.15, 0.20, 0.25)
h_figure <- 0.15 #h shown in Fig. A.2
year_agriculture <- 2003


# Data: mean series of each area (all pixels) and difference A - B ---------------------------
read_area <- function(area) {
  ndvi <- read_excel(paste0("data/processed/ndvi_bimonthly_area_", area, ".xlsx"), sheet = "ndvi")
  colMeans(as.matrix(ndvi[, -(1:2)]))
}
mean_series <- list(A = read_area("A"), B = read_area("B"))
series <- list(A = mean_series$A, B = mean_series$B, diff = mean_series$A - mean_series$B)
n <- length(mean_series$A)


# Harmonic matrix of order K -------------------------------------------------------------------
harmonics <- function(n, f, K) {
  tl <- 1:n
  X <- NULL
  for (k in 1:K) {
    X <- cbind(X, cos(2 * pi * tl * k / f))
    colnames(X)[ncol(X)] <- paste0("co", k)
    if (2 * k != f) { #sine is zero at the Nyquist frequency
      X <- cbind(X, sin(2 * pi * tl * k / f))
      colnames(X)[ncol(X)] <- paste0("si", k)
    }
  }
  X
}


# BFAST with adjustable K (replica of bfast::bfast, season = "harmonic") -----------------------
bfast_k <- function(Yt, h = 0.15, K = 3, max.iter = 1, level = 0.05, type = "OLS-MOSUM") {
  ti <- time(Yt); f <- frequency(Yt); nn <- length(Yt)
  X <- harmonics(nn, f, K)
  St <- stl(Yt, "periodic")$time.series[, "seasonal"]
  Vt.bp <- 0; Wt.bp <- 0; CheckTimeTt <- 1; CheckTimeSt <- 1; i <- 0
  output <- list()

  while ((!identical(CheckTimeTt, Vt.bp) | !identical(CheckTimeSt, Wt.bp)) & i < max.iter) {
    CheckTimeTt <- Vt.bp; CheckTimeSt <- Wt.bp

    #Trend
    Vt <- Yt - St
    p.Vt <- sctest(efp(Vt ~ ti, h = h, type = type))
    if (p.Vt$p.value <= level) {
      bp.Vt <- breakpoints(Vt ~ ti, h = h)
      nobp.Vt <- is.na(breakpoints(bp.Vt)[1])
    } else { nobp.Vt <- TRUE; bp.Vt <- NA }
    if (nobp.Vt) {
      Tt <- fitted(lm(Vt ~ ti)); Vt.bp <- 0
    } else {
      Tt <- fitted(lm(Vt ~ breakfactor(bp.Vt) / ti))
      Vt.bp <- confint(bp.Vt, het.err = FALSE)$confint[, 2]
    }
    Tt <- ts(as.numeric(Tt), start = start(Yt), frequency = f)

    #Seasonality (harmonic model of order K)
    Wt <- as.numeric(Yt - Tt)
    p.Wt <- sctest(efp(Wt ~ X, h = h, type = type))
    if (p.Wt$p.value <= level) {
      bp.Wt <- breakpoints(Wt ~ X, h = h)
      nobp.Wt <- is.na(breakpoints(bp.Wt)[1])
    } else { nobp.Wt <- TRUE; bp.Wt <- NA }
    if (nobp.Wt) {
      St <- fitted(lm(Wt ~ X)); Wt.bp <- 0
    } else {
      St <- fitted(lm(Wt ~ X %in% breakfactor(bp.Wt)))
      Wt.bp <- confint(bp.Wt, het.err = FALSE)$confint[, 2]
    }
    St <- ts(as.numeric(St), start = start(Yt), frequency = f)

    i <- i + 1
    output[[i]] <- list(Tt = Tt, St = St, Nt = Yt - Tt - St, Vt.bp = Vt.bp, Wt.bp = Wt.bp)
  }
  list(output = output, nobp = list(Vt = nobp.Vt, Wt = nobp.Wt), iterations = i)
}


# Metrics of a series --------------------------------------------------------------------------
series_metrics <- function(y, h, K, max.iter) {
  yt <- ts(y, start = c(2000, 1), frequency = freq)
  fit <- tryCatch(bfast_k(yt, h = h, K = K, max.iter = max.iter), error = function(e) NULL)
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
    summary = data.frame(iterations = fit$iterations,
                         n_trend_breaks = length(bp),
                         n_seasonal_breaks = if (fit$nobp$Wt) 0 else length(o$Wt.bp),
                         slope_first_segment = slope[1],
                         remainder_sd = sd(Nt)),
    breaks = if (length(bp) > 0) data.frame(year = time_yr[bp], magnitude = mag,
                                            slope_before = slope[-length(slope)], slope_after = slope[-1]) else NULL
  )
}


# 1. Validation: K = 3 must reproduce the original bfast() -------------------------------------
validation <- do.call(rbind, lapply(names(series), function(s) {
  yt <- ts(series[[s]], start = c(2000, 1), frequency = freq)
  orig <- bfast(yt, h = 0.15, season = "harmonic", max.iter = 1)
  rep3 <- bfast_k(yt, h = 0.15, K = 3, max.iter = 1)
  data.frame(series = s,
             breaks_bfast = paste(orig$output[[1]]$Vt.bp, collapse = " "),
             breaks_replica = paste(rep3$output[[1]]$Vt.bp, collapse = " "),
             max_diff_Tt = max(abs(orig$output[[1]]$Tt - rep3$output[[1]]$Tt)),
             max_diff_St = max(abs(orig$output[[1]]$St - rep3$output[[1]]$St)))
}))
cat("\n===== Validation (h = 0.15, max.iter = 1) =====\n"); print(validation)


# 2. Choice of K by BIC ------------------------------------------------------------------------
bic_k <- function(y) {
  t <- seq_along(y)
  sapply(k_values, function(K) BIC(lm(y ~ t + harmonics(length(y), freq, K))))
}
bic_series <- do.call(rbind, lapply(names(series), function(s) {
  b <- bic_k(series[[s]])
  data.frame(series = s, BIC_K1 = b[1], BIC_K2 = b[2], BIC_K3 = b[3], K_selected = which.min(b))
}))
cat("\n===== BIC by K (mean series) =====\n"); print(bic_series)


# 3. Mean series and difference ----------------------------------------------------------------
res_series <- list(); res_breaks <- list(); res_trend <- list()
cat("\n===== Mean series =====\n")
for (it in iter_values) for (K in k_values) for (h in h_values) for (s in names(series)) {
  r <- series_metrics(series[[s]], h, K, it)
  if (is.null(r)) { cat("Failed:", s, "K", K, "h", h, "iter", it, "\n"); next }
  q <- r$breaks
  cat(sprintf("iter %2d | K %d | h %.2f | %-4s: %s\n", it, K, h, s,
              if (is.null(q)) "none" else paste(sprintf("%.2f (mag %.3f; %.3f -> %.3f/yr)",
                                                        q$year, q$magnitude, q$slope_before, q$slope_after), collapse = " | ")))
  res_series[[length(res_series) + 1]] <- data.frame(max_iter = it, K = K, h = h, series = s, r$summary)
  if (!is.null(q)) res_breaks[[length(res_breaks) + 1]] <- data.frame(max_iter = it, K = K, h = h, series = s, q)
  if (h == h_figure) res_trend[[length(res_trend) + 1]] <- data.frame(max_iter = it, K = K, series = s, time_yr = r$time_yr, Tt = r$Tt)
}
res_series <- do.call(rbind, res_series)
res_breaks <- do.call(rbind, res_breaks)
res_trend <- do.call(rbind, res_trend)


# Figure A.2: trend (Tt) of the mean series by K and max.iter (h = 0.15) -----------------------
res_trend$series <- factor(res_trend$series, levels = c("A", "B", "diff"),
                           labels = c("Area A", "Area B", "Difference (A - B)"))
res_trend$max_iter <- paste("max.iter =", res_trend$max_iter)

#Values of K producing the same trend (maximum difference < tol_k NDVI) are drawn as a single curve
tol_k <- 0.005
label_k <- function(ks) {
  ks <- sort(ks)
  blocks <- split(ks, cumsum(c(1, diff(ks) != 1)))
  paste("K =", paste(sapply(blocks, function(b) if (length(b) == 1) b else paste0(b[1], "–", b[length(b)])), collapse = ", "))
}
groups_k <- do.call(rbind, lapply(split(res_trend, list(res_trend$series, res_trend$max_iter), drop = TRUE), function(d) {
  #Sequential grouping: each K joins the first group whose curve differs by less than tol_k
  ids <- list()
  for (k in sort(unique(d$K))) {
    tk <- d$Tt[d$K == k]
    j <- which(vapply(ids, function(g) max(abs(d$Tt[d$K == g[1]] - tk)) < tol_k, logical(1)))
    if (length(j) > 0) ids[[j[1]]] <- c(ids[[j[1]]], k) else ids[[length(ids) + 1]] <- k
  }
  ids <- ids[order(-lengths(ids))] #group with most values of K first
  do.call(rbind, lapply(seq_along(ids), function(k) {
    dk <- d[d$K == ids[[k]][1], ]
    data.frame(series = d$series[1], max_iter = d$max_iter[1], time_yr = dk$time_yr, Tt = dk$Tt,
               group = k, label = label_k(ids[[k]]))
  }))
}))
#Legend inside each panel, above the curves, to the right of the 2003 line
labels_k <- do.call(rbind, lapply(split(groups_k, list(groups_k$series, groups_k$max_iter), drop = TRUE), function(d) {
  r <- diff(range(d$Tt)); g_ <- sort(unique(d$group))
  data.frame(series = d$series[1], max_iter = d$max_iter[1], group = g_, label = d$label[match(g_, d$group)],
             time_yr = year_agriculture + 0.4, Tt = max(d$Tt) + r * 0.14 * rev(seq_along(g_)))
}))

g <- ggplot(groups_k, aes(time_yr, Tt, colour = factor(group))) +
  geom_vline(xintercept = year_agriculture, linetype = "dashed", colour = "grey40") +
  geom_line(data = groups_k[groups_k$group != 1, ], linewidth = 0.8) +
  #main group dashed and drawn on top, so it remains visible where it overlaps another curve
  geom_line(data = groups_k[groups_k$group == 1, ], linewidth = 0.8, linetype = "22") +
  geom_text(data = labels_k, aes(label = label), hjust = 0, size = text_pub(), family = "serif",
            fontface = "bold", show.legend = FALSE) +
  facet_grid(series ~ max_iter, scales = "free_y") +
  scale_colour_manual(values = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  labs(x = "Year", y = "NDVI trend component (Tt)") +
  theme_pub() + theme(legend.position = "none")
save_pub("FigA2", g, width_double, 170)


# Export results -------------------------------------------------------------------------------
write_xlsx(list(validation = validation, bic_series = bic_series, mean_series = res_series,
                breaks = if (is.null(res_breaks)) data.frame() else res_breaks),
           "data/processed/bfast_sensitivity_k.xlsx")
