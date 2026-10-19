#Manuscript: Section 2.2.2 item (iii); Results 3.2 (spatial pattern); Fig. 6 | run after 01 and before 06
#BFAST metrics for all pixels and their maps - Areas A and B ----------------------------------------
#Question: is the NDVI decline in 2000-2006 in Area A spatially concentrated (e.g. in the northern
#sector, with previous agricultural use - Silva, 2018)?
#
#Steps:
#1. Extraction of ALL pixels (same handling of missing values as 01_ndvi_extraction.R)
#2. BFAST on each pixel (default bfast(): harmonic order 3; h = 0.15; max.iter = 1)
#3. Maps of the metrics, with the soil sampling plots overlaid (visual interpretation)

#Control options --------------------------------------------------------------------------------
options(OutDec = ".")
rm(list = ls(all = TRUE))


# Packages -------------------------------------------------------------------------------------
for (p in c("terra", "zoo", "bfast", "readxl", "writexl", "ggplot2", "patchwork")) {
  if (!require(p, character.only = TRUE)) { install.packages(p, dependencies = TRUE); library(p, character.only = TRUE) }
}
source("R/00_publication_theme.R") #figure standard (font and print size)
source("R/00_bfast_functions.R")


# Parameters -----------------------------------------------------------------------------------
h <- 0.15
max_iter <- 1
ndvi_folders <- c(A = "data/raw/ndvi_area_A", B = "data/raw/ndvi_area_B")
plot_files <- c(A = "data/raw/plot_coordinates_area_A.xlsx", B = "data/raw/plot_coordinates_area_B.xlsx")


# 1 and 2. Extraction and BFAST for all pixels -------------------------------------------------
res <- list(); res_breaks <- list(); crs_utm <- NULL
for (a in names(ndvi_folders)) {
  d <- extract_all_pixels(ndvi_folders[a])
  crs_utm <- d$crs
  cat("Area", a, ":", nrow(d$ndvi), "pixels\n")
  t0 <- Sys.time()
  m <- lapply(seq_len(nrow(d$ndvi)), function(i) pixel_metrics(d$ndvi[i, ], h, max_iter, interp = d$interp[i, ], breaks = TRUE))
  ok <- !sapply(m, is.null)
  cat("  BFAST:", sum(ok), "ok |", sum(!ok), "failed |", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
  res[[a]] <- data.frame(area = a, d$xy[ok, ], do.call(rbind, lapply(m[ok], `[[`, "summary")))
  #All trend breakpoints (one row per breakpoint)
  res_breaks[[a]] <- do.call(rbind, lapply(which(ok), function(i)
    if (!is.null(m[[i]]$breaks)) data.frame(area = a, pixel = i, m[[i]]$breaks)))
}
pix <- do.call(rbind, res); rownames(pix) <- NULL
breaks_pix <- do.call(rbind, res_breaks); rownames(breaks_pix) <- NULL


# Soil sampling plots (lon/lat -> UTM 24S), overlaid on the maps -------------------------------
plots <- do.call(rbind, lapply(names(plot_files), function(a) {
  p <- as.data.frame(read_excel(plot_files[a]))
  v <- project(vect(p, geom = c("lon", "lat"), crs = "EPSG:4326"), crs_utm)
  data.frame(area = a, plot = p$plot, geom(v)[, c("x", "y")])
}))


# Summary per area -----------------------------------------------------------------------------
area_summary <- do.call(rbind, lapply(split(pix, pix$area), function(p) data.frame(
  area = p$area[1], pixels = nrow(p),
  pct_slope_2000_2006_negative = mean(p$slope_2000_2006 < 0) * 100,
  median_slope_2000_2006 = median(p$slope_2000_2006),
  pct_decline_2002_2006 = mean(p$decline_2002_2006) * 100,
  pct_decline_2011_2012 = mean(p$decline_2011_2012) * 100,
  mean_trend_breaks = mean(p$n_trend_breaks),
  median_remainder_sd = median(p$remainder_sd))))
cat("\n===== Summary (all pixels) =====\n"); print(area_summary, digits = 3)


# 3. Maps (Fig. 6A: Area A; Fig. 6B: Area B) ---------------------------------------------------
map_theme <- theme_pub() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", #horizontal legend: uses the width of each panel
        legend.key.width = unit(9, "mm"), legend.key.height = unit(2.5, "mm"),
        legend.margin = margin(0, 0, 0, 0))

metric_map <- function(a, v, title, scale) {
  p <- pix[pix$area == a, ]
  ggplot(p, aes(x, y)) +
    geom_raster(aes(fill = .data[[v]])) +
    geom_point(data = plots[plots$area == a, ], shape = 21, fill = "white", colour = "black", size = 1.2) +
    scale + coord_equal() + labs(title = title, x = NULL, y = NULL, fill = NULL) + map_theme
}
lim_slope <- max(abs(pix$slope_2000_2006))
lim_decline <- min(pix$largest_decline)
sc_slope <- scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0, limits = c(-lim_slope, lim_slope))
sc_decline <- scale_fill_gradient(low = "#B2182B", high = "#F7F7F7", limits = c(lim_decline, 0))
sc_remainder <- scale_fill_viridis_c(option = "magma", limits = range(pix$remainder_sd))
sc_ndvi <- scale_fill_gradient(low = "#F6E8C3", high = "#01665E", limits = range(pix$mean_ndvi))
sc_year <- scale_fill_viridis_c(option = "turbo", limits = c(2000, 2023), na.value = "grey85")

for (a in names(ndvi_folders)) {
  g <- (metric_map(a, "slope_2000_2006", "NDVI trend 2000-2006 (yr⁻¹)", sc_slope) |
          metric_map(a, "largest_decline", "Largest abrupt decline", sc_decline) |
          metric_map(a, "year_largest_decline", "Year of largest decline", sc_year)) /
       (metric_map(a, "remainder_sd", "Remainder SD", sc_remainder) |
          metric_map(a, "mean_ndvi", "Mean NDVI 2000-2022", sc_ndvi) |
          metric_map(a, "pct_interpolated", "Interpolated values (%)", scale_fill_viridis_c(limits = c(0, 20)))) +
    plot_annotation(title = paste("Area", a),
                    theme = theme(text = element_text(family = "serif"),
                                  plot.title = element_text(size = font_pub + 2, face = "bold")))
  save_pub(paste0("Fig6", a), g, width_double, 165)
}

#Year of trend breakpoints: % of pixels in each area with a breakpoint in each year (table exported below)
breaks_pix$year_int <- floor(breaks_pix$year)
breaks_pix$type <- ifelse(breaks_pix$magnitude < 0, "Decline", "Increase")
n_pix <- table(pix$area)
breaks_by_year <- aggregate(pixel ~ area + year_int + type, data = breaks_pix, FUN = function(z) length(unique(z)))
breaks_by_year$pct <- 100 * breaks_by_year$pixel / as.numeric(n_pix[breaks_by_year$area])

#Rasters of the metrics (multiband GeoTIFF)
for (a in names(ndvi_folders)) {
  p <- pix[pix$area == a, ]
  p$decline_2002_2006 <- as.numeric(p$decline_2002_2006); p$decline_2011_2012 <- as.numeric(p$decline_2011_2012)
  r <- rast(p[, setdiff(names(p), "area")], type = "xyz", crs = crs_utm)
  writeRaster(r, paste0("data/processed/bfast_pixel_metrics_area_", a, ".tif"), overwrite = TRUE)
}


# Export results -------------------------------------------------------------------------------
write_xlsx(list(pixels = pix, summary = area_summary, breaks = breaks_pix, breaks_by_year = breaks_by_year),
           "data/processed/bfast_pixel_metrics.xlsx")
