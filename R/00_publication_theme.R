#Manuscript: graphical standard of all figures | used by scripts 02 to 08 and 10 (requires the packages ggplot2, png and ragg)
#Figure standard for publication ------------------------------------------------------------------
#Figures are saved at their final print size, so that the font is not reduced by the journal:
#- width: 90 mm (single column), 140 mm (1.5 column) or 190 mm (full page) - Elsevier standard
#- text: 8-10 pt at print size (recommended minimum: 7 pt)
#- resolution: 600 dpi
#Usage: source("R/00_publication_theme.R") after loading ggplot2

font_pub <- 9 #base font size (pt)
width_single <- 90; width_mid <- 140; width_double <- 190 #widths (mm)
dpi_pub <- 600

#ggplot2 theme
theme_pub <- function(base = font_pub) {
  ggplot2::theme_bw(base_size = base, base_family = "serif") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base - 1, colour = "black"),
                   axis.title = ggplot2::element_text(size = base),
                   strip.text = ggplot2::element_text(size = base),
                   strip.background = ggplot2::element_rect(fill = "grey92"),
                   legend.text = ggplot2::element_text(size = base - 1),
                   legend.title = ggplot2::element_text(size = base),
                   plot.title = ggplot2::element_text(size = base + 1),
                   plot.caption = ggplot2::element_text(size = base - 1.5))
}

#Text size for geom_text (in mm) equivalent to 'pt' points
text_pub <- function(pt = font_pub - 1) pt / ggplot2::.pt

#Figures are saved in figures/ with the manuscript numbering (e.g. "Fig7"), in PNG and TIFF
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

#Copy of a PNG as TIFF (LZW compression; TIFF is the format preferred by Elsevier), with the same pixels
png_to_tiff <- function(png_file) {
  img <- png::readPNG(png_file)
  ragg::agg_tiff(sub("\\.png$", ".tif", png_file), width = dim(img)[2], height = dim(img)[1],
                 units = "px", res = dpi_pub, compression = "lzw", background = "white")
  grid::grid.raster(img, interpolate = FALSE)
  invisible(grDevices::dev.off())
}

#Save a ggplot2 figure at print size (mm): figures/<name>.png and .tif
save_pub <- function(name, g, width = width_double, height) {
  file <- file.path(fig_dir, paste0(name, ".png"))
  ggplot2::ggsave(file, g, width = width, height = height, units = "mm", dpi = dpi_pub)
  png_to_tiff(file)
}

#Open a png device for base R graphics (e.g. plot of bfast objects); close it with close_pub(name)
#Note: with 3 or more stacked panels (par(mfrow)), R reduces text to 66%; in these cases use
#pointsize = 13, which gives ~8.5 pt at print size
png_pub <- function(name, width = width_double, height, pointsize = font_pub) {
  png(filename = file.path(fig_dir, paste0(name, ".png")), width = width, height = height, units = "mm",
      res = dpi_pub, pointsize = pointsize, family = "serif")
}
close_pub <- function(name) {
  invisible(grDevices::dev.off())
  png_to_tiff(file.path(fig_dir, paste0(name, ".png")))
}
