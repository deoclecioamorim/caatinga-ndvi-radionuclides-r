#Manuscript: Appendix A - Supplementary data (Table A.1 and Figs. A.1-A.5) | run after scripts 02 to 08
#Builds the Word document of the supplementary material from the final figures in figures/

#Control options --------------------------------------------------------------------------------
options(OutDec = ".")
rm(list = ls(all = TRUE))


# Packages -------------------------------------------------------------------------------------
for (p in c("officer", "flextable", "png")) {
  if (!require(p, character.only = TRUE)) { install.packages(p, dependencies = TRUE); library(p, character.only = TRUE) }
}

max_width <- 16 / 2.54 #maximum figure width on the page (inches; 16 cm)
add_figure <- function(doc, file, caption) {
  dims <- dim(png::readPNG(file, native = TRUE)) #height x width (px)
  w <- max_width; h <- max_width * dims[1] / dims[2]
  doc <- body_add_img(doc, src = file, width = w, height = h, style = "centered")
  doc <- body_add_fpar(doc, fpar(ftext(caption$num, fp_text(bold = TRUE, font.family = "Times New Roman", font.size = 11)),
                                 ftext(paste0(" ", caption$text), fp_text(font.family = "Times New Roman", font.size = 11))))
  body_add_par(doc, "")
}
txt <- function(s, bold = FALSE, size = 11) ftext(s, fp_text(bold = bold, font.family = "Times New Roman", font.size = size))

doc <- read_docx()
doc <- body_add_fpar(doc, fpar(txt("Supplementary material", bold = TRUE, size = 14)))
doc <- body_add_fpar(doc, fpar(txt("Integrating remote sensing and natural radionuclides to assess environmental variability in a Brazilian desertification nucleus", size = 12)))
doc <- body_add_par(doc, "")


# Table A.1 ------------------------------------------------------------------------------------
#Confusion matrices of the LULC classification (2018 image), from "Appendix A- Supplementary data.docx"
#Rows: classified data; columns: reference data. In the original document, the overall accuracy (88.9% and 86.7%)
#was in the "Kappa Index" column; here overall accuracy and kappa are computed from the matrices.
conf_matrices <- list(A = matrix(c(29, 1, 0, 2, 28, 0, 0, 7, 23), 3, byrow = TRUE),
                      B = matrix(c(25, 5, 0, 3, 23, 4, 0, 0, 30), 3, byrow = TRUE))
classes <- c("C1", "C2", "C3")
tab_rows <- do.call(rbind, lapply(names(conf_matrices), function(a) {
  x <- conf_matrices[[a]]; n <- sum(x); po <- sum(diag(x)) / n; pe <- sum(rowSums(x) * colSums(x)) / n^2
  f <- function(v) formatC(v, format = "f", digits = 1)
  rbind(data.frame(Area = paste("Area", a), Class = classes, C1 = x[, 1], C2 = x[, 2], C3 = x[, 3], Total = rowSums(x),
                   UA = f(100 * diag(x) / rowSums(x))),
        data.frame(Area = paste("Area", a), Class = "Total", C1 = colSums(x)[1], C2 = colSums(x)[2], C3 = colSums(x)[3], Total = n, UA = ""),
        data.frame(Area = paste("Area", a), Class = "Producer’s accuracy (%)", C1 = f(100 * x[1, 1] / sum(x[, 1])),
                   C2 = f(100 * x[2, 2] / sum(x[, 2])), C3 = f(100 * x[3, 3] / sum(x[, 3])), Total = "", UA = ""),
        data.frame(Area = paste("Area", a), Class = "Overall accuracy (%)", C1 = f(100 * po), C2 = "", C3 = "", Total = "", UA = ""),
        data.frame(Area = paste("Area", a), Class = "Kappa", C1 = formatC((po - pe) / (1 - pe), format = "f", digits = 2), C2 = "", C3 = "", Total = "", UA = ""))
}))
tab_rows[] <- lapply(tab_rows, as.character)
ft1 <- flextable::flextable(tab_rows)
ft1 <- flextable::set_header_labels(ft1, UA = "User’s accuracy (%)")
ft1 <- flextable::merge_v(ft1, j = "Area")
ft1 <- flextable::font(ft1, fontname = "Times New Roman", part = "all")
ft1 <- flextable::fontsize(ft1, size = 10, part = "all")
ft1 <- flextable::bold(ft1, part = "header")
ft1 <- flextable::align(ft1, j = 3:7, align = "center", part = "all")
ft1 <- flextable::border_remove(ft1)
ft1 <- flextable::hline_top(ft1, part = "header", border = fp_border(width = 1))
ft1 <- flextable::hline_bottom(ft1, part = "header", border = fp_border(width = 1))
ft1 <- flextable::hline(ft1, i = 7, part = "body", border = fp_border(width = 0.5, color = "grey50"))
ft1 <- flextable::hline_bottom(ft1, part = "body", border = fp_border(width = 1))
ft1 <- flextable::autofit(ft1)
doc <- body_add_fpar(doc, fpar(txt("Table A.1. ", bold = TRUE), txt("Confusion matrices and accuracy metrics of the LULC classification for Areas A and B (2018 image). Rows: classified data; columns: reference data. C1: Dense Caatinga; C2: Sparse Caatinga; C3: Exposed soil.")))
doc <- flextable::body_add_flextable(doc, ft1)
doc <- body_add_par(doc, "")


# Figures A.1-A.5 ------------------------------------------------------------------------------
captions <- list(
  list(file = "figures/FigA1.png", num = "Fig. A.1.", text = "Sensitivity of the BFAST trend component of the mean bimonthly NDVI series (Area A, Area B and difference A − B) to the minimum segment size parameter h (0.05–0.45). Values of h producing identical trends are shown as a single line. Dashed vertical line: agricultural land use detected in 2003."),
  list(file = "figures/FigA2.png", num = "Fig. A.2.", text = "Sensitivity of the BFAST trend component of the mean bimonthly NDVI series to the order of the harmonic seasonal model (K = 1–3), with one and ten iterations (h = 0.15). Curves differing by less than 0.005 NDVI are shown as a single line."),
  list(file = "figures/FigA3.png", num = "Fig. A.3.", text = "Spatial correlograms (Moran’s I by distance class) of pixel-level BFAST metrics in Areas A and B. Dashed vertical line: minimum distance (90 m) between pixels used in the spatially thinned samples."),
  list(file = "figures/FigA4.png", num = "Fig. A.4.", text = "Distribution of Cliff’s delta (Area A vs. Area B) for pixel-level BFAST metrics across 1,000 spatially thinned samples (minimum distance of 90 m between pixels). Grey band: negligible effect (|δ| < 0.147)."),
  list(file = "figures/FigA5.png", num = "Fig. A.5.", text = "SHAP values of the natural radionuclides in the Random Forest model, expressed as contributions to the predicted probability of Area B. Colours indicate the activity concentration of each radionuclide (low to high)."))
for (cp in captions) doc <- add_figure(doc, cp$file, list(num = cp$num, text = cp$text))

print(doc, target = "Supplementary_material.docx")
cat("Created: Supplementary_material.docx\n")
