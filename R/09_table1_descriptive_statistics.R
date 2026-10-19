#Manuscript: Results 3.3; Table 1
#Descriptive statistics of natural radionuclide activity concentrations (Bq kg-1) by area
#Input: data/processed/soil_radionuclides.xlsx | Output: tables/Table1.docx and tables/Table1.xlsx

#Control options --------------------------------------------------------------------------------
options(OutDec = ".")
rm(list = ls(all = TRUE))


# Packages -------------------------------------------------------------------------------------
for (p in c("readxl", "writexl", "officer", "flextable")) {
  if (!require(p, character.only = TRUE)) { install.packages(p, dependencies = TRUE); library(p, character.only = TRUE) }
}


# Data -----------------------------------------------------------------------------------------
samples <- as.data.frame(read_excel("data/processed/soil_radionuclides.xlsx", sheet = "samples"))
radionuclides <- c(k40 = "40K", ra226 = "226Ra", ra228 = "228Ra", th228 = "228Th")

tab <- do.call(rbind, lapply(c("A", "B"), function(a) do.call(rbind, lapply(names(radionuclides), function(v) {
  y <- samples[samples$area == a, v]
  data.frame(Area = paste("Area", a), Radionuclide = radionuclides[[v]], n = length(y),
             Mean = mean(y), SD = sd(y), CV = 100 * sd(y) / mean(y),
             Min = min(y), Median = median(y), Max = max(y))
}))))
print(tab, digits = 4)


# Word table -----------------------------------------------------------------------------------
dir.create("tables", showWarnings = FALSE)
ft <- flextable(tab)
ft <- colformat_double(ft, j = c("Mean", "SD", "CV", "Min", "Median", "Max"), digits = 1)
ft <- merge_v(ft, j = "Area")
#Isotopes with the mass number as superscript
for (i in seq_len(nrow(tab))) {
  mass <- sub("^([0-9]+).*", "\\1", tab$Radionuclide[i]); el <- sub("^[0-9]+", "", tab$Radionuclide[i])
  ft <- compose(ft, i = i, j = "Radionuclide", value = as_paragraph(as_sup(mass), el))
}
ft <- set_header_labels(ft, CV = "CV (%)")
ft <- add_header_lines(ft, values = "Table 1. Descriptive statistics of natural radionuclide activity concentrations (Bq kg⁻¹) in soil samples (0–20 cm) from the study patches of Caatinga vegetation on Assunção Island, Cabrobó (PE, Brazil). SD: standard deviation; CV: coefficient of variation.")
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- fontsize(ft, size = 10, part = "all")
ft <- bold(ft, i = 2, part = "header")
ft <- align(ft, j = 3:9, align = "center", part = "all")
ft <- align(ft, i = 1, align = "left", part = "header")
ft <- border_remove(ft)
ft <- hline_top(ft, part = "body", border = officer::fp_border(width = 1))
ft <- hline(ft, i = 1, part = "header", border = officer::fp_border(width = 1))
ft <- hline_bottom(ft, part = "body", border = officer::fp_border(width = 1))
ft <- hline(ft, i = 4, part = "body", border = officer::fp_border(width = 0.5, color = "grey50"))
ft <- autofit(ft)

doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = "tables/Table1.docx")
write_xlsx(tab, "tables/Table1.xlsx")
cat("Created: tables/Table1.docx and tables/Table1.xlsx\n")
