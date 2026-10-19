#Manuscript: Section 2.3.3 (hierarchical model, Eq. 5); Results 3.3 (variability among and within plots); Fig. 7
#Variance partitioning of natural radionuclide activity concentrations: among plots x within plots
#Input: data/processed/soil_radionuclides.xlsx
#
#Data structure: 39 plots (21 in Area A, 18 in Area B), 3 subsamples per plot
#
#Model (one-way random-effects model):
#   y_ij = mu + a_i + e_ij,   a_i ~ N(0, s2_among),   e_ij ~ N(0, s2_within)
#where a_i is the effect of plot i and e_ij the variation among subsamples j within plot i.
#A subsample effect nested in plot is not included: each subsample has a single observation, so it
#would not be identifiable from the residual.
#Components are estimated by REML (lme4), separately for each area and radionuclide.

#Control options --------------------------------------------------------------------------------
options(OutDec = ".")
rm(list = ls(all = TRUE))


# Packages -------------------------------------------------------------------------------------
for (p in c("lme4", "readxl", "writexl", "ggplot2")) {
  if (!require(p, character.only = TRUE)) { install.packages(p, dependencies = TRUE); library(p, character.only = TRUE) }
}
source("R/00_publication_theme.R") #figure standard (font and print size)


# Data -----------------------------------------------------------------------------------------
radionuclides <- c(ra226 = "226Ra", k40 = "40K", ra228 = "228Ra", th228 = "228Th")
samples <- as.data.frame(read_excel("data/processed/soil_radionuclides.xlsx", sheet = "samples"))
samples$plot <- factor(samples$plot)

#Check the structure: number of plots and of subsamples per plot
n_sub <- aggregate(subsample ~ area + plot, data = samples, FUN = length)
cat("Plots per area:", tapply(n_sub$plot, n_sub$area, length),
    "| subsamples per plot:", paste(unique(n_sub$subsample), collapse = ", "), "\n")


# Model fitting and variance partitioning ------------------------------------------------------
partition <- do.call(rbind, lapply(c("A", "B"), function(a) {
  d <- droplevels(samples[samples$area == a, ])
  do.call(rbind, lapply(names(radionuclides), function(v) {
    m <- lmer(as.formula(paste(v, "~ 1 + (1 | plot)")), data = d, REML = TRUE)
    vc <- as.data.frame(VarCorr(m))
    s2_among <- vc$vcov[vc$grp == "plot"]
    s2_within <- vc$vcov[vc$grp == "Residual"]
    #95% confidence intervals of the standard deviations by profile likelihood
    ci <- tryCatch(suppressMessages(confint(m, parm = "theta_", method = "profile", signames = FALSE)),
                   error = function(e) matrix(NA, 2, 2))
    data.frame(area = a, radionuclide = v, n_plots = nlevels(d$plot), n_obs = nrow(d),
               mean = unname(fixef(m)),
               var_among_plots = s2_among, var_within_plots = s2_within,
               pct_among_plots = 100 * s2_among / (s2_among + s2_within),
               pct_within_plots = 100 * s2_within / (s2_among + s2_within),
               sd_among_ci95_lower = ci[1, 1], sd_among_ci95_upper = ci[1, 2],
               sd_within_ci95_lower = ci[2, 1], sd_within_ci95_upper = ci[2, 2],
               singular = isSingular(m))
  }))
}))
cat("\n===== Variance partitioning (REML; among x within plots) =====\n")
print(partition[, c("area", "radionuclide", "mean", "var_among_plots", "var_within_plots",
                    "pct_among_plots", "pct_within_plots", "singular")], digits = 3)


# Figure (Fig. 7) ------------------------------------------------------------------------------
long <- rbind(
  data.frame(partition[, c("area", "radionuclide")], component = "Among plots", pct = partition$pct_among_plots),
  data.frame(partition[, c("area", "radionuclide")], component = "Within plots", pct = partition$pct_within_plots))
long$radionuclide <- factor(long$radionuclide, levels = rev(names(radionuclides)))
long$component <- factor(long$component, levels = c("Within plots", "Among plots"))
long$area <- paste("Area", long$area)
radionuclide_labels <- c(ra226 = expression({}^{226}*Ra), k40 = expression({}^{40}*K),
                         ra228 = expression({}^{228}*Ra), th228 = expression({}^{228}*Th))

g <- ggplot(long, aes(pct, radionuclide, fill = component)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  #label colour mapped to the component (white on green, black on yellow)
  geom_text(aes(label = sprintf("%.0f", pct), colour = component), position = position_stack(vjust = 0.5),
            size = text_pub(), family = "serif", show.legend = FALSE) +
  facet_wrap(~ area, nrow = 1) +
  scale_y_discrete(labels = radionuclide_labels) +
  scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100.5), expand = c(0, 0)) +
  scale_fill_manual(values = c("Among plots" = "#21918c", "Within plots" = "#fde725"),
                    breaks = c("Among plots", "Within plots")) +
  scale_colour_manual(values = c("Among plots" = "white", "Within plots" = "black"), guide = "none") +
  labs(x = "Variance explained (%)", y = NULL, fill = NULL) +
  theme_pub() + theme(legend.position = "bottom", panel.grid = element_blank(),
                      panel.spacing = unit(6, "mm"), plot.margin = margin(4, 8, 4, 4))
save_pub("Fig7", g, width_double, 85)


# Export results -------------------------------------------------------------------------------
write_xlsx(list(partition = partition), "data/processed/radionuclide_variance_partitioning.xlsx")
