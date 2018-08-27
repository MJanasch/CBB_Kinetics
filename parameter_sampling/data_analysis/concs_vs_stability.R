#!/usr/bin/env Rscript

# Read infile names from command line
args = commandArgs(trailingOnly=T)
conc_file = args[1] # A metabolite concentration sets file
perc_file = args[2] # A percent stable steady state per concentration set file
meta_file = args[3] # A file with metabolomics concentrations and ratios
# conc_file = "/ssd/common/proj/Kinetic_Model/Metabolite_Sampling/Results/2018-03-07/a/all_metabolite_concentrations.tab"
# perc_file = "/ssd/common/proj/Kinetic_Model/SKM_Sampling/Results/2018-03-10/a/met_set_vs_percent_steady.tab"
# meta_file = "tools/maks/parameter_sampling/data_analysis/concentrations_and_cofactor_ratios_from_2016.tab"

library(reshape2)

# Load data
conc_data = read.table(conc_file, header=T, sep="\t")
perc_data = read.table(perc_file, header=F, sep="\t")
meta_data = read.table(meta_file, header=T, sep="\t")
colnames(perc_data) = c("set", "stable")

# Remove O2 and CO2 since they are fixed
conc_data = conc_data[,grep("O2", colnames(conc_data), invert=T)]
meta_data = meta_data[grep("O2", meta_data$metabolite, invert=T),]

# Calculate ATP/ADP and NADPH/NADP ratios
conc_data$ATP_ADP_ratio = conc_data$ATP / conc_data$ADP
conc_data$NADPH_NADP_ratio = conc_data$NADPH / conc_data$NADP
conc_data$ATP_NADPH_ratio = conc_data$ATP / conc_data$NADPH

# Calculate charge ratios
conc_data$energy_charge = conc_data$ATP / (conc_data$ATP + conc_data$ADP)
conc_data$redox_charge = conc_data$NADPH / (conc_data$NADPH + conc_data$NADP)

meta_charge = dcast(
  meta_data[grep("(^ATP$)|(^ADP$)", meta_data$metabolite, perl=T),],
  Organism + Study ~ metabolite, value.var="concentration"
)
meta_charge$energy_charge = meta_charge$ATP / (meta_charge$ATP + meta_charge$ADP)
meta_charge = melt(meta_charge[,c("Organism","Study","energy_charge")])
colnames(meta_charge) = c("Organism","Study","metabolite","concentration")
meta_charge$KEGG_ID = NA

meta_data = rbind(meta_data, meta_charge)

# Transform the data
conc_data = log10(conc_data)
meta_data$concentration = log10(meta_data$concentration)

# Use Synechocystis and Synechococcus data
meta_data = subset(meta_data, Organism %in% c("PCC 6803", "PCC 7942"))

# Split data into lower and upper stable steady state categories
perc_data = perc_data[order(perc_data$stable, decreasing=T),]

# Divide metabolite concentration sets into high and low % stable steady states

# By Deciles
deciles = quantile(perc_data$stable, prob = seq(0, 1, length = 11), type = 5)
d10 = deciles[10]
d1 = deciles[2]

perc_hi = subset(perc_data, stable > d10)
perc_lo = subset(perc_data, stable <= d1)

# By arbitrary defined cutoffs
# perc_hi = subset(perc_data, stable > 85)
# perc_lo = subset(perc_data, stable < 65)

# Subset concentrations to those in high and low % stable steady states groups
conc_hi = conc_data[perc_hi$set,]
conc_lo = conc_data[perc_lo$set,]

# Plot it
library(ggplot2)

c_lo_long = melt(conc_lo)
colnames(c_lo_long) = c("metabolite", "concentration")
c_lo_long$stability = rep("unstable", nrow(c_lo_long))
# c_lo_long$stability = rep("<65% stable steady states", nrow(c_lo_long))

c_hi_long = melt(conc_hi)
colnames(c_hi_long) = c("metabolite", "concentration")
c_hi_long$stability = rep("stable", nrow(c_hi_long))
# c_hi_long$stability = rep(">85% stable steady states", nrow(c_hi_long))

c_long = rbind(c_hi_long, c_lo_long)

meta_data$stability = "stable"
# meta_data$stability = ">85% stable steady states"


# Plot metabolite concentrations separately

gp = ggplot(
  c_long[grep("(charge)|(ratio)", c_long$metabolite, perl=T, invert=T),],
  aes(x=concentration, fill=stability)
)
gp = gp + geom_density(alpha=0.4)
gp = gp + geom_point(
  aes(x=concentration, colour=Organism, shape=Organism),
  subset(meta_data, metabolite %in% grep("(charge)|(ratio)", c_long$metabolite, perl=T, invert=T, value=T)),
  y=0, show.legend=F, colour="white", alpha=0.7, size=2, shape=19
)
gp = gp + geom_point(
  aes(x=concentration, colour=Organism, shape=Organism),
  subset(meta_data, metabolite %in% grep("(charge)|(ratio)", c_long$metabolite, perl=T, invert=T, value=T)),
  y=0, show.legend=F, fill="black", alpha=0.8, size=1, stroke=0.6
)
gp = gp + scale_shape_manual(values=c(6, 2))
gp = gp + scale_colour_manual(values=c("black", "red"))
gp = gp + theme_bw()
gp = gp + facet_wrap(~metabolite, ncol=6)
gp = gp + scale_fill_manual(values=c("#8073ac","#e08214"))
gp = gp + theme(strip.background = element_blank())

outfile = paste(dirname(perc_file), "concs_vs_stability.with_meta.pdf", sep="/")
ggsave(outfile, gp, width=400/25.4, height=225/25.4)

# Plot ratios separately
n_ratios = length(unique(grep("ratio", c_long$metabolite, value=T)))

gp = ggplot(
  c_long[grep("ratio", c_long$metabolite),],
  aes(x=concentration, fill=stability)
)
gp = gp + geom_density(alpha=0.4)
gp = gp + geom_point(
  aes(x=concentration, colour=Organism, shape=Organism),
  subset(meta_data, metabolite %in% grep("ratio", c_long$metabolite, value=T)),
  y=0, show.legend=F, colour="white", alpha=0.7, size=2, shape=19
)
gp = gp + geom_point(
  aes(x=concentration, colour=Organism, shape=Organism),
  subset(meta_data, metabolite %in% grep("ratio", c_long$metabolite, value=T)),
  y=0, show.legend=F, fill="black", alpha=0.8, size=1, stroke=0.6
)
gp = gp + scale_shape_manual(values=c(6, 2))
gp = gp + scale_colour_manual(values=c("black", "red"))
gp = gp + theme_bw()
gp = gp + facet_wrap(~metabolite, ncol=n_ratios)
gp = gp + scale_fill_manual(values=c("#8073ac","#e08214"))
gp = gp + theme(strip.background = element_blank())

outfile = paste(dirname(perc_file), "concs_vs_stability.ratios.pdf", sep="/")
ggsave(outfile, gp, width=(40+n_ratios*60)/25.4, height=60/25.4)


# Plot energy and redox charge separately
n_charges = length(unique(grep("charge", c_long$metabolite, value=T)))

gp = ggplot(
  c_long[grep("charge", c_long$metabolite),],
  aes(x=concentration, fill=stability)
)
gp = gp + geom_density(alpha=0.4)
gp = gp + geom_point(
  aes(x=concentration, colour=Organism, shape=Organism),
  subset(meta_data, metabolite %in% grep("charge", c_long$metabolite, value=T)),
  y=0, show.legend=F, colour="white", alpha=0.7, size=2, shape=19
)
gp = gp + geom_point(
  aes(x=concentration, colour=Organism, shape=Organism),
  subset(meta_data, metabolite %in% grep("charge", c_long$metabolite, value=T)),
  y=0, show.legend=F, fill="black", alpha=0.8, size=1, stroke=0.6
)
gp = gp + scale_shape_manual(values=c(6, 2))
gp = gp + scale_colour_manual(values=c("black", "red"))
gp = gp + theme_bw()
gp = gp + facet_wrap(~metabolite, ncol=n_charges)
gp = gp + scale_fill_manual(values=c("#8073ac","#e08214"))
gp = gp + theme(strip.background = element_blank())

outfile = paste(dirname(perc_file), "concs_vs_stability.charges.pdf", sep="/")
ggsave(outfile, gp, width=(40+n_charges*60)/25.4, height=60/25.4)
