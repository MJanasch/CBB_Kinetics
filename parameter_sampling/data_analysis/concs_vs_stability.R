#!/usr/bin/env Rscript

# Read infile names from command line
args = commandArgs(trailingOnly=T)
conc_file = args[1] # A metabolite concentration sets file
perc_file = args[2] # A percent stable steady state per concentration set file
meta_file = args[3] # A file with metabolomics concentrations and ratios

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

# Transform the data
conc_data = log10(conc_data)
meta_data$concentration = log10(meta_data$concentration)

# Use only Synechocystis 6803 data
meta_data = subset(meta_data, Organism == "PCC 6803")

# Split data into lower and upper stable steady state categories
perc_data = perc_data[order(perc_data$stable, decreasing=T),]

# Divide metabolite concentration sets into high and low % stable steady states

# By Deciles
# deciles = quantile(perc_data$stable, prob = seq(0, 1, length = 11), type = 5)
# d10 = deciles[10]
# d1 = deciles[2]
#
# perc_hi = subset(perc_data, stable > d10)
# perc_lo = subset(perc_data, stable <= d1)

# By arbitrary defined cutoffs
perc_hi = subset(perc_data, stable > 85)
perc_lo = subset(perc_data, stable < 65)

# Subset concentrations to those in high and low % stable steady states groups
conc_hi = conc_data[perc_hi$set,]
conc_lo = conc_data[perc_lo$set,]

# Plot it
library(ggplot2)
library(reshape2)

c_lo_long = melt(conc_lo)
colnames(c_lo_long) = c("metabolite", "concentration")
# c_lo_long$stability = rep("unstable", nrow(c_lo_long))
c_lo_long$stability = rep("<65% stable steady states", nrow(c_lo_long))

c_hi_long = melt(conc_hi)
colnames(c_hi_long) = c("metabolite", "concentration")
# c_hi_long$stability = rep("stable", nrow(c_hi_long))
c_hi_long$stability = rep(">85% stable steady states", nrow(c_hi_long))

c_long = rbind(c_hi_long, c_lo_long)

meta_data$stability = ">85% stable steady states"

gp = ggplot(c_long, aes(x=concentration, fill=stability))
gp = gp + geom_density(alpha=0.4)
gp = gp + geom_point(
  aes(x=concentration),
  subset(meta_data, metabolite %in% c_long$metabolite),
  y=0, show.legend=F, shape=25, fill="black", alpha=0.5
)
gp = gp + theme_bw()
gp = gp + facet_wrap(~metabolite, ncol=6)
gp = gp + scale_fill_manual(values=rev(c("#8073ac","#e08214")))
gp = gp + theme(strip.background = element_blank())

# outfile = "art/2018-05-09/concs_vs_stability.cutoff_test.pdf"
# ggsave(outfile, gp, width=400/25.4, height=225/25.4)

outfile = paste(dirname(perc_file), "concs_vs_stability.with_meta.pdf", sep="/")
ggsave(outfile, gp, width=400/25.4, height=225/25.4)
