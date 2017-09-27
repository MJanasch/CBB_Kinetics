#!/usr/bin/env Rscript

# Read infile names from command line
args = commandArgs(trailingOnly=T)
conc_file = args[1] # A metabolite concentration sets file
perc_file = args[2] # A percent stable steady state per concentration set file

# Load data
conc_data = read.table(conc_file, header=T, sep="\t")
perc_data = read.table(perc_file, header=F, sep="\t")
colnames(perc_data) = c("set", "stable")

# Remove O2 and CO2 since they are fixed
conc_data = conc_data[,grep("O2", colnames(conc_data), invert=T)]

# Calculate ATP/ADP and NADPH/NADP ratios
conc_data$ATP_ADP_ratio = conc_data$ATP / conc_data$ADP
conc_data$NADPH_NADP_ratio = conc_data$NADPH / conc_data$NADP

# Transform the data
conc_data = log10(conc_data)

# Split data into lower and upper stable steady state categories
perc_data = perc_data[order(perc_data$stable, decreasing=T),]

# Calculate deciles
deciles = quantile(perc_data$stable, prob = seq(0, 1, length = 11), type = 5)
d10 = deciles[10]
d1 = deciles[2]
perc_upper = subset(perc_data, stable > d10)
perc_lower = subset(perc_data, stable <= d1)
conc_upper = conc_data[perc_upper$set,]
conc_lower = conc_data[perc_lower$set,]

# Plot it
library(ggplot2)
library(reshape2)

c_lower_long = melt(conc_lower)
colnames(c_lower_long) = c("metabolite", "concentration")
c_lower_long$stability = rep("unstable", nrow(c_lower_long))

c_upper_long = melt(conc_upper)
colnames(c_upper_long) = c("metabolite", "concentration")
c_upper_long$stability = rep("stable", nrow(c_upper_long))

c_long = rbind(c_upper_long, c_lower_long)

gp = ggplot(c_long, aes(x=concentration, fill=stability))
gp = gp + geom_density(alpha=0.4)
gp = gp + theme_bw()
gp = gp + facet_wrap(~metabolite, ncol=6)
gp = gp + scale_fill_manual(values=c("#8073ac","#e08214"))


outfile = paste(dirname(perc_file), "concs_vs_stability.pdf", sep="/")
ggsave(outfile, gp, width=400/25.4, height=225/25.4)
