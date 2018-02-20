#!/usr/bin/env Rscript

### LOAD DATA ##################################################################

# Read infile names from command line
args = commandArgs(trailingOnly=T)
conc_file = args[1] # Metabolite concentration sets file
km_file = args[2] # Parameter sets file
rxn_file = args[3] # Reaction header file

# Testing
#conc_file = "/ssd/common/proj/Kinetic_Model/Metabolite_Sampling/Results/2018-02-13/a/all_metabolite_concentrations.tab"
#km_file = "/ssd/common/proj/Kinetic_Model/SKM_Sampling/Results/2018-02-16/a/concset_stability_parameters.tab.gz"
#rxn_file = "/ssd/common/proj/Kinetic_Model/SKM_Sampling/Results/2018-02-16/a/cbb_reaction_header.long.txt"

# Load data
library(data.table)
conc_data = as.data.frame(fread(conc_file, header=T, sep="\t"))
km_data = as.data.frame(fread(
  paste(c("gzip -dc ", km_file), collapse=""), header=T, sep="\t"
  ))
custom_rxn_labels = scan(rxn_file, character(), quote = "")

# Keep only Ka and Ki values
km_data = km_data[,c(1,2, grep("^Ka", colnames(km_data)), grep("^Ki", colnames(km_data)))]

# Add column with a set name
km_data$Parameter_set = as.numeric(rownames(km_data))

# Turn into long format
library(reshape2)
km_long = melt(km_data, id.vars=c("Conc_set", "Stable", "Parameter_set"))
conc_data$Conc_set = as.character(rownames(conc_data))
conc_long = melt(conc_data, id.vars=c("Conc_set"))

# Rename column names
colnames(conc_long)[(ncol(conc_long)-1):(ncol(conc_long))] = c("Metabolite", "Concentration")

# Add Metabolite and Reaction_n columns
km_long$variable_2 = gsub("K[ai]", "", as.character(km_long$variable), perl=T)
km_long$Metabolite = gsub("v[0-9]*", "", km_long$variable_2)
km_long$Reaction_n = gsub(".*v", "", km_long$variable_2)
km_long$Type = ifelse(grepl("^Ka", km_long$variable), "Ka", "Ki")

# Add Reaction column
km_long$Reaction = custom_rxn_labels[as.numeric(km_long$Reaction_n)]

# Clean up km_long
km_long = km_long[,c("Conc_set", "Parameter_set", "Stable", "value", "Metabolite", "Reaction", "Type")]

# Merge with concentration data
plot_data = merge(km_long, subset(conc_long, Metabolite %in% unique(km_long$Metabolite)))

plot_data$Conc_over_K = 0.001 * plot_data$Concentration / plot_data$value

plot_data$Stable = as.character(plot_data$Stable)

# Plot histograms of Km over concentration
library(ggplot2)

plot_data$Header = paste(
  paste(
    plot_data$Metabolite,
    plot_data$Reaction,
    sep=" / "),
  paste(
    rep("(", nrow(plot_data)),
    plot_data$Type,
    rep(")", nrow(plot_data)),
    sep = ""),
  sep=" ")

S = sample(1:nrow(plot_data), nrow(plot_data)*0.01)

gp = ggplot(plot_data, aes(x=Conc_over_K, fill=Stable))
gp = gp + geom_density(alpha=0.4)
gp = gp + theme_bw()
gp = gp + facet_wrap(~Header, ncol=5)
gp = gp + scale_fill_manual(values=c("#e08214","#8073ac"))
gp = gp + scale_x_log10()

outfile = paste(dirname(km_file), "conc_over_Kai.pdf", sep="/")

ggsave(outfile, gp, height=120/25.4, width=360/25.4)
