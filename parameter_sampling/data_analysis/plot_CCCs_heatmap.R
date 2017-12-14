#!/usr/bin/env Rscript

# Read infile names from command line
args = commandArgs(trailingOnly=T)
cccs_file = args[1] # An CCCs tab.gz file, MOD
met_list_file = args[2] # Met list file

# Load data
cccs_data = read.table(gzfile(cccs_file), header=T, sep="\t") # MOD

# The enzyme in the column influences the metabolite concentration in rows (CCC)

# Load custom labels for reactions
custom_met_labels = scan(met_list_file, character(), quote = "")

# Making CCCs that are close to zero become zero
cccs = cccs_data

cccs[,4:ncol(cccs)][abs(cccs[,4:ncol(cccs)]) < 1e-6] = 0 # Markus approves (I did? Okay :D)

# Plot the influence distributions
library(ggplot2)
library(reshape2)

cccs$Metabolite = custom_met_labels[cccs$Metabolite]
cccs_long = melt(cccs, id.vars=colnames(cccs[1:3]))

cccs_long$Value = abs(cccs_long$value)
cccs_long$Sign = ifelse(cccs_long$value > 0, "Positive", ifelse(cccs_long$value < 0, "Negative", "Zero"))

# Calculate Median and Median Absolute Deviation (MAD)
cccs_median = aggregate(value ~ variable + Metabolite, cccs_long, median)
colnames(cccs_median) = c("Effector", "Target", "Median_CCC")
cccs_mad = aggregate(value ~ variable + Metabolite, cccs_long, mad)
colnames(cccs_mad) = c("Effector", "Target", "MAD")

# Merge for plotting
cccs_med_mad = merge(cccs_median, cccs_mad)

cccs_med_mad$Effector = factor(as.character(cccs_med_mad$Effector), levels=levels(cccs_med_mad$Effector))
cccs_med_mad$Target = factor(as.character(cccs_med_mad$Target), levels=rev(levels(cccs_med_mad$Effector)))

library(scales)

gp = ggplot(cccs_med_mad, aes(x=Target, y=Effector, alpha=MAD, fill=Median_CCC))
gp = gp + geom_tile(colour="white", size=0.5)
gp = gp + theme_bw()
gp = gp + scale_fill_gradientn(colours=c("#fdae61","#ffffbf","#2c7bb6"), values=rescale(c(min(cccs_med_mad$Median_CCC),0,max(cccs_med_mad$Median_CCC)), c(0,1)))
gp = gp + scale_alpha_continuous(range=c(1,0.1))
gp = gp + theme(axis.text.x = element_text(angle = 90, hjust = 0))
gp = gp + scale_x_discrete(position="top")
gp = gp + theme(aspect.ratio=1)

outfile = paste(dirname(cccs_file), "CCCs_heatmap.pdf", sep="/")
ggsave(outfile, gp, height=150/25.4, width=180/25.4)
