#!/usr/bin/env Rscript

# Read infile names from command line
args = commandArgs(trailingOnly=T)
fccs_file = args[1] # An FCCs tab.gz file

# Load data
fccs_data = read.table(gzfile(fccs_file), header=T, sep="\t")

# The enzyme in the column influences the flux of reactions in rows (FCC)

# Load custom labels for reactions
custom_rxn_labels = scan("/tmp/cbb_reaction_header.long.txt", character(), quote = "")

# Making FCCs that are close to zero become zero
fccs = fccs_data

fccs[,4:ncol(fccs)][abs(fccs[,4:ncol(fccs)]) < 1e-6] = 0 # Markus approves

# Plot the influence distributions
library(ggplot2)
library(reshape2)

fccs$Reaction = custom_rxn_labels[fccs$Reaction]
fccs_long = melt(fccs, id.vars=colnames(fccs[1:3]))

fccs_long$Value = abs(fccs_long$value)
fccs_long$Sign = ifelse(fccs_long$value > 0, "Positive", ifelse(fccs_long$value < 0, "Negative", "Zero"))

# Calculate Median and Median Absolute Deviation (MAD)
fccs_median = aggregate(value ~ variable + Reaction, fccs_long, median)
colnames(fccs_median) = c("Effector", "Target", "Median_FCC")
fccs_mad = aggregate(value ~ variable + Reaction, fccs_long, mad)
colnames(fccs_mad) = c("Effector", "Target", "MAD")

# Merge for plotting
fccs_med_mad = merge(fccs_median, fccs_mad)

fccs_med_mad$Effector = factor(as.character(fccs_med_mad$Effector), levels=levels(fccs_med_mad$Effector))
fccs_med_mad$Target = factor(as.character(fccs_med_mad$Target), levels=rev(levels(fccs_med_mad$Effector)))

library(scales)

gp = ggplot(fccs_med_mad, aes(x=Target, y=Effector, alpha=MAD, fill=Median_FCC))
gp = gp + geom_tile(colour="white", size=0.5)
gp = gp + theme_bw()
gp = gp + scale_fill_gradientn(colours=c("#fdae61","#ffffbf","#2c7bb6"), values=rescale(c(min(fccs_med_mad$Median_FCC),0,max(fccs_med_mad$Median_FCC)), c(0,1)))
gp = gp + scale_alpha_continuous(range=c(1,0.1))
gp = gp + theme(axis.text.x = element_text(angle = 90, hjust = 0))
gp = gp + scale_x_discrete(position="top")
gp = gp + theme(aspect.ratio=1)

outfile = paste(dirname(fccs_file), "FCCs_heatmap.pdf", sep="/")
ggsave(outfile, gp, height=150/25.4, width=180/25.4)
