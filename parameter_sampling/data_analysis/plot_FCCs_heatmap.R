#!/usr/bin/env Rscript

# Read infile names from command line
args = commandArgs(trailingOnly=T)
fccs_file = args[1] # An FCCs tab.gz file
header_file = args[2] # Reaction header file

# Load data
library(data.table)
fccs_data = as.data.frame(fread(paste(c("gzip -dc ", fccs_file), collapse=""),
  header=T, sep="\t"
  ))

# The enzyme in the column influences the flux of reactions in rows (FCC)

# Load custom labels for reactions
custom_rxn_labels = scan(header_file, character(), quote = "")

# Making FCCs that are close to zero become zero
fccs = fccs_data

fccs[,4:ncol(fccs)][abs(fccs[,4:ncol(fccs)]) < 1e-6] = 0 # Markus approves

# Plot the influence distributions
library(ggplot2)
library(reshape2)

fccs$Reaction = custom_rxn_labels[fccs$Reaction]

# Calculate Median and Median Absolute Deviation (MAD)
fccs_med_mad_list = lapply(colnames(fccs)[4:ncol(fccs)], function(effector){
  # Calculate Median and MAD
  effector_median = aggregate(fccs[,effector], list(Target = fccs[,"Reaction"]), median)
  effector_mad = aggregate(fccs[,effector], list(Target = fccs[,"Reaction"]), mad)
  # The second column is the Median or the MAD
  colnames(effector_median)[2] = "Median_FCC"
  colnames(effector_mad)[2] = "MAD"
  # Add the Effector
  effector_median$Effector = effector
  effector_mad$Effector = effector
  #
  # Return the merged data frame
  merge(effector_median, effector_mad)
  })

fccs_med_mad = as.data.frame(rbindlist(fccs_med_mad_list))

# Specify order of Effectors and Targets for plotting
fccs_med_mad$Effector = factor(fccs_med_mad$Effector, levels=custom_rxn_labels)
fccs_med_mad$Target = factor(fccs_med_mad$Target, levels=rev(custom_rxn_labels))

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
