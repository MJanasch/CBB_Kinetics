#!/usr/bin/env Rscript

### LOAD DATA ##################################################################

# Read infile names from command line
args = commandArgs(trailingOnly=T)
conc_file = args[1] # Metabolite concentration sets file
km_file = args[2] # Parameter sets file
rxn_file = args[3] # Reaction header file

# Load data
library(data.table)
conc_data = as.data.frame(fread(conc_file, header=T, sep="\t"))
km_data = as.data.frame(fread(
  paste(c("gzip -dc ", km_file), collapse=""), header=T, sep="\t"
  ))
custom_rxn_labels = scan(rxn_file, character(), quote = "")

# Keep only Km values
km_data = km_data[,c(1,2, grep("^Km", colnames(km_data)), grep("K_PPool", colnames(km_data)))]

# Set up Km names translation
#Km_names = as.character(colnames(km_data)[3:ncol(km_data)])
#
#Km_names_translation = data.frame(variable=grep(".", Km_names, value=T, fixed=T))

#Km_names_translation$reaction = unlist(
#  lapply(
#    strsplit(as.character(Km_names_translation$variable), "v"),
#    "[[", 2
#    )
#  )

#reaction_original = data.frame(
#  reaction=c("6.1","5.1","8.1","18.1"),
#  reaction2=c("6","5","8","18")
#  )

#reaction_additional = data.frame(
#  reaction=c("6.1","5.1","8.1","18.1"),
#  reaction2=c("7","9","10","19")
#  )

#Km_names_translation_1 = merge(Km_names_translation, reaction_original)
#Km_names_translation_2 = merge(Km_names_translation, reaction_additional)

#Km_names_translation_1$proper_name = unlist(lapply(strsplit(as.character(Km_names_translation_1$variable), ".", fixed=T), "[[", 1))
#Km_names_translation_2$proper_name = paste(unlist(lapply(strsplit(as.character(Km_names_translation_1$variable), "v"), "[[", 1)), Km_names_translation_2$reaction2, sep="v")

# Promiscuous reactions only use the second Km value
# Therefore, replace values in first columns
#for (i in 1:nrow(Km_names_translation_1)){
#  data_name = as.character(Km_names_translation_1[i,"variable"])
#  proper_name = as.character(Km_names_translation_1[i,"proper_name"])
#  km_data[,proper_name] = km_data[,data_name]
#}

# Then rename second reactions
#for (i in 1:nrow(Km_names_translation_2)){
#  bad_name = Km_names_translation_2[i,"variable"]
#  proper_name = Km_names_translation_2[i,"proper_name"]
#  colnames(km_data)[grep(bad_name, colnames(km_data))] = proper_name
#}

# Add column with a set name
km_data$Parameter_set = as.numeric(rownames(km_data))

# Turn into long format
library(reshape2)
km_long = melt(km_data, id.vars=c("Conc_set", "Stable", "Parameter_set"))
conc_data$Conc_set = as.character(rownames(conc_data))
conc_long = melt(conc_data, id.vars=c("Conc_set"))

# Rename column names
colnames(km_long)[(ncol(km_long)-1):(ncol(km_long))] = c("km_name", "K_m")
colnames(conc_long)[(ncol(conc_long)-1):(ncol(conc_long))] = c("Metabolite", "Concentration")

# Add Metabolite and Reaction_n columns
km_long$km_name_2 = gsub("Km", "", as.character(km_long$km_name), fixed=T)
km_long$Metabolite = gsub("v[0-9]*", "", km_long$km_name_2)
km_long$Reaction_n = gsub(".*v", "", km_long$km_name_2)

# Add Reaction column
km_long$Reaction = custom_rxn_labels[as.numeric(km_long$Reaction_n)]

# Clean up km_long
km_long = km_long[,c("Conc_set", "Parameter_set", "Stable", "K_m", "Metabolite", "Reaction")]

# Merge with concentration data
library(dplyr)

# Change variable type for inner_join
conc_long$Conc_set = as.numeric(conc_long$Conc_set)
conc_long$Metabolite = as.character(conc_long$Metabolite)

# Edit PPool metabolite
km_long$Metabolite[km_long$Metabolite == "K_PPool"] = "PPool"

plot_data = inner_join(km_long, conc_long)

plot_data$Conc_over_Km = 0.001 * plot_data$Concentration / plot_data$K_m

plot_data$Stable = as.character(plot_data$Stable)

plot_data$Header = paste(plot_data$Metabolite, plot_data$Reaction, sep=" / ")

# Split off PPool data
plot_data_ppool = subset(plot_data, Metabolite == "PPool")
plot_data = subset(plot_data, Metabolite != "PPool")

# Plot distributions of Km over concentration
library(ggplot2)

gp = ggplot(plot_data_ppool, aes(x=Conc_over_Km, fill=Stable))
gp = gp + geom_density(alpha=0.4)
gp = gp + theme_bw()
gp = gp + facet_wrap(~Header, ncol=10)
gp = gp + scale_fill_manual(values=c("#e08214","#8073ac"))
gp = gp + scale_x_log10()

outfile = paste(dirname(km_file), "conc_over_K_PPool.pdf", sep="/")

ggsave(outfile, gp, height=60/25.4, width=80/25.4)

gp = ggplot(plot_data, aes(x=Conc_over_Km, fill=Stable))
gp = gp + geom_density(alpha=0.4)
gp = gp + theme_bw()
gp = gp + facet_wrap(~Header, ncol=10)
gp = gp + scale_fill_manual(values=c("#e08214","#8073ac"))
gp = gp + scale_x_log10()

outfile = paste(dirname(km_file), "conc_over_Km.pdf", sep="/")

ggsave(outfile, gp, height=480/25.4, width=640/25.4)
