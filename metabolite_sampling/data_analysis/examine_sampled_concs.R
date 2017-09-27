#!/usr/bin/env Rscript

# Load command line arguments
args = commandArgs(trailingOnly=T)
infile = args[1] # A metabolite concentration sets file

# Load data
conc_data = read.table(infile, header=T, sep="\t")

# Remove O2 and CO2 since they are fixed
conc_data = conc_data[,grep("O2", colnames(conc_data), invert=T)]

# Log-transform and scale the data
conc_data = scale(log(conc_data))

# Perform PCA
pcaA = princomp(conc_data)
pcaB = prcomp(conc_data)

# Create plotting dataframes
plot_sets = as.data.frame(pcaB$x)
plot_cpds = as.data.frame(pcaB$rotation)
plot_cpds$Compound = rownames(plot_cpds)

# Plot it
library(ggplot2)
library(ggrepel)

gp = ggplot(plot_sets, aes(x=PC1, y=PC2))
gp = gp + geom_point(size=0.1)
gp = gp + theme_bw()

gp = ggplot(plot_cpds, aes(x=PC1, y=PC2, colour=PC3, label=Compound))
gp = gp + geom_point(aes(size=PC3))
gp = gp + scale_size(range=c(1,3))
gp = gp + scale_colour_gradient2(low="#d0d1e6", mid="#3690c0", high="#014636")
gp = gp + geom_text_repel(force=3, size=4)
gp = gp + theme_bw()

outfile = paste(sub("\\.tab$", "", infile), ".pca.pdf", sep="")

ggsave(outfile, gp, height=15/2.54, width=18/2.54)

# Without ATP, ADP, NADP and NADPH
conc_data_red = conc_data[,grep("ATP|AD", colnames(conc_data), perl=T, invert=T)]

# Perform PCA
pcaA = princomp(conc_data_red)
pcaB = prcomp(conc_data_red)

# Create plotting dataframes
plot_sets = as.data.frame(pcaB$x)
plot_cpds = as.data.frame(pcaB$rotation)
plot_cpds$Compound = rownames(plot_cpds)

gp = ggplot(plot_cpds, aes(x=PC1, y=PC2, colour=PC3, label=Compound))
gp = gp + geom_point(aes(size=PC3))
gp = gp + scale_size(range=c(1,3))
gp = gp + scale_colour_gradient2(low="#d0d1e6", mid="#3690c0", high="#014636")
gp = gp + geom_text_repel(force=3, size=4)
gp = gp + theme_bw()

outfile = paste(sub("\\.tab$", "", infile), ".pca_nocofactor.pdf", sep="")

ggsave(outfile, gp, height=15/2.54, width=18/2.54)


################################################################################
# Plot all vs all concentrations

# Load data
conc_data = read.table(infile, header=T, sep="\t")

# Remove O2 and CO2 since they are fixed
conc_data = conc_data[,grep("O2", colnames(conc_data), invert=T)]

# Transform the data
conc_data = log10(conc_data)

# Plot it
library(GGally)
library(ggplot2)

outfile = paste(sub("\\.tab$", "", infile), ".all_vs_all.png", sep="")

png(outfile, height = 3000, width = 3000)
g <- ggpairs(
  conc_data,
  lower = list(continuous = wrap("points", alpha = 0.3, size=0.1), combo = wrap("dot", alpha = 0.4, size=0.2) )
)
print(g)
dev.off()
