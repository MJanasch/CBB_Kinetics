#!/usr/bin/env Rscript

# Read infile names from command line
args = commandArgs(trailingOnly=T)
conc_file = args[1] # A metabolite concentration sets file
perc_file = args[2] # A percent stable steady state per concentration set file
stdg_file = args[3] # A table with reactions and standard delta G's

# Load data
conc_data = read.table(conc_file, header=T, sep="\t")
perc_data = read.table(perc_file, header=F, sep="\t")
colnames(perc_data) = c("set", "stable")
stdg_data = read.table(stdg_file, header=T, sep="\t")

# Split data into lower and upper stable steady state categories
perc_data = perc_data[order(perc_data$stable, decreasing=T),]

# Calculate deciles
deciles = quantile(perc_data$stable, prob = seq(0, 1, length = 11), type = 5)
d10 = deciles[10]
d1 = deciles[2]
perc_upper = subset(perc_data, stable > d10)
perc_lower = subset(perc_data, stable <= d1)

# Transform concentrations from mM to M
conc_data = conc_data / 1000

# Clean up equations
stdg_data$Eq_kimo = gsub(">|<", "", stdg_data$Eq_kimo, perl=T)

# Calculate the dfGs for all concentration sets
RT = 8.31e-3*303.15

df_per_conc = do.call(cbind, lapply(stdg_data$Reaction, function(reaction) {
  equation = strsplit(subset(stdg_data, Reaction == reaction)[,"Eq_kimo"], " = ", fixed=T)[[1]]
  eq_left = strsplit(equation[1], " + ", fixed=T)[[1]]
  eq_right = strsplit(equation[2], " + ", fixed=T)[[1]]
  drg = subset(stdg_data, Reaction == reaction)$dG
  for (elem in eq_left) {
    elem = strsplit(elem, " ")[[1]]
    if (length(elem) > 1) {
      stoich = as.numeric(elem[1])
      met = paste(elem[2], collapse="")
    } else {
      stoich = 1
      met = paste(elem, collapse="")
    }
    # Calculate delta G
    drg = drg - RT*stoich*log(conc_data[,met])
  }
  for (elem in eq_right) {
    elem = strsplit(elem, " ")[[1]]
    if (length(elem) > 1) {
      stoich = as.numeric(elem[1])
      met = paste(elem[2], collapse="")
    } else {
      stoich = 1
      met = paste(elem, collapse="")
    }
    # Calculate delta G
    drg = drg + RT*stoich*log(conc_data[,met])
  }
  drg = data.frame(Reaction = drg*-1)
  colnames(drg) = reaction
  drg
}))

df_upper = df_per_conc[perc_upper$set,]
df_lower = df_per_conc[perc_lower$set,]

# Plot it
library(ggplot2)
library(reshape2)

df_lower_long = melt(df_lower)
colnames(df_lower_long) = c("Reaction", "DF")
df_lower_long$stability = rep("unstable", nrow(df_lower_long))

df_upper_long = melt(df_upper)
colnames(df_upper_long) = c("Reaction", "DF")
df_upper_long$stability = rep("stable", nrow(df_upper_long))

df_long = rbind(df_upper_long, df_lower_long)

gp = ggplot(df_long, aes(x=DF, fill=stability))
gp = gp + geom_density(alpha=0.4)
gp = gp + geom_vline(xintercept=13.1, alpha=0.4)
gp = gp + theme_bw()
gp = gp + facet_wrap(~Reaction, ncol=6)
gp = gp + scale_fill_manual(values=c("#8073ac","#e08214"))
library(scales)
mysqrt_trans <- function() {
  trans_new("mysqrt",
            transform = base::sqrt,
            inverse = function(x) ifelse(x<0, 0, x^2),
            domain = c(0, Inf))
}
gp <- gp + scale_x_continuous(trans = "mysqrt", breaks=c(0,5,13.1,25,50,100))


outfile = paste(dirname(perc_file), "deltaG_vs_stability.pdf", sep="/")
ggsave(outfile, gp, width=400/25.4, height=225/25.4)
