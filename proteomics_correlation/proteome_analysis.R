# Define infiles
enzyme_file = "/home/johannes/proj/kimo/tools/maks/proteomics_correlation/CBB_Enzymes.tab"
proteomics_file = "/home/johannes/proj/kimo/data/2017-06-26/proteomics_20170410_long.csv"

# Load data
enzymes = read.table(enzyme_file, stringsAsFactors=F, header=T, sep="\t")
proteome = read.table(proteomics_file, stringsAsFactors=F, header=T, sep=",")

# Subset proteome to Calvin cycle
colnames(enzymes)[3] = "protein"
prot_cbb = merge(enzymes, proteome)

# Plot it
prot_cbb$Label = paste(prot_cbb$Abbrev, prot_cbb$protein, sep=" - ")
prot_cbb$Upper = prot_cbb$rel_intensity + prot_cbb$rel_intensity * prot_cbb$CV
prot_cbb$Lower = prot_cbb$rel_intensity - prot_cbb$rel_intensity * prot_cbb$CV

library(ggplot2)

gp = ggplot(prot_cbb, aes(x=growthrate, y=rel_intensity))
gp = gp + geom_line()
gp = gp + geom_errorbar(aes(ymin=Lower, ymax=Upper))
gp = gp + geom_point()
gp = gp + theme_bw()
gp = gp + facet_wrap(~Gene_Name + Label)
gp = gp + scale_y_continuous(trans="log2",
          breaks=c(1/8, 1/4, 1/2, 1, 2, 3),
          labels=c("1/8", "1/4", "1/2", "1", "2", "3"))
gp = gp + ylab("Relative intensity (error bars = CV)") + xlab("Growth rate")

ggsave("/home/johannes/proj/kimo/art/2017-06-26/CBB_enzymes_vs_growth_rate.pdf",
       gp, width=420/25.4, height=297/25.4)

# Perform hierarchical clustering
library(reshape2)
library(MKmisc)
prot_cbb_wide = dcast(prot_cbb[,c("Label", "rel_intensity", "growthrate")], Label ~ growthrate, value.var="rel_intensity")
clustering = hclust(corDist(as.matrix(prot_cbb_wide[,2:ncol(prot_cbb_wide)])))
prot_cbb_wide$Trend = c("Down","Unclear","Up")[cutree(clustering, 3)]

# Merge with full Calvin cycle dataset
prot_cbb = merge(prot_cbb, prot_cbb_wide[,c("Label", "Trend")])

# Plot again
gp = ggplot(prot_cbb, aes(x=growthrate, y=rel_intensity, colour=Trend))
gp = gp + geom_line()
gp = gp + geom_errorbar(aes(ymin=Lower, ymax=Upper))
gp = gp + geom_point()
gp = gp + theme_bw()
gp = gp + facet_wrap(~Gene_Name + Label)
gp = gp + scale_y_continuous(trans="log2",
          breaks=c(1/8, 1/4, 1/2, 1, 2, 3),
          labels=c("1/8", "1/4", "1/2", "1", "2", "3"))
gp = gp + ylab("Relative intensity (error bars = CV)") + xlab("Growth rate")

ggsave("/home/johannes/proj/kimo/art/2017-06-26/CBB_enzymes_vs_growth_rate.clustered.pdf",
       gp, width=420/25.4, height=297/25.4)

# Plot versus light intensity
gp = ggplot(prot_cbb, aes(x=light, y=rel_intensity, colour=Trend))
gp = gp + geom_line()
gp = gp + geom_errorbar(aes(ymin=Lower, ymax=Upper))
gp = gp + geom_point()
gp = gp + theme_bw()
gp = gp + facet_wrap(~Gene_Name + Label)
gp = gp + scale_y_continuous(trans="log2",
          breaks=c(1/8, 1/4, 1/2, 1, 2, 3),
          labels=c("1/8", "1/4", "1/2", "1", "2", "3"))
gp = gp + ylab("Relative intensity (error bars = CV)") + xlab("Light intensity")

ggsave("/home/johannes/proj/kimo/art/2017-06-26/CBB_enzymes_vs_light.clustered.pdf",
       gp, width=420/25.4, height=297/25.4)

# Make a plot of the clustering
clustering$labels = prot_cbb_wide$Label
pdf("/home/johannes/proj/kimo/art/2017-06-26/CBB_enzymes_vs_light.clustered.pdf", width=200/25.4, height=200/25.4)
plot(clustering)
dev.off()
