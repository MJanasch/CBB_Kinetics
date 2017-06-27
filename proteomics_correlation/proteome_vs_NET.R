# Define infiles
enzyme_file = "/home/johannes/proj/kimo/tools/maks/proteomics_correlation/CBB_Enzymes.tab"
proteomics_file = "/home/johannes/proj/kimo/data/2017-06-26/proteomics_20170410_long.csv"
reaction_gene_file = "/home/johannes/proj/kimo/tools/maks/proteomics_correlation/data/syn_Knoop.reaction_gene.long.tab"
NET_class_file = "/home/johannes/proj/kimo/tools/maks/proteomics_correlation/data/syn_Knoop.reaction_NETclass.A-F.tab"

# Load data
enzymes = read.table(enzyme_file, stringsAsFactors=F, header=T, sep="\t")

proteome = read.table(proteomics_file, stringsAsFactors=F, header=T, sep=",")

rxn_gene = read.table(reaction_gene_file, stringsAsFactors=F, header=F, sep="\t")
colnames(rxn_gene) = c("Reaction", "protein")

NETclass = read.table(NET_class_file, stringsAsFactors=F, header=F, sep="\t")
colnames(NETclass) = c("Reaction", "NETclass")

# Subset proteome to Calvin cycle
colnames(enzymes)[3] = "protein"
prot_cbb = merge(enzymes, proteome)

# Subset to relevant columns
prot_cbb = prot_cbb[,c("protein", "Abbrev", "light", "growthrate", "CV", "rel_intensity", "Gene_Name")]

# Add NET and Reaction annotations
prot_cbb = merge(merge(rxn_gene, NETclass), prot_cbb, all.y=T)

# Remove 1000 µE setting
prot_cbb = subset(prot_cbb, light < 1000)

# Label the genes
prot_cbb$Label = paste(prot_cbb$Abbrev, prot_cbb$protein, prot_cbb$Reaction, sep=" - ")

# Perform hierarchical clustering
library(reshape2)
library(MKmisc)
prot_cbb_wide = dcast(prot_cbb[,c("Label", "rel_intensity", "growthrate")], Label ~ growthrate, value.var="rel_intensity")
clustering = hclust(corDist(as.matrix(prot_cbb_wide[,2:ncol(prot_cbb_wide)])))
prot_cbb_wide$Trend = c("Down","Unclear", "Unclear", "Up", "Unclear", "Unclear")[cutree(clustering, 6)]

# Merge with full Calvin cycle dataset
prot_cbb = merge(prot_cbb, prot_cbb_wide[,c("Label", "Trend")])

# Calculate upper and lower bound values based on CV
prot_cbb$Upper = prot_cbb$rel_intensity + prot_cbb$rel_intensity * prot_cbb$CV
prot_cbb$Lower = prot_cbb$rel_intensity - prot_cbb$rel_intensity * prot_cbb$CV

# Plot it
library(ggplot2)

prot_cbb$NETclass[is.na(prot_cbb$NETclass)] = "Not in NET"

gp = ggplot(prot_cbb, aes(x=growthrate, y=rel_intensity, colour=Trend))
gp = gp + geom_line(aes(linetype=NETclass))
gp = gp + geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.005)
gp = gp + geom_point()
gp = gp + theme_bw()
gp = gp + facet_wrap(~Gene_Name + Label, ncol=5)
gp = gp + scale_y_continuous(trans="log2",
          breaks=c(1/8, 1/4, 1/2, 1, 2, 3),
          labels=c("1/8", "1/4", "1/2", "1", "2", "3"))
gp = gp + scale_linetype_manual(values=c(1,2,3))
gp = gp + scale_colour_manual(values=c("#8073ac", "#bdbdbd", "#e08214"))
gp = gp + ylab("Relative intensity (error bars = CV)") + xlab("Growth rate")

ggsave("/home/johannes/proj/kimo/art/2017-06-27/CBB_enzymes_vs_NET.clustered.pdf",
       gp, width=420/25.4, height=300/25.4)

# Make a plot of the clustering
clustering$labels = prot_cbb_wide$Label
pdf("/home/johannes/proj/kimo/art/2017-06-27/CBB_enzyme_proteomics.clustering.2.pdf", width=200/25.4, height=200/25.4)
plot(clustering)
dev.off()
