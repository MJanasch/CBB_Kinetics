#!/usr/bin/env Rscript

### LOAD DATA ##################################################################

# Read infile names from command line
args = commandArgs(trailingOnly=T)
fccs_file = args[1] # An FCCs tab.gz file
header_file = args[2] # Reaction header file

# Load data
library(data.table)
fccs = as.data.frame(fread(paste(c("gzip -dc ", fccs_file), collapse=""),
  header=T, sep="\t"
  ))

# The enzyme in the column (effector) influences
# the flux of reactions in rows (target)

# Load custom labels for reactions
custom_rxn_labels = scan(header_file, character(), quote = "")

# Set FCCs that are close to zero to zero
fccs[,4:ncol(fccs)][abs(fccs[,4:ncol(fccs)]) < 1e-6] = 0 # Markus approves

# Replace reaction numbers with reaction names
fccs$Reaction = custom_rxn_labels[fccs$Reaction]

# Load pvclust clustering library
library(pvclust)

# Store future rownames of fccs
fccs_rownames = fccs$Reaction

# Remove unnecessary columns of fccs and make matrix
fccs = as.matrix(fccs[,4:ncol(fccs)])
rownames(fccs) = fccs_rownames

# Transpose and transform fccs for further analysis
fccs = asinh(t(fccs))
# Transposed; effector now in row and target in column

### EFFECTOR PATTERNS (CONTROL) ################################################

# Reduce to 10% of original data
sampled_eff_cols = sample(1:ncol(fccs), ceiling(ncol(fccs)*0.1))

# Scale data in preparation of PCA
fccs_eff_scaled = scale(fccs[,sampled_eff_cols])
# Scale applied to columns, i.e. normalized by scaling per target reaction

# Perform PCA
fccs_eff_pca = prcomp(fccs_eff_scaled[,is.finite(colSums(fccs_eff_scaled))])

# Perform clustering on rotated and reduced data
eff_clst = pvclust(t(fccs_eff_pca$x), nboot=10000, parallel=T)

# Plot clustering
out_eff = paste(dirname(fccs_file), "effector_patterns_pvclust.pdf", sep="/")
pdf(out_eff, height=240/25.4, width=360/25.4)
plot(eff_clst)
pvrect(eff_clst, alpha=0.95)
dev.off()

### TARGET PATTERNS (BEING CONTROLLED) #########################################

# Fold over data frame to bind vectors of target patterns
fccs_chunks = split(1:ncol(fccs), ceiling(seq_along(1:ncol(fccs))/nrow(fccs)))
fccs_tar = t(as.matrix(rbindlist(
 lapply(fccs_chunks, function(x){
     as.data.frame(fccs[,unlist(x)])
   })
 )))
colnames(fccs_tar) = rep(rownames(fccs_tar), ncol(fccs_tar)/nrow(fccs_tar))

# Reduce to 10% of original data
sampled_tar_cols = sample(1:ncol(fccs_tar), ceiling(ncol(fccs_tar)*0.1))

# Scale data in preparation of PCA
fccs_tar_scaled = scale(fccs_tar[,sampled_tar_cols])
# Scale applied to columns, i.e. normalized by scaling per target reaction

# Perform PCA
fccs_tar_pca = prcomp(fccs_tar_scaled[,is.finite(colSums(fccs_tar_scaled))])

# Perform clustering on rotated and reduced data
tar_clst = pvclust(t(fccs_tar_pca$x), nboot=10000, parallel=T)

# Plot clustering
out_tar = paste(dirname(fccs_file), "target_patterns_pvclust.pdf", sep="/")
pdf(out_tar, height=240/25.4, width=360/25.4)
plot(tar_clst)
pvrect(tar_clst, alpha=0.95)
dev.off()
