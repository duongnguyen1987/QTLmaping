#===========================================================
# Install and load package and function for QTL mapping
#===========================================================

install.packages("qtl")
library(qtl)

source("https://rqtl.org/dotfunc.R")

#============================
#  Set working directory
#============================

setwd("C:/Users/NGU205/OneDrive - CSIRO/Documents/CSIRO/Oat Crown Rust resistance project/DArTseq/DArTSeq 2024/OrderAppendix_1_DO24-8998/ProvenaA_GS7B")


Map_QTL <-read.cross("csvr", file="ProvenaxGS7_ASMap_1e-19_OCR_2023.csv", 
                     genotypes=c("AA","BB"),
                     na.strings = "-",crosstype="riself", estimate.map=FALSE)
# Jitter the marker positions on the map
jittermap(Map_QTL)

# Print summary statistics of the genetic map
summaryMap(Map_QTL)
summary(Map_QTL)

# Plot phenotype data for the first phenotype column
plot.pheno(Map_QTL, pheno.col=1)

# Plot the missing data pattern
plot.missing(Map_QTL)

# Plot the genetic map
plot.map(Map_QTL, horizontal=FALSE, main='Genetic Map of Provena x GS7 (ASMap Final)')

# Plot a heatmap of the LOD scores
heatMap(Map_QTL, what = "lod", lmax = 12)

#============================
#  Simple interval mapping 
#============================


# Perform simple interval mapping for the first phenotype column
out.pheno.em <- scanone(Map_QTL, pheno.col=c(1), method="em")

# Plot the LOD scores for the entire genome
plot(out.pheno.em, lodcolumn=c(1), lty=1, col=c("blue"))

# zooming in on a particular chromosome

# Plot the LOD scores for specific chromosomes
plot(out.pheno.em, lodcolumn=c(1), chr=c("2A"), lty=1, col=c("blue"))
plot(out.pheno.em, lodcolumn=c(1), chr=c("4A"), lty=1, col=c("blue"))
plot(out.pheno.em, lodcolumn=c(1), chr=c("4D"), lty=1, col=c("blue"))
plot(out.pheno.em, lodcolumn=c(1), chr=c("7A"), lty=1, col=c("blue"))

# Get confidence intervals for the QTL on specific chromosomes
lodint(out.pheno.em, chr="2A")
lodint(out.pheno.em, chr="4A")
lodint(out.pheno.em, chr="4D")
lodint(out.pheno.em, chr="7A")

# Print a summary of the simple interval mapping results
summary(out.pheno.em)

# Generate a cross-tabulation of genotypes for two markers
geno.crosstab(Map_QTL, "SNP-459779164_4A", "SNP-408202291_4D", eliminate.zeros=TRUE)

# Generate a genotype table for chromosome 4A and filter by p-value
gt <- geno.table(Map_QTL, chr="4A")
gt[gt$P.value < 0.001, ]

#====================================================================================================
#  Use scanone to perform a permutation test to get a genome-wide LOD significance threshold
#====================================================================================================

# Perform a permutation test for simple interval mapping
operm.em <- scanone(Map_QTL, pheno.col=c(1), method="em", n.perm=1000)

# Plot the permutation test results
plot(operm.em, lodcolumn=c(1), lty=1, col=c("blue"))

# Print a summary of the permutation test results
summary(operm.em, alpha=0.05)

# Combine the permutation results with the mapping results and print the summary
summary(out.em, perms=operm.em, alpha=0.05, pvalues=TRUE)

#===============================================================================
#  Create summary table including effect and lod infomration for each chromosome
#===============================================================================

# Perform a genome scan for the effect of the markers on the phenotype
effects <- effectscan(sim.geno(Map_QTL))

# Get a summary of the LOD scores
Lod <- summary(out.pheno.em)

# Convert row names of 'effects' to a new column 'Marker'
effects$Marker <- rownames(effects)

# Match markers and merge effect values into the summary table
summary_with_effects <- merge(Lod, effects[, c("Marker", "a")], by.x = "row.names", by.y = "Marker", all.x = TRUE)

# Print the updated summary table with effect values
print(summary_with_effects)

#================================================
#  Use CIM to perform composite interval mapping
#================================================

# Perform composite interval mapping for the first phenotype column
out_pheno_cim.em <- cim(Map_QTL, pheno.col = (1), method=c("em"), map.function=c("kosambi"))

# Plot the composite interval mapping results
plot(out_pheno_cim.em, lodcolumn=c(1), lty=1, show.marker.names=FALSE, col=c("red"))

# Get confidence intervals for the QTL on specific chromosomes
lodint(out_pheno_cim.em, chr="2A")
lodint(out_pheno_cim.em, chr="4A")
lodint(out_pheno_cim.em, chr="4D")
lodint(out_pheno_cim.em, chr="7A")

# Plot the composite interval mapping results for specific chromosomes
plot(out_pheno_cim.em, chr="4A", show.marker.names=FALSE, col=c("red"))
plot(out_pheno_cim.em, chr="2A", show.marker.names=FALSE, col=c("red"))
plot(out_pheno_cim.em, chr="4D", show.marker.names=FALSE, col=c("red"))
plot(out_pheno_cim.em, chr="7A", show.marker.names=FALSE, col=c("red"))

# Get the marker covariance
attr(out_pheno_cim.em, "marker.covar.pos")

# Get a summary of the composite interval mapping results
Lod <- summary(out_pheno_cim.em)

# Perform a genome scan for the effect of the markers on the phenotype
effects <- effectscan(sim.geno(Map_QTL))

# Convert row names of 'effects' to a new column 'Marker'
effects$Marker <- rownames(effects)

# Match markers and merge effect values into the summary table
summary_with_effects <- merge(Lod, effects[, c("Marker", "a")], by.x = "row.names", by.y = "Marker", all.x = TRUE)

# Print the updated summary table with effect values
print(summary_with_effects)

# Plot the results of both simple interval mapping (SIM) and composite interval mapping (CIM)
plot(out.pheno.em, out_pheno_cim.em, show.marker.names=FALSE, col=c("blue", "red"))

# Add threshold lines to the plot
abline(h=2.92, col="blue", lty=2) # SIM threshold line

# Add a legend to the plot
legend("topright", legend=c("Composite Interval Mapping (CIM)", "Simple Interval Mapping (SIM)", "LoD Threshold"), 
       col=c("red", "blue", "blue"), 
       lty=c(1, 1, 2), cex=0.8)

# Add another legend to the plot
legend("topright", legend=c("SIM", "CIM"), col=c("red", "blue"), lty=1, cex=0.8)

# Plot phenotype values for a specific marker
plotyield <- plotPXG(Map_QTL, marker = "SNP_458139651_4A", pheno.col = (1))

# Plot the effect of a specific marker on the phenotype
plotyieldeffect <- effectplot(Map_QTL, mname1 = "SNP_458139651_4A", pheno.col = (1))
