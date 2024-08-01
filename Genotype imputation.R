
#=========================================================================
# Install and load package  ABHgenotypeR for marker imputation
#=========================================================================

install.packages("ABHgenotypeR") 
library(ABHgenotypeR)

#============================
#  Set working directory
#============================

setwd("C:/Users/NGU205/OneDrive - CSIRO/Documents/CSIRO/Oat Crown Rust resistance project/DArTseq/DArTSeq 2024/OrderAppendix_1_DO24-8998/ProvenaA_GS7B")


#==========================================================================
# After convert vcf file  to ABH format in TASSEL data was loaded back in R 
# using ABHgenotypeR to impute the missing
#==========================================================================

# Data was loaded
genotypes <- readABHgenotypes("Provena_A_GS7_B.csv", nameA = "Provena", nameB = "GS7")

# Genotypes can be plotted by ABHgenotypeR package

plotGenos(genotypes) # Plot all chromosomes
plotGenos(genotypes, chromToPlot = 1) # Plot chromosome 1

# Imputation of missing genotypes
# Basically, if the genotypes left and right of a stretch of missing data are identical the genotypes are filled in

postImpGenotypes <- imputeByFlanks(genotypes)

plotGenos(postImpGenotypes) # Plot all chromosomes after imputation

# Error corrections

ErrCorr1Genotypes <- correctUndercalledHets(postImpGenotypes, maxHapLength = 3)

plotGenos(ErrCorr1Genotypes)

ErrCorr2Genotypes <- correctStretches(ErrCorr1Genotypes, maxHapLength = 3)

plotGenos(ErrCorr2Genotypes)

# Check allele frequency and marker density after imputation and correction steps

plotMarkerDensity(genos = ErrCorr2Genotypes)
plotAlleleFreq(genos = ErrCorr2Genotypes)


# Writed file to csv format after imputation and correction steps

writeABHgenotypes(ErrCorr2Genotypes, outfile = "Provena_A_GS7_B_imputed.csv")
