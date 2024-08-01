#=========================================
# Install and load Load DArTR package
#=========================================

install.packages("dartR") 
library(dartR)


#============================
#  Set working directory
#============================

setwd("C:/Users/NGU205/OneDrive - CSIRO/Documents/CSIRO/Oat Crown Rust resistance project/DArTseq/DArTSeq 2024/OrderAppendix_1_DO24-8998/ProvenaA_GS7B")

#====================================================================================
# Load DArTSeq data into R using dartR package and conducted initial filtering
#====================================================================================


gl <- gl.read.dart(filename= "Report_DO24-8998_SNP_mapping_2.csv",
                   ind.metafile ="Group_pop.csv", recalc = TRUE) 

#=========================================================================
# Check number of markers, individuals, and population in the data set 
#=========================================================================

nLoc(gl) #number of loci
locNames(gl) #list of loci
nInd(gl) #number of individuals (specimens or samples)
indNames(gl) #list of individuals
pop(gl)     #list of population assignments for each individual

#=======================================================================================
# Get the data for parental lines and RILs population of derived from Provena and GS7
#========================================================================================

# Extract data

gl_ProvenaxGS7 <- gl.keep.pop(gl, pop.list=c('ProvenaxGS7', 'GS7',"Provena"))

# Save the data for the extracted population into r data. 
save(gl_ProvenaxGS7, file="gl_ProvenaxGS7_SNP.rdata")

# The data could be loaded in using load("gl_ProvenaxGS7_SNP.rdata")

#=========================================================================
# Filtering steps of the data
#=========================================================================

gl_filter<-gl.filter.allna(gl_ProvenaxGS7, by.pop = FALSE, recalc = FALSE, verbose = 3) # Filter individuals and markers that are scored as all missing (NA)

gl_filter <- gl.filter.monomorphs(gl_filter, verbose = 3) #  Remove monomorphic loci 

gl_filter <- gl.filter.callrate(gl_filter, method='loc', threshold=0.5, verbose=3) # Remove loci based on Call Rate, threshold = 0.5

gl_filter <- gl.filter.maf(gl_filter, threshold=0.01, verbose=3) #  Remove loci with MAF < 0.01 over all the dataset

gl_filter <-gl.filter.heterozygosity(gl_filter, t.upper = 0.95, t.lower = 0, verbose = 3) # Remove individuals with heterozygosity > 0.95

#=============================================================
# Convert greenlight file to vcf file to load into TASSEL 
# for convering genotypic data format into ABH format
#=============================================================


gl2vcf(gl_filter, plink_path = getwd(),outfile = "gl_ProvenaxGS7_SNP", outpath = getwd(), snp_pos = "ChromPosSnp_Oat_OT3098_v2",  snp_chr = "Chrom_Oat_OT3098_v2")
