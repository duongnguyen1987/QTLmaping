#===================================================
# Install and load package for map construction
#==================================================

install.packages("ASMap")
install.packages("qtl")

library(qtl)
library(ASMap)

#============================
#  Set working directory
#============================

setwd("C:/Users/NGU205/OneDrive - CSIRO/Documents/CSIRO/Oat Crown Rust resistance project/DArTseq/DArTSeq 2024/OrderAppendix_1_DO24-8998/ProvenaA_GS7B")



#================================
# Data manipulation and clean up 
#================================

# Load data in R using read.cross()

Provena_A_GS7_B <-read.cross("csvr", file="Provena_A_GS7_B_imputed_csvr.csv", 
                           genotypes=c("A","B"),
                           na.strings = c("N", "H"),crosstype="riself", estimate.map=FALSE)

#=============================================================
# Remove marker and individual with high missing data
#=============================================================

summary(Provena_A_GS7_B)

#Look at the pattern of Missing data, through the function plotMissing()

plotMissing(Provena_A_GS7_B)	#The results indicate individuals with missing data (horizontal lines), as well as  markers with  missing data (vertical lines).

#To plot individuals vs number of succesful genotyped marker

plot(ntyped(Provena_A_GS7_B), ylab="No. typed markers", xlab = "Individuals", main="No. genotypes by individual")

#To plot marker succesful genotyped in number of individuals 

plot(ntyped(Provena_A_GS7_B, "mar"), ylab="No. typed individuals", xlab = "Markers", main="No. genotypes by marker")

# To omit the individuals with lots of missing genotype data (>3000), we may use the function subset.cross(),as follows

Provena_A_GS7_B <- subset(Provena_A_GS7_B, ind=(ntyped(Provena_A_GS7_B)>3000))

# To omit the markers with lots of missing data, we first need to identify the names of the markers. 

# assign markers

nt.bymar <- ntyped(Provena_A_GS7_B, "mar")

# set threshold of markers have genotyping rate of less than 80% of individuals

todrop <- names(nt.bymar[nt.bymar < 80])

# drop.markers() is used to drop the markers that were matching the threshold

Provena_A_GS7_B_1 <- drop.markers(Provena_A_GS7_B, todrop)


#===========================================================
# Remove duplicated markers and individuals 
#===========================================================

# Compare genotypes across pairs of individuals and store the result in a matrix
cg <- comparegeno(Provena_A_GS7_B_1)

# Create a histogram of the lower triangular part of the comparison matrix (excluding the diagonal)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")

# Add a rug plot to the histogram to show individual data points along the x-axis
rug(cg[lower.tri(cg)])

#Identify pairs with well over 95% matching genotypes

Provena_A_GS7_B_1_duplicate <- genClones(Provena_A_GS7_B_1, tol = 0.95, id="Genotype")

Provena_A_GS7_B_1_duplicate$cgd

# Look for duplicated markers 

dup <- findDupMarkers(Provena_A_GS7_B_1, exact.only=FALSE)

Provena_A_GS7_B_1 <- drop.markers(Provena_A_GS7_B_1, dup)
summary(Provena_A_GS7_B_1)


#===========================================================
# Look for markers with distorted segregation patterns
#===========================================================

# We use the function geno.table to inspect the segregation patterns

gt <- geno.table(Provena_A_GS7_B_1)

gt[gt$P.value < 0.05/totmar(Provena_A_GS7_B_1),]

todrop <- rownames(gt[gt$P.value < 1e-10,])

Provena_A_GS7_B_2 <- drop.markers(Provena_A_GS7_B_1, todrop)


#=======================================
# Generate linkage map 
#=======================================

# Calculate a genetic map

# try different thresholds of pvalue for map creation

mstMC <- mstmap.cross(Provena_A_GS7_B_2,  id = "Genotype", bychr = FALSE, suffix = "numeric", dist.fun = "kosambi",
                      objective.fun = "COUNT", p.value = 1e-19, noMap.dist = 30, 
                      noMap.size = 2,miss.thresh = 1, mvest.bc = FALSE, traces=FALSE,
                      as.cross=TRUE)

# create a list of the common chromosome name for each linkage groups

chromosomeL <- list()
for (i in 1:length(mymap))
{
  lgnames <-names(sort(summary(as.factor(sub("^.*[_-]",'',names(mymap[[i]])))),decreasing=TRUE))
  lgnames <- lgnames[lgnames!='']### remove missing names
  if (lgnames[1]%in%names(chromosomeL)) ### append list
    chromosomeL[lgnames[1]] <- list(c(chromosomeL[[lgnames[1]]],names(mymap)[i]))
  else ### add a new chromosome to the list
    chromosomeL[lgnames[1]] <- list(names(mymap)[i])
}

# print the list of chromosomes and the linkage groups belonging to these chromosomes

print(chromosomeL)

# try to merge linkage groups that come from the same chromosome

for (chro in sort(names(chromosomeL))) ### loop through the sorted chromosome names
{
  lgs <- chromosomeL[chro]
  if (length(lgs[[1]]) > 1) ### try to merge linkage groups if there are several
  {
    print(chro)
    mstMClinked <- mergeCross(mstMClinked,merge=lgs,gap=10) ### merge
    mstMClinked <- mstmap(mstMClinked,p.value=1e-2,bychr=TRUE,chr=chro) ### separate linkage groups only if they are still not linked using a very permissive p-value of 1e-2
  }
  else ### rename the linkage group by chromosome
    names(mstMClinked$geno)[names(mstMClinked$geno)==lgs] <- chro
}

# Extract the map information from the cross object
mymap <- pull.map(mstMClinked)

# Plot the map horizontally with a title
plot(mymap, horizontal = TRUE, main = 'Provena x GS7 - LG to Chromosome map')

# Summarize the map object
summary(mymap)

# Summarize the linked cross object
summary(mstMClinked)

# Get the number of markers per linkage group
lg_marker_counts <- nmar(mstMClinked)

# Identify linkage groups with fewer than 7 markers
small_lgs <- names(lg_marker_counts[lg_marker_counts < 7])

# Subset the cross object to remove linkage groups with fewer than 7 markers
mstMClinked_1 <- subset(mstMClinked, chr = !(names(lg_marker_counts) %in% small_lgs))

# Add jitter to the map to avoid overplotting markers
jittermap(mstMClinked_1)

# Plot the filtered genetic map
plot.map(mstMClinked_1, horizontal = TRUE, main = 'Genetic Map of Provena x GS7 (filter LG < 7)')

# Get the number of markers per linkage group in the filtered map
nmar(mstMClinked_1)

# Get the chromosome lengths
chrlen(mstMClinked_1)

# Extract the chromosome lengths from the filtered cross object
chr_lengths <- chrlen(mstMClinked_1)

# Identify chromosomes with non-zero lengths
non_zero_chrs <- names(chr_lengths[chr_lengths > 0])

# Subset the cross object to keep only chromosomes with non-zero lengths
mstMClinked_2 <- subset(mstMClinked_1, chr = non_zero_chrs)

# Add jitter to the map again after filtering
jittermap(mstMClinked_2)

# Get the number of markers per linkage group in the further filtered map
nmar(mstMClinked_2)

# Get the chromosome lengths again after further filtering
chrlen(mstMClinked_2)

# Plot the genetic map after filtering out small linkage groups and zero-length chromosomes
plot.map(mstMClinked_2, horizontal = TRUE, main = 'Provena x GS7 (filtered LG < 7 & length > 0)')

# Extract the map information from the further filtered cross object
mymap <- pull.map(mstMClinked_2)

# Function to generate the merge list dynamically
generate_merge_list <- function(map) {
  chromosome_names <- names(map)
  prefixes <- unique(gsub("\\..*", "", chromosome_names))
  merge_list <- list()
  
  for (prefix in prefixes) {
    matching_chromosomes <- chromosome_names[grep(paste0("^", prefix, "\\."), chromosome_names)]
    if (length(matching_chromosomes) > 1) {
      merge_list[[prefix]] <- matching_chromosomes
    }
  }
  
  return(merge_list)
}

# Generate the merge list
merge_list <- generate_merge_list(mymap)

# Print the merge list to verify
print(merge_list)

# Apply the mergeCross function with the dynamically generated merge list
mymap_combined <- mergeCross(mstMClinked_2, merge = merge_list)

# Check the summary of the new combined cross object
summary(mymap_combined)

# Extract the map from the combined cross object
mymap <- pull.map(mymap_combined)

# Plot the final genetic map vertically
plot.map(mymap, horizontal = FALSE, main = 'Genetic Map of Provena x GS7 (ASMap final)')

# Example renaming for specific chromosomes
names(mymap)[names(mymap) == "1A.1"] <- "1A"
names(mymap)[names(mymap) == "1C.2"] <- "1C"
names(mymap)[names(mymap) == "2C.2"] <- "2C"
names(mymap)[names(mymap) == "2D.1"] <- "2D"
names(mymap)[names(mymap) == "4D.1"] <- "4D"
names(mymap)[names(mymap) == "6D.1"] <- "6D"

# Plot the final genetic map after renaming chromosomes
plot.map(mymap, horizontal = FALSE, main = 'Genetic Map of Provena x GS7 (ASMap final)')

# Check the summary of the final modified cross object
summary(mymap)

# Write the final combined cross object to a CSV file
write.cross(mymap_combined, format = c("csvr"), filestem = "ProvenaxGS7_ASMap_1e-19")