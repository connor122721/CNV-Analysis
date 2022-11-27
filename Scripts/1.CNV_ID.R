# Identify CNVs from Bams
# 9.5.2022
# ijob -c 10 --mem=50G -p standard -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(cn.mops)
library(data.table)
library(foreach)

# Working directory
setwd("/project/berglandlab/connor/cnvs")

# Bam files
BAMFiles1 <- list.files("/project/berglandlab/connor/BACKUP_scratch/all_bam/bam_dorthe_fin", pattern="RG.bam$", full.names=TRUE)
BAMFiles2 <- list.files("/project/berglandlab/connor/BACKUP_scratch/all_bam/Euro_bams", pattern=".bam$", full.names=TRUE)
BAMFiles3 <- list.files("/project/berglandlab/connor/BACKUP_scratch/all_bam/final_bam", pattern=".bam$", full.names=TRUE)
BAMFiles <- c(BAMFiles1, BAMFiles2, BAMFiles3)

# Metadata
fin <- data.frame(read.csv("../metadata/samples.fin.9.8.22.csv"))

# Merge with metadata
dt <- data.table(merge(data.frame(Sample=sub(pattern = "_$", 
                            replacement =  "", 
                            x = IRanges::reverse(data.table::tstrsplit(IRanges::reverse(data.table::tstrsplit(BAMFiles, 
                                                              "finalmap")[[1]]), "/")[[1]])),
                 File=BAMFiles), fin, by="Sample"))

# Chromosomes
chroms <- scan("/project/berglandlab/connor/metadata/goodChrom.txt", what="")
conts <- unique(dt[Species =="Daphnia pulex"]$cont)

# Register cores
doParallel::registerDoParallel(cores = 10)

# Go through each species
dt.out <- foreach(i=1:length(conts), .combine = "rbind", .errorhandling = "remove") %do% {
  
  # Bam file
  print(paste(conts[i], i, sep=" "))
  bams <- c(dt[cont==conts[i]]$File)
  
  # Get read counts
  bamDataRanges <- getReadCountsFromBAM(bams, 
                                        refSeqNames=chroms,
                                        parallel = 10, 
                                        WL=10000)
  
  # Identify CNVs
  resCNMOPS <- cn.mops(bamDataRanges, 
                       useMedian = T, 
                       parallel = 10)
  
  # Number of CNVs
  resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
  
  # Summarize
  segm <- as.data.frame(segmentation(resCNMOPS))
  CNVs <- as.data.frame(cnvs(resCNMOPS))
  CNVregions <- as.data.frame(cnvr(resCNMOPS))
  
  # Output
  write.csv(segm, paste("/project/berglandlab/connor/cnvs/", 
                        conts[i],"segmentation.csv", sep=""))
  
  write.csv(CNVs, paste("/project/berglandlab/connor/cnvs/", 
                        conts[i], "cnvs.csv", sep=""))
  
  write.csv(CNVregions, paste("/project/berglandlab/connor/cnvs/", 
                              conts[i],"cnvr.csv", sep=""))
  
  # Finish
  return(CNVs)
}
