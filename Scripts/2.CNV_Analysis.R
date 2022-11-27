# CNV analysis plotting
# 11.27.2022
# ijob -c 1 --mem=50G -p largemem -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(viridis)
library(doParallel)
library(readxl)
library(cowplot)
library(adegenet)
library(seqinr)

# Working directory
setwd("/project/berglandlab/connor/")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Executable in command line
out <- c("new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")

# Read in total SNPs
tot_snps <- data.table(fread("metadata/snps_new"))

# Load GDS
genofile <- seqOpen(out)

# Samples in GDS
samps <- seqGetData(genofile, var.name = "sample.id")
fin <- fin[Sample %in% samps]

# Read in gene annotations
panth <- data.table(read_excel("../daphnia_ref/Daphnia_annotation_PANTHER.xls"))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

# Gene analyses
pro <- read.fasta("/project/berglandlab/Karen/genomefiles/Daphnia.proteins.aed.0.6.fasta", seqtype="AA")
pro <- data.table(gene=getName(pro), AA.length=getLength(pro))
pro <- pro[,splice:=tstrsplit(gene, "-")[[2]]]

# Read in gtf file
gene.gtf <- data.table(fread("../daphnia_ref/Daphnia.aed.0.6.gtf"))
colnames(gene.gtf)[1:5] <- c("chrom", "file", "sec", "start", "stop")

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}
gene.gtf$gene <- unlist(str_remove_all(lapply(gene.gtf$V9, extract_attributes, "transcript_id"), pattern = ";"))

# SNP metadata - TSPs
tot <- data.table(readRDS(file="/project/berglandlab/connor/data/classified_snps_filt.rds"))

# Samples
compi=c("Daphnia.pulex.NorthAmerica",
        "Daphnia.pulex.Europe")
samp1 <- as.vector(fin[cont %in% as.character(compi[1])]$Sample)
samp2 <- as.vector(fin[cont %in% as.character(compi[2])]$Sample)

# Read in cnv info
files <- system("ls -f -R /project/berglandlab/connor/cnvs/*cnvs.csv", intern = TRUE)
listi <- lapply(files, fread)
setattr(listi, 'names', list.files(path = "/project/berglandlab/connor/cnvs", pattern = "cnvs.csv"))

# Bind list 
dt <- data.table(rbindlist(listi, use.names = T, idcol = T))

# Count number of unique CNVs per continent
dt2 <- data.table(dt %>% 
         group_by(.id) %>% 
         summarize(num.cnv = length(.id),
                   mean= mean(mean, na.rm = T)))

dt %>% 
  ggplot(., aes(x=(start+end)/2,
                y=mean,
                color=CN,
                group=sampleName
                )) +
  geom_smooth(se=F, linetype=2) +
  geom_point() +
  geom_hline(yintercept=0) +
  theme_bw() +
  labs(x = "Position", 
       color = "",
       y = "Normalized Read Depth") +
  theme(strip.text = element_text(face="bold.italic", size=20),
        legend.text = element_text(size=20), 
        legend.position = "bottom",
        legend.title = element_text(face="bold", size=20),
        axis.text.x = element_text(face="bold", size=20),
        axis.text.y = element_text(face="bold", size=20),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.title = element_text(face="bold", size=20))

# Overlap with exons
toti <- data.table(tot[classified=="shared_poly"] %>% select(start=position, stop=position, chrom, variant.id))
cnvi <- dt %>% select(start, stop=end, chrom=seqnames, CN)
setkey(cnvi, chrom, start, stop)

# Calculate exon length
laps <- na.omit(foverlaps(toti, cnvi, type="within") %>% select(-c("i.stop")))
colnames(laps)[c(2:3, 5)] <- c("cnv.start", "cnv.stop", "position")

# Proportion of TSPs within CNVs by category
cnv.tsp <- data.table(laps %>% 
  group_by(CN) %>% 
  summarize(n=length(unique(variant.id))) %>% 
  mutate(n.stan=n/length(unique(toti$variant.id))))

cnv.plot <- {cnv.tsp %>% 
  ggplot(aes(x=n, 
             y=CN)) +
  geom_col(fill="purple", color="black") +
  theme_bw() +
  labs(x = "Number of TSPs",
       y = "CNV Category") +
  theme(strip.text = element_text(face="bold.italic", size=18),
        legend.text = element_text(face="bold", size=18),
        legend.title = element_text(face="bold", size=18),
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.title = element_text(face="bold", size=20))}
  
pdf("figures/cnv.num.tsp.plot.pdf")
cnv.plot
dev.off()
  
# Filtered bed file
bed="data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.merge50.bed"




