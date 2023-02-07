# CNV analysis plotting & annotation 
# 12.2.2022
# ijob -c 1 --mem=50G -p largemem -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(cowplot)
library(GenomicRanges)
library(GenomicFeatures)
library(patchwork)
library(cn.mops)

# Working directory
setwd("/project/berglandlab/connor/")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# SNP metadata - TSPs
tot <- data.table(readRDS(file="/project/berglandlab/connor/data/classified_snps_filt.rds"))

# Samples
compi=c("Daphnia.pulex.NorthAmerica",
        "Daphnia.pulex.Europe")
samp1 <- as.vector(fin[cont %in% as.character(compi[1])]$Sample)
samp2 <- as.vector(fin[cont %in% as.character(compi[2])]$Sample)

# Read in CN.mops
resCNMOPS <- readRDS(file="/project/berglandlab/connor/cnvs/Daphnia.pulex.Europe.resCNMOPS.rds")

# Read in cnv info
files <- system("ls -f -R /project/berglandlab/connor/cnvs/*cnvs.csv", intern = TRUE)
listi <- lapply(files, fread)
setattr(listi, 'names', list.files(path = "/project/berglandlab/connor/cnvs", pattern = "cnvs.csv"))

# Bind list & make unique CNV
dt <- data.table(rbindlist(listi, use.names = T, idcol = T) %>% 
                   mutate(cnv=paste(seqnames, start, end, CN, sep="_")))

# Get frequency of CNVs
dt.nofilt <- data.table(dt %>% 
             group_by(.id, cnv) %>% 
             mutate(freq.cnv=length(cnv)))

# Output bed files
#write.table(dt.filt %>% select(seqnames, start, end),
#  file = "cnvs/CNV_all_genes.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# Count number of unique CNVs per continent
dt.nofilt <- data.table(dt.nofilt %>% 
                        group_by(.id, cnv) %>% 
                        summarize(CN = unique(CN),
                                  freq.cnv = unique(freq.cnv)) %>% 
                        mutate(num.samps = case_when(
                          .id %like% "Europe" ~ length(unique(fin[cont=="Daphnia.pulex.Europe"]$Sample)),
                          .id %like% "America" ~ length(unique(fin[cont=="Daphnia.pulex.NorthAmerica"]$Sample)))) %>% 
                        mutate(pop.cnv = freq.cnv/num.samps))

# Count number of unique CNVs per continent
dt.nofilt.tot <- data.table(dt.nofilt %>% 
                    group_by(.id, CN) %>% 
                    summarize(num.cnv = length(unique(cnv))))

# Frequency of CNVs
cnv.hist <- {dt.nofilt %>% 
    ggplot(., aes(x=freq.cnv,
                  fill=CN)) +
    geom_histogram() +
    facet_wrap(~.id) +
    theme_bw() +
    labs(x = "Counts of CNV", 
         fill = "",
         y = "Density") +
    theme(strip.text = element_text(face="bold.italic", size=20),
          legend.text = element_text(size=20), 
          legend.position = "bottom",
          legend.title = element_text(face="bold", size=20),
          axis.text.x = element_text(face="bold", size=20),
          axis.text.y = element_text(face="bold", size=20),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))}

case.cnv <- cnvr(resCNMOPS)
target <- makeGRangesFromDataFrame(dt[cnv %in% 
                dt.nofilt[.id %like% "Euro"][freq.cnv>20]$cnv], 
                         keep.extra.columns = T)

# Overlap 
olaps <- data.frame(findOverlaps(case.cnv, target))

# Forloop plotting
foreach(i=1:length(unique(olaps$queryHits))) %do% {

  print(i)
  quer=unique(olaps$queryHits)[i]
  
  pdf(paste("figures/cnv.validation.cnvr", quer, ".pdf", sep=""))
  plot(resCNMOPS, 
       which=quer, 
       toFile=T, 
       margin=c(5,5))
  dev.off()

}

# Overlap with Toll-like receptor 2 gene
toll = GRanges(seqnames='Scaffold_9200_HRSCAF_10757', ranges=IRanges(6676234, 6686097))
toll.lap <- data.frame(findOverlaps(case.cnv, toll))

# Read in gene annotations
library(readxl)
panth <- data.table(read_excel("/project/berglandlab/daphnia_ref/Daphnia_annotation_PANTHER.xls"))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

# GTF
gene.gtf <- data.table(fread("/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf"))
colnames(gene.gtf)[1:5] <- c("chrom", "file", "sec", "start", "stop")

# Get gene info
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}
gene.gtf$gene <- unlist(str_remove_all(lapply(gene.gtf$V9, extract_attributes, "transcript_id"), pattern = ";"))
gene.gtf <- data.table(gene.gtf %>% 
                mutate(gene.id=tstrsplit(gene, "-")[[1]],
                       splice=tstrsplit(gene, "-")[[2]]))

# Example CNV
exam=c(260000:420000)

# Gene structure plot
gene.struc <- {gene.gtf[chrom=="Scaffold_9200_HRSCAF_10757"][start %in% 6676234:6686097][!sec=="transcript"] %>% 
    ggplot(aes(x=start/1000, y=gene)) +
    geom_linerange(aes(xmin=start/1000, 
                       xmax=stop/1000,
                       color=splice), 
                   size=7) +
    theme_classic() +
    labs(x="Position (Kbp)",
         y="Gene",
         color="Splice") +
    scale_color_brewer(palette = "Set2") +
    theme(strip.text = element_text(face="bold", size=18),
          title = element_text(face="bold.italic", size=20),
          legend.position = "bottom",
          legend.title = element_text(face="bold", size=20),
          legend.text = element_text(face="bold", size=18),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_blank(),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))}

panth[qseqid %in% gene.gtf[chrom=="Scaffold_9200_HRSCAF_10757"][start %in% 6676234:6686097][sec=="transcript"]$gene]

dt.nofilt[cnv %like% "Scaffold_9200_HRSCAF_10757_310001_350000"]

pdf("figures/cnv.genestructure.pdf", width=10, height=8)
gene.struc
dev.off()
