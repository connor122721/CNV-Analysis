# CNV analysis plotting & annotation 
# 11.28.2022
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

# Working directory
setwd("/project/berglandlab/connor/")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Samples
compi=c("Daphnia.pulex.NorthAmerica",
        "Daphnia.pulex.Europe")
samp1 <- as.vector(fin[cont %in% as.character(compi[1])]$Sample)
samp2 <- as.vector(fin[cont %in% as.character(compi[2])]$Sample)

# Read in cnv info
files <- system("ls -f -R /project/berglandlab/connor/cnvs/*cnvs.csv", intern = TRUE)
listi <- lapply(files, fread)
setattr(listi, 'names', list.files(path = "/project/berglandlab/connor/cnvs", pattern = "cnvs.csv"))

# Bind list & make unique CNV
dt <- data.table(rbindlist(listi, use.names = T, idcol = T) %>% 
                   mutate(cnv=paste(seqnames, start, end, CN, sep="_")))

# Get frequency of CNVs
dt.filt <- data.table(dt %>% 
             group_by(.id, cnv) %>% 
             mutate(freq.cnv=length(cnv)))

### Filter ###

# Filter bed file 
bed <- fread("data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.merge50.bed", header = T)
colnames(bed) <- c("seqnames", "start", "end")

# Filter bed
bed.i <- makeGRangesFromDataFrame(bed, keep.extra.columns = T)

# Annotating CNVs w/metadata
cnv.reg <- makeGRangesFromDataFrame(dt.filt, keep.extra.columns = T)

# Filter out sites
cnv.reg.result <- setdiff(x=cnv.reg, y=bed.i, ignore.strand=T)

# Make reverse mapping
revmap <- findOverlaps(cnv.reg.result, cnv.reg, select="arbitrary")
mcols(cnv.reg.result) <- mcols(cnv.reg)[revmap, ,drop=FALSE]

# Convert
dt.filt.fin <- data.table(data.frame(cnv.reg.result))

# Count number of unique CNVs per continent
dt.hist <- data.table(dt.filt.fin %>% 
                    group_by(.id, cnv) %>% 
                    summarize(CN = unique(CN),
                              freq.cnv = unique(freq.cnv)) %>% 
                    mutate(num.samps = case_when(
                      .id %like% "Europe" ~ length(unique(fin[cont=="Daphnia.pulex.Europe"]$Sample)),
                      .id %like% "America" ~ length(unique(fin[cont=="Daphnia.pulex.NorthAmerica"]$Sample)))) %>% 
                    mutate(pop.cnv = freq.cnv/num.samps))

# Count number of unique CNVs per continent
dt2 <- data.table(dt.filt.fin %>% 
         group_by(.id, CN) %>% 
         summarize(num.cnv = length(unique(cnv))))

# Frequency of CNVs
cnv.hist <- {dt.hist %>% 
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

cnv.hist2 <- {dt.hist %>% 
    ggplot(., aes(x=pop.cnv,
                  fill=CN)) +
    geom_histogram() +
    facet_wrap(~.id, scales="free") +
    theme_bw() +
    labs(x = "Population frequency of CNV", 
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

# Regions in genome
cnv.region <- {dt.filt.fin %>% 
  ggplot(., aes(x=(start+end)/2,
                y=median,
                color=CN,
                group=sampleName
                )) +
  geom_point() +
  geom_hline(yintercept=0) +
  facet_wrap(~.id) +
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
        axis.title = element_text(face="bold", size=20))}

### Annotating CNVs ###

# Europe 
dt.filt.euro <- dt.filt.fin[.id %like% "Europe"]

# NAm.
dt.filt.nam <- dt.filt.fin[.id %like% "America"]

# Add gff
txdb <- makeTxDbFromGFF("../daphnia_ref/Daphnia.aed.0.6.gtf")

# Annotate and output gene list
annotateCnvs <- function(cnv, txdb) {

  cnv=dt.filt.euro; txdb=txdb
  
  # Annotating CNVs w/metadata
  cnv.reg <- makeGRangesFromDataFrame(cnv, keep.extra.columns = T)
  
  # Gene annotations
  genes <- transcripts(txdb)
  exonic <- exonicParts(txdb)
  
  # Overlap with genes
  olaps <- findOverlaps(cnv.reg, genes)
  
  # Annotate genes in overlap 
  long_annotated <- cnv.reg[queryHits(olaps)]
  long_annotated$gene <- genes[subjectHits(olaps)]$tx_name
  
  # Extract genes
  olaps.dt <- data.table(data.frame(long_annotated))
  
  # Output genes
  unique(olaps.dt$gene)
  
}

# Europe genes
genes.euro <- annotateCnvs(cnv=dt.filt.euro, txdb=txdb)

# Nam. genes
genes.nam <- annotateCnvs(cnv=dt.filt.nam, txdb=txdb)

# QC
table(genes.nam %in% genes.euro)

# Write gene lists
write.table(genes.euro, file = "cnvs/CNV_all_genes_filtprivate_euro", 
            row.names = F, col.names = F, quote = F)

write.table(genes.nam, file = "cnvs/CNV_all_genes_filtprivate_nam", 
            row.names = F, col.names = F, quote = F)

# Output plots
cnv.mega <- (cnv.hist / cnv.hist2 | cnv.region) +
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')

pdf("figures/cnv.freq.hist.mega.pdf", width = 24, height = 12)
cnv.mega
dev.off()

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



