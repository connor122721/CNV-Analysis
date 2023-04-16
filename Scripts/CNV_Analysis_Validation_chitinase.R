# CNV analysis plotting & annotation 
# 4.9.2022
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

case.cnv <- cnvr(resCNMOPS)

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

# GRANGES GTF
gene.gtf2 <- GRanges(seqnames=gene.gtf$chrom,
                     sec=gene.gtf$sec,
                     gene.id=gene.gtf$gene.id,
                     gene=gene.gtf$gene,
                     ranges=IRanges(gene.gtf$start, 
                                    gene.gtf$stop))

# Find chitinases
chit <- unique(panth[bio_func %like% "chitin"]$qseqid)

# Find gram-negative binding proteins
gram <- unique(panth[bio_func %like% "Gram-negative"]$qseqid)

# Make bed file of genes
chitin <- gene.gtf[gene %in% chit]
chitin1 <- chitin %>% dplyr::select(chrom, start, stop, gene)

# Output bed
write.table(chitin1, file = "control_cnvs/chitin.genes",
            quote = F, row.names = F, col.names = F, sep='\t')

# Overlapped
pos <- data.table(fread("control_cnvs/sort.unique.chitin.cnv.hits"))
colnames(pos) <- c("seqnames", "start", "end")

# Find individuals w/ CNVs
pos.dt <- data.table(dt %>% right_join(pos[seqnames=="Scaffold_9201_HRSCAF_10758"]))
pos.dt[.id %like% "Europe"]
# April_2017_DBunk_13_finalmap_mdup.bam; Scaffold_9201_HRSCAF_10758_70001_200000_CN1




# Find 
gene.lap <- findOverlaps(case.cnv, gene.gtf2)
hits <- data.frame(chrom = seqnames(case.cnv)[queryHits(gene.lap)],
                   ranges(case.cnv)[queryHits(gene.lap)])

genes.list <- lapply(gene.gtf[chrom %in% unique(hits)], FUN = )




# Toll like gene locations 
toll <- gene.gtf[gene %in% panth[bio_func%like%"toll"]$qseqid][sec=="transcript"]

# Overlap with Toll-like receptor 2 gene
toll.cand <- GRanges(seqnames=toll$chrom, 
               ranges=IRanges(toll$start, toll$stop))

# TEST
#toll.cand <- GRanges(seqnames="Scaffold_9200_HRSCAF_10757", ranges=IRanges(1,1000000))

# Find 
toll.lap <- findOverlaps(case.cnv, toll.cand)
hits <- data.frame(chrom = seqnames(case.cnv)[queryHits(toll.lap)],
                   ranges(case.cnv)[queryHits(toll.lap)])

genes.list <- lapply(gene.gtf[chrom %in% unique(hits)], FUN = )

# Gene structure plot
gene.struc <- {gene.gtf[chrom %in% unique(hits$chrom)][start %in% 
                        c(unique(hits$start):unique(hits$end))][!sec=="transcript"] %>% 
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
