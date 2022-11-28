# CNVs - GO enrichment
# 11.28.2022
# ijob -c 1 --mem=50G -p largemem -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(cowplot)
library(readxl)
library(viridis)
library(clusterProfiler)
library(enrichplot)
library(seqinr)

# Working directory
setwd("/project/berglandlab/connor/")

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

# Read in gene annotations
panth <- data.table(read_excel("../daphnia_ref/Daphnia_annotation_PANTHER.xls"))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

# Metadata
fin <- data.table(read.csv("metadata/samples.fin.9.8.22.csv"))

# Gene analyses
pro <- read.fasta("/project/berglandlab/Karen/genomefiles/Daphnia.proteins.aed.0.6.fasta", seqtype="AA")
pro <- data.table(gene=getName(pro), AA.length=getLength(pro))
pro <- pro[,splice:=tstrsplit(gene, "-")[[2]]]

# Read in Gene data
gene.dt.euro <- data.table(read.delim("cnvs/CNV_all_genes_filtprivate_euro", header = F))
gene.dt.nam <- data.table(read.delim("cnvs/CNV_all_genes_filtprivate_nam", header = F))

# Read in gmatrix
gmap <- readRDS(file = "/project/berglandlab/connor/data/gmap.daphnia.pulex.rds")

### Europe ###

# GO analysis
gene.cand <- unique(gene.dt.euro$V1)

# Universe
univ <- unique(gene.gtf$gene)

# Enrichment analysis
go.test <- enricher(gene = gene.cand,
         pvalueCutoff = 0.05,
         pAdjustMethod = "fdr",
         minGSSize = 10,
         maxGSSize = 500,
         qvalueCutoff = 0.01,
         universe = univ, 
         TERM2GENE = gmap)

# Find term description for each GO term
go.test2 <- data.table(go.test@result) %>% filter(p.adjust < 0.05 & qvalue < 0.01)
pp.euro <- foreach(i=1:length(go.test2$ID), .combine = "rbind", .errorhandling = "remove") %do% {
  tmp.df = go.test2 %>%
    filter(ID == go.test2$ID[i])
  
  message(paste(i, go.test2$ID[i], sep = " | "))
  
  go.test2[ID==go.test2$ID[i]] %>% 
    mutate(term=go2term(ID)$Term) %>% 
    mutate(x = eval(parse(text=GeneRatio)),
           y = eval(parse(text=BgRatio))) %>% 
    mutate(odd = x/y)
}

### NAM ###

# GO analysis
gene.cand <- unique(gene.dt.nam$V1)

# Universe
univ <- unique(gene.gtf$gene)

# Enrichment analysis
go.test <- enricher(gene = gene.cand,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "fdr",
                    minGSSize = 10,
                    maxGSSize = 500,
                    qvalueCutoff = 0.01,
                    universe = univ, 
                    TERM2GENE = gmap)

# Find term description for each GO term
go.test2 <- data.table(go.test@result) %>% filter(p.adjust < 0.05 & qvalue < 0.01)
pp.nam <- foreach(i=1:length(go.test2$ID), .combine = "rbind", .errorhandling = "remove") %do% {
  tmp.df = go.test2 %>%
    filter(ID == go.test2$ID[i])
  
  message(paste(i, go.test2$ID[i], sep = " | "))
  
  go.test2[ID==go.test2$ID[i]] %>% 
    mutate(term=go2term(ID)$Term) %>% 
    mutate(x = eval(parse(text=GeneRatio)),
           y = eval(parse(text=BgRatio))) %>% 
    mutate(odd = x/y)
}

# Plot GO term 
plot.go.euro <- {pp.euro %>% 
    filter(p.adjust<0.05) %>% 
    arrange(desc(p.adjust)) %>% 
    ggplot(aes(x = log2(odd),
               y = reorder(term, log2(odd)),
               size = Count,
               color = log10(p.adjust))) +
    geom_point() +
    labs(x=expression(bold(paste("log"[2], "(Enrichment)"))), 
         y="",
         color=expression(bold(paste("FDR log"[10], "(", italic(p),"-value)"))),
         size="Candidate gene number",
         title=expression(bold(paste("All CNVs", bold(italic(" D. pulex ")), "Europe")))) +
    theme_bw() +
    scale_color_viridis(option = "plasma") +
    theme(title = element_text(face="bold", size=18),
          legend.text = element_text(face="bold", size=18),
          legend.title = element_text(face="bold", size=18),
          legend.background = element_blank(),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18))}

plot.go.nam <- {pp.nam %>% 
    filter(p.adjust<0.05) %>% 
    arrange(desc(p.adjust)) %>% 
    ggplot(aes(x = log2(odd),
               y = reorder(term, log2(odd)),
               size = Count,
               color = log10(p.adjust))) +
    geom_point() +
    labs(x=expression(bold(paste("log"[2], "(Enrichment)"))), 
         y="",
         color=expression(bold(paste("FDR log"[10], "(", italic(p),"-value)"))),
         size="Candidate gene number",
         title=expression(bold(paste("All CNVs", bold(italic(" D. pulex ")), "North America")))) +
    theme_bw() +
    scale_color_viridis(option = "plasma") +
    theme(title = element_text(face="bold", size=18),
          legend.text = element_text(face="bold", size=18),
          legend.title = element_text(face="bold", size=18),
          legend.background = element_blank(),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18))}

# Output for Revigo
write.table(pp.euro %>% 
              select(ID, p.adjust), 
            file = "cnvs/clusterprofiler.euro.private.txt",
            quote = F, row.names = F, col.names = F, sep = "\t")

write.table(pp.nam %>% 
              select(ID, p.adjust), 
            file = "cnvs/clusterprofiler.nam.private.txt",
            quote = F, row.names = F, col.names = F, sep = "\t")

pdf("figures/goterm_CNVs_euro_filtprivate.pdf", width=16, height=10)
plot.go.euro
dev.off()

pdf("figures/goterm_CNVs_nam_filtprivate.pdf", width=20, height=10)
plot.go.nam
dev.off()
