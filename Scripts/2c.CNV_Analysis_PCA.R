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

### CNVR ###

# Read in cnv info
files1 <- system("ls -f -R /project/berglandlab/connor/cnvs/*cnvr.csv", intern = TRUE)
listi1 <- lapply(files1, fread)
setattr(listi1, 'names', list.files(path = "/project/berglandlab/connor/cnvs", pattern = "cnvr.csv"))

# Bind list & make unique CNV
dt.cnvr <- data.table(fread("/project/berglandlab/connor/cnvs/Daphnia.pulex.Europecnvr.csv"))

dt.t <- data.table(t(dt.cnvr %>% dplyr::select(!colnames(dt.cnvr)[1:6]))) %>% 
                     mutate(samp=colnames(dt.cnvr)[7:dim(dt.cnvr)[2]])

# Wide to long
rows.name <- dt.cnvr %>% 
  dplyr::select(colnames(dt.cnvr)[2:4]) %>% 
  summarize(name=paste(seqnames, start, end, sep="_"))
colnames(dt.t) <- c(rows.name$name, "samp")

# Make CNVs numeric
dt.t <- data.table(dt.t %>% 
                     pivot_longer(cols=c(colnames(dt.t)[1:dim(dt.t)[2]-1]),
                                  names_to="name") %>% 
                     mutate(value=as.numeric(str_replace(value, pattern = "CN", replacement = ""))) %>% 
                     pivot_wider(names_from = "name"))

# PCA of CNVs
pca <- prcomp(dt.t %>% dplyr::select(!samp))

pca.summ <- summary(pca)

# Merge by Sample
pca1 <- data.table(data.table(sample=dt.t$samp, PC1=pca$x[,1],
                  PC2=pca$x[,2], PC3=pca$x[,3],
                  PC4=pca$x[,4], PC5=pca$x[,5],
                  PC6=pca$x[,6], PC7=pca$x[,7], 
                  PC8=pca$x[,8]) %>% 
  mutate(Sample=gsub("\\RG", "",
                gsub("\\.bam", "",
                gsub("\\mdup.bam", "", 
                gsub("\\_finalmap_*", "",sample))))))

# Remove missing samples
pca1 <- data.table(pca1[Sample %in% fin$Sample])

# Merge PCA and metadata
pca1 <- data.table(merge(pca1, fin, by="Sample") %>% 
           mutate(Species=case_when(cont=="Daphnia.pulexcaria.NorthAmerica"~
                                          "D.pulex nam x D.pulicaria nam",
                                    Species=="Daphnia pulex" ~ "D.pulex",
                                    Species=="Daphnia pulicaria" ~ "D.pulicaria",
                                    TRUE ~ Species)))

# PC 1/2 Plot
library(ggforce)

# Remove Unknown SCs
pca1 <- pca1[!SC=="OO"][!pondID==""]

pca12 <- {pca1 %>%
  ggplot(., 
         aes(x=PC1, 
             y=PC2, 
             fill=pondID), 
             shape=) +
  geom_point(size=6, alpha=0.8) +
  geom_mark_ellipse() +
  theme_minimal() + 
  labs(x=paste("PC1 (", round(pca.summ$importance[2,1], digits=3)*100, " %)", sep=""),
       y=paste("PC2 (", round(pca.summ$importance[2,2], digits=3)*100, " %)", sep="")) +
  #scale_fill_continuous(name = "Species complex") +
  theme(strip.text = element_text(face="bold.italic", size=16),
        legend.text = element_text(size=18, face="bold.italic"),
        legend.title = element_text(face="bold", size=20),
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.title = element_text(face="bold", size=20))}

# PCA 1/2
pdf("figures/pca.cnv.euro.pond.pdf")
pca12
dev.off()

# SKree plot
pdf("figures/pca.skree.cnv.pdf")
plot(pca)
dev.off()
