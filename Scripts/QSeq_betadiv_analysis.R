# load libraries
library(vegan)
library(phyloseq)
library(QSeq)
library(ggplot2)
library(viridis)


# set working directory to github directory
# setwd("~/Documents/GitHub/QSeq_Model")

# import datasets (GLU, FSP, TFW) ####
GLU<-readRDS("Data/GLU_comdat.RDS")
#sample_data(GLU)$effort<-sample_sums(GLU)/sample_data(GLU)$QPCR_16s
#sample_data(GLU)$seqdepth<-sample_sums(GLU)

# move copy into project directory
FSP<-readRDS("Data/FSP_comdat.rds")

# import TFW, prepare meta data
TFW<-readRDS("Data/TFW_comdat.RDS")
meta<-read.csv("Data/TFW_metadata.csv")
rownames(meta)<-meta$ID
sample_data(TFW)<-sample_data(meta)
sample_names(TFW)<-paste("Sample_", sample_names(TFW), sep="")

# function to run beta diveristy method:
# ps = phyloseq object
# col = list of sample metadata column names that will be used for permanova:
# type = "cov" or "factor"
# strata = blocking factor
# method = type of distance matrix to implement
run.adonis<-function(ps, col, type, method, strata=NULL){
  ps2<-ps
  # remove taxa and samples that do not have any data:
  if(any(taxa_sums(ps2)==0)){ps2<-prune_taxa(taxa_sums(ps2)>0, ps2)}
  if(any(sample_sums(ps2)==0)){ps2<-prune_samples(sample_sums(ps2)>0, ps2)}
  print(ps2)
  # make array for output
  #out<-data.frame("Factors"=col, "R2"=rep(NA, length(col)), "P_val"=rep(NA, length(col)))

  # for each model, run permanova, extract R SQ values, put in order, output data frame with one column for each model design
  out<-NULL
  
  
  
  if(type=="cov"){
    out$P<-matrix(ncol=nrow(col), nrow=ncol(col))
    out$R<-matrix(ncol=nrow(col), nrow=ncol(col))
    print(col[1,])
    print(col)
    print(dim(out$P))
    rownames(out$P)<-col[1,]
    rownames(out$R)<-col[1,]
    otu<-as.data.frame(as.matrix(otu_table(ps2)))
    # make sure otu table is in the correct orientation
    if(taxa_are_rows(ps2)){otu <- t(otu)}
    for(i in 1:nrow(col)){
      # subset to samples with metadata
      # ps2<-prune_samples(!is.na(sample_data(ps)[[i]]), ps)
      # output otu table
      
      # run permanova
      #perm<-adonis(otu~as.numeric(as.character(sample_data(ps2)[[i]])), 
      perm<-adonis(otu~as.numeric(as.character(sample_data(ps2)[[col[i,1]]]))+as.numeric(as.character(sample_data(ps2)[[col[i,2]]]))+as.numeric(as.character(sample_data(ps2)[[col[i,3]]]))+as.numeric(as.character(sample_data(ps2)[[col[i,4]]])),
                   #strata=sample_data(ps2)[[strata]], 
                   method="jaccard")
      # extract P value and R squared
      mtch<-match(col[i,], rownames(out$R))
      print(col[i,])
      print(rownames(out$R))
      print(mtch)
      print(perm$aov.tab$R2)
      out$R[,i]<-perm$aov.tab$R2[1:ncol(col)][mtch]
      out$P[,i]<-perm$aov.tab$`Pr(>F)`[1:ncol(col)][mtch]
    }
  }
  # ignore factor for now
  if(type=="factor"){
    for(i in col){
      # output otu table
      otu<-as.data.frame(as.matrix(otu_table(ps2)))
      print(dim(otu))
      # make sure otu table is in the correct orientation
      if(taxa_are_rows(ps2)){otu <- as.data.frame(t(as.matrix(otu_table(ps2))))}
      # run permanova
      out<-adonis2(as.matrix(otu)~as.factor(as.character(sample_data(ps2)[[col]])),
                   #strata=sample_data(ps2)[[strata]], 
                   method=method)
      # extract P value and R squared
      
    }
  }
  
  out
}

# total abundance relative standard deviation: ####

sd(sample_data(GLU)$QPCR_16s)/mean(sample_data(GLU)$QPCR_16s) # 0.75
sd(sample_data(FSP)$Bac_QPCR)/mean(sample_data(FSP)$Bac_QPCR) # 1.28
sd(sample_data(TFW)$QPCR_16S)/mean(sample_data(TFW)$QPCR_16S) # 0.98

# Qseq transformation ####
Q.FSP<-QSeq(FSP, "Bac_QPCR")
Q.GLU<-QSeq(GLU, "QPCR_16s")
Q.TFW<-QSeq(TFW, "QPCR_16S")

# run permanovas: data for Table 3 ####
# FSP permanova
# relative abundance:
run.adonis(transform_sample_counts(FSP, function(x) x/sum(x)), col="Treatment", type="factor", method="jaccard") 
run.adonis(transform_sample_counts(FSP, function(x) x/sum(x)), col="Treatment", type="factor", method="bray")
run.adonis(transform_sample_counts(FSP, function(x) x/sum(x)), col="Depth", type="factor", method="jaccard")
run.adonis(transform_sample_counts(FSP, function(x) x/sum(x)), col="Depth", type="factor", method="bray")

# aitchison distance
run.adonis(FSP, col="Treatment", type="factor", method="robust.aitchison")
run.adonis(transform_sample_counts(FSP, function(x) x+1), col="Treatment", type="factor", method="aitchison")
run.adonis(FSP, col="Depth", type="factor", method="robust.aitchison")
run.adonis(transform_sample_counts(FSP, function(x) x+1), col="Depth", type="factor", method="aitchison")

# QSeq
run.adonis(Q.FSP, col="Treatment", type="factor", method="jaccard")
run.adonis(Q.FSP, col="Treatment", type="factor", method="bray")
run.adonis(Q.FSP, col="Depth", type="factor", method="jaccard")
run.adonis(Q.FSP, col="Depth", type="factor", method="bray")

# GLUSEEN
# relative abundance:
run.adonis(transform_sample_counts(GLU, function(x) x/sum(x)), col="Cities", type="factor", method="jaccard")
run.adonis(transform_sample_counts(GLU, function(x) x/sum(x)), col="Cities", type="factor", method="bray")
run.adonis(transform_sample_counts(GLU, function(x) x/sum(x)), col="Codes", type="factor", method="jaccard")
run.adonis(transform_sample_counts(GLU, function(x) x/sum(x)), col="Codes", type="factor", method="bray")

# aitchison:
run.adonis(GLU, col="Cities", type="factor", method="robust.aitchison")
run.adonis(transform_sample_counts(GLU, function(x) x+1), col="Cities", type="factor", method="aitchison")
run.adonis(GLU, col="Codes", type="factor", method="robust.aitchison")
run.adonis(transform_sample_counts(GLU, function(x) x+1), col="Codes", type="factor", method="aitchison")

# QSeq:
run.adonis(Q.GLU, col="Cities", type="factor", method="jaccard")
run.adonis(Q.GLU, col="Cities", type="factor", method="bray")
run.adonis(Q.GLU, col="Codes", type="factor", method="jaccard")
run.adonis(Q.GLU, col="Codes", type="factor", method="bray")

# TFW
# relative abundance:
run.adonis(transform_sample_counts(TFW, function(x) x/sum(x)), col="SITE", type="factor", method="jaccard")
run.adonis(transform_sample_counts(TFW, function(x) x/sum(x)), col="SITE", type="factor", method="bray")
run.adonis(transform_sample_counts(TFW, function(x) x/sum(x)), col="PLANT", type="factor", method="jaccard")
run.adonis(transform_sample_counts(TFW, function(x) x/sum(x)), col="PLANT", type="factor", method="bray")

# aitchison
run.adonis(TFW, col="SITE", type="factor", method="robust.aitchison")
run.adonis(transform_sample_counts(TFW, function(x) x+1), col="SITE", type="factor", method="aitchison")
run.adonis(TFW, col="PLANT", type="factor", method="robust.aitchison")
run.adonis(transform_sample_counts(TFW, function(x) x+1), col="PLANT", type="factor", method="aitchison")

# QSeq:
run.adonis(Q.TFW, col="SITE", type="factor", method="jaccard")
run.adonis(Q.TFW, col="SITE", type="factor", method="bray")
run.adonis(Q.TFW, col="PLANT", type="factor", method="jaccard")
run.adonis(Q.TFW, col="PLANT", type="factor", method="bray")

# Figure 3####
# make data frame of abundance data across datasets:
density.plot<-data.frame("dataset"=c(rep("TFW", nsamples(TFW)), rep("GLUSEEN", nsamples(GLU)), rep("FSP", nsamples(FSP))), "Abundance"=c(as.numeric(sample_data(TFW)$QPCR_16S), as.numeric(sample_data(GLU)$QPCR_16s),as.numeric(sample_data(FSP)$Bac_QPCR)))

# plot:
ggplot(density.plot, aes(x = log(Abundance), fill = dataset)) + geom_density(alpha = 0.5) + theme_bw() + scale_fill_viridis_d()

ggplot(density.plot, aes(x = log(Abundance), fill = dataset)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C"))

# statistically evaluate homogeneity of variances:
leveneTest(Abundance ~ dataset, data = density.plot)

