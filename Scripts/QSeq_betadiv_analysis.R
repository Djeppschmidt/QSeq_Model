# load libraries
library(vegan)
library(phyloseq)
library(QSeq)

# set working directory to github directory
# setwd("~/Documents/GitHub/QSeq_Model")
# import datasets (GLU, FSP, TFW)
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

# FSP<-readRDS("SnakeMake_pipeline/DataFiles/projects/GLU/comdat.RDS")

# for beta diversity analysis:
# Do permanova for experimental categories, and each environmental variable individually
# make barchart of % var explained with and without QSeq

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

# check variance of total abundance for each:

sd(sample_data(GLU)$QPCR_16s)/mean(sample_data(GLU)$QPCR_16s)
sd(sample_data(FSP)$Bac_QPCR)/mean(sample_data(FSP)$Bac_QPCR)
sd(sample_data(TFW)$QPCR_16S)/mean(sample_data(TFW)$QPCR_16S)


# output otu table
totu<-as.data.frame(t(as.matrix(otu_table(Q.FSP))))


# make sure otu table is in the correct orientation

# run permanova
tout<-adonis2(as.matrix(totu)~as.factor(as.character(sample_data(FSP)$Treatment)),
             #strata=sample_data(ps2)[[strata]], 
             method="jaccard")
adonis2(as.matrix(totu)~as.factor(as.character(sample_data(FSP)$Treatment)),
              #strata=sample_data(ps2)[[strata]], 
              method="bray")
# not published:
# make list of factors in metadata for each dataset
# cov.list<-c("pH.H2O", "C_org", "NH4_N", "NO3_N")
# FSP.cov<-c("pH", "C_percent", "Nh4_ugPerg", "No3_ugPerg")
# TFW.cov<-c("PH", "TOTALC", "TOTALN", "WATER")

# make sure no NA exist in sample data !!!!! Also prepare filtered and aggregated data
GLU<-prune_samples(!is.na(sample_data(GLU)$pH.H2O) | !is.na(sample_data(GLU)$C_org) | !is.na(sample_data(GLU)$NH4_N) | !is.na(sample_data(GLU)$NO3_N), GLU)
GLU.t<-tax_glom(GLU, taxrank="Genus")
GLU.f<-filter_taxa(GLU, function(x) mean(x)>1, prune=T)
GLU.tf<-filter_taxa(GLU.t, function(x) mean(x)>1, prune=T)

FSP<-prune_samples(!is.na(sample_data(FSP)$pH) | !is.na(sample_data(FSP)$C_percent) | !is.na(sample_data(FSP)$Nh4_ugPerg) | !is.na(sample_data(FSP)$No3_ugPerg), FSP)
FSP.t<-tax_glom(FSP, taxrank="Genus")
FSP.f<-filter_taxa(FSP, function(x) mean(x)>1, prune=T)
FSP.tf<-filter_taxa(FSP.t, function(x) mean(x)>1, prune=T)

TFW<-prune_samples(!is.na(sample_data(TFW)$PH) | !is.na(sample_data(TFW)$TOTALC) | !is.na(sample_data(TFW)$TOTALN) | !is.na(sample_data(TFW)$WATER), TFW)
TFW.t<-tax_glom(TFW, taxrank="Genus")
TFW.f<-filter_taxa(TFW, function(x) mean(x)>1, prune=T)
TFW.tf<-filter_taxa(TFW.t, function(x) mean(x)>1, prune=T)

# Qseq transformation
Q.FSP<-QSeq(FSP, "Bac_QPCR")
Q.GLU<-QSeq(GLU, "QPCR_16s")
Q.TFW<-QSeq(TFW, "QPCR_16S")

# run permanova with different normalization and distance matrix:
# FSP permanova ####
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

# GLUSEEN permanovas ####
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

# TFW permanova ####
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

#QSeq:
run.adonis(Q.TFW, col="SITE", type="factor", method="jaccard")
run.adonis(Q.TFW, col="SITE", type="factor", method="bray")
run.adonis(Q.TFW, col="PLANT", type="factor", method="jaccard")
run.adonis(Q.TFW, col="PLANT", type="factor", method="bray")

# Figure 3####

# get distribution of 16S QPCR measure of abundance across datasets
density.plot<-data.frame("dataset"=c(rep("TFW", nsamples(TFW)), rep("GLUSEEN", nsamples(GLU)), rep("FSP", nsamples(FSP))), "Abundance"=c(as.numeric(sample_data(TFW)$QPCR_16S), as.numeric(sample_data(GLU)$QPCR_16s),as.numeric(sample_data(FSP)$Bac_QPCR)))

# plot:
ggplot(density.plot, aes(x = log(Abundance), fill = dataset)) + geom_density(alpha = 0.5) + theme_bw() + scale_fill_viridis_d() #+ facet_wrap(~dataset)

#plot(density(log(sample_data(GLU)$QPCR_16s)))
#lines(density(log(sample_data(FSP)$Bac_QPCR)))
#lines(density(log(sample_data(TFW)$QPCR_16S)))

# statistically evaluate homogeneity of variances:
leveneTest(Abundance ~ dataset, data = density.plot)



# remove all these if not referenced for analysis:
# data exploration (not in paper) ####


# rarefaction function
#ps.rarefy<-function(ps, sample){
#  ps2<-ps
#  otu<-as.data.frame(as.matrix(otu_table(ps)))
#  # make sure otu table is in the correct orientation
#  if(taxa_are_rows(ps2)){otu <- t(otu)}
#  otu<-rrarefy(otu, sample)
#  print(dim(otu))
#  otu_table(ps2)<-otu_table(otu, taxa_are_rows = F)
#  ps2
#}

# prepare data frame of all combinations of environmental factors for each dataset:
#combinations<-gtools::permutations(4, 4, cov.list)
#FSP.comb<-gtools::permutations(4, 4, FSP.cov)
#TFW.comb<-gtools::permutations(4, 4, TFW.cov)

# Run permuted model permanova test
# look at how each level linteracts
#       rarefaction to different depths * aggregation * filter taxa * (aggregate + filter)


# gene = name of column in metadata for QPCR
# x = phyloseq object
# env.factor = data.frame with all combinations / orders of environmental variables to be tested
# title = title for plot
# rare.val = value for average rarefaction depth
compare.permanova<-function(x, gene,env.factors, title, rare.val){
  Qx<-QSeq(x, abundance=gene)
  rx<-ps.rarefy(x, round(rare.val*(sample_data(x)[[gene]]/mean(sample_data(x)[[gene]]))))
  
  # run permanova
  permanova<-run.adonis(x, col=env.factors, type = "cov")
  Q.permanova<-run.adonis(Qx, col=env.factors, type = "cov")
  r.permanova<-run.adonis(rx, col=env.factors, type = "cov")
  
  # Chi square analysis:
  chi1<-data.frame("value"=c(Q.permanova$R[1,]>permanova$R[1,],Q.permanova$R[1,]<permanova$R[1,]),"category"=c(rep("QSeq", length(permanova$R[1,])), rep("None", length(permanova$R[1,]))))
  chi2<-data.frame("value"=c(Q.permanova$R[2,]>permanova$R[2,],Q.permanova$R[2,]<permanova$R[2,]),"category"=c(rep("QSeq", length(permanova$R[2,])), rep("None", length(permanova$R[2,]))))
  chi3<-data.frame("value"=c(Q.permanova$R[3,]>permanova$R[3,],Q.permanova$R[3,]<permanova$R[3,]),"category"=c(rep("QSeq", length(permanova$R[3,])), rep("None", length(permanova$R[3,]))))
  chi4<-data.frame("value"=c(Q.permanova$R[4,]>permanova$R[4,],Q.permanova$R[4,]<permanova$R[4,]),"category"=c(rep("QSeq", length(permanova$R[4,])), rep("None", length(permanova$R[4,]))))
  print(chi1)
  print(rownames(permanova$R)[1])
  print(chisq.test(chi1$value, chi1$category))
  print(rownames(permanova$R)[2])
  print(chisq.test(chi2$value, chi2$category))
  print(rownames(permanova$R)[3])
  print(chisq.test(chi3$value,chi3$category))
  print(rownames(permanova$R)[4])
  print(chisq.test(chi4$value,chi4$category))
  
  # compile runs for plotting
  pdf<-data.frame("Normalization"=c(rep("none", 4), rep("QSeq",4), rep("rarefaction", 4)),"Factor"=c(rep(rownames(permanova$R),3)), "mean"=c(rep(NA, 12)), "SD"=c(rep(NA,12)))
  pdf$mean<-c(rowSums(permanova$R)/ncol(permanova$R),rowSums(Q.permanova$R)/ncol(Q.permanova$R),rowSums(r.permanova$R)/ncol(r.permanova$R))
  pdf$SD<-c(apply(permanova$R, 1, sd),apply(Q.permanova$R, 1, sd),apply(r.permanova$R, 1, sd))
  
  print(ggplot(pdf, aes(x=Factor, y=mean, fill=Normalization))+
    geom_bar(stat="identity", position = position_dodge())+
    geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), position = position_dodge(0.9), width=0.4)+
    scale_fill_viridis_d()+
    ggtitle(title)+
    theme_bw())
  
  return(pdf)
}

# gene = name of column in metadata for QPCR
# x = phyloseq object
# factor = name of column of factor to be tested
# title = title for plot
# rare.val = value for average rarefaction depth
compare.permanova2<-function(x, gene, factor, title, rare.val){
  Qx<-QSeq(x, abundance=gene)
  
  # run permanova
  permanova<-run.adonis(x, col=factor, type = "factor")
  Q.permanova<-run.adonis(Qx, col=factor, type = "factor")
  
  # compile runs for plotting
  pdf<-data.frame("Normalization"=c(rep("none", 4), rep("QSeq",4), rep("rarefaction", 4)),"Factor"=c(rep(rownames(permanova$R),3)), "mean"=c(rep(NA, 12)), "SD"=c(rep(NA,12)))
  pdf$mean<-c(rowSums(permanova$R)/ncol(permanova$R),rowSums(Q.permanova$R)/ncol(Q.permanova$R),rowSums(r.permanova$R)/ncol(r.permanova$R))
  pdf$SD<-c(apply(permanova$R, 1, sd),apply(Q.permanova$R, 1, sd),apply(r.permanova$R, 1, sd))
  
  print(ggplot(pdf, aes(x=Factor, y=mean, fill=Normalization))+
          geom_bar(stat="identity", position = position_dodge())+
          geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), position = position_dodge(0.9), width=0.4)+
          scale_fill_viridis_d()+
          ggtitle(title)+
          theme_bw())
  
  return(pdf)
}
# run data ####
library(ggplot2)
# fsp
perm.fsp.50k<-compare.permanova(FSP, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP 50K", rare.val=50000)
perm.fsp.25k<-compare.permanova(FSP, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP 25K", rare.val=25000)
perm.fsp.10k<-compare.permanova(FSP, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP 10K", rare.val=10000)

perm.fsp.t50K<-compare.permanova(FSP.t, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP t50K", rare.val=50000)
perm.fsp.t25K<-compare.permanova(FSP.t, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP t25K", rare.val=25000)
perm.fsp.t10K<-compare.permanova(FSP.t, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP t10K", rare.val=10000)

perm.fsp.f50K<-compare.permanova(FSP.f, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP f50K", rare.val=50000)
perm.fsp.f25K<-compare.permanova(FSP.f, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP f25K", rare.val=25000)
perm.fsp.f10K<-compare.permanova(FSP.f, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP f10K", rare.val=10000)

perm.fsp.tf50K<-compare.permanova(FSP.tf, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP ft50K", rare.val=50000)
perm.fsp.tf25K<-compare.permanova(FSP.tf, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP ft25K", rare.val=25000)
perm.fsp.tf10K<-compare.permanova(FSP.tf, gene="Bac_QPCR" , env.factors=FSP.comb, title="FSP ft10K", rare.val=10000)


# GLUSEEN
perm.glu.50k<-compare.permanova(GLU, gene="QPCR_16s" , env.factors=combinations, title="GLU 50K", rare.val=50000)
perm.glu.25k<-compare.permanova(GLU, gene="QPCR_16s" , env.factors=combinations, title="GLU 25K", rare.val=25000)
perm.glu.10k<-compare.permanova(GLU, gene="QPCR_16s" , env.factors=combinations, title="GLU 10K", rare.val=10000)



perm.glu.t50K<-compare.permanova(GLU.t, gene="QPCR_16s" , env.factors=combinations, title="GLU t50K", rare.val=50000)
perm.glu.t25K<-compare.permanova(GLU.t, gene="QPCR_16s" , env.factors=combinations, title="GLU t25K", rare.val=25000)
perm.glu.t10K<-compare.permanova(GLU.t, gene="QPCR_16s" , env.factors=combinations, title="GLU t10K", rare.val=10000)


perm.glu.f50K<-compare.permanova(GLU.f, gene="QPCR_16s" , env.factors=combinations, title="GLU f50K", rare.val=50000)
perm.glu.f25K<-compare.permanova(GLU.f, gene="QPCR_16s" , env.factors=combinations, title="GLU f25K", rare.val=25000)
perm.glu.f10K<-compare.permanova(GLU.f, gene="QPCR_16s" , env.factors=combinations, title="GLU f10K", rare.val=10000)

perm.glu.tf50K<-compare.permanova(GLU.tf, gene="QPCR_16s" , env.factors=combinations, title="GLU ft50K", rare.val=50000)
perm.glu.tf25K<-compare.permanova(GLU.tf, gene="QPCR_16s" , env.factors=combinations, title="GLU ft25K", rare.val=25000)
perm.glu.tf10K<-compare.permanova(GLU.tf, gene="QPCR_16s" , env.factors=combinations, title="GLU ft10K", rare.val=10000)


# TFW

perm.TFW.1.5<-compare.permanova(TFW, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW 1.5K", rare.val=1500)
perm.TFW.1<-compare.permanova(TFW, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW 1K", rare.val=1000)
perm.TFW.0.5<-compare.permanova(TFW, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW 0.5K", rare.val=500)

perm.TFW.t1.5<-compare.permanova(TFW.t, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW t1.5", rare.val=1500)
perm.TFW.t1<-compare.permanova(TFW.t, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW t1", rare.val=1000)
perm.TFW.t0.5<-compare.permanova(TFW.t, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW t0.5", rare.val=500)

perm.TFW.f1.5<-compare.permanova(TFW.f, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW f1.5", rare.val=1500)
perm.TFW.f1<-compare.permanova(TFW.f, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW f1", rare.val=1000)
perm.TFW.f0.5<-compare.permanova(TFW.f, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW f0.5", rare.val=500)

perm.TFW.tf1.5<-compare.permanova(TFW.tf, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW tf1.5", rare.val=1500)
perm.TFW.tf1<-compare.permanova(TFW.tf, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW tf1", rare.val=1000)
perm.TFW.tf0.5<-compare.permanova(TFW.tf, gene="QPCR_16S" , env.factors=TFW.comb, title="TFW tf0.5", rare.val=500)


# combine plots for each dataset to make line plots for each factor:
# FSP

FSP.data<-list(perm.fsp.50k, perm.fsp.25k, perm.fsp.10k, perm.fsp.t50K, perm.fsp.t25K, perm.fsp.t10K, perm.fsp.f50K, perm.fsp.f25K, perm.fsp.f10K, perm.fsp.tf50K, perm.fsp.tf25K, perm.fsp.tf10K)

FSP.covlist<-NULL
for(i in FSP.cov){
  FSP.covlist[[i]]<-do.call(rbind, lapply(FSP.data, function(x) x[x$Factor==i,]))
  FSP.covlist[[i]]$Rarefaction<-c(rep(c(rep(50000, 3), rep(25000,3), rep(10000, 3)),4))
  FSP.covlist[[i]]$Filter<-c(rep("none", 9), rep("filtered",9), rep("none", 9), rep("filtered", 9))
  FSP.covlist[[i]]$AggregateTax<-c(rep("none", 18), rep("genus", 18))
}

for(i in FSP.cov){
  print(ggplot(FSP.covlist[[i]], aes(x=as.factor(Rarefaction), y=mean, color=Normalization))+
    geom_point(position = position_dodge(width = 0.9))+
    geom_errorbar(aes(ymin=mean-SD,ymax=mean+SD), position = "dodge")+
    facet_wrap(~Filter+AggregateTax)+
    theme_bw()+
    ggtitle(unique(FSP.covlist[[i]]$Factor))+
    ylab("Variance Explained")+
    xlab("Rarefaction Depth")+
    scale_fill_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")))
}

# GLU
GLU.data<-list(perm.glu.50k, perm.glu.25k, perm.glu.10k, perm.glu.t50K, perm.glu.t25K, perm.glu.t10K, perm.glu.f50K, perm.glu.f25K, perm.glu.f10K, perm.glu.tf50K, perm.glu.tf25K, perm.glu.tf10K)

GLU.covlist<-NULL
for(i in cov.list){
  GLU.covlist[[i]]<-do.call(rbind, lapply(GLU.data, function(x) x[x$Factor==i,]))
  GLU.covlist[[i]]$Rarefaction<-c(rep(c(rep(50000, 3), rep(25000,3), rep(10000, 3)),4))
  GLU.covlist[[i]]$Filter<-c(rep("none", 9), rep("filtered",9), rep("none", 9), rep("filtered", 9))
  GLU.covlist[[i]]$AggregateTax<-c(rep("none", 18), rep("genus", 18))
}

for(i in cov.list){
  print(ggplot(GLU.covlist[[i]], aes(x=as.factor(Rarefaction), y=mean, color=Normalization))+
          geom_point(position = position_dodge(width = 0.9))+
          geom_errorbar(aes(ymin=mean-SD,ymax=mean+SD), position = "dodge")+
          facet_wrap(~Filter+AggregateTax)+
          theme_bw()+
          ggtitle(unique(GLU.covlist[[i]]$Factor))+
          ylab("Variance Explained")+
          xlab("Rarefaction Depth")+
          scale_fill_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
          scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")))
}

# TFW
TFW.data<-list(perm.TFW.1.5, perm.TFW.1, perm.TFW.0.5, perm.TFW.t1.5, perm.TFW.t1, perm.TFW.t0.5, perm.TFW.f1.5, perm.TFW.f1, perm.TFW.f0.5, perm.TFW.tf1.5, perm.TFW.tf1, perm.TFW.tf0.5)

TFW.covlist<-NULL
for(i in TFW.cov){
  TFW.covlist[[i]]<-do.call(rbind, lapply(TFW.data, function(x) x[x$Factor==i,]))
  TFW.covlist[[i]]$Rarefaction<-c(rep(c(rep(50000, 3), rep(25000,3), rep(10000, 3)),4))
  TFW.covlist[[i]]$Filter<-c(rep("none", 9), rep("filtered",9), rep("none", 9), rep("filtered", 9))
  TFW.covlist[[i]]$AggregateTax<-c(rep("none", 18), rep("genus", 18))
}

for(i in TFW.cov){
  print(ggplot(TFW.covlist[[i]], aes(x=as.factor(Rarefaction), y=mean, color=Normalization))+
          geom_point(position = position_dodge(width = 0.9))+
          geom_errorbar(aes(ymin=mean-SD,ymax=mean+SD), position = "dodge")+
          facet_wrap(~Filter+AggregateTax)+
          theme_bw()+
          ggtitle(unique(TFW.covlist[[i]]$Factor))+
          ylab("Variance Explained")+
          xlab("Rarefaction Depth")+
          scale_fill_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
          scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")))
}

# make array of factors: Effects of rarefaction, filtering, aggregating 

# exploration ####
run.adonis(GLU.r, cov.list, strata="Cities", type="cov")
run.adonis(tax_glom(GLU.r, taxrank = "Genus"), cov.list, strata="Cities", type="cov")
run.adonis(Q.GLU.r, cov.list, strata="Cities", type="cov")

# aggregated to genus level
run.adonis(tax_glom(GLU.r, taxrank = "Genus"), cov.list, strata="Cities", type="cov") # rarefied
run.adonis(tax_glom(Q.GLU.r, taxrank = "Genus"), cov.list, strata="Cities", type="cov") # Qseq + rarefied
run.adonis(tax_glom(GLU, taxrank = "Genus"), cov.list, strata="Cities", type="cov") # raw
run.adonis(tax_glom(Q.GLU, taxrank = "Genus"), cov.list, strata="Cities", type="cov") # Qseq

# by asv
run.adonis(GLU.r, cov.list, strata="Cities", type="cov") # rarefied
run.adonis(Q.GLU.r, cov.list, strata="Cities", type="cov") # Qseq + rarefied
run.adonis(GLU, cov.list, strata="Cities", type="cov") # raw
run.adonis(Q.GLU, cov.list, strata="Cities", type="cov") # Qseq

fac.list<-c("Codes")
run.adonis(GLU.r, fac.list, strata="Cities", type="factor") # rarefied
run.adonis(Q.GLU.r, fac.list, strata="Cities", type="factor") # Qseq + rarefied
run.adonis(GLU, fac.list, strata="Cities", type="factor") # raw
run.adonis(Q.GLU, fac.list, strata="Cities", type="factor") # Qseq

run.adonis(tax_glom(GLU, taxrank = "Genus"), fac.list, strata="Cities", type="factor") # raw
run.adonis(tax_glom(Q.GLU, taxrank = "Genus"), fac.list, strata="Cities", type="factor") # Qseq

# comment out strata:
sample_data(GLU.r)$QPCR_16s
run.adonis(tax_glom(GLU.r, taxrank = "Genus"), fac.list, strata=NULL, type="factor")
run.adonis(tax_glom(Q.GLU.r, taxrank = "Genus"), fac.list, strata=NULL, type="factor")

#run.adonis(tax_glom(transform_sample_counts(GLU.r, function(x) x/sum(x)), taxrank = "Genus"), cov.list, strata=NULL, type="cov")
run.adonis(tax_glom(GLU, taxrank = "Genus"), fac.list, strata=NULL, type="factor")
run.adonis(tax_glom(Q.GLU, taxrank = "Genus"), fac.list, strata=NULL, type="factor")

#glu.otu<-as.data.frame(as.matrix(otu_table(GLU)))
#adonis(glu.adonis~as.numeric(as.character(sample_data(GLU)$C_org)))


# calculate percent of taxa that have different response in QSeq vs normal analysis for each dataset