# function to initialize data
# make list of matrices, with and without an overlying gradient
# 
# Prepare workspace ####
library(ALDEx2)
library(phyloseq)
library(corrplot)
library(ANCOMBC)
library(ggplot2)

# will need to install FastSpar independently !!!!

# functions for analysis ####

# implements FastSpar on community data from a phyloseq object
# FastSpar must be installed. You will need to know where it is installed to direct the function to implement it.
#' subsample community
#' @param ps phyloseq object (required)
#' @param name name of file to save for input to FastSpar (required)
#' @param threads number of cores to implement FastSpar (required)
#' @param path path to local installation of FastSpar. (required) e.g. ("Users/me/opt/anaconda3/bin/")
#' @export
#' @examples
#' run.FastSpar()
run.FastSpar<-function(ps, name, threads=10, path){
  otu<-as.data.frame(as.matrix(otu_table(ps)))
  otu<-data.frame("OTU ID"=paste("Taxon", c(1:nrow(otu)), sep=""), otu)
  colnames<-colnames(otu)
  otu<-rbind(colnames, otu)
  otu[1,1]<-"#OTU ID"
  utils::write.table(otu, file=name, quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
  
  # save path to local installation of R 
  bin<-Sys.getenv("PATH")
  
  # set path to local installation of FastSpar to implement run
  Sys.setenv(PATH=path)
  
  # run FastSpar correlations with multithreading.
  system(paste0("fastspar --otu_table ", name, " --correlation median_correlation.tsv --covariance median_covariance.tsv --threads ", threads)) 
  dir.create("bootstrap_counts/")
  system(paste0("fastspar_bootstrap --otu_table ", name, " --number 1000 --prefix bootstrap_counts/bootstrap"))
  dir.create("bootstrap_correlation/")
  
  index<-c(0:999)
  index<-as.character(index)
  parallel::mclapply(index, function(x) system(paste0("echo y | fastspar --otu_table ", "bootstrap_counts/bootstrap_", x, ".tsv ", "--correlation bootstrap_correlation/cor_", x, ".tsv ", "--covariance bootstrap_correlation/cov_", x, ".tsv ", "-i 5")),mc.preschedule=F,mc.cores=threads,mc.cleanup = TRUE)
  
  system(paste0("fastspar_pvalues --otu_table ", name, " --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_ --permutations 1000 --outfile pvalues.tsv"))
  
  # return PATH to local installation of R so packages keep working
  Sys.setenv(PATH=bin)
  
  # import results of FastSpar to R
  out<-NULL
  out$cor<-as.matrix(read.table("median_correlation.tsv", header = FALSE, row.names = 1, sep = "\t"))
  colnames(out$cor)<-rownames(out$cor)
  out$pval<-as.matrix(read.table("pvalues.tsv", header = FALSE, row.names = 1, sep = "\t"))
  colnames(out$pval)<-rownames(out$pval)
  out
  
}

# implements aldex2 on community data from a phyloseq object
# 
#' subsample community
#' @param ps phyloseq object (required)
#' @export
#' @examples
#' run.FastSparl()
runALDEx<-function(ps){
  require(ALDEx2)
  x <- aldex.clr(as.data.frame(as.matrix(otu_table(ps))), sample_data(ps)$Group, mc.samples=128, denom="all", verbose=FALSE)
  out<-NULL
  out$effect<-aldex.effect(x, CI=TRUE, verbose=FALSE)
  out$test<-aldex.kw(x)
  rownames(out$test)<-taxa_names(ps)
  out
}

# take basis model from above, and subsample to simulate sequencing
# functions to model sequencing
#' subsample community
#' @param comm1 phyloseq object
#' @param sample vector specifying sampling depth
#' @param b depth
#' @param c variation
#' @keywords reference community model microbiome
#' @export
#' @examples
#' sample.model()
sample.model<-function(comm1, sample, b, c){
  if(any(sample_sums(comm1)>sample)){
    while(any(sample_sums(comm1)<sample)){
      sample<-round(rnorm(10, b, c))
      sample
    }}
  #print(sample)
  a<-make.table(comm1, sample)
  a
}

#' subsample community
#' @param comm1 phyloseq object
#' @param sample vector of arbitrary sampling depth set by set.seqDepth()
#' @keywords reference community model microbiome
#' @export
#' @examples
#' make.table()
make.table<-function(comm1, sample){
  comm2<-comm1
  otu<-as.data.frame(as.matrix(otu_table(comm1)))
  otu[]<-0
  otu<-otu[order(rownames(otu)),]
  otu_table(comm2)<-otu_table(otu, taxa_are_rows=TRUE)
  comm<-as.data.frame(as.matrix(otu_table(comm1)))
  m<-as.data.frame(t(table(sample(rownames(comm), sample[1], replace=T, prob=comm[,1]/sum(comm[,1])))))
  m<-m[,colnames(m)!="Var1"]
  colnames(m)[colnames(m)=="Freq"]<-paste("Sample", 1, sep="")
  #print("step one")
  m2<-as.data.frame(t(table(sample(rownames(comm),sample[2], replace=T, prob=comm[,2]/sum(comm[,2])))))
  m<-merge(m, m2, by="Var2", all=T)
  #print("step two")
  m<-m[,colnames(m)!="Var1"]
  colnames(m)[colnames(m)=="Freq"]<-paste("Sample", 2, sep="")
  for(i in 3:ncol(comm)){
    a<-as.data.frame(t(table(sample(rownames(comm), sample[i], replace=T, prob=comm[,i]/sum(comm[,i])))))
    a<-a[,colnames(a)!="Var1"]
    m<-merge(m,a, by="Var2", all=T)
    m<-m[,colnames(m)!="Var1"]
    colnames(m)[colnames(m)=="Freq"]<-paste("Sample", i, sep="")
    
  }
  rownames(m)<-m$Var2
  m<-m[,-1]
  m[is.na(m)]<-0
  m<-m[order(rownames(m)),]
  otu_table(comm1)<-otu_table(m, taxa_are_rows=T)
  comm3<-merge_phyloseq(comm1, comm2)
  m2<-as.data.frame(as.matrix(otu_table(comm3)))
  m2<-m2[order(rownames(m2)),]
  print(rownames(m2))
  print(rownames(otu))
  if(!identical(colnames(m2),colnames(otu))){
    stop("error: colnames do not match in rarefaction")
  }
  if(!identical(rownames(m2), rownames(otu))){
    stop("error: rownames do not match in rarefaction")
  }
  otu_table(comm3)<-otu_table(m2, taxa_are_rows=T)
  comm3
}

#' wrapper for Qscale - workhorse for QSeq transformation
#' @param ps phyloseq object
#' @param type select "B" if gradient, "A" if no underlying gradient
#' @keywords reference community model microbiome
#' @export
#' @examples
#' QScale()
QScale<-function(ps, type){
  
  if(type=="B"){
    scale<-as.numeric(as.character(sample_data(ps)$R2.Gradient)) # for gradient datasets
    print(scale)
    out<-Qscale(ps, scale)
  }
  
  if(type=="A"){
    scale<-as.numeric(as.character(sample_data(ps)$R1.Gradient)) # for nongradient datsets
    print(scale)
    out<-Qscale(ps, scale)
  }
  out
}

#' Qscale - workhorse for QSeq transformation
#' @param ps phyloseq object
#' @param scale vector of abundances to scale sample data
#' @keywords reference community model microbiome
#' @export
#' @examples
#' Qscale()
Qscale<-function(ps, scale){
  # extract and transform otu table from phyloseq
  scaled<-data.frame(mapply(`*`, data.frame(as.matrix(otu_table(transform_sample_counts(ps, function(x) x/sum(x))))), scale))
  # add names to matrix
  names<-rownames(data.frame(as.matrix(otu_table(ps))))
  print(names)
  rownames(scaled)<-names
  # update and return phyloseq object with QSeq otu table
  scaled<-round(scaled)
  p2<-ps
  otu_table(p2)<- otu_table(scaled, taxa_are_rows=T)
  p2
  
}

#' initializes the model space by creating base
#' @param nreps number of iterations for the simulation
#' @keywords 
#' @export
#' @examples
#' initialize.data()
initialize.data<-function(nreps){
  out<-as.list(1:nreps)
  names(out)<-paste("Rep", c(1:nreps), sep="")
  for(i in 1:length(out)){
    out[[i]]<-matrix(nrow=50, ncol=10)
    out[[i]]<-matrix(round(rnorm(500,25,10)), nrow=50, ncol=10)
    out[[i]][out[[i]]<0]<-0
    rownames(out[[i]])<-paste("Taxon", c(1:50), sep="")
    colnames(out[[i]])<-paste("Sample", c(1:10), sep="")
  }
  
  out
} # delete maybe?


#' initializes the model space by creating n matrices
#' @param nreps number of iterations for the simulation
#' @param dist vector of 50 means for skewness
#' @keywords 
#' @export
#' @examples
#' initialize.data()
initialize.skew<-function(nreps, dist){
  out<-as.list(1:nreps)
  names(out)<-paste("Rep", c(1:nreps), sep="")
  for(i in 1:length(out)){
    out[[i]]<-matrix(nrow=50, ncol=10)
    out[[i]]<-t(matrix(round(sapply(dist, function(x) rnorm(10, x, 10))), nrow=10, ncol=50))
    out[[i]][out[[i]]<0]<-0
    rownames(out[[i]])<-paste("Taxon", c(1:50), sep="")
    colnames(out[[i]])<-paste("Sample", c(1:10), sep="")
  }
  
  out
}


#' initialize simulation with skewness and sparseness
#' @param nreps number of times to run the simulation
#' @param dist vector of mean values for taxa
#' @param sparse degree of arbitrary sparsity (total number of observations to replace with zero)
#' @keywords 
#' @export
#' @examples
#' initialize.sparseskew()
initialize.sparseskew<-function(nreps, dist, sparse){
  out<-as.list(1:nreps)
  names(out)<-paste("Rep", c(1:nreps), sep="")
  for(i in 1:length(out)){
    out[[i]]<-matrix(nrow=50, ncol=10)
    out[[i]]<-t(matrix(round(sapply(dist, function(x) rnorm(10, x, .1*x))), nrow=10, ncol=50))
    out[[i]][out[[i]]<0]<-0
    out[[i]][sample(c(1:500), sparse)]<-0
    rownames(out[[i]])<-paste("Taxon", c(1:50), sep="")
    colnames(out[[i]])<-paste("Sample", c(1:10), sep="")
  }
  
  out
}


# uses output list from initialize.sparseskew() 
#' initialize simulation with skewness and sparseness
#' @param nreps number of times to run the simulation
#' @param dist vector of mean values for taxa
#' @param sparse degree of arbitrary sparsity (total number of observations to replace with zero)
#' @keywords 
#' @export
#' @examples
#' initialize.sparseskew()
impose.gradient<-function(init.list){
  
  gradient<-seq(1,30, by=3)
  out<-init.list
  for(j in 1:length(init.list)){
    for(i in 1:length(gradient)){
    out[[j]][,i]<-init.list[[j]][,i]*gradient[i]
    
    }
  }
  out
}

#' create gradient, give some flexibility to impose a variety of gradients:
#' @param init.list Reference matrices for simulation
#' @param grad vector sample total abundances to apply
#' @keywords 
#' @export
#' @examples
#' impose.gradient2()
impose.gradient2<-function(init.list, grad){
  
  out<-init.list
  for(j in 1:length(init.list)){
    gradient<-sort(round(c(rnorm(10, 10, grad))))
    if(any(gradient<1)){
      while(any(gradient<1)){gradient<-sort(round(c(rnorm(10, 10, grad))))}
    }
    for(i in 1:length(gradient)){
      out[[j]][,i]<-init.list[[j]][,i]*gradient[i]
      
    }
  }
  out
}
#' implements the core analysis of the simulation
#' @param ref reference with no gradient (base)
#' @param gradientref reference with gradient
#' @param error QSeq error rate as a percentage of total abudance(how accurate is QPCR?)
#' @keywords 
#' @export
#' @examples
#' analysis.guts()
analysis.guts<-function(ref, gradientref, error){
  # generate metadata
  samdat<-data.frame("Sample"=paste("Sample", c(1:10), sep=""), "Gradient" = seq(1,30, by=3), "Group" = c(1,1,1,1,1,2,2,2,2,2))
  rownames(samdat)<-samdat$Sample

  # impose an error penalty on QSeq (measurement error from QPCR)
  samdat$R1.Gradient<-colSums(ref)+rnorm(10,0,error*colSums(ref))
  if(any(samdat$R1.Gradient<1)){
  while(any(samdat$R1.Gradient<1)){
    samdat$R1.Gradient<-samdat$R1.Gradient+1
  }}
  samdat$R2.Gradient<-colSums(gradientref)+rnorm(10,0,error*colSums(gradientref))
  if(any(samdat$R2.Gradient<1)){
  while(any(samdat$R2.Gradient<1)){
    samdat$R2.Gradient<-samdat$R2.Gradient+1
  }}

  #print(c(samdat$R1.Gradient, samdat$R2.Gradient))
  
  # create reference phyloseq objects
  ps.ref<-phyloseq(otu_table(ref, taxa_are_rows=TRUE), sample_data(samdat))
  ps.ref.gradient<-phyloseq(otu_table(gradientref, taxa_are_rows=TRUE), sample_data(samdat))
  
  # create sample phyloseq objects
  print("start sample model")
  sample<-round(rnorm(10, 800, 30)) # check to make sure samples are reasonable
  S.cmap<-sample.model(ps.ref, sample, 800, 30)
  print("end sample model")
  otu_table(S.cmap)<-otu_table(as.matrix(otu_table(S.cmap))[rownames(ref),], taxa_are_rows = TRUE)
  S.Gradient<-sample.model(ps.ref.gradient, sample, 800, 30)
  otu_table(S.Gradient)<-otu_table(as.matrix(otu_table(S.Gradient))[rownames(gradientref),], taxa_are_rows = TRUE)
  
  out<-NULL
  
  Q.cmap<-QScale(S.cmap, type = "A") # add function ........ []
  Q.cmapGradient<-QScale(S.Gradient, type = "B")  # add function ........ []
  print("models finished")
  # get summary stats for accuracy of correlations [given assumptions] finished []
  # maybe add in fastspar on Qseq values
  sparcor.cmat<-run.FastSpar(S.cmap, name="cmat", threads=8)  # add function ........ []
  sparcor.gradient<-run.FastSpar(S.Gradient, name="cmat_gradient", threads=8)
  
  ref.cor<-cor(t(ref))
  ref.gradient.cor<-cor(t(gradientref))
  
  out$coraccuracy<-matrix(nrow=6,ncol=2)
  colnames(out$coraccuracy)<-c("Reference", "Reference.Gradient")
  rownames(out$coraccuracy)<-c("Subsample", "Subsample.Gradient", "FastSpar", "FastSpar.Gradient", "QSeqCor", "QSeqCor.Gradient")
  
  # double check all these value inputs ........... []
  
  out$coraccuracy[1,1]<-cor(as.vector(ref.cor), as.vector(cor(t(as.matrix(otu_table(S.cmap))[rownames(ref)]))))
  out$coraccuracy[2,1]<-cor(as.vector(ref.cor), as.vector(cor(t(as.matrix(otu_table(S.Gradient))[rownames(ref)]))))
  out$coraccuracy[3,1]<-cor(as.vector(ref.cor), as.vector(sparcor.cmat$cor))
  out$coraccuracy[4,1]<-cor(as.vector(ref.cor), as.vector(sparcor.gradient$cor))
  out$coraccuracy[5,1]<-cor(as.vector(cor(t(as.matrix(otu_table(Q.cmap))[rownames(ref)]))), as.vector(ref.cor))
  out$coraccuracy[6,1]<-cor(as.vector(cor(t(as.matrix(otu_table(Q.cmapGradient))[rownames(ref)]))), as.vector(ref.cor))
  out$coraccuracy[1,2]<-cor(as.vector(ref.gradient.cor),as.vector(cor(t(as.matrix(otu_table(S.cmap))[rownames(gradientref)]))))
  out$coraccuracy[2,2]<-cor(as.vector(ref.gradient.cor), as.vector(cor(t(as.matrix(otu_table(S.Gradient))[rownames(gradientref)]))))
  out$coraccuracy[3,2]<-cor(as.vector(ref.gradient.cor),  as.vector(sparcor.cmat$cor))
  out$coraccuracy[4,2]<-cor(as.vector(ref.gradient.cor),  as.vector(sparcor.gradient$cor))
  out$coraccuracy[5,2]<-cor(as.vector(cor(t(as.matrix(otu_table(Q.cmap))[rownames(ref)]))), as.vector(ref.gradient.cor))
  out$coraccuracy[6,2]<-cor(as.vector(cor(t(as.matrix(otu_table(Q.cmapGradient))[rownames(gradientref)]))), as.vector(ref.gradient.cor))
  
  # run differential abundances
  
  diff<-matrix(nrow=nrow(ref), ncol=26)
  colnames(diff)<-c("Coef.ref", "pval.ref", "Coef.grad.ref", "pval.grad.ref", "Coef.aldex","pval.aldex","Coef.aldex.grad","pval.aldex.grad","Coef.aldex.Q","pval.aldex.Q","Coef.aldex.gQ", "pval.aldex.gQ","Coef.ancombc", "pval.ancombc","Coef.ancombc.Q","pval.ancombc.Q","Coef.ancombc.grad", "pval.ancombc.grad","Coef.ancombc.gQ", "pval.ancombc.gQ","Coef.QS","pval.QS","Coef.QS.grad", "pval.QS.grad", "Coeff.ref.ancomb", "pval.ref.ancomb")
  QS.cmat<-as.matrix(as.data.frame(as.matrix(otu_table(Q.cmap))))
  QS.cmat.grad<-as.matrix(as.data.frame(as.matrix(otu_table(Q.cmapGradient))))
  
  for(i in 1:nrow(ref)){
    # get reference glm
    complm<-NULL
    complm<-glm(ref[i,]~samdat$Group)
    summlm<-summary(aov(complm))
    # get reference glm for gradient
    complm2<-NULL
    complm2<-glm(gradientref[i,]~samdat$Group)
    summlm2<-summary(aov(complm2))
    # run qseq glm
    QSlm<-NULL
    QSlm<-glm(QS.cmat[i,]~samdat$Group)
    QSlm1<-summary(aov(QSlm))
    # run qseq glm on gradient data
    QSlmG<-NULL
    QSlmG<-glm(QS.cmat.grad[i,]~samdat$Group)
    QSlmG1<-summary(aov(QSlmG))
    
    diff[i,1]<-complm$coefficients[2]
    diff[i,3]<-complm2$coefficients[2]
    diff[i,2]<-summlm[[1]]$`Pr(>F)`[1]
    diff[i,4]<-summlm2[[1]]$`Pr(>F)`[1]
    
    diff[i,21]<-QSlm$coefficients[2]
    diff[i,23]<-QSlmG$coefficients[2]
    diff[i,22]<-QSlm1[[1]]$`Pr(>F)`[1]
    diff[i,24]<-QSlmG1[[1]]$`Pr(>F)`[1]
  }
  
  # turn effects into logical for direction
  diff[diff[,21]<0,21]<--1
  diff[diff[,23]<0,23]<--1
  diff[diff[,21]>0,21]<-1
  diff[diff[,23]>0,23]<-1
  
  diff[diff[,1]<0,1]<--1
  diff[diff[,3]<0,3]<--1
  diff[diff[,1]>0,1]<-1
  diff[diff[,3]>0,3]<-1
  
  print(S.cmap)
  print(S.Gradient)
  
  aldex1<-runALDEx(S.cmap) # add this function .... []
  aldex1q<-runALDEx(Q.cmap)
  aldexgrad<-runALDEx(S.Gradient)
  aldexgradQ<-runALDEx(Q.cmapGradient)
  
  diff[,6]<-as.vector(aldex1$test$glm.ep)
  diff[,8]<-as.vector(aldex1q$test$glm.ep)
  diff[,10]<-as.vector(aldexgrad$test$glm.ep)
  diff[,12]<-as.vector(aldexgradQ$test$glm.ep)
  
  diff[,5]<-as.vector(aldex1$effect$effect)
  diff[,7]<-as.vector(aldex1q$effect$effect)
  diff[,9]<-as.vector(aldexgrad$effect$effect)
  diff[,11]<-as.vector(aldexgradQ$effect$effect)
  
  diff[diff[,5]>0,5]<-1
  diff[diff[,5]<0,5]<--1
  diff[diff[,7]>0,7]<-1
  diff[diff[,7]<0,7]<--1
  diff[diff[,9]>0,9]<-1
  diff[diff[,9]<0,9]<--1
  diff[diff[,11]>0,11]<-1
  diff[diff[,11]<0,11]<--1
  
  
  ancomb.1<-ancombc(phyloseq = S.cmap, formula = "Group", # add this function ........ []
                    p_adj_method = "none", zero_cut = 0.90, lib_cut = 0, 
                    group = "Group", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
  
  ancomb.gradient<-ancombc(phyloseq = S.Gradient, formula = "Group", 
                           p_adj_method = "none", zero_cut = 0.90, lib_cut = 0, 
                           group = "Group", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                           max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
  
  ancomb.1Q<-ancombc(phyloseq = Q.cmap, formula = "Group", 
                     p_adj_method = "none", zero_cut = 0.90, lib_cut = 0, 
                     group = "Group", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                     max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
  
  ancomb.gradientQ<-ancombc(phyloseq = Q.cmapGradient, formula = "Group", 
                            p_adj_method = "none", zero_cut = 0.90, lib_cut = 0, 
                            group = "Group", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                            max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
  
  diff[,13]<-as.vector(ancomb.1$res$beta$Group)  # make sure beta is the right thing to extract ... []
  diff[,15]<-as.vector(ancomb.gradient$res$beta$Group)
  diff[,17]<-as.vector(ancomb.1Q$res$beta$Group)
  diff[,19]<-as.vector(ancomb.gradientQ$res$beta$Group)
  
  diff[diff[,13]<0,13]<--1
  diff[diff[,15]<0,15]<--1
  diff[diff[,13]>0,13]<-1
  diff[diff[,15]>0,15]<-1
  
  diff[diff[,17]<0,17]<--1
  diff[diff[,19]<0,19]<--1
  diff[diff[,17]>0,17]<-1
  diff[diff[,19]>0,19]<-1
  
  diff[,14]<-as.vector(ancomb.1$res$p_val$Group)
  diff[,16]<-as.vector(ancomb.gradient$res$p_val$Group)
  diff[,18]<-as.vector(ancomb.1Q$res$p_val$Group)
  diff[,20]<-as.vector(ancomb.gradientQ$res$p_val$Group)
  
  # define true positive for difference in relative abundance  ...  []
  ancomb<-ancombc(phyloseq = S.cmap, formula = "Group", # add this function ........ []
                    p_adj_method = "none", zero_cut = 0.90, lib_cut = 0, 
                    group = "Group", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
  diff[,25]<-as.vector(ancomb$res$beta$Group)# effect direction
  diff[,26]<-as.vector(ancomb$res$p_val$Group) # p val
  #View(diff)
  # get summary stats for differential abundance testing
  
  out$effectaccuracy<-matrix(nrow=10,ncol=8)
  colnames(out$effectaccuracy)<-c("CTruecall.p", "CTruecall.n", "CFalsecall.p", "CFalsecall.n","GTruecall.p", "GTruecall.n", "GFalsecall.p", "GFalsecall.n")
  rownames(out$effectaccuracy)<-c("Aldex", "Aldex.gradient", "Aldex.Q", "Aldex.gradientQ", "ancombc", "ancombc.gradient", "ancombc.Q", "ancombc.gradientQ", "QSeq", "QSeq.gradient")
  
  #colnames(diff)<-c("Coef.ref", "pval.ref", "Coef.grad.ref", "pval.grad.ref", "Coef.aldex","pval.aldex","Coef.aldex.grad","pval.aldex.grad","Coef.aldex.Q","pval.aldex.Q","Coef.aldex.gQ", "pval.aldex.gQ","Coef.ancombc", "pval.ancombc","Coef.ancombc.Q","pval.ancombc.Q","Coef.ancombc.grad", "pval.ancombc.grad","Coef.ancombc.gQ", "pval.ancombc.gQ","Coef.QS","pval.QS","Coef.QS.grad", "pval.QS.grad")
  
  # compositional analysis
  out$effectaccuracy[1,1]<-sum(diff[,1]==diff[,5]&diff[,2]<0.05&diff[,6]<0.05)
  out$effectaccuracy[1,2]<-sum(diff[,2]>0.05&diff[,6]>0.05)
  out$effectaccuracy[1,3]<-sum(diff[,2]>0.05&diff[,6]<0.05) # rewrite this statement diff[,26]>0.05&diff[,6]<0.05
  out$effectaccuracy[1,4]<-sum(diff[,2]<0.05&diff[,6]>0.05)# coefficient doesn't matter if no detection is made
  
  out$effectaccuracy[2,1]<-sum(diff[,1]==diff[,7]&diff[,2]<0.05&diff[,8]<0.05)
  out$effectaccuracy[2,2]<-sum(diff[,2]>0.05&diff[,8]>0.05)
  out$effectaccuracy[2,3]<-sum(diff[,2]>0.05&diff[,8]<0.05)
  out$effectaccuracy[2,4]<-sum(diff[,2]<0.05&diff[,8]>0.05)
  
  out$effectaccuracy[3,1]<-sum(diff[,1]==diff[,9]&diff[,2]<0.05&diff[,10]<0.05)
  out$effectaccuracy[3,2]<-sum(diff[,2]>0.05&diff[,10]>0.05)
  out$effectaccuracy[3,3]<-sum(diff[,2]>0.05&diff[,10]<0.05)
  out$effectaccuracy[3,4]<-sum(diff[,2]<0.05&diff[,10]>0.05)
  
  out$effectaccuracy[4,1]<-sum(diff[,1]==diff[,11]&diff[,2]<0.05&diff[,12]<0.05)
  out$effectaccuracy[4,2]<-sum(diff[,2]>0.05&diff[,12]>0.05)
  out$effectaccuracy[4,3]<-sum(diff[,2]>0.05&diff[,12]<0.05)
  out$effectaccuracy[4,4]<-sum(diff[,2]<0.05&diff[,12]>0.05)
  
  out$effectaccuracy[5,1]<-sum(diff[,1]==diff[,13]&diff[,2]<0.05&diff[,14]<0.05)
  out$effectaccuracy[5,2]<-sum(diff[,2]>0.05&diff[,14]>0.05)
  out$effectaccuracy[5,3]<-sum(diff[,2]>0.05&diff[,14]<0.05)
  out$effectaccuracy[5,4]<-sum(diff[,2]<0.05&diff[,14]>0.05)
  
  out$effectaccuracy[6,1]<-sum(diff[,1]==diff[,15]&diff[,2]<0.05&diff[,16]<0.05)
  out$effectaccuracy[6,2]<-sum(diff[,2]>0.05&diff[,16]>0.05)
  out$effectaccuracy[6,3]<-sum(diff[,2]>0.05&diff[,16]<0.05)
  out$effectaccuracy[6,4]<-sum(diff[,2]<0.05&diff[,16]>0.05)
  
  out$effectaccuracy[7,1]<-sum(diff[,1]==diff[,17]&diff[,2]<0.05&diff[,18]<0.05)
  out$effectaccuracy[7,2]<-sum(diff[,2]>0.05&diff[,18]>0.05)
  out$effectaccuracy[7,3]<-sum(diff[,2]>0.05&diff[,18]<0.05)
  out$effectaccuracy[7,4]<-sum(diff[,2]<0.05&diff[,18]>0.05)
  
  out$effectaccuracy[8,1]<-sum(diff[,1]==diff[,19]&diff[,2]<0.05&diff[,20]<0.05)
  out$effectaccuracy[8,2]<-sum(diff[,2]>0.05&diff[,20]>0.05)
  out$effectaccuracy[8,3]<-sum(diff[,2]>0.05&diff[,20]<0.05)
  out$effectaccuracy[8,4]<-sum(diff[,2]<0.05&diff[,20]>0.05)
  
  out$effectaccuracy[9,1]<-sum(diff[,1]==diff[,21]&diff[,2]<0.05&diff[,22]<0.05)
  out$effectaccuracy[9,2]<-sum(diff[,2]>0.05&diff[,22]>0.05)
  out$effectaccuracy[9,3]<-sum(diff[,2]>0.05&diff[,22]<0.05)
  out$effectaccuracy[9,4]<-sum(diff[,2]<0.05&diff[,22]>0.05)
  
  out$effectaccuracy[10,1]<-sum(diff[,1]==diff[,23]&diff[,2]<0.05&diff[,24]<0.05)
  out$effectaccuracy[10,2]<-sum(diff[,2]>0.05&diff[,24]>0.05)
  out$effectaccuracy[10,3]<-sum(diff[,2]>0.05&diff[,24]<0.05)
  out$effectaccuracy[10,4]<-sum(diff[,2]<0.05&diff[,24]>0.05)
  
 # non compositional approach:
  out$effectaccuracy[1,5]<-sum(diff[,3]==diff[,5]&diff[,4]<0.05&diff[,6]<0.05)
  out$effectaccuracy[1,6]<-sum(diff[,4]>0.05&diff[,6]>0.05)
  out$effectaccuracy[1,7]<-sum(diff[,4]>0.05&diff[,6]<0.05)
  out$effectaccuracy[1,8]<-sum(diff[,4]<0.05&diff[,6]>0.05)
  
  out$effectaccuracy[2,5]<-sum(diff[,3]==diff[,7]&diff[,4]<0.05&diff[,8]<0.05)
  out$effectaccuracy[2,6]<-sum(diff[,4]>0.05&diff[,8]>0.05)
  out$effectaccuracy[2,7]<-sum(diff[,4]>0.05&diff[,8]<0.05)
  out$effectaccuracy[2,8]<-sum(diff[,4]<0.05&diff[,8]>0.05)
  
  out$effectaccuracy[3,5]<-sum(diff[,3]==diff[,9]&diff[,4]<0.05&diff[,10]<0.05)
  out$effectaccuracy[3,6]<-sum(diff[,4]>0.05&diff[,10]>0.05)
  out$effectaccuracy[3,7]<-sum(diff[,4]>0.05&diff[,10]<0.05)
  out$effectaccuracy[3,8]<-sum(diff[,4]<0.05&diff[,10]>0.05)
  
  out$effectaccuracy[4,5]<-sum(diff[,3]==diff[,11]&diff[,4]<0.05&diff[,12]<0.05)
  out$effectaccuracy[4,6]<-sum(diff[,4]>0.05&diff[,12]>0.05)
  out$effectaccuracy[4,7]<-sum(diff[,4]>0.05&diff[,12]<0.05)
  out$effectaccuracy[4,8]<-sum(diff[,4]<0.05&diff[,12]>0.05)
  
  out$effectaccuracy[5,5]<-sum(diff[,3]==diff[,13]&diff[,4]<0.05&diff[,14]<0.05)
  out$effectaccuracy[5,6]<-sum(diff[,4]>0.05&diff[,14]>0.05)
  out$effectaccuracy[5,7]<-sum(diff[,4]>0.05&diff[,14]<0.05)
  out$effectaccuracy[5,8]<-sum(diff[,4]<0.05&diff[,14]>0.05)
  
  out$effectaccuracy[6,5]<-sum(diff[,3]==diff[,15]&diff[,4]<0.05&diff[,16]<0.05)
  out$effectaccuracy[6,6]<-sum(diff[,4]>0.05&diff[,16]>0.05)
  out$effectaccuracy[6,7]<-sum(diff[,4]>0.05&diff[,16]<0.05)
  out$effectaccuracy[6,8]<-sum(diff[,4]<0.05&diff[,16]>0.05)
  
  out$effectaccuracy[7,5]<-sum(diff[,3]==diff[,17]&diff[,4]<0.05&diff[,18]<0.05)
  out$effectaccuracy[7,6]<-sum(diff[,4]>0.05&diff[,18]>0.05)
  out$effectaccuracy[7,7]<-sum(diff[,4]>0.05&diff[,18]<0.05)
  out$effectaccuracy[7,8]<-sum(diff[,4]<0.05&diff[,18]>0.05)
  
  out$effectaccuracy[8,5]<-sum(diff[,3]==diff[,19]&diff[,4]<0.05&diff[,20]<0.05)
  out$effectaccuracy[8,6]<-sum(diff[,4]>0.05&diff[,20]>0.05)
  out$effectaccuracy[8,7]<-sum(diff[,4]>0.05&diff[,20]<0.05)
  out$effectaccuracy[8,8]<-sum(diff[,4]<0.05&diff[,20]>0.05)
  
  out$effectaccuracy[9,5]<-sum(diff[,3]==diff[,21]&diff[,4]<0.05&diff[,22]<0.05)
  out$effectaccuracy[9,6]<-sum(diff[,4]>0.05&diff[,22]>0.05)
  out$effectaccuracy[9,7]<-sum(diff[,4]>0.05&diff[,22]<0.05)
  out$effectaccuracy[9,8]<-sum(diff[,4]<0.05&diff[,22]>0.05)
  
  out$effectaccuracy[10,5]<-sum(diff[,3]==diff[,23]&diff[,4]<0.05&diff[,24]<0.05)
  out$effectaccuracy[10,6]<-sum(diff[,4]>0.05&diff[,24]>0.05)
  out$effectaccuracy[10,7]<-sum(diff[,4]>0.05&diff[,24]<0.05)
  out$effectaccuracy[10,8]<-sum(diff[,4]<0.05&diff[,24]>0.05)
  
  out$diftab<-diff
  
  # return summary data
  out
}

#' wrapper to create the data for the paper / run the analysis
#' @param Reflist list of reference matrices with no gradient (base)
#' @param GReflist List of reference matrices with gradient
#' @param error QSeq error rate as a percentage of total abudance(how accurate is QPCR?)
#' @keywords 
#' @export
#' @examples
#' impose.gradient2()
run.analysis<-function(Reflist, GReflist, error=0.1){
  require(phyloseq)
  require(ANCOMBC)
  require(ALDEx2)
  
  out<-as.list(names(Reflist))
  names(out)<-names(Reflist)
  
  for(i in 1:length(Reflist)){
    
    out[[i]]<-analysis.guts(Reflist[[i]],GReflist[[i]], error=error)
    
  }
# run / summarize analysis ####

  # correlation between analysis and reflist of taxon correlations for conditions
  
  out
}


# take basis model from above, and subsample to simulate sequencing
# functions to model sequencing
#' subsample community
#' @param comm1 phyloseq object
#' @param sample vector specifying sampling depth
#' @param b depth
#' @param v variation
#' @keywords reference community model microbiome
#' @export
#' @examples
#' sample.model()
sample.model<-function(comm1, sample, b, v){
  if(any(sample_sums(comm1)>sample)){
    while(any(sample_sums(comm1)<sample)){
      sample<-round(rnorm(10, b, v))
      sample
    }}
  #print(sample)
  a<-make.table(comm1, sample)
  a
}


#' compile summary statistics on taxon correlations using each method
#' @param x output object from run.analysis
#' @keywords 
#' @export
#' @examples
#' compile.summary.table.correlations()
compile.summary.table.correlations<-function(x){
  outtab<-matrix(ncol=4, nrow=3)
  colnames(outtab)<-c("Reference Mean", "Reference SE", "Ref.Gradient Mean", "Ref.Gradient SE")
  rownames(outtab)<-c("Subsample", "FastSpar", "QSeq")
  subsample.ref<-c(rep(NA,length(x)))
  subsample.gradient<-c(rep(NA,length(x)))
  FastSpar.ref<-c(rep(NA,length(x)))
  FastSpar.gradient<-c(rep(NA,length(x)))
  QSeq<-c(rep(NA,length(x)))
  QSeq.ref<-c(rep(NA,length(x)))
  
  for(i in 1:length(x)){
    subsample.ref[i]<-x[[i]]$coraccuracy[1,1]
    subsample.gradient[i]<-x[[i]]$coraccuracy[2,2]
    FastSpar.ref[i]<-x[[i]]$coraccuracy[3,1]
    FastSpar.gradient[i]<-x[[i]]$coraccuracy[4,2]
    QSeq[i]<-x[[i]]$coraccuracy[5,1]
    QSeq.ref[i]<-x[[i]]$coraccuracy[6,2]
  }
  
  valdf<-data.frame(subsample.ref, subsample.gradient, FastSpar.ref, FastSpar.gradient, QSeq, QSeq.ref)
  
  outtab[1,1]<-mean(subsample.ref)
  outtab[1,2]<-sd(subsample.ref)/sqrt(length(subsample.ref))
  outtab[1,3]<-mean(subsample.gradient)
  outtab[1,4]<-sd(subsample.gradient)/sqrt(length(subsample.gradient))
  outtab[2,1]<-mean(FastSpar.ref)
  outtab[2,2]<-sd(FastSpar.ref)/sqrt(length(FastSpar.ref))
  outtab[2,3]<-mean(FastSpar.gradient)
  outtab[2,4]<-sd(FastSpar.gradient)/sqrt(length(FastSpar.gradient))
  outtab[3,1]<-mean(QSeq)
  outtab[3,2]<-sd(QSeq)/sqrt(length(QSeq))
  outtab[3,3]<-mean(QSeq.ref)
  outtab[3,4]<-sd(QSeq.ref)/sqrt(length(QSeq.ref))
  
  out<-NULL
  out$cortab<-outtab
  out$df<-valdf
  out
}

#' make plots from correlation tables
#' @param cortab.list list of correlation tables from multiple run.analysis
#' @keywords 
#' @export
#' @examples
#' make.plot.gradients()
make.plot.gradients<-function(cortab.list){
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  df<-ldply(cortab.list, cbind)
  df$level<-c(rep(names(cortab.list)[1], 3),rep(names(cortab.list)[2], 3),rep(names(cortab.list)[3], 3),rep(names(cortab.list)[4], 3),rep(names(cortab.list)[5], 3))
  df$TRT<-c(rep(c("Subsample", "FastSpar", "QSeq")))
  #df1<-melt(df[,c(1,2)], id=c("Reference Mean", "Reference SE"))
  #df2<-melt(df[,c(3,4)], id=c("Ref.Gradient Mean", "Ref.Gradient SE"))
  #df1<-df[,c(2,3,6,7)]
  #print(df1)
  #df2<-df[,c(3,4,5)]       
  gplt1<-ggplot(df, aes(x=level, y=`Reference Mean`, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=`Reference Mean`-`Reference SE`, ymax=`Reference Mean`+`Reference SE`), alpha=.1)+
    ggtitle("No Gradient") +
    ylab("Mean Correlation") + 
    xlab("Percent Variation in Total Abundance") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  gplt2<-ggplot(df, aes(x=level, y=`Ref.Gradient Mean`, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=`Ref.Gradient Mean`-`Ref.Gradient SE`, ymax=`Ref.Gradient Mean`+`Ref.Gradient SE`, fill=TRT), alpha=.1) +
    ggtitle("Gradient") +
    ylab("Mean Correlation") + 
    xlab("Percent Variation in Total Abundance") +
    scale_fill_manual(breaks=c("Subsample", "FastSpar", "QSeq"), values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  out<-NULL
  out$NG<-gplt1
  out$GD<-gplt2
  out
}

#' make plots from differential abundance analysis
#' @param effect.list list of differential abundance effects from multiple runs of run.analysis
#' @keywords 
#' @export
#' @examples
#' make.plot.gradients()
make.plot2<-function(effect.list){
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  df<-ldply(effect.list, cbind)
  df$level<-c(rep(names(effect.list)[1], 3),rep(names(effect.list)[2], 3),rep(names(effect.list)[3], 3),rep(names(effect.list)[4], 3),rep(names(effect.list)[5], 3))
  df$TRT<-c(rep(c("ALDEx", "ancomBC", "QSeq")))
  
  TP<-ggplot(df, aes(x=level, y=TP.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=TP.Mean-TP.SE, ymax=TP.Mean+TP.SE), alpha=.1)+
    ggtitle("TRUE Positive, No Gradient") +
    ylab("Mean True Positive") + 
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#9D5A6C", "#E47250")) +
    scale_colour_manual(values = c("#5A4A6F", "#9D5A6C", "#E47250")) +
    theme_bw()
  TN<-ggplot(df, aes(x=level, y=TN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=TN.Mean-TN.SE, ymax=TN.Mean+TN.SE), alpha=.1)+
    ggtitle("TRUE Negative, No Gradient") +
    ylab("Mean True Negative") + 
    xlab("Percent Variation in Total Abundance") +
    scale_fill_manual(values=c("#5A4A6F", "#9D5A6C", "#E47250")) +
    scale_colour_manual(values = c("#5A4A6F", "#9D5A6C", "#E47250")) +
    theme_bw()
  
  GTP<-ggplot(df, aes(x=level, y=GTP.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GTP.Mean-GTP.SE, ymax=GTP.Mean+GTP.SE), alpha=.1)+
    ggtitle("TRUE Positive, Gradient") +
    ylab("Average Proportion True Positive Detection") + 
    xlab("Percent Variation in Total Abundance") +
    scale_fill_manual(values=c("#5A4A6F", "#9D5A6C", "#E47250")) +
    scale_colour_manual(values = c("#5A4A6F", "#9D5A6C", "#E47250")) +
    theme_bw()
  GTN<-ggplot(df, aes(x=level, y=GTN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GTN.Mean-GTN.SE, ymax=GTN.Mean+GTN.SE), alpha=.1)+
    ggtitle("TRUE Negative, Gradient") +
    ylab("Mean True Negative") + 
    xlab("Percent Variation in Total Abundance") +
    scale_fill_manual(values=c("#5A4A6F", "#9D5A6C", "#E47250")) +
    scale_colour_manual(values = c("#5A4A6F", "#9D5A6C", "#E47250")) +
    theme_bw()
  
  FP<-ggplot(df, aes(x=level, y=FP.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=FP.Mean-FP.SE, ymax=FP.Mean+FP.SE), alpha=.1)+
    ggtitle("FALSE Positive, No Gradient") +
    ylab("Mean FALSE Positive") + 
    xlab("Percent Variation in Total Abundance") +
    scale_fill_manual(values=c("#5A4A6F", "#9D5A6C", "#E47250")) +
    scale_colour_manual(values = c("#5A4A6F", "#9D5A6C", "#E47250")) +
    theme_bw()
  FN<-ggplot(df, aes(x=level, y=FN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=FN.Mean-FN.SE, ymax=FN.Mean+FN.SE), alpha=.1)+
    ggtitle("FALSE Negative, No Gradient") +
    ylab("Mean FALSE Negative") + 
    xlab("Percent Variation in Total Abundance") +
    scale_fill_manual(values=c("#5A4A6F", "#9D5A6C", "#E47250")) +
    scale_colour_manual(values = c("#5A4A6F", "#9D5A6C", "#E47250")) +
    theme_bw()
  
  GFP<-ggplot(df, aes(x=level, y=GFP.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GFP.Mean-GFP.SE, ymax=GFP.Mean+GFP.SE), alpha=.1)+
    ggtitle("FALSE Positive, Gradient") +
    ylab("Mean FALSE Positive") + 
    xlab("Percent Variation in Total Abundance") +
    scale_fill_manual(values=c("#5A4A6F", "#9D5A6C", "#E47250")) +
    scale_colour_manual(values = c("#5A4A6F", "#9D5A6C", "#E47250")) +
    theme_bw()
  GFN<-ggplot(df, aes(x=level, y=GFN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GFN.Mean-GFN.SE, ymax=GFN.Mean+GFN.SE), alpha=.1)+
    ggtitle("FALSE Negative, Gradient") +
    ylab("Mean FALSE Negative") + 
    xlab("Percent Variation in Total Abundance") +
    scale_fill_manual(values=c("#5A4A6F", "#9D5A6C", "#E47250")) +
    scale_colour_manual(values = c("#5A4A6F", "#9D5A6C", "#E47250")) +
    theme_bw()
  
  
  out<-NULL
  out$TP<-TP
  out$TN<-TN
  out$FP<-FP
  out$FN<-FN
  
  out$GTP<-GTP
  out$GTN<-GTN
  out$GFP<-GFP
  out$GFN<-GFN
  out
}

#' run statistics across simulations to see which differential abundance methods are equivalent
#' @param x list outputs from run.analysis
#' @keywords 
#' @export
#' @examples
#' make.plot.gradients()
diffabundStats<-function(x){
  df<-x$df[,c(5,13,21)]
  df<-melt(df)
  print(summary(aov(value~variable, data=df)))
  print(TukeyHSD(aov(value~variable, data=df)))
}

#' run statistics across simulations to see which correlation methods are equivalent
#' @param cortab list outputs from run.analysis
#' @keywords 
#' @export
#' @examples
#' make.plot.gradients()
correlationStats<-function(cortab){
  require(reshape2)
  df1<-melt(cortab$df[,c(1,3,5)])
  df2<-melt(cortab$df[,c(2,4,6)])
  out<-NULL
  out$anova1<-with(df1,summary(aov(value~variable)))
  out$tukey1<-with(df1,TukeyHSD(aov(value~variable)))
  out$anova2<-with(df2,summary(aov(value~variable)))
  out$tukey2<-with(df2,TukeyHSD(aov(value~variable)))
  out
}

#' summarize true/false predictions from differential abundance analysis
#' @param x list outputs from run.analysis
#' @keywords 
#' @export
#' @examples
#' getSummary.NG()
getSummary.NG<-function(x){
  outtab<-matrix(ncol=16, nrow=3)
  colnames(outtab)<-c("TP.Mean", "TP.SE", "FP.Mean", "FP.SE", "TN.Mean", "TN.SE", "FN.Mean", "FN.SE","GTP.Mean","GTP.SE", "GFP.Mean", "GFP.SE","GTN.Mean","GTN.SE", "GFN.Mean", "GFN.SE")#
  rownames(outtab)<-c("ALDEx", "ancomBC", "QSeq")
  ALDEx.TP<-c(rep(NA,length(x)))
  ALDEx.TN<-c(rep(NA,length(x)))
  ALDEx.FP<-c(rep(NA,length(x)))
  ALDEx.FN<-c(rep(NA,length(x)))
  
  ALDEx.GTP<-c(rep(NA,length(x)))
  ALDEx.GTN<-c(rep(NA,length(x)))
  ALDEx.GFP<-c(rep(NA,length(x)))
  ALDEx.GFN<-c(rep(NA,length(x)))
  
  ancom.TP<-c(rep(NA,length(x)))
  ancom.TN<-c(rep(NA,length(x)))
  ancom.FP<-c(rep(NA,length(x)))
  ancom.FN<-c(rep(NA,length(x)))
  
  ancom.GTP<-c(rep(NA,length(x)))
  ancom.GTN<-c(rep(NA,length(x)))
  ancom.GFP<-c(rep(NA,length(x)))
  ancom.GFN<-c(rep(NA,length(x)))
  
  QSeq.TP<-c(rep(NA,length(x)))
  QSeq.TN<-c(rep(NA,length(x)))
  QSeq.FP<-c(rep(NA,length(x)))
  QSeq.FN<-c(rep(NA,length(x)))
  
  QSeq.GTP<-c(rep(NA,length(x)))
  QSeq.GTN<-c(rep(NA,length(x)))
  QSeq.GFP<-c(rep(NA,length(x)))
  QSeq.GFN<-c(rep(NA,length(x)))
  
  for(i in 1:length(x)){
    ALDEx.TP[i]<-x[[i]]$effectaccuracy[1,1]
    ALDEx.TN[i]<-x[[i]]$effectaccuracy[1,2]
    ALDEx.FP[i]<-x[[i]]$effectaccuracy[1,3]
    ALDEx.FN[i]<-x[[i]]$effectaccuracy[1,4]
    
    ALDEx.GTP[i]<-x[[i]]$effectaccuracy[1,5]/sum(x[[i]]$diftab[,4]<0.05)
    ALDEx.GTP[i][is.na(ALDEx.GTP[i])]<-0
    ALDEx.GTN[i]<-x[[i]]$effectaccuracy[1,6]
    ALDEx.GFP[i]<-x[[i]]$effectaccuracy[1,7]/sum(x[[i]]$diftab[,4]>0.05)
    ALDEx.GFN[i]<-x[[i]]$effectaccuracy[1,8]
    
    ancom.TP[i]<-x[[i]]$effectaccuracy[5,1]
    ancom.TN[i]<-x[[i]]$effectaccuracy[5,2]
    ancom.FP[i]<-x[[i]]$effectaccuracy[5,3]
    ancom.FN[i]<-x[[i]]$effectaccuracy[5,4]
    
    ancom.GTP[i]<-x[[i]]$effectaccuracy[5,5]/sum(x[[i]]$diftab[,4]<0.05)
    ancom.GTP[i][is.na(ancom.GTP[i])]<-0
    ancom.GTN[i]<-x[[i]]$effectaccuracy[5,6]
    ancom.GFP[i]<-x[[i]]$effectaccuracy[5,7]/sum(x[[i]]$diftab[,4]>0.05)
    ancom.GFN[i]<-x[[i]]$effectaccuracy[5,8]
    
    QSeq.TP[i]<-x[[i]]$effectaccuracy[9,1]
    QSeq.TN[i]<-x[[i]]$effectaccuracy[9,2]
    QSeq.FP[i]<-x[[i]]$effectaccuracy[9,3]
    QSeq.FN[i]<-x[[i]]$effectaccuracy[9,4]
    
    QSeq.GTP[i]<-x[[i]]$effectaccuracy[10,5]/sum(x[[i]]$diftab[,4]<0.05)
    QSeq.GTP[i][is.na(QSeq.GTP[i])]<-0
    QSeq.GTN[i]<-x[[i]]$effectaccuracy[10,6]
    QSeq.GFP[i]<-x[[i]]$effectaccuracy[10,7]/sum(x[[i]]$diftab[,4]>0.05)
    QSeq.GFN[i]<-x[[i]]$effectaccuracy[10,8]
    
  }
  
  valdf<-data.frame(ALDEx.TP, ALDEx.TN, ALDEx.FP, ALDEx.FN, ALDEx.GTP,ALDEx.GTN,ALDEx.GFP,ALDEx.GFN,ancom.TP,ancom.TN,ancom.FP,ancom.FN,ancom.GTP,ancom.GTN,ancom.GFP,ancom.GFN,QSeq.TP,QSeq.TN,QSeq.FP,QSeq.FN, QSeq.GTP,QSeq.GTN,QSeq.GFP,QSeq.GFN)
  
  outtab[1,1]<-mean(ALDEx.TP)
  outtab[1,2]<-sd(ALDEx.TP)/sqrt(length(ALDEx.TP))
  outtab[1,3]<-mean(ALDEx.FP)
  outtab[1,4]<-sd(ALDEx.FP)/sqrt(length(ALDEx.FP))
  
  outtab[1,5]<-mean(ALDEx.TN)
  outtab[1,6]<-sd(ALDEx.TN)/sqrt(length(ALDEx.TN))
  outtab[1,7]<-mean(ALDEx.FN)
  outtab[1,8]<-sd(ALDEx.FN)/sqrt(length(ALDEx.FN))
  
  outtab[1,9]<-mean(ALDEx.GTP)
  outtab[1,10]<-sd(ALDEx.GTP)/sqrt(length(ALDEx.GTP))
  outtab[1,11]<-mean(ALDEx.GFP)
  outtab[1,12]<-sd(ALDEx.GFP)/sqrt(length(ALDEx.GFP))
  
  outtab[1,13]<-mean(ALDEx.GTN)
  outtab[1,14]<-sd(ALDEx.GTN)/sqrt(length(ALDEx.GTN))
  outtab[1,15]<-mean(ALDEx.GFN)
  outtab[1,16]<-sd(ALDEx.GFN)/sqrt(length(ALDEx.GFN))
  
  outtab[2,1]<-mean(ancom.TP)
  outtab[2,2]<-sd(ancom.TP)/sqrt(length(ancom.TP))
  outtab[2,3]<-mean(ancom.FP)
  outtab[2,4]<-sd(ancom.FP)/sqrt(length(ancom.FP))
  
  outtab[2,5]<-mean(ancom.TN)
  outtab[2,6]<-sd(ancom.TN)/sqrt(length(ancom.TN))
  outtab[2,7]<-mean(ancom.FN)
  outtab[2,8]<-sd(ancom.FN)/sqrt(length(ancom.FN))
  
  outtab[2,9]<-mean(ancom.GTP)
  outtab[2,10]<-sd(ancom.GTP)/sqrt(length(ancom.GTP))
  outtab[2,11]<-mean(ancom.GFP)
  outtab[2,12]<-sd(ancom.GFP)/sqrt(length(ancom.GFP))
  
  outtab[2,13]<-mean(ancom.GTN)
  outtab[2,14]<-sd(ancom.GTN)/sqrt(length(ancom.GTN))
  outtab[2,15]<-mean(ancom.GFN)
  outtab[2,16]<-sd(ancom.GFN)/sqrt(length(ancom.GFN))
  
  outtab[3,1]<-mean(QSeq.TP)
  outtab[3,2]<-sd(QSeq.TP)/sqrt(length(QSeq.TP))
  outtab[3,3]<-mean(QSeq.FP)
  outtab[3,4]<-sd(QSeq.FP)/sqrt(length(QSeq.FP))
  
  outtab[3,5]<-mean(QSeq.TN)
  outtab[3,6]<-sd(QSeq.TN)/sqrt(length(QSeq.TN))
  outtab[3,7]<-mean(QSeq.FN)
  outtab[3,8]<-sd(QSeq.FN)/sqrt(length(QSeq.FN))
  
  outtab[3,9]<-mean(QSeq.GTP)
  outtab[3,10]<-sd(QSeq.GTP)/sqrt(length(QSeq.GTP))
  outtab[3,11]<-mean(QSeq.GFP)
  outtab[3,12]<-sd(QSeq.GFP)/sqrt(length(QSeq.GFP))
  
  outtab[3,13]<-mean(QSeq.GTN)
  outtab[3,14]<-sd(QSeq.GTN)/sqrt(length(QSeq.GTN))
  outtab[3,15]<-mean(QSeq.GFN)
  outtab[3,16]<-sd(QSeq.GFN)/sqrt(length(QSeq.GFN))
  
  out<-NULL
  out$accuracytab<-outtab
  out$df<-valdf
  out
}

#' combine halves of two different correlation matrices
#' @param a matrix 1
#' @param b matrix 2
#' @keywords 
#' @export
#' @examples
#' combineMatrix()
combineMatrix<-function(a, b){
  new <- matrix(NA, nrow = nrow(a), ncol = ncol(a))
  new[upper.tri(new)] <- a[upper.tri(a)]
  new[lower.tri(new)] <- b[lower.tri(b)]
  new[is.na(new)]<-0
  rownames(new)<-rownames(a)
  colnames(new)<-colnames(a)
  new
}

#' make plots for total abundance estimation error
#' @param aeffect.list list of outputs from getSummary.NG
#' @keywords 
#' @export
#' @examples
#' make.plot3
make.plot3<-function(effect.list){
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  df<-ldply(effect.list, cbind)
  df$level<-c(rep(names(effect.list)[1], 3),rep(names(effect.list)[2], 3),rep(names(effect.list)[3], 3),rep(names(effect.list)[4], 3),rep(names(effect.list)[5], 3),rep(names(effect.list)[6], 3),rep(names(effect.list)[7], 3),rep(names(effect.list)[8], 3))
  df$TRT<-c(rep(c("ALDEx", "ancomBC", "QSeq")))
  #df1<-melt(df[,c(1,2)], id=c("Reference Mean", "Reference SE"))
  #df2<-melt(df[,c(3,4)], id=c("Ref.Gradient Mean", "Ref.Gradient SE"))
  #df1<-df[,c(2,3,6,7)]
  #print(df1)
  #df2<-df[,c(3,4,5)]       
  TP<-ggplot(df, aes(x=level, y=TP.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=TP.Mean-TP.SE, ymax=TP.Mean+TP.SE), alpha=.1)+
    ggtitle("TRUE Positive, No Gradient") +
    ylab("Mean True Positive") + 
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  TN<-ggplot(df, aes(x=level, y=TN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=TN.Mean-TN.SE, ymax=TN.Mean+TN.SE), alpha=.1)+
    ggtitle("TRUE Negative, No Gradient") +
    ylab("Mean True Negative") + 
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  
  GTP<-ggplot(df, aes(x=level, y=GTP.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GTP.Mean-GTP.SE, ymax=GTP.Mean+GTP.SE), alpha=.1)+
    ggtitle("TRUE Positive, Gradient") +
    ylab("Mean True Positive") + 
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  GTN<-ggplot(df, aes(x=level, y=GTN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GTN.Mean-GTN.SE, ymax=GTN.Mean+GTN.SE), alpha=.1)+
    ggtitle("TRUE Negative, Gradient") +
    ylab("Mean True Negative") + 
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  
  FP<-ggplot(df, aes(x=level, y=FP.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=FP.Mean-FP.SE, ymax=FP.Mean+FP.SE), alpha=.1)+
    ggtitle("FALSE Positive, No Gradient") +
    ylab("Mean FALSE Positive") + 
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  FN<-ggplot(df, aes(x=level, y=FN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=FN.Mean-FN.SE, ymax=FN.Mean+FN.SE), alpha=.1)+
    ggtitle("FALSE Negative, No Gradient") +
    ylab("Mean FALSE Negative") + 
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  
  GFP<-ggplot(df, aes(x=level, y=GFP.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GFP.Mean-GFP.SE, ymax=GFP.Mean+GFP.SE), alpha=.1)+
    ggtitle("FALSE Positive, Gradient") +
    ylab("Mean FALSE Positive") + 
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  GFN<-ggplot(df, aes(x=level, y=GFN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GFN.Mean-GFN.SE, ymax=GFN.Mean+GFN.SE), alpha=.1)+
    ggtitle("FALSE Negative, Gradient") +
    ylab("Mean FALSE Negative") + 
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  
  
  out<-NULL
  out$TP<-TP
  out$TN<-TN
  out$FP<-FP
  out$FN<-FN
  
  out$GTP<-GTP
  out$GTN<-GTN
  out$GFP<-GFP
  out$GFN<-GFN
  out
}

# Analysis ####

# import data from paper. To run on new data, skip this step:

# set working directory:
# setwd("~/Documents/GitHub/QSeq_Model/")

# init.dksp05
init.sksp05<-readRDS("Data/initsksp05.rds")
model.out.Gradient1<-readRDS("Data/model_gradient1.rds")
model.out.Gradient10<-readRDS("Data/model_gradient10.rds")
model.out.Gradient20<-readRDS("Data/model_gradient20.rds")
model.out.Gradient30<-readRDS("Data/model_gradient30.rds")
model.out.Gradient40<-readRDS("Data/model_gradient40.rds")

model.out.er0.01<-readRDS("Data/model_error01.RDS")
model.out.er0.02<-readRDS("Data/model_error02.RDS")
model.out.er0.03<-readRDS("Data/model_error03.RDS")
model.out.er0.05<-readRDS("Data/model_error05.RDS")
model.out.er0.1<-readRDS("Data/model_error10.RDS")
model.out.er0.2<-readRDS("Data/model_error20.RDS")
model.out.er0.3<-readRDS("Data/model_error30.RDS")
model.out.er0.5<-readRDS("Data/model_error50.RDS")

# lognormal distribution (supplemental figure ...):

# evaluate increasing QSeq measurement error (supplemental figure 4):
effectTab0.1<-getSummary.NG(model.out.er0.1)
effectTab0.2<-getSummary.NG(model.out.er0.2)
effectTab0.3<-getSummary.NG(model.out.er0.3)
effectTab0.5<-getSummary.NG(model.out.er0.5)

effectTab0.01<-getSummary.NG(model.out.er0.01)
effectTab0.02<-getSummary.NG(model.out.er0.02)
effectTab0.03<-getSummary.NG(model.out.er0.03)
effectTab0.05<-getSummary.NG(model.out.er0.05)

effectlist<-list("0.1"=effectTab0.1$accuracytab, "0.2"=effectTab0.2$accuracytab,"0.3"=effectTab0.3$accuracytab,"0.5"=effectTab0.5$accuracytab,"0.01"=effectTab0.01$accuracytab,"0.02"=effectTab0.02$accuracytab,"0.03"=effectTab0.03$accuracytab,"0.05"=effectTab0.05$accuracytab)

effectplots<-make.plot3(effectlist)

effectplots$TP # supplemental figure 4A ####
effectplots$FP
effectplots$TN
effectplots$FN

effectplots$GTP # supplemental figure 4B ####
effectplots$GFP
effectplots$GTN
effectplots$GFN

# Generate datasets ####
# skip if using data from paper

# establish gradient
grad1<-c(rep(1, 10)) # no variation
grad2<-sort(round(c(rnorm(10, 10, 1)))) # 10 % variation
grad3<-sort(round(c(rnorm(10, 10, 2)))) # 20 % variation
grad4<-sort(round(c(rnorm(10, 10, 3)))) # 30 % variation
grad5<-sort(round(c(rnorm(10, 10, 4)))) # 40 % variation

# initialize base community
init.sksp05<-initialize.sparseskew(nreps=100, dist = c(rlnorm(50, log(50), log(2.2))), sparse=25)

# impose gradients at different variances
gradient.test1<-impose.gradient2(init.sksp05, grad=0)
gradient.test10<-impose.gradient2(init.sksp05, grad=1)
gradient.test20<-impose.gradient2(init.sksp05, grad=2)
gradient.test30<-impose.gradient2(init.sksp05, grad=3)
gradient.test40<-impose.gradient2(init.sksp05, grad=4)

# analysis (start here if using data from publication) ####

# run diff abundance test at each level
model.out.Gradient1<-run.analysis(init.sksp05, gradient.test1, error = 0.05)
model.out.Gradient10<-run.analysis(init.sksp05, gradient.test10, error = 0.05)
model.out.Gradient20<-run.analysis(init.sksp05, gradient.test20, error = 0.05)
model.out.Gradient30<-run.analysis(init.sksp05, gradient.test30, error = 0.05)
model.out.Gradient40<-run.analysis(init.sksp05, gradient.test40, error = 0.05)

# aggregate summary stats 
effect.Gradient1<-getSummary.NG(model.out.Gradient1)
effect.Gradient10<-getSummary.NG(model.out.Gradient10)
effect.Gradient20<-getSummary.NG(model.out.Gradient20)
effect.Gradient30<-getSummary.NG(model.out.Gradient30)
effect.Gradient40<-getSummary.NG(model.out.Gradient40)

# plots/stats for paper

# gradient plot

Mload<-data.frame("Variation"=c(rep(1,10), rep(10,10), rep(20,10),rep(30,10),rep(40,10)), "Treatment"= c(c(rep(1,5), rep(2,5)),c(rep(1,5), rep(2,5)),c(rep(1,5), rep(2,5)),c(rep(1,5), rep(2,5)),c(rep(1,5), rep(2,5))), "value"=c(apply(gradient.test1$Rep1, MARGIN=2, sum),apply(gradient.test10$Rep1, MARGIN=2, sum),apply(gradient.test20$Rep1, MARGIN=2, sum),apply(gradient.test30$Rep1, MARGIN=2, sum),apply(gradient.test40$Rep1, MARGIN=2, sum))) 

Mload %>% group_by(Variation, Treatment) %>% summarize(mean=mean(value), SE=sd(value)/sqrt(5)) ->mload

mload$Treatment<-as.factor(mload$Treatment)

p <- ggplot(mload, aes(x=Variation, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=2,
                position=position_dodge(8))

p + theme_bw() + scale_fill_viridis_d(alpha=0.7) 

# Figure 2C ####
refcor0<-cor(t(gradient.test1$Rep1))
refcor40<-cor(t(gradient.test40$Rep1))
combined.refcor<-combineMatrix(refcor40, refcor0)
corrplot::corrplot(combined.refcor,tl.cex = 0.5)

# plot correlations
cortab.grad1<-compile.summary.table.correlations(model.out.Gradient1)
cortab.grad20<-compile.summary.table.correlations(model.out.Gradient20)
cortab.grad30<-compile.summary.table.correlations(model.out.Gradient30)
cortab.grad40<-compile.summary.table.correlations(model.out.Gradient40)
cortab.grad10<-compile.summary.table.correlations(model.out.Gradient10)

# cor stats ####
correlationStats(cortab.grad1)
correlationStats(cortab.grad10)
correlationStats(cortab.grad20)
correlationStats(cortab.grad30)
correlationStats(cortab.grad40)

# plot correlation
gradients.corlist<-list("1"=cortab.grad1$cortab,"10"=cortab.grad10$cortab,"20"=cortab.grad20$cortab,"30"=cortab.grad30$cortab,"40"=cortab.grad40$cortab)

grad.trtplots<-make.plot.gradients(gradients.corlist)
grad.trtplots$NG
grad.trtplots$GD # figure 2D ####

# plot differential abundance accuracy
effectlist.grad<-list("0"=effect.Gradient1$accuracytab,"10"=effect.Gradient10$accuracytab, "20"=effect.Gradient20$accuracytab, "30"=effect.Gradient30$accuracytab, "40"=effect.Gradient40$accuracytab)

effectplots.grad<-make.plot2(effectlist.grad)

effectplots.grad$TP
effectplots.grad$FP
effectplots.grad$TN
effectplots.grad$FN

effectplots.grad$GTP # plot figure 2 A ####
effectplots.grad$GFP # plot figure 2 B ####
effectplots.grad$GTN
effectplots.grad$GFN

# differential abundance stats ####

diffabundStats(effect.Gradient1)
diffabundStats(effect.Gradient10)
diffabundStats(effect.Gradient20)
diffabundStats(effect.Gradient30)
diffabundStats(effect.Gradient40)
