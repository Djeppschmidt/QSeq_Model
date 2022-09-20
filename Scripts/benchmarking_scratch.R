# other patterns ####

# initialize datasets
testinitialize<-initialize.data(nreps=100)
test.grad<-impose.gradient(testinitialize)
saveRDS(testinitialize,"~/Desktop/PhD/Compositionality/initialize.rds")
saveRDS(test.grad,"~/Desktop/PhD/Compositionality/test_gradient.rds")

# make gradient of skewness
init.skew1<-initialize.skew(nreps=100, dist=c(rlnorm(50, log(50), log(1)))) # sd 1
init.skew1.1<-initialize.skew(nreps=100, dist=c(rlnorm(50, log(50), log(1.1))))
init.skew1.2<-initialize.skew(nreps=100, dist=c(rlnorm(50, log(50), log(1.2)))) # sd 1.2
init.skew1.5<-initialize.skew(nreps=100, dist=c(rlnorm(50, log(50), log(1.5))))
init.skew1.8<-initialize.skew(nreps=100, dist=c(rlnorm(50, log(50), log(1.8)))) # sd 1.8
#init.skew2.0<-initialize.skew(nreps=100, dist=c(rlnorm(50, log(50), log(2.0))))
init.skew2.2<-initialize.skew(nreps=100, dist=c(rlnorm(50, log(50), log(2.2)))) # sd 2.2
init.skew2.5<-initialize.skew(nreps=100, dist=c(rlnorm(50, log(50), log(2.5))))
init.skew2.8<-initialize.skew(nreps=100, dist=c(rlnorm(50, log(50), log(2.8)))) # sd 2.8

skew1.grad<-impose.gradient(init.skew1)
saveRDS(init.skew1,"~/Desktop/PhD/Compositionality/initskew1.rds")
saveRDS(skew1.grad,"~/Desktop/PhD/Compositionality/skewgradient1.rds")

skew1.2.grad<-impose.gradient(init.skew1.2)
saveRDS(init.skew1.02,"~/Desktop/PhD/Compositionality/initskew1.2.rds")
saveRDS(skew1.02.grad,"~/Desktop/PhD/Compositionality/skewgradient1.2.rds")

skew1.8.grad<-impose.gradient(init.skew1.8)
saveRDS(init.skew1.8,"~/Desktop/PhD/Compositionality/initskew1.8.rds")
saveRDS(skew1.8.grad,"~/Desktop/PhD/Compositionality/skewgradient1.8.rds")

skew2.2.grad<-impose.gradient(init.skew2.2)
saveRDS(init.skew2.2,"~/Desktop/PhD/Compositionality/initskew2.2.rds")
saveRDS(skew2.2.grad,"~/Desktop/PhD/Compositionality/skewgradient2.2.rds")

skew2.8.grad<-impose.gradient(init.skew2.8)
saveRDS(init.skew2.8,"~/Desktop/PhD/Compositionality/initskew2.8.rds")
saveRDS(skew2.8.grad,"~/Desktop/PhD/Compositionality/skewgradient2.8.rds")

skew1.5.grad<-impose.gradient(init.skew1.5)
saveRDS(init.skew1.5,"~/Desktop/PhD/Compositionality/initskew1.5.rds")
saveRDS(skew1.5.grad,"~/Desktop/PhD/Compositionality/skewgradient1.5.rds")

skew1.1.grad<-impose.gradient(init.skew1.1)
saveRDS(init.skew1.1,"~/Desktop/PhD/Compositionality/initskew1.1.rds")
saveRDS(skew1.1.grad,"~/Desktop/PhD/Compositionality/skewgradient1.1.rds")

skew2.5.grad<-impose.gradient(init.skew2.5)
saveRDS(init.skew2.5,"~/Desktop/PhD/Compositionality/initskew2.5.rds")
saveRDS(skew2.5.grad,"~/Desktop/PhD/Compositionality/skewgradient2.5.rds")

# make gradient of sparsity on top of skewness
# there are 500 total; randomly replaces n=sparse values with 0
init.sksp01<-initialize.sparseskew(nreps=100, dist = c(rlnorm(50, log(50), log(2.2))), sparse=5) # 1 %
init.sksp05<-initialize.sparseskew(nreps=100, dist = c(rlnorm(50, log(50), log(2.2))), sparse=25) # 5 %
init.sksp10<-initialize.sparseskew(nreps=100, dist = c(rlnorm(50, log(50), log(2.2))), sparse=50) # 10 %
init.sksp20<-initialize.sparseskew(nreps=100, dist = c(rlnorm(50, log(50), log(2.2))), sparse=100) # 20 % NOT RUN: ERROR

saveRDS(init.sksp01, "~/Desktop/PhD/Compositionality/initsksp01.rds")
saveRDS(init.sksp05, "~/Desktop/PhD/Compositionality/initsksp05.rds")
saveRDS(init.sksp10, "~/Desktop/PhD/Compositionality/initsksp10.rds")
saveRDS(init.sksp20, "~/Desktop/PhD/Compositionality/initsksp20.rds")

sksp01.grad<-impose.gradient(init.sksp01)
sksp05.grad<-impose.gradient(init.sksp05)
sksp10.grad<-impose.gradient(init.sksp10)
sksp20.grad<-impose.gradient(init.sksp20)

saveRDS(sksp01.grad, "~/Desktop/PhD/Compositionality/sksp01grad.rds")
saveRDS(sksp05.grad, "~/Desktop/PhD/Compositionality/sksp05grad.rds")
saveRDS(sksp10.grad, "~/Desktop/PhD/Compositionality/sksp10grad.rds")
saveRDS(sksp20.grad, "~/Desktop/PhD/Compositionality/sksp20grad.rds")

sksp01.grad<-impose.gradient(init.sksp01)
saveRDS(init.sksp01,"~/Desktop/PhD/Compositionality/initsksp01.rds")
saveRDS(sksp01.grad,"~/Desktop/PhD/Compositionality/initsksp01_gradient.rds")

t1<-Sys.time()
model.out.er0.1<-run.analysis(testinitialize, test.grad, error=0.1)
saveRDS(model.out.er0.1, "~/Desktop/PhD/Compositionality/model_error10.RDS")
Sys.time()-t1

# make gradient of gradients
# need to be ordered so that differential abundance division later on is meaningful!!
grad1<-c(rep(1, 10)) # no variation
grad2<-sort(round(c(rnorm(10, 10, 1)))) # 10 % variation
grad3<-sort(round(c(rnorm(10, 10, 2)))) # 20 % variation
grad4<-sort(round(c(rnorm(10, 10, 3)))) # 30 % variation
grad5<-sort(round(c(rnorm(10, 10, 4)))) # 40 % variation

gradient.test1<-impose.gradient2(init.sksp05, grad=0)
gradient.test10<-impose.gradient2(init.sksp05, grad=1)
gradient.test20<-impose.gradient2(init.sksp05, grad=2)
gradient.test30<-impose.gradient2(init.sksp05, grad=3)
gradient.test40<-impose.gradient2(init.sksp05, grad=4)

# run model: ####

# test measurement (QPCR) error on perfectly even dataset
model.out.er0.01<-run.analysis(testinitialize, test.grad, error=0.01)
saveRDS(model.out.er0.01, "~/Desktop/PhD/Compositionality/model_error01.RDS") # run
model.out.er0.02<-run.analysis(testinitialize, test.grad, error=0.02)
saveRDS(model.out.er0.02, "~/Desktop/PhD/Compositionality/model_error02.RDS") # run
model.out.er0.03<-run.analysis(testinitialize, test.grad, error=0.03)
saveRDS(model.out.er0.03, "~/Desktop/PhD/Compositionality/model_error03.RDS")
model.out.er0.05<-run.analysis(testinitialize, test.grad, error=0.05)
saveRDS(model.out.er0.05, "~/Desktop/PhD/Compositionality/model_error05.RDS")
model.out.er0.1<-run.analysis(testinitialize, test.grad, error=0.1)
saveRDS(model.out.er0.1, "~/Desktop/PhD/Compositionality/model_error10.RDS")
model.out.er0.2<-run.analysis(testinitialize, test.grad, error=0.2)
saveRDS(model.out.er0.2, "~/Desktop/PhD/Compositionality/model_error20.RDS") # run
model.out.er0.3<-run.analysis(testinitialize, test.grad, error=0.3)
saveRDS(model.out.er0.3, "~/Desktop/PhD/Compositionality/model_error30.RDS") # run
model.out.er0.5<-run.analysis(testinitialize, test.grad, error=0.5)
saveRDS(model.out.er0.5, "~/Desktop/PhD/Compositionality/model_error50.RDS") # notrun

# test effect of increasing skewness
model.out.sk1<-run.analysis(init.skew1, skew1.grad, error=0.05)
saveRDS(model.out.sk1, "~/Desktop/PhD/Compositionality/model_sk1.RDS") # run
model.out.sk1.2<-run.analysis(init.skew1.2, skew1.2.grad, error=0.05)
saveRDS(model.out.sk1.2, "~/Desktop/PhD/Compositionality/model_sk102.RDS")
model.out.sk1.8<-run.analysis(init.skew1.8, skew1.8.grad, error=0.05)
saveRDS(model.out.sk1.8, "~/Desktop/PhD/Compositionality/model_sk18.RDS")
model.out.sk2.2<-run.analysis(init.skew2.2, skew2.2.grad, error=0.05)
saveRDS(model.out.sk2.2, "~/Desktop/PhD/Compositionality/model_sk22.RDS")
model.out.sk2.8<-run.analysis(init.skew2.8, skew2.8.grad, error=0.05)
saveRDS(model.out.sk2.8, "~/Desktop/PhD/Compositionality/model_sk2.8.RDS")
model.out.sk1.5<-run.analysis(init.skew1.5, skew1.5.grad, error=0.05)
saveRDS(model.out.sk1.5, "~/Desktop/PhD/Compositionality/model_sk1_5.RDS")
model.out.sk2.5<-run.analysis(init.skew2.5, skew2.5.grad, error=0.05)
saveRDS(model.out.sk2.5, "~/Desktop/PhD/Compositionality/model_sk2_5.RDS")

# test effect of increasing sparsity

model.out.spsk1<-run.analysis(init.sksp01, sksp01.grad, error = 0.05)
saveRDS(model.out.spsk1, "~/Desktop/PhD/Compositionality/model_spsk1.rds")
model.out.spsk5<-run.analysis(init.sksp05, sksp05.grad, error = 0.05)
saveRDS(model.out.spsk5, "~/Desktop/PhD/Compositionality/model_spsk5.rds")
model.out.spsk10<-run.analysis(init.sksp10, sksp10.grad, error = 0.05)
saveRDS(model.out.spsk10, "~/Desktop/PhD/Compositionality/model_spsk10.rds")
model.out.spsk20<-run.analysis(init.sksp20, sksp20.grad, error = 0.05)
saveRDS(model.out.spsk20, "~/Desktop/PhD/Compositionality/model_spsk20.rds")

# test gradient degree

model.out.Gradient1<-run.analysis(init.sksp05, gradient.test1, error = 0.05)
saveRDS(model.out.Gradient1, "~/Desktop/PhD/Compositionality/model_gradient1.rds")
model.out.Gradient10<-run.analysis(init.sksp05, gradient.test10, error = 0.05)
saveRDS(model.out.Gradient10, "~/Desktop/PhD/Compositionality/model_gradient10.rds")
model.out.Gradient20<-run.analysis(init.sksp05, gradient.test20, error = 0.05)
saveRDS(model.out.Gradient20, "~/Desktop/PhD/Compositionality/model_gradient20.rds")
model.out.Gradient30<-run.analysis(init.sksp05, gradient.test30, error = 0.05)
saveRDS(model.out.Gradient30, "~/Desktop/PhD/Compositionality/model_gradient30.rds")
model.out.Gradient40<-run.analysis(init.sksp05, gradient.test40, error = 0.05)
saveRDS(model.out.Gradient40, "~/Desktop/PhD/Compositionality/model_gradient40.rds")

init.sksp05<-readRDS("~/Desktop/PhD/Compositionality/initsksp05.rds")

#model.out.Gradient1<-readRDS("~/Desktop/PhD/Compositionality/model_gradient1.rds")
#model.out.Gradient10<-readRDS("~/Desktop/PhD/Compositionality/model_gradient10.rds")
#model.out.Gradient20<-readRDS("~/Desktop/PhD/Compositionality/model_gradient20.rds")
#model.out.Gradient30<-readRDS("~/Desktop/PhD/Compositionality/model_gradient30.rds")
#model.out.Gradient40<-readRDS("~/Desktop/PhD/Compositionality/model_gradient40.rds")

# error on gradient (not necessary ...)

model.out.Gradient1e10<-run.analysis(init.sksp05, gradient.test1, error = 0.1)
saveRDS(model.out.Gradient1e10, "~/Desktop/PhD/Compositionality/model_gradient1e10.rds")
model.out.Gradient10e10<-run.analysis(init.sksp05, gradient.test10, error = 0.1)
saveRDS(model.out.Gradient10e10, "~/Desktop/PhD/Compositionality/model_gradient10e10.rds")
model.out.Gradient20e10<-run.analysis(init.sksp05, gradient.test20, error = 0.1)
saveRDS(model.out.Gradient20e10, "~/Desktop/PhD/Compositionality/model_gradient20e10.rds")
model.out.Gradient30e10<-run.analysis(init.sksp05, gradient.test30, error = 0.1)
saveRDS(model.out.Gradient30e10, "~/Desktop/PhD/Compositionality/model_gradient30e10.rds")
model.out.Gradient40e10<-run.analysis(init.sksp05, gradient.test40, error = 0.1)
saveRDS(model.out.Gradient40e10, "~/Desktop/PhD/Compositionality/model_gradient40e10.rds")

model.out.Gradient1e20<-run.analysis(init.sksp05, gradient.test1, error = 0.2)
saveRDS(model.out.Gradient1e20, "~/Desktop/PhD/Compositionality/model_gradient1e20.rds")
model.out.Gradient10e20<-run.analysis(init.sksp05, gradient.test10, error = 0.2)
saveRDS(model.out.Gradient10e20, "~/Desktop/PhD/Compositionality/model_gradient10e20.rds")
model.out.Gradient20e20<-run.analysis(init.sksp05, gradient.test20, error = 0.2)
saveRDS(model.out.Gradient20e20, "~/Desktop/PhD/Compositionality/model_gradient20e20.rds")
model.out.Gradient30e20<-run.analysis(init.sksp05, gradient.test30, error = 0.2)
saveRDS(model.out.Gradient30e20, "~/Desktop/PhD/Compositionality/model_gradient30e20.rds")
model.out.Gradient40e20<-run.analysis(init.sksp05, gradient.test40, error = 0.2)
saveRDS(model.out.Gradient40e20, "~/Desktop/PhD/Compositionality/model_gradient40e20.rds")

# read in datasets ####

# skewed datasets
model.out.sk1<-readRDS("~/Desktop/PhD/Compositionality/model_sk1.RDS") # run
model.out.sk1.2<-readRDS("~/Desktop/PhD/Compositionality/model_sk102.RDS")
model.out.sk1.8<-readRDS("~/Desktop/PhD/Compositionality/model_sk18.RDS")
model.out.sk2.2<-readRDS("~/Desktop/PhD/Compositionality/model_sk22.RDS")
model.out.sk2.8<-readRDS("~/Desktop/PhD/Compositionality/model_sk2.8.RDS")
model.out.sk1.5<-readRDS("~/Desktop/PhD/Compositionality/model_sk1_5.RDS")
model.out.sk2.5<-readRDS("~/Desktop/PhD/Compositionality/model_sk2_5.RDS")

model.out.spsk1<-readRDS("~/Desktop/PhD/Compositionality/model_spsk1.rds")

# gradient datasets
model.out.Gradient1<-readRDS("~/Desktop/PhD/Compositionality/model_gradient1.rds")
model.out.Gradient10<-readRDS("~/Desktop/PhD/Compositionality/model_gradient10.rds")
model.out.Gradient20<-readRDS("~/Desktop/PhD/Compositionality/model_gradient20.rds")
model.out.Gradient30<-readRDS("~/Desktop/PhD/Compositionality/model_gradient30.rds")
model.out.Gradient40<-readRDS("~/Desktop/PhD/Compositionality/model_gradient40.rds")

# gradients with increased error:
model.out.Gradient1e10<-readRDS("~/Desktop/PhD/Compositionality/model_gradient1e10.rds")
model.out.Gradient10e10<-readRDS("~/Desktop/PhD/Compositionality/model_gradient10e10.rds")
model.out.Gradient20e10<-readRDS("~/Desktop/PhD/Compositionality/model_gradient20e10.rds")
model.out.Gradient30e10<-readRDS("~/Desktop/PhD/Compositionality/model_gradient30e10.rds")
model.out.Gradient40e10<-readRDS("~/Desktop/PhD/Compositionality/model_gradient40e10.rds")

model.out.Gradient1e20<-readRDS("~/Desktop/PhD/Compositionality/model_gradient1e20.rds") # error, doesn't exist
model.out.Gradient10e20<-readRDS("~/Desktop/PhD/Compositionality/model_gradient10e20.rds")
model.out.Gradient20e20<-readRDS("~/Desktop/PhD/Compositionality/model_gradient20e20.rds")
model.out.Gradient30e20<-readRDS("~/Desktop/PhD/Compositionality/model_gradient30e20.rds")
model.out.Gradient40e20<-readRDS("~/Desktop/PhD/Compositionality/model_gradient40e20.rds")

# plotting functions ####
# x = output of run.analysis
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

# outputs for paper ####
# re-import data if necessary

model.out.er0.01<-readRDS("~/Desktop/PhD/Compositionality/model_error01.RDS") # run
model.out.er0.02<-readRDS("~/Desktop/PhD/Compositionality/model_error02.RDS") # run
model.out.er0.03<-readRDS("~/Desktop/PhD/Compositionality/model_error03.RDS")
model.out.er0.05<-readRDS("~/Desktop/PhD/Compositionality/model_error05.RDS")
model.out.er0.1<-readRDS("~/Desktop/PhD/Compositionality/model_error10.RDS")
model.out.er0.2<-readRDS("~/Desktop/PhD/Compositionality/model_error20.RDS") # run
model.out.er0.3<-readRDS("~/Desktop/PhD/Compositionality/model_error30.RDS") # run
model.out.er0.5<-readRDS("~/Desktop/PhD/Compositionality/model_error50.RDS") #

# corrplot to demonstrate principles
refcor0<-cor(t(gradient.test1$Rep1))
refcor40<-cor(t(gradient.test40$Rep1))

# make combined correlation for Figure 2
combined.refcor<-combineMatrix(refcor40, refcor0)
corrplot::corrplot(combined.refcor,tl.cex = 0.5)

# Boxplot taxon cor correlations

cortab.1<-compile.summary.table.correlations(model.out.er0.1)
cortab.2<-compile.summary.table.correlations(model.out.er0.2)
cortab.3<-compile.summary.table.correlations(model.out.er0.3)
cortab.5<-compile.summary.table.correlations(model.out.er0.5)
cortab.01<-compile.summary.table.correlations(model.out.er0.01)
cortab.02<-compile.summary.table.correlations(model.out.er0.02)
cortab.03<-compile.summary.table.correlations(model.out.er0.03)
cortab.05<-compile.summary.table.correlations(model.out.er0.05)


cortab.1$cortab
cortab.2$cortab
cortab.1$df

# cortab for skewness:
cortab.sk1<-compile.summary.table.correlations(model.out.sk1)
#cortab.sk1.1<-compile.summary.table.correlations(model.out.sk1.1)
cortab.sk1.2<-compile.summary.table.correlations(model.out.sk1.2)
cortab.sk1.5<-compile.summary.table.correlations(model.out.sk1.5)
cortab.sk1.8<-compile.summary.table.correlations(model.out.sk1.8)
cortab.sk2.2<-compile.summary.table.correlations(model.out.sk2.2)
cortab.sk2.5<-compile.summary.table.correlations(model.out.sk2.5)
cortab.sk2.8<-compile.summary.table.correlations(model.out.sk2.8)

# skewness with sparsity
cortab.spsk1<-compile.summary.table.correlations(model.out.spsk1)
cortab.spsk5<-compile.summary.table.correlations(model.out.spsk5)
cortab.spsk10<-compile.summary.table.correlations(model.out.spsk10)

# gradients

cortab.grad1<-compile.summary.table.correlations(model.out.Gradient1)
cortab.grad20<-compile.summary.table.correlations(model.out.Gradient20)
cortab.grad30<-compile.summary.table.correlations(model.out.Gradient30)
cortab.grad40<-compile.summary.table.correlations(model.out.Gradient40)
cortab.grad10<-compile.summary.table.correlations(model.out.Gradient10)

# gradients with error

cortab.grad1e10<-compile.summary.table.correlations(model.out.Gradient1e10)
cortab.grad20e10<-compile.summary.table.correlations(model.out.Gradient20e10)
cortab.grad30e10<-compile.summary.table.correlations(model.out.Gradient30e10)
cortab.grad40e10<-compile.summary.table.correlations(model.out.Gradient40e10)
cortab.grad10e10<-compile.summary.table.correlations(model.out.Gradient10e10)

cortab.grad1e20<-compile.summary.table.correlations(model.out.Gradient1e20)
cortab.grad20e20<-compile.summary.table.correlations(model.out.Gradient20e20)
cortab.grad30e20<-compile.summary.table.correlations(model.out.Gradient30e20)
cortab.grad40e20<-compile.summary.table.correlations(model.out.Gradient40e20)
cortab.grad10e20<-compile.summary.table.correlations(model.out.Gradient10e20)



corstats.1<-correlationStats(cortab.1)
corstats.2<-correlationStats(cortab.2)
corstats.3<-correlationStats(cortab.3)
corstats.5<-correlationStats(cortab.5)

corstats.01<-correlationStats(cortab.01)
corstats.02<-correlationStats(cortab.02)
corstats.03<-correlationStats(cortab.03)
corstats.05<-correlationStats(cortab.05)

# corelation stats for skewness
corstats.sk1<-correlationStats(cortab.sk1)
corstats.sk1.2<-correlationStats(cortab.sk1.2)
corstats.sk1.5<-correlationStats(cortab.sk1.5)
corstats.sk1.8<-correlationStats(cortab.sk1.8)
corstats.sk2.2<-correlationStats(cortab.sk2.2)
corstats.sk2.5<-correlationStats(cortab.sk2.5)
corstats.sk2.8<-correlationStats(cortab.sk2.8)

# stats for skewness + sparsity
corstats.spsk1<-correlationStats(cortab.spsk1)
corstats.spsk5<-correlationStats(cortab.spsk5)
corstats.spsk10<-correlationStats(cortab.spsk10)

# stats for gradients 

corstats.grad1<-correlationStats(cortab.grad1)
corstats.grad10<-correlationStats(cortab.grad10)
corstats.grad20<-correlationStats(cortab.grad20)
corstats.grad30<-correlationStats(cortab.grad30)
corstats.grad40<-correlationStats(cortab.grad40)

corstats.grad1e10<-correlationStats(cortab.grad1e10)
corstats.grad10e10<-correlationStats(cortab.grad10e10)
corstats.grad20e10<-correlationStats(cortab.grad20e10)
corstats.grad30e10<-correlationStats(cortab.grad30e10)
corstats.grad40e10<-correlationStats(cortab.grad40e10)

corstats.grad1e20<-correlationStats(cortab.grad1e20)
corstats.grad10e20<-correlationStats(cortab.grad10e20)
corstats.grad20e20<-correlationStats(cortab.grad20e20)
corstats.grad30e20<-correlationStats(cortab.grad30e20)
corstats.grad40e20<-correlationStats(cortab.grad40e20)

# cortab.list = list of cortab$cortab
make.plot<-function(cortab.list){
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  df<-ldply(cortab.list, cbind)
  df$level<-c(rep(names(cortab.list)[1], 3),rep(names(cortab.list)[2], 3),rep(names(cortab.list)[3], 3),rep(names(cortab.list)[4], 3),rep(names(cortab.list)[5], 3),rep(names(cortab.list)[6], 3),rep(names(cortab.list)[7], 3),rep(names(cortab.list)[8], 3))
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
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  gplt2<-ggplot(df, aes(x=level, y=`Ref.Gradient Mean`, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=`Ref.Gradient Mean`-`Ref.Gradient SE`, ymax=`Ref.Gradient Mean`+`Ref.Gradient SE`, fill=TRT), alpha=.1) +
    ggtitle("Gradient") +
    scale_fill_manual(breaks=c("Subsample", "FastSpar", "QSeq"), values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  out<-NULL
  out$NG<-gplt1
  out$GD<-gplt2
  out
}
library(plyr)
corlist<-list("0.1"=cortab.1$cortab, "0.2"=cortab.2$cortab,"0.3"=cortab.3$cortab,"0.5"=cortab.5$cortab,"0.01"=cortab.01$cortab,"0.02"=cortab.02$cortab,"0.03"=cortab.03$cortab,"0.05"=cortab.05$cortab)

trtplots<-make.plot(corlist)
trtplots$NG
trtplots$GD

# make plot for skew
make.plot.sk<-function(cortab.list){
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  df<-ldply(cortab.list, cbind)
  df$level<-c(rep(names(cortab.list)[1], 3),rep(names(cortab.list)[2], 3),rep(names(cortab.list)[3], 3),rep(names(cortab.list)[4], 3),rep(names(cortab.list)[5], 3),rep(names(cortab.list)[6], 3),rep(names(cortab.list)[7], 3))
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
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  gplt2<-ggplot(df, aes(x=level, y=`Ref.Gradient Mean`, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=`Ref.Gradient Mean`-`Ref.Gradient SE`, ymax=`Ref.Gradient Mean`+`Ref.Gradient SE`, fill=TRT), alpha=.1) +
    ggtitle("Gradient") +
    scale_fill_manual(breaks=c("Subsample", "FastSpar", "QSeq"), values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  out<-NULL
  out$NG<-gplt1
  out$GD<-gplt2
  out
}

sk.corlist<-list("1.0"=cortab.sk1$cortab,"1.2"=cortab.sk1.2$cortab,"1.5"=cortab.sk1.5$cortab,"1.8"=cortab.sk1.8$cortab,"2.2"=cortab.sk2.2$cortab,"2.5"=cortab.sk2.5$cortab,"2.8"=cortab.sk2.8$cortab)

sk.trtplots<-make.plot.sk(sk.corlist)
sk.trtplots$NG
sk.trtplots$GD

# make plot for sparsity and skew

make.plot.spsk<-function(cortab.list){
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  df<-ldply(cortab.list, cbind)
  df$level<-c(rep(names(cortab.list)[1], 3),rep(names(cortab.list)[2], 3),rep(names(cortab.list)[3], 3))
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
    xlab("Total Abundance Estimation Error (SD)") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  gplt2<-ggplot(df, aes(x=level, y=`Ref.Gradient Mean`, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=`Ref.Gradient Mean`-`Ref.Gradient SE`, ymax=`Ref.Gradient Mean`+`Ref.Gradient SE`, fill=TRT), alpha=.1) +
    ggtitle("Gradient") +
    scale_fill_manual(breaks=c("Subsample", "FastSpar", "QSeq"), values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  out<-NULL
  out$NG<-gplt1
  out$GD<-gplt2
  out
}

spsk.corlist<-list("1.0"=cortab.spsk1$cortab,"5"=cortab.spsk5$cortab,"10"=cortab.spsk10$cortab)

spsk.trtplots<-make.plot.spsk(spsk.corlist)
spsk.trtplots$NG
spsk.trtplots$GD


# gradients

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

gradients.corlist<-list("1"=cortab.grad1$cortab,"10"=cortab.grad10$cortab,"20"=cortab.grad20$cortab,"30"=cortab.grad30$cortab,"40"=cortab.grad40$cortab)

grad.trtplots<-make.plot.gradients(gradients.corlist)
grad.trtplots$NG
grad.trtplots$GD

# gradients with estimation error variation!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

make.plot.gradients<-function(cortab.list){
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  df<-as.data.frame(do.call(rbind, cortab.list))
  df$level<-c(rep(0, 3),rep(10, 3),rep(20, 3),rep(30, 3),rep(40, 3),rep(c(0,10,20,30,40), 2))
  df$TRT<-c(rep(c("Subsample", "FastSpar", "QSeq e-5"), 5), rep("QSeq e-10", 5), rep("QSeq e-20",5))
  df$TRT<-factor(df$TRT, levels = c("Subsample", "FastSpar", "QSeq e-5","QSeq e-10","QSeq e-20"))    
  gplt1<-ggplot(df, aes(x=level, y=`Reference Mean`, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=`Reference Mean`-`Reference SE`, ymax=`Reference Mean`+`Reference SE`), alpha=.1)+
    ggtitle("No Gradient") +
    ylab("Mean Correlation") + 
    xlab("Variance in Total Abundance") +
    scale_fill_manual(values=c("#808080","#696969","#FF0000", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#808080","#696969","#FF0000", "#E47250", "#9D5A6C")) +
    theme_bw()
  gplt2<-ggplot(df, aes(x=level, y=`Ref.Gradient Mean`, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=`Ref.Gradient Mean`-`Ref.Gradient SE`, ymax=`Ref.Gradient Mean`+`Ref.Gradient SE`, fill=TRT), alpha=.1) +
    ggtitle("Gradient") +
    ylab("Mean Correlation") + 
    xlab("Percent Variation in Total Abundance") +
    scale_fill_manual(breaks=c("Subsample", "FastSpar", "QSeq"), values=c("#808080","#696969","#FF0000", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#808080","#696969","#FF0000", "#E47250", "#9D5A6C")) +
    theme_bw()
  out<-NULL
  out$NG<-gplt1
  out$GD<-gplt2
  out
}

gradients.corlist<-list("1e5"=cortab.grad1$cortab,"10e5"=cortab.grad10$cortab,"20e5"=cortab.grad20$cortab,"30e5"=cortab.grad30$cortab,"40e5"=cortab.grad40$cortab,"1e10"=cortab.grad1e10$cortab[3,],"10e10"=cortab.grad10e10$cortab[3,],"20e10"=cortab.grad20e10$cortab[3,],"30e10"=cortab.grad30e10$cortab[3,],"40e10"=cortab.grad40e10$cortab[3,],"1e20"=cortab.grad1e20$cortab[3,],"10e20"=cortab.grad10e20$cortab[3,],"20e20"=cortab.grad20e20$cortab[3,],"30e20"=cortab.grad30e20$cortab[3,],"40e20"=cortab.grad40e20$cortab[3,])

grad.trtplots<-make.plot.gradients(gradients.corlist)
grad.trtplots$NG
grad.trtplots$GD

# plots of T/F detection rates with increasing error
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

# extract summaries and make list input for summary tables and figures
# estimation error
effectTab0.1<-getSummary.NG(model.out.er0.1)
effectTab0.2<-getSummary.NG(model.out.er0.2)
effectTab0.3<-getSummary.NG(model.out.er0.3)
effectTab0.5<-getSummary.NG(model.out.er0.5)

effectTab0.01<-getSummary.NG(model.out.er0.01)
effectTab0.02<-getSummary.NG(model.out.er0.02)
effectTab0.03<-getSummary.NG(model.out.er0.03)
effectTab0.05<-getSummary.NG(model.out.er0.05)

effectlist<-list("0.1"=effectTab0.1$accuracytab, "0.2"=effectTab0.2$accuracytab,"0.3"=effectTab0.3$accuracytab,"0.5"=effectTab0.5$accuracytab,"0.01"=effectTab0.01$accuracytab,"0.02"=effectTab0.02$accuracytab,"0.03"=effectTab0.03$accuracytab,"0.05"=effectTab0.05$accuracytab)

# skewness
effectsk1<-getSummary.NG(model.sk1)
effectsk11<-getSummary.NG(model.sk11)
effectsk15<-getSummary.NG(model.sk15)
effectsk18<-getSummary.NG(model.sk18)
effectsk22<-getSummary.NG(model.sk22)
effectsk25<-getSummary.NG(model.sk25)
effectsk28<-getSummary.NG(model.sk28)


effectlist.sk<-list("1.0"=effectsk1$accuracytab, "1.1"=effectsk11$accuracytab,"1.5"=effectsk15$accuracytab,"1.8"=effectsk18$accuracytab,"2.2"=effectsk22$accuracytab,"2.5"=effectsk25$accuracytab,"2.8"=effectsk28$accuracytab)

# sparsity on skewness

effectsksp1<-getSummary.NG(model.out.spsk1)
effectsksp5<-getSummary.NG(model.out.spsk5)
effectsksp10<-getSummary.NG(model.out.spsk10)
effectlist.sksp<-list("1"=effectsksp1$accuracytab, "5"=effectsksp5$accuracytab,"10"=effectsksp10$accuracytab)

# increasing Gradients

effect.Gradient1<-getSummary.NG(model.out.Gradient1)
#effect.Gradient1$df[is.na(effect.Gradient1$df)]<-0

effect.Gradient10<-getSummary.NG(model.out.Gradient10)
effect.Gradient20<-getSummary.NG(model.out.Gradient20)
effect.Gradient30<-getSummary.NG(model.out.Gradient30)
effect.Gradient40<-getSummary.NG(model.out.Gradient40)

effect.Gradient1e10<-getSummary.NG(model.out.Gradient1e10)
effect.Gradient10e10<-getSummary.NG(model.out.Gradient10e10)
effect.Gradient20e10<-getSummary.NG(model.out.Gradient20e10)
effect.Gradient30e10<-getSummary.NG(model.out.Gradient30e10)
effect.Gradient40e10<-getSummary.NG(model.out.Gradient40e10)

effect.Gradient1e20<-getSummary.NG(model.out.Gradient1e20)
effect.Gradient10e20<-getSummary.NG(model.out.Gradient10e20)
effect.Gradient20e20<-getSummary.NG(model.out.Gradient20e20)
effect.Gradient30e20<-getSummary.NG(model.out.Gradient30e20)
effect.Gradient40e20<-getSummary.NG(model.out.Gradient40e20)

# stats for differential abundance: ####

diffabundStats<-function(x){
  df<-x$df[,c(5,13,21)]
  df<-melt(df)
  print(summary(aov(value~variable, data=df)))
  print(TukeyHSD(aov(value~variable, data=df)))
}
diffabundStats(effect.Gradient1)
diffabundStats(effect.Gradient10)
diffabundStats(effect.Gradient20)
diffabundStats(effect.Gradient30)
diffabundStats(effect.Gradient40)
# plot true positives / false positives (figure 1 in paper)
make.plot2<-function(effect.list){
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  df<-ldply(effect.list, cbind)
  df$level<-c(rep(names(effect.list)[1], 3),rep(names(effect.list)[2], 3),rep(names(effect.list)[3], 3),rep(names(effect.list)[4], 3),rep(names(effect.list)[5], 3))
  df$TRT<-c(rep(c("ALDEx", "ancomBC", "QSeq")))
  #df$GTP.Mean<-df$GTP.Mean/c(rep(mean(unlist(lapply(model.out.Gradient1, function(x) sum(x$diftab[,4]<0.05)))),3), rep(mean(unlist(lapply(model.out.Gradient10, function(x) sum(x$diftab[,4]<0.05)))),3),rep(mean(unlist(lapply(model.out.Gradient20, function(x) sum(x$diftab[,4]<0.05)))),3),rep(mean(unlist(lapply(model.out.Gradient30, function(x) sum(x$diftab[,4]<0.05)))),3),rep(mean(unlist(lapply(model.out.Gradient40, function(x) sum(x$diftab[,4]<0.05)))),3))
  
  #df$GTP.SE<-df$GTP.SE/c(rep(mean(unlist(lapply(model.out.Gradient1, function(x) sum(x$diftab[,4]<0.05)))),3),
  #                           rep(mean(unlist(lapply(model.out.Gradient10, function(x) sum(x$diftab[,4]<0.05)))),3),
  #                           rep(mean(unlist(lapply(model.out.Gradient20, function(x) sum(x$diftab[,4]<0.05)))),3),
  #                           rep(mean(unlist(lapply(model.out.Gradient30, function(x) sum(x$diftab[,4]<0.05)))),3),
  #                           rep(mean(unlist(lapply(model.out.Gradient40, function(x) sum(x$diftab[,4]<0.05)))),3))
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

effectlist.grad<-list("0"=effect.Gradient1$accuracytab,"10"=effect.Gradient10$accuracytab, "20"=effect.Gradient20$accuracytab, "30"=effect.Gradient30$accuracytab, "40"=effect.Gradient40$accuracytab)

effectplots.grad<-make.plot2(effectlist.grad)

effectplots.grad$TP
effectplots.grad$FP
effectplots.grad$TN
effectplots.grad$FN

effectplots.grad$GTP # plot figure 2 A
effectplots.grad$GFP # plot figure 2 B
effectplots.grad$GTN
effectplots.grad$GFN
# plot true positives / false positives
make.plot2<-function(effect.list){
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

effectplots<-make.plot2(effectlist)

effectplots$TP
effectplots$FP
effectplots$TN
effectplots$FN

effectplots$GTP
effectplots$GFP
effectplots$GTN
effectplots$GFN

effects.skew<-make.plot2(effectlist.sk)

make.plot2.spsk<-function(effect.list){
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  df<-ldply(effect.list, cbind)
  df$level<-factor(c(rep(names(effect.list)[1], 3),rep(names(effect.list)[2], 3),rep(names(effect.list)[3], 3)), levels=c(names(effect.list)))
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
    xlab("percent sparsity") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  TN<-ggplot(df, aes(x=level, y=TN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=TN.Mean-TN.SE, ymax=TN.Mean+TN.SE), alpha=.1)+
    ggtitle("TRUE Negative, No Gradient") +
    ylab("Mean True Negative") + 
    xlab("percent sparsity") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  
  GTP<-ggplot(df, aes(x=level, y=GTP.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GTP.Mean-GTP.SE, ymax=GTP.Mean+GTP.SE), alpha=.1)+
    ggtitle("TRUE Positive, Gradient") +
    ylab("Mean True Positive") + 
    xlab("percent sparsity") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  GTN<-ggplot(df, aes(x=level, y=GTN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GTN.Mean-GTN.SE, ymax=GTN.Mean+GTN.SE), alpha=.1)+
    ggtitle("TRUE Negative, Gradient") +
    ylab("Mean True Negative") + 
    xlab("percent sparsity") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  
  FP<-ggplot(df, aes(x=level, y=FP.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=FP.Mean-FP.SE, ymax=FP.Mean+FP.SE), alpha=.1)+
    ggtitle("FALSE Positive, No Gradient") +
    ylab("Mean FALSE Positive") + 
    xlab("percent sparsity") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  FN<-ggplot(df, aes(x=level, y=FN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=FN.Mean-FN.SE, ymax=FN.Mean+FN.SE), alpha=.1)+
    ggtitle("FALSE Negative, No Gradient") +
    ylab("Mean FALSE Negative") + 
    xlab("percent sparsity") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  
  GFP<-ggplot(df, aes(x=level, y=GFP.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GFP.Mean-GFP.SE, ymax=GFP.Mean+GFP.SE), alpha=.1)+
    ggtitle("FALSE Positive, Gradient") +
    ylab("Mean FALSE Positive") + 
    xlab("percent sparsity") +
    scale_fill_manual(values=c("#5A4A6F", "#E47250", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250", "#9D5A6C")) +
    theme_bw()
  GFN<-ggplot(df, aes(x=level, y=GFN.Mean, group=TRT))+
    geom_line(aes(color=TRT)) +
    geom_point(aes(color=TRT)) + 
    geom_ribbon(aes(ymin=GFN.Mean-GFN.SE, ymax=GFN.Mean+GFN.SE), alpha=.1)+
    ggtitle("FALSE Negative, Gradient") +
    ylab("Mean FALSE Negative") + 
    xlab("percent sparsity") +
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
effects.sksp<-make.plot2.spsk(effectlist.sksp)
effects.sksp$TP
effects.sksp$FP
effects.sksp$GTP
effects.sksp$GFP



# true effects:

plot.references<-function(ref.list){
  require(ggplot2)
  #require()
  GTD<-unlist(lapply(ref.list, function(x) sum(x$diftab[,4]<0.05)))
  #print(GTD)
  GND<-unlist(lapply(ref.list, function(x) sum(x$diftab[,4]>0.05)))
  
  TD<-unlist(lapply(ref.list, function(x) sum(x$diftab[,2]<0.05)))
  ND<-unlist(lapply(ref.list, function(x) sum(x$diftab[,2]>0.05)))
  
  #df<-data.frame("Value"=c(mean(TD), mean(ND), mean(GTD), mean(GND)),"SE"=c(sd(TD)/sqrt(100), sd(ND)/sqrt(100), sd(GTD)/sqrt(100), sd(GND)/sqrt(100)) ,"Gradient"=c(rep("none",2), rep("Gradient", 2)), "Effect"=c("TRUE","FALSE","TRUE","FALSE"))
  df<-data.frame("Value"=c(TD, ND, GTD, GND) ,"Category"=factor(c(rep("Compositional T",100), rep("Compositional F", 100), rep("Gradient T", 100), rep("Gradient F", 100)), levels=c("Compositional T", "Compositional F", "Gradient T", "Gradient F")))
  
  p<-ggplot(df, aes(x=Category, y=Value))+
    geom_boxplot(aes(color=Category)) +
    #geom_point(aes(color=Gradient)) + 
    #geom_ribbon(aes(ymin=Value-SE, ymax=Value+SE), alpha=.1)+
    ggtitle("TRUE Positive and Negative Rates") +
    ylab("Effect (Positive vs Negative)") + 
    xlab("Number of Taxa (total is 50)") +
    scale_fill_manual(values=c("#E47250", "#9D5A6C","#336699","#5A4A6F")) +
    scale_colour_manual(values = c("#E47250", "#9D5A6C","#336699","#5A4A6F")) +
    theme_bw()
  
  p
  
}

plotreference<-plot.references(model.out.Gradient40)
plotreference

plot.references(model.out.Gradient10)
plot.references(model.out.Gradient20)
plot.references(model.out.Gradient30)
plot.references(model.out.Gradient40)



# debug workspace ####
testinitialize<-initialize.data(nreps=100)
test.grad<-impose.gradient(testinitialize)
testwrapper<-analysis.guts(testinitialize[[1]],test.grad[[1]], error=0.1)

testwrapper$effectaccuracy
testwrapper$coraccuracy
testwrapper$diftab
table(testwrapper$diftab[,4]<0.05) # test gives 44 T, 6 F


colSums(testinitialize[[1]])+rnorm(10,0,0.1*colSums(testinitialize[[1]]))

testcortab<-compile.summary.table.correlations(model.out.er0.1)
testcortab$cortab

# confidence index

cdf<-data.frame("reads"=c(5000,7500,11000,15000,18000,21500,24500,28000), "cutoff"=c(0.093,0.082,0.06,0.054,0.048,0.043
                                                                                     ,0.038,0.038))
library(aomisc)
library(drc)

#model <- drm(cutoff~reads, fct = DRC.expoDecay(), data=cdf)
#model <- drm(cutoff~reads, fct = DRC.asymReg(), data=cdf)
confidence.model <- drm(cutoff~reads, fct = DRC.powerCurve(), data=cdf)
#model <- drm(cutoff~reads, fct = DRC.logCurve(), data=cdf)

plot(model)

# proof of concept
pred<-data.frame("reads"=c(2000,10000,100000,5000000))
predict(model, newdata = pred)

# differential abundance testing with and without filter criteria

