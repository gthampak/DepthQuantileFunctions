## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- eval=FALSE, include=FALSE,results='hide'------------------------------------------------------------------------------------------------
## source("00. datasets.R")
## source("01. dqf-outlier.R")
## source("02. dqf-subset.R")
## source("03. extract-dqf-subset.R")
## source("10. dqf-transforms.R")
## source("11. dqf-metrics.R")


## ---- eval=FALSE, include=FALSE,results='hide'------------------------------------------------------------------------------------------------
## knitr::purl("40. dqf-diagnostics.Rmd")


## ---------------------------------------------------------------------------------------------------------------------------------------------
dqf.diagnostics <- function(dqf,labels=NULL){
  if(is.null(labels)){
    labels <- rep(1,nrow(dqf))
  }
  
  dqf.maxed <- scale.dqf.max(dqf)
  dqf.summed <- scale.dqf.sum(dqf)
  
  par(mfrow=c(1,3))
  plot.dqf(dqf,labels=labels,main="DQF")
  plot.dqf(dqf.maxed,labels=labels,main="DQF Scaled to Function Max")
  plot.dqf(dqf.summed,labels=labels,main="DQF Scaled to Function Sum")
  
  plot.dqf(dqf,labels=labels,main="DQF - with Bounds")
  draw.mean.bounds(dqf.mean(dqf),dqf.sd(dqf))
  plot.dqf(dqf.maxed,labels=labels,main="DQF Scaled to Function Max")
  draw.mean.bounds(dqf.mean(dqf.maxed),dqf.sd(dqf.maxed))
  plot.dqf(dqf.summed,labels=labels,main="DQF Scaled to Function Sum")
  draw.mean.bounds(dqf.mean(dqf.summed),dqf.sd(dqf.summed))
  
  # par(mfrow=c(1,3))
  plot.dqf.2norm(dqf,n.sd=2,main="2-norm (DQF,Mean(DQF))")
  plot.dqf.2norm(dqf.maxed,n.sd=2,main="(Maxed)")
  plot.dqf.2norm(dqf.summed,n.sd=2,main="(Summed)")
  
  # par(mfrow=c(1,3))
  plot.dqf.supnorm(dqf,n.sd=2,main="sup-norm (DQF,Mean(DQF))")
  plot.dqf.supnorm(dqf.maxed,n.sd=2,main="(Maxed)")
  plot.dqf.supnorm(dqf.summed,n.sd=2,main="(Summed)")
  
  # par(mfrow=c(1,3))
  plot.prop.outside.bounds(dqf,threshold=.2,n.sd=2,main="Prop. of dqf Outisde Bounds")
  plot.prop.outside.bounds(dqf.maxed,threshold=.2,n.sd=2,main="(Maxed)")
  plot.prop.outside.bounds(dqf.maxed,threshold=.2,n.sd=2,main="(Summed)")
  
  # par(mfrow=c(1,3))
  plot.zscore.dqf(dqf,labels,main="dqf z-scores")
  plot.zscore.dqf(dqf.maxed,labels,main="(Maxed)")
  plot.zscore.dqf(dqf.summed,labels,main="(Summed)")
  
  # par(mfrow=c(1,3))
  plot.mean.zscores(dqf,labels,main="mean z-score")
  plot.mean.zscores(dqf.maxed,labels,"(maxed)")
  plot.mean.zscores(dqf.summed,labels,"(summed)")
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
dqf.c.diagnostics <- function(dqf,labels=NULL){
  if(is.null(labels)){
    labels <- rep(1,nrow(dqf))
  }
  
  dqf.maxed <- scale.dqf.max(dqf)
  dqf.summed <- scale.dqf.sum(dqf)
  
  par(mfrow=c(2,3))
  
  plot.dqf(dqf,labels=labels,main="DQF - with Bounds")
  draw.mean.bounds(dqf.mean(dqf),dqf.sd(dqf))
  plot.dqf(dqf.summed,labels=labels,main="DQF Scaled to Function Sum")
  draw.mean.bounds(dqf.mean(dqf.summed),dqf.sd(dqf.summed))
  
  plot.zscore.dqf(dqf,labels,main="dqf z-scores")
  plot.zscore.dqf(dqf.summed,labels,main="(Summed)")
  
  plot.mean.zscores(dqf,labels,main="mean z-score")
  plot.mean.zscores(dqf.summed,labels,main="(summed)")
  
}

