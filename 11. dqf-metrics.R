## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- eval=FALSE, include=FALSE,results='hide'------------------------------------------------------------------------------------------------
## source("00. datasets.R")
## source("01. dqf-outlier.R")
## source("02. dqf-subset.R")
## source("03. extract-dqf-subset.R")
## source("10. dqf-transforms.R")


## ---- eval=FALSE, include=FALSE,results='hide'------------------------------------------------------------------------------------------------
## knitr::purl("11. dqf-metrics.Rmd")


## ---------------------------------------------------------------------------------------------------------------------------------------------
dqf.mean <- function(dqf){
  means <- c()
  
  for(i in 1:length(dqf[1,])){
    means <- c(means, mean(dqf[,i]))
  }
  
  return(means)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
dqf.sd <- function(dqf){
  sds <- c()
  
  for(i in 1:length(dqf[1,])){
    sds <- c(sds, sd(dqf[,i]))
  }
  
  return(sds)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
dqf.upperbound <- function(mean,sd,n.sd=2){
  return(mean+n.sd*sd)
}

dqf.lowerbound <- function(mean,sd,n.sd=2){
  return(mean-n.sd*sd)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
draw.mean.bounds <- function(mean, sd, n.sd=2){
  x <- seq(.01,1,.01)
  lines(x, mean, col='purple',lwd=2.0)
  lines(x,dqf.lowerbound(mean,sd,n.sd),col='purple',lwd=1.0)
  lines(x,dqf.upperbound(mean,sd,n.sd),col='purple',lwd=1.0)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
func.2norm <- function(dqf1,dqf2){
  sqrt(sum((dqf1-dqf2)^2))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
dqf.2norm <- function(dqf){
  dqf.2norm <- c()
  
  mean <- dqf.mean(dqf)
  sd <- dqf.sd(dqf)
  
  for(i in 1:length(dqf[,1])){
    norm <- func.2norm(mean,dqf[i,])
    dqf.2norm <- c(dqf.2norm, norm)
  }
  
  return(dqf.2norm)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
plot.dqf.2norm <- function(dqf, n.sd=2,main=""){
  dqf.2norm <- c()
  labels <- c()
  
  mean <- dqf.mean(dqf)
  sd <- dqf.sd(dqf)
  
  bound.dqf.2norm <- func.2norm(dqf.upperbound(mean,sd,n.sd),mean)
  
  for(i in 1:length(dqf[,1])){
    norm <- func.2norm(mean,dqf[i,])
    dqf.2norm <- c(dqf.2norm, norm)
    if(norm < bound.dqf.2norm) labels[i] <- 1
    else labels[i] <- 2
  }
  
  plot(dqf.2norm,col=labels,main=main)
  abline(h=bound.dqf.2norm)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
func.supnorm <- function(dqf1,dqf2){
  sqrt(max((dqf1-dqf2)^2))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
dqf.supnorm <- function(dqf){
  dqf.supnorm <- c()
  
  mean <- dqf.mean(dqf)
  sd <- dqf.sd(dqf)
  
  for(i in 1:length(dqf[,1])){
    norm <- func.supnorm(mean,dqf[i,])
    dqf.2norm <- c(dqf.supnorm, norm)
  }
  
  return(dqf.supnorm)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
plot.dqf.supnorm <- function(dqf, n.sd=2,main=""){
  dqf.supnorm <- c()
  labels <- c()
  
  mean <- dqf.mean(dqf)
  sd <- dqf.sd(dqf)
  
  bound.dqf.supnorm <- func.supnorm(dqf.upperbound(mean,sd,n.sd),mean)
  
  for(i in 1:length(dqf[,1])){
    norm <- func.supnorm(mean,dqf[i,])
    dqf.supnorm <- c(dqf.supnorm, norm)
    if(norm < bound.dqf.supnorm) labels[i] <- 1
    else labels[i] <- 2
  }
  
  plot(dqf.supnorm,col=labels,main=main)
  abline(h=bound.dqf.supnorm)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
prop.outside.bounds <- function(dqf, n.sd=2){
  n <- length(dqf[,1])
  count <- rep(0,n)
  
  mean <- dqf.mean(dqf)
  sd <- dqf.sd(dqf)
  
  lower.bound <- dqf.lowerbound(mean,sd,n.sd)
  upper.bound <- dqf.upperbound(mean,sd,n.sd)
  
  for(i in 1:n){
    for(j in 1:length(dqf[i,])){
      if(dqf[i,j] < lower.bound[j] | dqf[i,j] > upper.bound[j]){
        count[i] <- count[i]+1
      }
    }
  }
  
  return(count/n)
  
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
plot.prop.outside.bounds <- function(dqf,threshold=.2,n.sd=2,main=""){
  n <- length(dqf[,1])
  count <- rep(0,n)
  labels <- rep(1,n)
  
  mean <- dqf.mean(dqf)
  sd <- dqf.sd(dqf)
  
  lower.bound <- dqf.lowerbound(mean,sd,n.sd)
  upper.bound <- dqf.upperbound(mean,sd,n.sd)
  
  for(i in 1:n){
    for(j in 1:length(dqf[i,])){
      if(dqf[i,j] < lower.bound[j] | dqf[i,j] > upper.bound[j])count[i] <- count[i]+1
    }
    count[i] <- count[i]/n
    if(count[i] > threshold) labels[i] <- 2  
  }
  
  plot(count,col=labels,main=main)
  abline(h=threshold)
  
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
dqf.zscore <- function(dqf){
  return(abs(scale(dqf)))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
dqf.mean.zscore <- function(dqf){
  return(rowSums(dqf.zscore(dqf))/length(dqf[,1]))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
plot.zscore.dqf <- function(dqf,labels=NULL,main=""){
  scaled <- abs(scale(dqf))
  
  if(is.null(labels)) labels <- rep(1,length(dqf[,1]))
  
  x <- seq(.01,1,.01)
  plot(x,scaled[1,],t='l',col=labels[1],ylim=c(min(scaled,na.rm=TRUE),max(scaled,na.rm=TRUE)),main=main,ylab="z-score")
  for(i in 2:length(scaled[,1])){
    lines(x,scaled[i,],col=labels[i])
  }
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
plot.mean.zscores <- function(dqf,labels=NULL,main=""){
  mean.zscores <- rowSums(dqf.zscore(dqf),na.rm=TRUE)/length(dqf[1,])
  if(is.null(labels)) labels <- rep(1,length(dqf[,1]))
  
  plot(mean.zscores,col=labels,ylim=c(min(mean.zscores,na.rm=TRUE),max(mean.zscores,na.rm=TRUE)),main=main)
}

