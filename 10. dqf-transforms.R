## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- eval=FALSE, include=FALSE,results='hide'------------------------------------------------------------------------------------------------
## source("00. datasets.R")
## source("01. dqf-outlier.R")
## source("02. dqf-subset.R")
## source("03. extract-dqf-subset.R")


## ---- eval=FALSE, include=FALSE,results='hide'------------------------------------------------------------------------------------------------
## knitr::purl("10. dqf-transforms.Rmd")


## ---------------------------------------------------------------------------------------------------------------------------------------------
scale.dqf.max <- function(dqf){
  n.functions <- length(dqf[,1])
  ret <- dqf
  
  for(i in 1:n.functions){
    func.max <- max(dqf[i,])
    for(j in 1:length(dqf[i,])){
      ret[i,j] <- dqf[i,j]/func.max
    }
  }
  
  return(ret)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
scale.dqf.sum <- function(dqf){
  n.functions <- length(dqf[,1])
  ret <- dqf
  
  for(i in 1:n.functions){
    func.sum <- sum(dqf[i,])
    for(j in 1:length(dqf[i,])){
      ret[i,j] <- dqf[i,j]/func.sum
    }
  }
  
  return(ret)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
scale.dqf.globalmax <- function(dqf){
  n.functions <- length(dqf[,1])
  globalmax <- max(dqf)
  
  ret <- dqf
  
  for(i in 1:n.functions){
    for(j in 1:length(dqf[i,])){
      ret[i,j] <- dqf[i,j]/globalmax
    }
  }
  
  return(ret)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
scale.dqf.globalsum <- function(dqf){
  n.functions <- length(dqf[,1])
  globalsum <- sum(dqf)
  
  ret <- dqf
  
  for(i in 1:n.functions){
    
    for(j in 1:length(dqf[i,])){
      ret[i,j] <- dqf[i,j]/globalsum
    }
  }
  
  return(ret)
}

