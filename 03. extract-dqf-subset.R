## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- eval=FALSE, include=FALSE,results='hide'------------------------------------------------------------------------------------------------
## source("00. datasets.R")
## source("01. dqf-outlier.R")
## source("02. dqf-subset.R")


## ---- eval=FALSE, include=FALSE,results='hide'------------------------------------------------------------------------------------------------
## knitr::purl("03. extract-dqf-subset.Rmd")


## ---------------------------------------------------------------------------------------------------------------------------------------------
extract.dqfs <- function(dqf.s,subset){

  pairs <- dqf.s$ret.pairs
  pair.indices <- which(pairs[,1] %in% subset & pairs[,2] %in% subset)
  
  ret.pairs <- pairs[pair.indices,]
  ret.goods <- dqf.s$ret.goods
  
  depthity <- rep(0,nrow(ret.goods[[1]]))
  qfs <- matrix(0, nrow=length(pair.indices), ncol=100)
  
  for (i in 1:length(pair.indices)) {
    
    goods <- ret.goods[[pair.indices[i]]]
    
    for (c in 1:nrow(ret.goods[[1]])) {
      # 100
      good <- goods[c,][subset]
      depthity[c] <- min(c(sum(good == -1), sum(good == 1)))
    }
    
    qfs[i, ] <- quantile(depthity, seq(0, 1, length = 100), na.rm = TRUE)
    
  }
  
  
  n.obs <- length(subset)
  
  dqf <- matrix(0, n.obs, 100)
  
  if(length(ret.pairs)==2){ # avoid incorrect dimensions
    for (i in 1:n.obs) {
      obs <- subset[i]
      dqf[i,] <- apply(qfs, 2, mean, na.rm = TRUE)
    }
  } 
  else{
    for (i in 1:n.obs) {
      obs <- subset[i]
      dqf[i,] <- apply(qfs[which(ret.pairs[, 1] == obs | ret.pairs[, 2] == obs), ], 2, mean, na.rm = TRUE)
    }
  }
  
 
  return(dqf) 
}

