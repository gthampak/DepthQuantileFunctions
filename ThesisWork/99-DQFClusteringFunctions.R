library(glue)

sd.w <- function(x, k) {
  # computes the windsorized standard deviation, called by dqf.outlier
  # Inputs:
  ## k: (integer) number of observations at each extreme to alter
  ## x: (vector) numeric data values
  k <- floor(k)
  if (k == 0) {  #corresponds to non-robust adaptive DQF
    return(sd(x))
  } else {
    x <- sort(x)
    n <- length(x)
    x[1:k] <- x[k]
    x[(n-k+1):n] <- x[n-k+1]
    return(sd(x))
  }
}

subsamp.dqf <- function(n.obs, subsample) {
  # called by dqf.outlier, computes random subset of pairs
  pairs <- c()
  subsample <- floor(subsample/2)*2
  for (i in 1:n.obs) {
    for (j in (i+1):(i+subsample/2)) {
      pairs <- rbind(pairs, c(i,j*(j <= n.obs) + (j-n.obs)*(j > n.obs)))
    }
  }
  return(pairs)
}

plot.dqf <- function(dqf,labels=NULL,xlab='Quantiles',ylab='Depth',main=''){
  x <- seq(.01,1,.01)
  
  if(nrow(dqf)==1){ # only one function
    plot(x,dqf,t='l',xlab=xlab,ylab=ylab,main=main)
  }
  else{
    n.functions <- nrow(dqf)
    if(is.null(labels)) labels <- rep(1,n.functions) 
    
    plot(x,dqf[1,],t='l',ylim=c(0,max(dqf)),col=labels[1],xlab=xlab,ylab=ylab,main=main)
    for(i in 2:n.functions){
      lines(x,dqf[i,],col=labels[i])
    }
  }
}

show <- function(length,s){
  labels <- rep(1,length)
  labels[s] <- 2
  return(labels)
}

dqf.subset <- function(data = NULL, gram.mat = NULL, g.scale=2, angle=c(45), kernel="linear", p1=1, p2=0, n.splits=100, subsample=50, z.scale=TRUE, k.w=3, adaptive=TRUE, G="norm") {
  # kernelized version of depthity
  # 
  # inputs: 
  #   data (matrix or data frame) - a data matrix of explanatory variables
  #   kernel - of form "linear", "rbf" or "poly", or a user defined function 
  #   g.scale (scalar) - scales the base distribution G
  #   angle (numeric vector of length 3)- angles of cone from midline, must live between 0 and 90
  #   p1 - first parameter for kernel
  #   p2 - second parameter for kernel
  #   n.splits (integer) - the number of split points at which the DQF is computed
  #   subsample (integer)- the number of random pairs for each observation
  #   z.scale (logical) - should the data be z-scaled first
  #   k.w (integer) - the number of points altered in the windsorized standard deviation
  #   adaptive - if TRUE, uses windsorized standard deviation to scale base distribution
  #   G - base distribution:   "norm" or "unif"
  ##
  # Output:
  #   angle  - vector of angles used, same as inputted
  #   dqf1, dqf2, dqf3 - matrices of depth quantile functions, rows are observations
  if (G=="norm") {
    param1 <- 0; param2 <- 1 }
  if (G=="unif") {
    param1 <- -1; param2 <- 1 }
  if (is.null(data) & is.null(gram.mat)) 
    stop("Either a data set or Gram matrix must be provided")
  if (min(angle) <= 0 | max(angle) >= 90) 
    stop("Angles must be between 0 and 90")
  if (is.null(data))
    n.obs <- nrow(gram.mat)
  if (is.null(gram.mat))
    n.obs <- nrow(data)
  
  
  #scram <- sample(n.obs)
  #unscram <- c()
  #for(i in 1:length(scram)){
  #  unscram <- c(unscram,which(scram==i))
  #}
  
  pairs <- subsamp.dqf(n.obs, subsample)
  if (is.null(gram.mat)) {  
    if (z.scale==TRUE)
      data <- apply(data, 2, scale) #z-scale data
    if (is.function(kernel)==TRUE)
      kern <- kernel
    if (kernel == "linear") {
      kern <- function(x,y)
        return(sum(x*y))
    }
    if (kernel == "rbf") {
      kern <- function(x,y)
        return(exp(-sum((x-y)^2)/p1))
    }             
    if (kernel == "poly") {
      kern <- function(x,y)
        return((sum(x*y)+p2)^p1)
    }
    #data <- data[scram,]
    gram <- matrix(0,n.obs, n.obs)
    for (i in 1:n.obs) {
      for (j in i:n.obs) {
        gram[i,j] <- kern(data[i,], data[j,])
        gram[j,i] <- gram[i,j]
      }
    }
  } else {
    if (diff(dim(gram.mat))!=0)
      stop("Gram matrix must be square")
    if (isSymmetric(gram.mat) == FALSE)
      stop("Gram matrix must be symmetric")
    gram <- gram.mat
    #gram <- gram[scram,]
    #gram <- gram[,scram]
  }
  splits <- get(paste("q", G, sep=""))((1:n.splits)/(n.splits+1),param1, param2) * g.scale
  depthity1 <- rep(0,length(splits))
  norm.k2 <- error.k <- k.to.mid <- rep(0, n.obs)
  dep1 <- matrix(0, nrow=nrow(pairs), ncol=n.splits)
  qfs1 <- matrix(0, nrow=nrow(pairs), ncol=100)
  
  ret.pairs <- matrix(0,nrow=nrow(pairs),ncol=2)
  k.to.mids <- matrix(ncol = n.obs, nrow = nrow(pairs))
  
  ret.gs <- matrix(ncol = n.obs, nrow = n.splits)  # skeleton for goods data.frame
  ret.goods <- list()
  
  # pairs.goods <- list(pairs=list(c(0,0)),goods=list(data.frame(matrix(ncol = n.obs, nrow = 100)))) # create list of pairs and their corresponding 100 (num.splits) good vectors
  
  for (i.subs in 1:nrow(pairs)) {
    i <- pairs[i.subs,1];  j <- pairs[i.subs,2]
    
    # pairs.goods$pairs <- append(pairs.goods$pairs,list(c(scram[i],scram[j]))) # record pairs of points (original indices)
    
    #ret.pairs[i.subs,] <- c(scram[i],scram[j])
    ret.pairs[i.subs,] <- c(i,j)
    
    for (k in 1:n.obs) {
      norm.k2[k] <- gram[k,k] + 1/4*(gram[i,i]+gram[j,j]) + 1/2*gram[i,j]-gram[k,i]-gram[k,j]
      k.to.mid[k] <- (gram[k,i]-gram[k,j]+1/2*(gram[j,j]-gram[i,i]))/sqrt(gram[i,i]+gram[j,j]-2*gram[i,j]) # return this
      error.k[k] <- sqrt(abs(norm.k2[k] - k.to.mid[k]^2))
    }
    
    k.to.mids[i.subs,] <- k.to.mid
    
    for (c in 1:length(splits)) {
      good <- rep(1, n.obs)
      s <- splits[c] * (sd.w(k.to.mid, k.w)*adaptive + (adaptive==FALSE)) 
      good[k.to.mid/s > 1] <- 0  #points on other side of cone tip removed
      d.to.tip <- abs(k.to.mid - s)
      good1 <- good * (abs(atan(error.k / d.to.tip)) < (angle[1]/360*2*pi))  #points outside of cone removed
      good1 <- good1 * (1 - 2*(sign(k.to.mid)==sign(s)))  #which side of midpoint are they on
      ret.gs[c,] <- good1
      depthity1[c] <- min(c(sum(good1==-1), sum(good1==1)))
    }
    ret.goods[[i.subs]] <- ret.gs
    qfs1[i.subs,] <- quantile(depthity1, seq(0,1,length=100), na.rm=TRUE)
  }
  dqf1 <- matrix(0,n.obs, 100)
  for (i in 1:n.obs) {
    dqf1[i,] <- apply(qfs1[which(pairs[,1]==i | pairs[,2]==i),],2,mean, na.rm=TRUE)
  }
  
  # reorder
  # dqf1 <- dqf1[unscram,]  #map back to original indicies
  
  return(list(dqf1=dqf1,ret.pairs=ret.pairs,ret.goods=ret.goods,k.to.mids=k.to.mids,splits=splits))
}

extract.dqf <- function(dqf,subset){
  n.subset <- length(subset)
  depthity1 <- rep(0,100)
  
  indices <- which(dqf$ret.pairs[,1] %in% subset & dqf$ret.pairs[,2] %in% subset)
  pairs <- dqf$ret.pairs[indices,]
  
  qfs1 <- matrix(0, nrow=length(indices), ncol=100)
  
  qfs.count <- 1
  for(j in indices){
    gs <- dqf$ret.goods[[j]][,subset]
    for(i in 1:nrow(gs)){
      g <- gs[i,]
      depthity1[i] <- min(c(sum(g==-1), sum(g==1)))
    }
    qfs1[qfs.count,] <- quantile(depthity1, seq(0,1,length=100), na.rm=TRUE)
    qfs.count <- qfs.count+1
  }
  
  dqf1 <- matrix(0,n.subset, 100)
  for (i in 1:n.subset) {
    dqf1[i,] <- apply(qfs1[which(pairs[,1]==i | pairs[,2]==i),],2,mean,na.rm=TRUE)
  }
  
  return(dqf1)
}

extract.depths <- function(dqf,subset){
  n.subset <- length(subset)
  depthity1 <- rep(0,100)
  
  indices <- which(dqf$ret.pairs[,1] %in% subset & dqf$ret.pairs[,2] %in% subset)
  subset.pairs <- dqf$ret.pairs[indices,]
  
  depths <- matrix(0, nrow=length(indices), ncol=100)
  
  depths.count <- 1
  for(j in indices){
    gs <- dqf$ret.goods[[j]][,subset]
    for(i in 1:nrow(gs)){
      g <- gs[i,]
      depthity1[i] <- min(c(sum(g==-1), sum(g==1)))
    }
    depths[depths.count,] <- depthity1
    depths.count <- depths.count+1
  }
  
  return(list(depths=depths,subset.pairs=subset.pairs))
}

calculate.qfs <- function(depths,k2mid,splits){
  qfs <- matrix(0,nrow=nrow(depths),ncol=length(splits))
  
  for(i in 1:nrow(depths)){
    qf <- rep(0,100)
    d <- depths[i,]
    k.sd <- sd.w(k2mid[i,],3)
    
    q.prev <- 1
    min.s <- 1; max.s <- -1
    for(j in 0:max(d)){
      s <- splits[which(d==j)]
      min.s <- min(c(s,min.s)); max.s <- max(c(s,max.s))
      q <- ceiling((pnorm(max.s,sd=k.sd)-pnorm(min.s,sd=k.sd))*100)
      
      if(q>q.prev & q <= 100){
        qf[q.prev:q] <- j
        q.prev <- q+1
      }
      
    }
    qfs[i,] <- qf
    
  }
  return(qfs)
}

calculate.dqfs <- function(qfs,pairs,subset,n.obs){
  
  if(length(pairs)==2) return(qfs)
  
  dqfs <- matrix(0,n.obs, 100)
  for (i in 1:n.obs) {
    dqfs[i,] <- apply(qfs[which(pairs[,1]==i | pairs[,2]==i),],2,mean, na.rm=TRUE)
  }
  return(dqfs[subset,])
}

row.col <- function(dataframe,index){
  n.row <- nrow(dataframe); n.col <- ncol(dataframe)
  
  row <- index%%n.row
  if(row==0){ row <- n.row; col <- index/n.row}
  else{ col <- ceiling(index/n.row) }
  
  return(list(row=row,col=col))
  
}

initial.cluster <- function(data,n.clusters){
  n.obs <- nrow(data)
  distance_mat <- dist(data, method = 'euclidean')
  Hierar_cl <- hclust(distance_mat, method = 'single')
  clusters <- cutree(Hierar_cl, k = n.clusters)
  
  # calculate minimum distance between clusters
  
  dist.mat <- as.matrix(distance_mat)
  clus.indices <- list()
  for(i in 1:n.clusters){
    clus.indices[[i]] <- which(clusters==i)
  }
  
  # calculate min distance from cluster i to j by:
  # calculate min of distances between point i_1 to point js
  # calculate min of those
  
  inter.dists <- matrix(0,nrow=n.clusters,ncol=n.clusters) # distance between closest points of clusters
  closest.pts <- matrix(0,nrow=n.clusters,ncol=n.clusters) # M(i,j) is index of closest point in cluster j to cluster i
  for(i in 1:n.clusters){
    for(j in i:n.clusters){
      if(i==j){inter.dists[i,j] <- Inf}
      else{
        inter.dists[i,j] <- inter.dists[j,i] <- min(dist.mat[clus.indices[[i]],clus.indices[[j]]])
        
        min.index <- which(dist.mat==inter.dists[i,j])[1]
        rc <- row.col(dist.mat,min.index)
        row <- rc$row; col <- rc$col
        
        if(row %in% clus.indices[[i]]){ closest.pts[j,i] <- row; closest.pts[i,j] <- col }
        else{ closest.pts[i,j] <- row; closest.pts[j,i] <- col }
        
      }
    }
  }
  
  return(list(clusters=clusters,inter.dists=inter.dists,closest.pts=closest.pts))
}

max.dist <- function(inter.dists){
  id <- inter.dists
  for(i in 1:length(id)){ if(id[i] == Inf) id[i] <- 0 }
  return(max(id))
}

combine.clusters <- function(combined.clusters,c1,c2){
  combined.clusters[[c1]] <- sort(c(combined.clusters[[c1]],combined.clusters[[c2]]))
  for(c in combined.clusters[[c1]]) combined.clusters[[c]] <- combined.clusters[[c1]]
  return(combined.clusters)
}

all.indices <- function(clusters, cluster.group){
  # to ensure function works with both single and vector cluster.group inputs
  v <- c()
  v <- c(v, cluster.group)
  
  indices <- c()
  for(c in v) indices <- c(indices,which(clusters == c))
  
  return(indices)
}

closest.pt.clusters <- function(inter.dists, cg1, cg2){
  n.cg1 <- length(c(cg1))
  n.cg2 <- length(c(cg2))
  m <- min(inter.dists[cg1,cg2])
  min.index <- which(inter.dists==m)[1]
  
  rc <- row.col(inter.dists,min.index)
  row <- rc$row; col <- rc$col
  
  return(list(c1=row,c2=col))
  
}

compile.clusters <- function(clusters, combined.clusters){
  final.clusters <- clusters
  for(i in 1:length(combined.clusters)){
    c <- combined.clusters[[i]][1]
    indices <- c()
    for(j in combined.clusters[[i]]){
      indices <- c(indices,which(final.clusters==j))
    }
    final.clusters[indices] <- c
  }
  
  return(final.clusters)
}

combine.prompt <- function(row,col,inter.dists,combined.clusters){
  yes <- '"Yes"'; no <- '"No"'; pp <- '"Postpone"'
  print(glue('Clusters {cluster.string(combined.clusters[[row]])} and {cluster.string(combined.clusters[[col]])} are the closest by Euclidean distance (dist = {inter.dists[row,col]}).'))
  print(glue('Type "{as.character(combined.clusters[[row]][1])}" for dqfs of points in cluster "{cluster.string(combined.clusters[[row]])}" and closest point from cluster "{cluster.string(combined.clusters[[col]])}".'))
  print(glue('Type "{as.character(combined.clusters[[col]][1])}" for dqfs of points in cluster "{cluster.string(combined.clusters[[col]])}" and closest point from cluster "{cluster.string(combined.clusters[[row]])}".'))
  print(glue('Enter {yes} to combine clusters {cluster.string(combined.clusters[[row]])} and {cluster.string(combined.clusters[[col]])}?'))
  print(glue('Enter {no} to move on to next clusters.'))
  print(glue('Enter {pp} to consider combining these clusters later.'))
}

cluster.string <- function(cluster){
  
  ret <- '{'
  
  for(i in 1:length(cluster)){
    if(i==length(cluster)){ret <- paste(ret,as.character(cluster[i]),sep='')}
    else{ret <- paste(ret, as.character(cluster[i]),' ',sep='')}
  }
  
  ret <- paste(ret,'}',sep='')
  
  return(ret)
}

dqf.clustering <- function(data = NULL,dqf.s=NULL,n.clusters=20, gram.mat = NULL, g.scale=2, angle=c(45), kernel="linear", p1=1, p2=0, n.splits=100, subsample=50, z.scale=TRUE, k.w=3, adaptive=TRUE, G="norm"){
  
  # data <- scale(data)
  
  # set n.obs - to number of observations
  if (is.null(data)) n.obs <- nrow(gram.mat)
  if (is.null(gram.mat)) n.obs <- nrow(data)
  
  subsample <- n.obs # calculate all sub-samples
  
  # initial clustering
  ic <- initial.cluster(data,n.clusters)
  clusters <- ic$clusters
  inter.dists <- ic$inter.dists
  closest.pts <- ic$closest.pts
  
  
  combined.clusters <- list()
  for(i in 1:n.clusters) combined.clusters[[i]] <- c(i)
  
  # calculate dqf.subset
  if(is.null(dqf.s)){
    print("Calculating `dqf.subset`...")
    dqf.s <- dqf.subset(data = data, gram.mat = gram.mat, g.scale=g.scale, angle=c(45), kernel=kernel, p1=p1, p2=p2, n.splits=n.splits, subsample=subsample, z.scale=z.scale, k.w=k.w, adaptive=adaptive, G=G)
    print("`dqf.subset` calculations complete.")
  }
  
  print("We will now begin combining clusters.")
  print("You may type exit at any time and the function will return `dqf.subset` and current clusters.")
  print("Press Enter to continue")
  s = readline()
  while(s != 'exit'){
    move.on <- FALSE
    
    if(min(inter.dists)==Inf){
      final.clusters <- compile.clusters(clusters, combined.clusters)
      return(list(dqf.s=dqf.s,final.clusters=final.clusters))
    }
    mi <- which.min(inter.dists) # mi is short for min.index
    rc <- row.col(inter.dists,mi)
    row <- rc$row; col <- rc$col
    
    # Prepare dqfs
    cpc <- closest.pt.clusters(inter.dists,combined.clusters[[row]],combined.clusters[[col]])
    pt2.index <- closest.pts[cpc$c1,cpc$c2] # index of closest point in cluster group 2 to cg1
    pt1.index <- closest.pts[cpc$c2,cpc$c1] # index of closest point in cluster group 1 to cg2
    
    subset1 <- c(all.indices(clusters,combined.clusters[[row]]),pt2.index)
    ed1 <- extract.depths(dqf.s,subset1)
    depths1 <- ed1$depths; subset1.pairs <- ed1$subset.pairs
    qfs1 <- calculate.qfs(depths=depths1,dqf.s$k.to.mids,dqf.s$splits)
    dqfs1 <- calculate.dqfs(qfs1,subset1.pairs,subset1,n.obs)
    labels1 <- rep(1,length(subset1)); labels1[which(subset1==pt2.index)] <- 2
    
    subset2 <- c(all.indices(clusters,combined.clusters[[col]]),pt1.index)
    ed2 <- extract.depths(dqf.s,subset2)
    depths2 <- ed2$depths; subset2.pairs <- ed2$subset.pairs
    qfs2 <- calculate.qfs(depths=depths2,dqf.s$k.to.mids,dqf.s$splits)
    dqfs2 <- calculate.dqfs(qfs2,subset2.pairs,subset2,n.obs)
    labels2 <- rep(1,length(subset2)); labels2[which(subset2==pt1.index)] <- 2
    
    combine.prompt(row,col,inter.dists,combined.clusters)
    
    s <- readline()
    while(!move.on){
      if(s == "Yes"){
        print(glue("Clusters {cluster.string(combined.clusters[[row]])} and {cluster.string(combined.clusters[[col]])} were combined."))
        combined.clusters <- combine.clusters(combined.clusters,row,col)
        inter.dists[combined.clusters[[row]],combined.clusters[[col]]] <- Inf # don't need to consider combining them anymore
        inter.dists[combined.clusters[[col]],combined.clusters[[row]]] <- Inf
        move.on <- TRUE
      }else if(s == "No"){
        # print("Are you sure? If yes, combinning clusters {cluster.string(combined.clusters[[row]])} and {cluster.string(combined.clusters[[col]])} will not be revisited.")
        # s <- readline()
        inter.dists[combined.clusters[[row]],combined.clusters[[col]]] <- Inf
        inter.dists[combined.clusters[[col]],combined.clusters[[row]]] <- Inf
        move.on <- TRUE
      }else if(s == "Later"){
        print("Later")
        s <- readline()
      }else if(s == row){
        plot.dqf(dqfs1,labels1)
        combine.prompt(row,col,inter.dists,combined.clusters)
        s <- readline()
      }else if(s == col){
        plot.dqf(dqfs2,labels2)
        combine.prompt(row,col,inter.dists,combined.clusters)
        s <- readline()
      } else if(s == 'exit'){
        move.on = TRUE
      } else{
        print("Invalid input")
        combine.prompt(row,col,inter.dists,combined.clusters)
        s <- readline()
      }
    }
    
  }
  
  final.clusters <- compile.clusters(clusters, combined.clusters)
  
  return(list(dqf.s=dqf.s,final.clusters=final.clusters))
  
}


