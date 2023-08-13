## ----setup, include=FALSE--------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- eval=FALSE, include=FALSE,results='hide'-----------------------------------
## source("00. datasets.R")
## source("01. dqf-outlier.R")
## source("02. dqf-subset.R")
## source("03. extract-dqf-subset.R")


## ---- eval=FALSE, include=FALSE,results='hide'-----------------------------------
## knitr::purl("04. dqf-clustering.Rmd")


## ---- eval=FALSE, include=FALSE,results='hide'-----------------------------------
## install.packages("glue")
## require(glue)


## --------------------------------------------------------------------------------
row.col <- function(dataframe,index){
  n.row <- nrow(dataframe); n.col <- ncol(dataframe)
  
  row <- index%%n.row
  if(row==0){ row <- n.row; col <- index/n.row}
  else{ col <- ceiling(index/n.row) }
  
  return(list(row=row,col=col))
  
}


## --------------------------------------------------------------------------------
initial.cluster <- function(data,n.clusters,cluster.method,p1){
  n.obs <- nrow(data)
  
  distance_mat <- dist(data, method = "euclidean")
  clusters <- NULL
  
  tryCatch(expr = {
    if (cluster.method == "HC") {
      
      if(is.null(p1)){ p1 <- "single"}
      
      Hierar_cl <- hclust(distance_mat, method = p1) # p1
      clusters <- cutree(Hierar_cl, k = n.clusters)
    }
    else if (cluster.method == "KM") {
      set.seed(47)
      
      if(is.null(p1)){ p1 <- 40 }
      
      kmeans.re <- kmeans(data, centers = n.clusters, nstart = p1)
      clusters <- kmeans.re$cluster
      print(p1)
    }
    else {
      return("error.clustermethod")
    }
  },
  error = function(e) { return("error.parameter")},
  warning = function(w) { return("error.parameter")},
  finally = {}
  )
  
  if(is.null(clusters)){ return("error.parameter") }
  
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


## --------------------------------------------------------------------------------
calculate.inter.dists <- function(data, clusters) {
  
  n.clusters <- length(unique(clusters))
  
  distance_mat <- dist(data, method = "euclidean")
  
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


## --------------------------------------------------------------------------------
max.dist <- function(inter.dists){
  id <- inter.dists
  for(i in 1:length(id)){ if(id[i] == Inf) id[i] <- 0 }
  return(max(id))
}


## --------------------------------------------------------------------------------
combine.clusters <- function(combined.clusters,c1,c2){
  combined.clusters[[c1]] <- sort(c(combined.clusters[[c1]],combined.clusters[[c2]]))
  for(c in combined.clusters[[c1]]) combined.clusters[[c]] <- combined.clusters[[c1]]
  return(combined.clusters)
}


## --------------------------------------------------------------------------------
all.indices <- function(clusters, cluster.group){
  # to ensure function works with both single and vector cluster.group inputs
  v <- c()
  v <- c(v, cluster.group)
  
  indices <- c()
  for(c in v) indices <- c(indices,which(clusters == c))
  
  return(indices)
}


## --------------------------------------------------------------------------------
closest.pt.clusters <- function(inter.dists, cg1, cg2){
  n.cg1 <- length(c(cg1))
  n.cg2 <- length(c(cg2))
  m <- min(inter.dists[cg1,cg2])
  min.index <- which(inter.dists==m)[1]
  
  rc <- row.col(inter.dists,min.index)
  row <- rc$row; col <- rc$col
  
  return(list(c1=row,c2=col))
  
}


## --------------------------------------------------------------------------------
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


## --------------------------------------------------------------------------------
organize.clusters <- function(clusters){
  min.pointer <- 1
  max.pointer <- max(clusters)
  while(max(clusters) != length(unique(clusters))){
    while(min.pointer %in% unique(clusters)) min.pointer <- min.pointer+1
    clusters[which(clusters==max(clusters))] <- min.pointer
  }
  
  return(clusters)
}


## --------------------------------------------------------------------------------
combine.prompt <- function(row,col,inter.dists,combined.clusters){
  yes <- '"Yes"'; no <- '"No"'; pp <- '"pp"'; d <- '"d"'; b <- '"b"'
  print(glue('Clusters {cluster.string(combined.clusters[[row]])} and {cluster.string(combined.clusters[[col]])} are the closest by Euclidean distance (dist = {inter.dists[row,col]}).'))
  print(glue('Type "{as.character(combined.clusters[[row]][1])}" for dqfs of points in cluster "{cluster.string(combined.clusters[[row]])}" and closest point from cluster "{cluster.string(combined.clusters[[col]])}".'))
  print(glue('Type "{as.character(combined.clusters[[col]][1])}" for dqfs of points in cluster "{cluster.string(combined.clusters[[col]])}" and closest point from cluster "{cluster.string(combined.clusters[[row]])}".'))
  print(glue('Type {b} (for both) to view both DQF plots.'))
  print(glue('Type {d} (for diagnostics) to view full dqf diagnostics for displayed data.'))
  print(glue('Enter {yes} to combine clusters {cluster.string(combined.clusters[[row]])} and {cluster.string(combined.clusters[[col]])}.'))
  print(glue('Enter {no} to move on to next clusters.'))
  print(glue('Enter {pp} (for postpone) to consider combining these clusters later.'))
}


## --------------------------------------------------------------------------------
cluster.string <- function(cluster){
  
  ret <- '{'
  
  for(i in 1:length(cluster)){
    if(i==length(cluster)){ret <- paste(ret,as.character(cluster[i]),sep='')}
    else{ret <- paste(ret, as.character(cluster[i]),' ',sep='')}
  }
  
  ret <- paste(ret,'}',sep='')
  
  return(ret)
}


## --------------------------------------------------------------------------------
dqf.clustering <- function(data = NULL,dqf.s=NULL,initial.clusters=NULL,n.clusters=NULL,cluster.method=NULL,cluster.param=NULL, gram.mat = NULL, g.scale=2, angle=c(45), kernel="linear", p1=1, p2=0, n.splits=100, subsample=50, z.scale=TRUE, k.w=3, adaptive=TRUE, G="norm"){
  
  # data <- scale(data)
  
  # set n.obs - to number of observations
  if (is.null(data)) n.obs <- nrow(gram.mat)
  if (is.null(gram.mat)) n.obs <- nrow(data)
  
  # calculate all sub-samples
  subsample <- n.obs

  # clustering
  # check if user provided initial clusters
  if(is.null(initial.clusters)){
    # perform initial clustering
    print("Starting clusters not inputted. Performing initial clustering.")
    print("---")
    
    if(is.null(cluster.method)){
      print("Starting clusters not inputted. Performing single linkage hierarchical clustering.")
      print("---")
      cluster.method <- "HC"
      cluster.param <- "single"
    }
    
    if(is.null(n.clusters)){
      print("Starting number of clusters not inputted. Defaulting to 10 initial clusters")
      print("---")
      n.clusters <- 10
    }
    
    print("Performing initial clustering.")
    
    if(is.null(cluster.method)){
      cluster.method <- "HC"
      cluster.param <- "single"
    }
    
    ic <- initial.cluster(data,n.clusters,cluster.method=cluster.method,p1=cluster.param)
    
    if(is.character(ic)){
      if(ic == "error.parameter"){
        print("Invalid input for cluster.param")
        return()
      }
      else if(ic == "error.clustermethod"){
        print("Invalid input for cluster.method")
        return()
      }
    }
    else{
      clusters <- ic$clusters
      inter.dists <- ic$inter.dists
      closest.pts <- ic$closest.pts
    }
    
  }
  else{
    print("Starting clusters inputted. Calculating inter-cluster distances.")
    print("---")
    ic <- calculate.inter.dists(data, clusters)
    clusters <- ic$clusters
    inter.dists <- ic$inter.dists
    closest.pts <- ic$closest.pts
  }
  

  # set up clusters data structure for combining clusters
  # contains lists of combined clusters, indexed by initial clusters
  combined.clusters <- list()
  for(i in 1:n.clusters) combined.clusters[[i]] <- c(i)
  
  # calculate dqf.subset
  if(is.null(dqf.s)){
    print("Calculating `dqf.subset`...")
    dqf.s <- dqf.subset(data = data, gram.mat = gram.mat, g.scale=g.scale, angle=angle, kernel=kernel, p1=p1, p2=p2, n.splits=n.splits, subsample=subsample, z.scale=z.scale, k.w=k.w, adaptive=adaptive, G=G)
    print("`dqf.subset` calculations complete.")
  }
  else{
    print("dqf.s object provided.")
    print("---")
  }
  
  print("We will now begin combining clusters.")
  print("You may type exit at any time and the function will return `dqf.subset` and current clusters.")
  print("Press Enter to continue")
  s = readline()
  while(s != 'exit'){
    move.on <- FALSE
    
    if(min(inter.dists)==Inf){
      final.clusters <- compile.clusters(clusters, combined.clusters)
      
      print("Clustering Complete. dqf.s object and final.clusters returned.")
      
      return(list(dqf.s=dqf.s,final.clusters=final.clusters))
    }
    
    max.dist <- max.dist(inter.dists)
    mi <- which.min(inter.dists) # mi is short for min.index
    rc <- row.col(inter.dists,mi)
    row <- rc$row; col <- rc$col
    
    # Prepare dqfs
    cpc <- closest.pt.clusters(inter.dists,combined.clusters[[row]],combined.clusters[[col]])
    pt2.index <- closest.pts[cpc$c1,cpc$c2] # index of closest point in cluster group 2 to cg1
    pt1.index <- closest.pts[cpc$c2,cpc$c1] # index of closest point in cluster group 1 to cg2
    
    subset1 <- c(all.indices(clusters,combined.clusters[[row]]),pt2.index)
    dqfs1 <- extract.dqfs(dqf.s,subset1)
    labels1 <- rep(1,length(subset1)); labels1[which(subset1==pt2.index)] <- 2
    
    subset2 <- c(all.indices(clusters,combined.clusters[[col]]),pt1.index)
    dqfs2 <- extract.dqfs(dqf.s,subset2)
    labels2 <- rep(1,length(subset2)); labels2[which(subset2==pt1.index)] <- 2
    
    diag.num <- row
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
        inter.dists[combined.clusters[[row]],combined.clusters[[col]]] <- Inf
        inter.dists[combined.clusters[[col]],combined.clusters[[row]]] <- Inf
        move.on <- TRUE
      }else if(s == "pp"){
        inter.dists[combined.clusters[[row]],combined.clusters[[col]]] <- inter.dists[combined.clusters[[row]],combined.clusters[[col]]] + max.dist
        inter.dists[combined.clusters[[col]],combined.clusters[[row]]] <- inter.dists[combined.clusters[[col]],combined.clusters[[row]]] + max.dist
        move.on <- TRUE
      }else if(s %in% combined.clusters[[row]]){
        par(mfrow=c(1,1))
        plot.dqf(dqfs1,labels1)
        combine.prompt(row,col,inter.dists,combined.clusters)
        s <- readline()
        if(s %in% c(combined.clusters[[row]],combined.clusters[[col]])) diag.num <- s
      }else if(s %in% combined.clusters[[col]]){
        par(mfrow=c(1,1))
        plot.dqf(dqfs2,labels2)
        combine.prompt(row,col,inter.dists,combined.clusters)
        s <- readline()
        if(s %in% c(combined.clusters[[row]],combined.clusters[[col]])) diag.num <- s
      }else if(s == 'd'){
        if(diag.num %in% combined.clusters[[row]]){
          dqf.c.diagnostics(dqfs1,labels1)
          combine.prompt(row,col,inter.dists,combined.clusters)
          s <- readline()
          if(s %in% c(combined.clusters[[row]],combined.clusters[[col]])) diag.num <- s
        } else{
          dqf.c.diagnostics(dqfs2,labels2)
          combine.prompt(row,col,inter.dists,combined.clusters)
          s <- readline()
          if(s %in% c(combined.clusters[[row]],combined.clusters[[col]])) diag.num <- s
        }
        par(mfrow=c(1,1))
      }
      else if(s == 'exit'){
        move.on = TRUE
      }else{
        print("Invalid input")
        combine.prompt(row,col,inter.dists,combined.clusters)
        s <- readline()
      }
    }
    
  }
  
  final.clusters <- compile.clusters(clusters, combined.clusters)
  
  print("Process terminated. dqf.s object and final.clusters returned.")
  
  return(list(dqf.s=dqf.s,final.clusters=final.clusters))
  
}

