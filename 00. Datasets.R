## ----setup, include=FALSE, eval=FALSE---------------------------------------------------------------------------------------------------------
## knitr::opts_chunk$set(echo = TRUE)


## ---- eval=FALSE, include=FALSE,results='hide'------------------------------------------------------------------------------------------------
## knitr::purl("00. datasets.Rmd")


## ---------------------------------------------------------------------------------------------------------------------------------------------
# one outlier
data.one.outlier <- function(n.data=50, n.outliers=3, dist.outliers=8){
  # outputs a standardized 2 dimensional standard normal distribution with outliers
  
  # Inputs:
  #   n.data: (integer) number of data points in normals
  #   n.outliers: (integer) number of outliers
  #   dist.outliers: distance of the center of outlier distributions
  # Outputs:
  #   data: (dataframe) dataset with n.data points in the normal and n.outliers outliers 
  #   labels: (vector) 1s for indices of data points, 2s for indices of outliers 
  
  x <- rnorm(n.data)
  y <- rnorm(n.data)
  o1 <- rnorm(n.outliers, dist.outliers)
  o2 <- rnorm(n.outliers, dist.outliers)
  data <- cbind(x,y)
  outliers <- cbind(o1,o2)
  data <- data.frame(rbind(data,outliers))
  data <- scale(data)
  labels <- c(rep(1,n.data),rep(2,n.outliers))
  
  return(list(data=data,labels=labels))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
# bimodal with outliers
data.bimodal.outlier <- function() {
  x <- rnorm(50)
  y <- rnorm(50)
  data <- cbind(x,y)
  
  x <- rnorm(3,6)
  y <- rnorm(3,6)
  d <- cbind(x,y)
  data <- rbind(data,d)
  
  x <- rnorm(50, 12)
  y <- rnorm(50, 12)
  d <- cbind(x, y)
  data <- rbind(data,d)

  data <- scale(data)
  labels <- c(rep(1,50),rep(2,3),rep(3,50))
  
  return(list(data=data,labels=labels))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
data.n.modal <- function(n.obs=150,n.modal=3){
  n <- n.obs/n.modal
  
  data <- c()
  labels <- c()
  
  for(i in 1:n.modal){
    x <- rnorm(n,5*i)
    y <- rnorm(n,5*i)
    d <- cbind(x,y)
    data <- rbind(data,d)
    labels <- c(labels,rep(i,n))
  }
  
  return(list(data=data,labels=labels))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
# half moon
data.halfmoon <- function(){
  x <- seq(-20,20,.3)
  x1 <- x+20
  y <- (x-2)^2 + rnorm(length(x),0,30)
  y2 <- -(x+2)^2 + 600 + rnorm(length(x),0,30)
  data1 <- cbind(x,y)
  data2 <- cbind(x1,y2)
  m.data <- rbind(data1,data2)
  m.labels <- c(rep(1,length(data1[,1])),rep(2,length(data2[,1])))
  m.data <- scale(m.data)
  
  return(list(data=m.data,labels=m.labels))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
data.three.circles <- function(){
  x <- c(runif(20,-1,1),seq(-1,-.9,.03),seq(.9,1,.03))
  y <- sqrt(1-x^2) + rnorm(length(x),0,.1)
  y1.2 <- -sqrt(1-x^2) + rnorm(length(x),0,.1)
  c.data <- rbind(cbind(x,y),cbind(x,y1.2))
  x2 <- c(runif(40,-2,2),seq(-2,-1.8,.03),seq(1.8,2,.03))
  y2.1 <- sqrt(4-x2^2) + rnorm(length(x2),0,.1)
  y2.2 <- -sqrt(4-x2^2) + rnorm(length(x2),0,.1)
  c.data <- rbind(c.data,rbind(cbind(x2,y2.1),cbind(x2,y2.2)))
  x3 <- c(runif(60,-3,3),seq(-3,-2.7,.03),seq(2.7,3,.03))
  y3.1 <- sqrt(9-x3^2) + rnorm(length(x3),0,.1)
  y3.2 <- -sqrt(9-x3^2) + rnorm(length(x3),0,.1)
  c.data <- rbind(c.data,rbind(cbind(x3,y3.1),cbind(x3,y3.2)))
  c.labels <- c(rep(1,2*length(x)),rep(2,2*length(x2)),rep(3,2*length(x3)))
  c.data <- scale(c.data)
  
  return(list(data=c.data,labels=c.labels))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
data.double.helix <- function(n.nodes=150,noise=.05,n.spirals=1.5){
  z <- seq(0,1,1/n.nodes)[1:n.nodes]*2*pi*n.spirals
  z1 <- z + rnorm(n.nodes,0,noise)
  x <- sin(z1) + rnorm(n.nodes,0,noise)
  y <- cos(z1) + rnorm(n.nodes,0,noise)
  data <- cbind(x,y,z1)
  
  offset <- pi
  z2 <- z + rnorm(n.nodes,0,noise)
  x <- sin(z2 + offset) + rnorm(n.nodes,0,noise)
  y <- cos(z2 + offset) + rnorm(n.nodes,0,noise)
  d <- cbind(x,y,z2)
  data <- rbind(data,d)
  data <- scale(data)
  
  labels <- c(rep(1,n.nodes),rep(2,n.nodes))
  
  return(list(data=data,labels=labels))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------
generate.filament <- function(n=150, step_size=0.05, bool_func=TRUE) {
  # bool_func: generate function by restricting theta
  x <- c(); y <- c()
  x[1] <- runif(1)
  y[1] <- runif(1)
  theta <- runif(1, 0, 360)
  if (bool_func) {
    theta <- runif(1, -30, 30)
  }
  
  data.x <- c()
  data.y <- c()
  
  for (i in 2:n) {
    # idea 1: generate data on axis orthogonal to direction
    # idea 2: generate data in bivariate normal around point
    
    perp_theta <- theta + 90
    # should implement random number of points
    # num_points <- runif(1, 0, 2)
    dist <- rnorm(1, 0, 0.08)
    data.x[i-1] <- x[i-1] + dist*cos(perp_theta*pi/180)
    data.y[i-1] <- y[i-1] + dist*sin(perp_theta*pi/180)
    
    
    x[i] <- x[i-1] + step_size*cos(theta*pi/180)
    y[i] <- y[i-1] + step_size*sin(theta*pi/180)
    
    if (bool_func) {
      old_theta <- theta
      theta <- theta + rnorm(1, 0, 15)
      while (theta > 90 | theta < -90) {
        theta <- old_theta + rnorm(1, 0, 15)
      }
    }
    else { theta <- theta + rnorm(1, 0, 15) }
    # could implement some idea of inertia to prevent sharp switchbacks/encourage smoothness?
  }

  x <- (x-min(x))/(max(x)-min(x)) 
  y <- (y-min(y))/(max(y)-min(y))
  data.x <- (data.x-min(data.x))/(max(data.x)-min(data.x))
  data.y <- (data.y-min(data.y))/(max(data.y)-min(data.y))
  #plot(x, y, type='l')
  #points(data.x, data.y)
  return(cbind(data.x, data.y))
}

