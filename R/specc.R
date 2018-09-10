specc_p<-function(x, centers, kernel = "rbfdot", kpar = "automatic", 
 nystrom.red = FALSE, nystrom.sample = dim(x)[1]/6, iterations = 200, 
 mod.sample =  0.75, na.action = na.omit)
{
  print("specc starts")
  print(Sys.time())
  x <- na.action(x)                 
  rown <- rownames(x)               
  x <- as.matrix(x)                 
  m <- nrow(x)                      
  if (missing(centers))
    stop("centers must be a number or a matrix")
  if (length(centers) == 1) {       
    nc <-  centers
    if (m < centers)                
      stop("more cluster centers than data points.")
  }
  else
    nc <- dim(centers)[2]           
  
  
  if(is.character(kpar)) {          
    kpar <- match.arg(kpar,c("automatic","local"))
    
    if(kpar == "automatic")
    {
      if (nystrom.red == TRUE)
        sam <- sample(1:m, floor(mod.sample*nystrom.sample))
      else
        sam <- sample(1:m, floor(mod.sample*m)) 
      
      
      sx <- unique(x[sam,])       
      ns <- dim(sx)[1]            
      dota <- rowSums(sx*sx)/2    
      ktmp <- crossprod(t(sx))     
      for (i in 1:ns)
        ktmp[i,]<- 2*(-ktmp[i,] + dota + rep(dota[i], ns))
      
      
      ## fix numerical prob.
      ktmp[ktmp<0] <- 0
      ktmp <- sqrt(ktmp)          
      print("distance matrix completes")
      print(Sys.time())
      
      kmax <- max(ktmp)
      kmin <- min(ktmp + diag(rep(Inf,dim(ktmp)[1])))
      kmea <- mean(ktmp)
      lsmin <- log2(kmin)
      lsmax <- log2(kmax)
      midmax <- min(c(2*kmea, kmax))
      midmin <- max(c(kmea/2,kmin))
      rtmp <- c(seq(midmin,0.9*kmea,0.05*kmea), seq(kmea,midmax,0.08*kmea))
      if ((lsmax - (Re(log2(midmax))+0.5)) < 0.5) step <- (lsmax - (Re(log2(midmax))+0.5))
      else step <- 0.5
      if (((Re(log2(midmin))-0.5)-lsmin) < 0.5 ) stepm <-  ((Re(log2(midmin))-0.5) - lsmin)
      else stepm <- 0.5
      
      tmpsig <- c(2^(seq(lsmin,(Re(log2(midmin))-0.5), stepm)), rtmp, 2^(seq(Re(log2(midmax))+0.5, lsmax,step)))
      diss <- matrix(rep(Inf,length(tmpsig)*nc),ncol=nc)
      
      print("length of tmpsig: ")
      print(length(tmpsig))
      print("for loop starts")
      print(Sys.time())
      
      cl = makeCluster(20, outfile="")
      registerDoParallel(cl)
     
      list_diss=foreach (i = 1:length(tmpsig))%dopar%
      # for (i in 1:length(tmpsig)) 
      {
        ka <- exp((-(ktmp^2))/(2*(tmpsig[i]^2)))
        diag(ka) <- 0
        
        d <- 1/sqrt(rowSums(ka))
        
        if(!any(d==Inf) && !any(is.na(d))&& (max(d)[1]-min(d)[1] < 10^4))
        {
          l <- d * ka %*% diag(d)
          xi <- eigen(l,symmetric=TRUE)$vectors[,1:nc]
          yi <- xi/sqrt(rowSums(xi^2))
          res <- kmeans(yi, centers, iterations)
          # diss[i,] <- res$withinss
          # print(i)
          # print(res$withinss)
          return(res$withinss)
        }
      }
      stopCluster(cl)
      print("end of for loop")
      print(Sys.time())
      # print("diss: ")
      # print(diss)
     
      for(i in 1:length(list_diss))
      {
        diss[i,]=list_diss[[i]]
      }
      # print("diss")
      # print(diss)
     
      ms <- which.min(rowSums(diss))
      kernel <- rbfdot((tmpsig[ms]^(-2))/2)
      
      ## Compute Affinity Matrix
      if (nystrom.red == FALSE)
        km <- kernelMatrix(kernel, x)
    }
    if (kpar=="local")
    {
      if (nystrom.red == TRUE)
        stop ("Local Scaling not supported for nystrom reduction.")
      s <- rep(0,m)
      dota <- rowSums(x*x)/2
      dis <- crossprod(t(x))
      for (i in 1:m)
        dis[i,]<- 2*(-dis[i,] + dota + rep(dota[i],m))
      
      ## fix numerical prob.
      dis[dis < 0] <- 0
      
      for (i in 1:m)
        s[i] <- median(sort(sqrt(dis[i,]))[1:5])
      
      ## Compute Affinity Matrix
      km <- exp(-dis / s%*%t(s))
      kernel <- "Localy scaled RBF kernel"
      
      
    }
  }
  else 
  {
    if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
    if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")
    
    ## Compute Affinity Matrix
    if (nystrom.red == FALSE)
      km <- kernelMatrix(kernel, x)
  }
  
  
  
  if (nystrom.red == TRUE){
    
    n <- floor(nystrom.sample)
    ind <- sample(1:m, m)
    x <- x[ind,]
    
    tmps <- sort(ind, index.return = TRUE)
    reind <- tmps$ix
    A <- kernelMatrix(kernel, x[1:n,])
    B <- kernelMatrix(kernel, x[-(1:n),], x[1:n,])
    d1 <- colSums(rbind(A,B))
    d2 <- rowSums(B) + drop(matrix(colSums(B),1) %*% .ginv(A)%*%t(B))
    dhat <- sqrt(1/c(d1,d2))
    
    A <- A * (dhat[1:n] %*% t(dhat[1:n]))
    B <- B * (dhat[(n+1):m] %*% t(dhat[1:n]))
    
    Asi <- .sqrtm(.ginv(A))
    Q <- A + Asi %*% crossprod(B) %*% Asi
    tmpres <- svd(Q)
    U <- tmpres$u
    L <- tmpres$d
    V <- rbind(A,B) %*% Asi %*% U %*% .ginv(sqrt(diag(L)))
    yi <- matrix(0,m,nc)
    
    ## for(i in 2:(nc +1))
    ##   yi[,i-1] <- V[,i]/V[,1]
    
    for(i in 1:nc) ## specc
      yi[,i] <- V[,i]/sqrt(sum(V[,i]^2))
    
    res <- kmeans(yi[reind,], centers, iterations)
    
  }
  else{
    if(is(kernel)[1] == "rbfkernel")
      diag(km) <- 0
    
    d <- 1/sqrt(rowSums(km))
    l <- d * km %*% diag(d)
    xi <- eigen(l)$vectors[,1:nc]
    yi <- xi/sqrt(rowSums(xi^2))
    res <- kmeans(yi, centers, iterations)
  }
  
  cent <- matrix(unlist(lapply(1:nc,ll<- function(l){colMeans(x[which(res$cluster==l), ,drop=FALSE])})),ncol=dim(x)[2], byrow=TRUE)
  
  withss <- unlist(lapply(1:nc,ll<- function(l){sum((x[which(res$cluster==l),, drop=FALSE] - cent[l,])^2)}))
  names(res$cluster) <- rown
  print("end of specc")
  print(Sys.time())
  return(new("specc", .Data=res$cluster, size = res$size, centers=cent, withinss=withss, kernelf= kernel))
  
}

