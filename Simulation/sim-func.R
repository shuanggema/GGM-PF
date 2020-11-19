#####################################################################################
# Codes for simulation studies: 
# The functions for generating simulated data & evaluating performances,
# which are used to support numerical simulation studies.
#####################################################################################


#####################################################################################
# Functions for generation of simulated data
# This section includes four functions: 
# generate.data()  tridiag.cor()  Power.law.network()  nearest.neighbor.network()
#####################################################################################
generate.data = function(N,Mu0.list,Theta0.list,Sigma0.list){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: generate.data
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the simulated data.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R packages: mvtnorm
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ N: K0 * 1 vector, the sample sizes of subgroups.
  ## @ Mu0.list: a list including K0 mean vectors (p * 1).
  ## @ Theta0.list: a list including K0 precision matrices (p * p).
  ## @ Sigma0.list: a list including K0 correlation matrices (p * p).
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list "whole.data" including:
  ## @ L0: n * 1 vector, the subgroup labels to which each sample belongs.
  ## @ Mu0: K0 * p matrix, K0 mean vectors.
  ## @ Theta0: K0 * p * p array, K0 precision matrices.
  ## @ data: n * p matrix, the design matrix.
  ## @ n_all: int, the total sample size.
  ## @ K0: int, the true number of subgroups.
  ## ---------------------------------------------------------------------------------------------------------------
  
  K0 = length(Mu0.list)
  p = length(Mu0.list[[1]])
  Mu0=matrix(0,K0,p);L0=NULL;Theta0=array(0, dim = c(p, p, K0));data=NULL
  for (k in 1:K0) {
    Mu0[k,] <- Mu0.list[[k]]
    L0 <- c(L0,rep(k,N[k]))
  }
  for (k in 1:K0) {
    Theta0[,,k] <- as.matrix(Theta0.list[[k]])
  }
  for (k in 1:K0) {
    data <- rbind(data,mvrnorm(N[k],Mu0[k,],Sigma0.list[[k]]))
  }
  n_all = dim(data)[1]
  whole.data=list()
  whole.data$L0=L0
  whole.data$Mu0=Mu0
  whole.data$Theta0=Theta0
  whole.data$data=data
  whole.data$n_all=n_all
  whole.data$K0=K0
  return(whole.data)
}

tridiag.cor = function(p, s=4, rho=0.3){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: tridiag.cor
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the tri-diagonal precision matrices.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding data: No
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ p: Dimensions of the precision matrix.
  ## @ s: The number of sub-networks.
  ## @ rho: The value of non-zero elements on non-diagonal elements.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ A: The precision matrix.
  ## ---------------------------------------------------------------------------------------------------------------
  
  m=p/s
  if(p%%s != 0){
    print("warning! Matrix dimensions cannot be rounded by sub-matrix dimensions.")
  }
  sig.G = matrix(0,m,m)
  for (j in 1:m) {
    ncol0 = c(j-1,j,j+1)
    sig = c(rho, 1, rho)
    ncol = ncol0[which(ncol0>0 & ncol0<=m)]
    sig = sig[which(ncol0>0 & ncol0<=m)]
    sig.G[j,ncol] = sig
  }
  submatrix=list()
  for (ss in 1:s) {
    submatrix[[ss]] = sig.G
  }
  A=submatrix[[1]]
  if(s > 1){
    for (ss in 2:s) {
      A=bdiag(A,submatrix[[ss]])
    }
  }
  return(A)
}

Power.law.network = function(p,s=10,umin=0.1,umax=0.4,I2=0,I3=0){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: Power.law.network
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the s-block power-law precision matrices.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding data: No
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ p: Dimensions of the precision matrix.
  ## @ s: The number of sub-networks.
  ## @ umin: The lower bound of non-zero elements on non-diagonal elements.
  ## @ umax: The upper bound of non-zero elements on non-diagonal elements.
  ## @ I2: The replacement blocks for the precision matrix of the second subgroup.
  ## @ I3: The replacement blocks for the precision matrix of the third subgroup.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ A list including The precision matrices of three subgroups.
  ## ---------------------------------------------------------------------------------------------------------------
  
  pp=p/s
  if(p%%s != 0){
    print("warning! Matrix dimensions cannot be rounded by sub-matrix dimensions.")
  }
  submatrix=list()
  for (ss in 1:s) {
    g = ba.game(pp, m=2,directed = F)
    Eg = as.data.frame(get.edgelist(g))
    subi = diag(1,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }
    for (i in 1:pp) {
      subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+0.1
    }
    submatrix[[ss]]=subi
  }
  submatrix2=submatrix
  submatrix3=submatrix
  if(length(I2)>1 | I2[1]!=0){
    for (ii in 1:length(I2)) {
      iii=I2[ii]
      g = ba.game(pp, m=2,directed = F)
      Eg = as.data.frame(get.edgelist(g))
      subi = diag(1,pp)
      for (q in 1:dim(Eg)[1]) {
        i=Eg[q,1];j=Eg[q,2]
        ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
        subi[i,j]=ij;subi[j,i]=ij
      }
      for (i in 1:pp) {
        subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+0.1
      }
      submatrix2[[iii]] = subi
    }
  }
  if(length(I3)>1 | I3[1]!=0){
    for (ii in 1:length(I3)) {
      iii=I3[ii]
      g = ba.game(pp, m=2,directed = F)
      Eg = as.data.frame(get.edgelist(g))
      subi = diag(1,pp)
      for (q in 1:dim(Eg)[1]) {
        i=Eg[q,1];j=Eg[q,2]
        ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
        subi[i,j]=ij;subi[j,i]=ij
      }
      for (i in 1:pp) {
        subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+0.1
      }
      submatrix3[[iii]] = subi
    }
  }
  A=submatrix[[1]]
  for (ss in 2:s) {
    A=bdiag(A,submatrix[[ss]])
  }
  A = as.matrix(A)
  A2=submatrix2[[1]]
  for (ss in 2:s) {
    A2=bdiag(A2,submatrix2[[ss]])
  }
  A2 = as.matrix(A2)
  A3=submatrix3[[1]]
  for (ss in 2:s) {
    A3=bdiag(A3,submatrix3[[ss]])
  }
  A3 = as.matrix(A3)
  return(list(A1=A,A2=A2,A3=A3))
}

nearest.neighbor.network = function(p,s=10,near=3,umin=0.1,umax=0.4,I2=0,I3=0){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: nearest.neighbor.network
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the s-block nearest neighbor precision matrices.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding data: No
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ p: Dimensions of the precision matrix.
  ## @ s: The number of sub-networks.
  ## @ near: The parameter of nearest neighbors, which controls the degree of sparsity.
  ## @ umin: The lower bound of non-zero elements on non-diagonal elements.
  ## @ umax: The upper bound of non-zero elements on non-diagonal elements.
  ## @ I2: The replacement blocks for the precision matrix of the second subgroup.
  ## @ I3: The replacement blocks for the precision matrix of the third subgroup.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ A list including The precision matrices of three subgroups.
  ## ---------------------------------------------------------------------------------------------------------------
  
  pp=p/s
  if(p%%s != 0){
    print("warning! Matrix dimensions cannot be rounded by sub-matrix dimensions.")
  }
  submatrix=list()
  for (ss in 1:s) {
    subi = diag(1,pp,pp)
    point.matrix = matrix(runif(2*pp,0,1),pp,2)
    num.nearest = matrix(0,pp,near)
    for (j in 1:pp) {
      distancej = apply((t(point.matrix[-j,]) - point.matrix[j,])^2,2,sum)
      corj = setdiff(1:p,j)[order(distancej)[1:near]]
      num.nearest[j,] = corj
    }
    for (oi in 1:dim(num.nearest)[1]) {
      for (oj in 1:dim(num.nearest)[2]) {
        i=oi;j=num.nearest[oi,oj]
        ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
        subi[i,j]=ij;subi[j,i]=ij
      }
    }
    for (i in 1:pp) {
      subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+0.1
    }
    submatrix[[ss]]=subi
  }
  submatrix2=submatrix
  submatrix3=submatrix
  if(length(I2)>1 | I2[1]!=0){
    for (ii in 1:length(I2)) {
      iii=I2[ii]
      subi = diag(1,pp,pp)
      point.matrix = matrix(runif(2*pp,0,1),pp,2)
      num.nearest = matrix(0,pp,near)
      for (j in 1:pp) {
        distancej = apply((t(point.matrix[-j,]) - point.matrix[j,])^2,2,sum)
        corj = setdiff(1:p,j)[order(distancej)[1:near]]
        num.nearest[j,] = corj
      }
      for (oi in 1:dim(num.nearest)[1]) {
        for (oj in 1:dim(num.nearest)[2]) {
          i=oi;j=num.nearest[oi,oj]
          ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
          subi[i,j]=ij;subi[j,i]=ij
        }
      }
      for (i in 1:pp) {
        subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+0.1
      }
      submatrix2[[iii]] = subi
    }
  }
  if(length(I3)>1 | I3[1]!=0){
    for (ii in 1:length(I3)) {
      iii=I3[ii]
      subi = diag(1,pp,pp)
      point.matrix = matrix(runif(2*pp,0,1),pp,2)
      num.nearest = matrix(0,pp,near)
      for (j in 1:pp) {
        distancej = apply((t(point.matrix[-j,]) - point.matrix[j,])^2,2,sum)
        corj = setdiff(1:p,j)[order(distancej)[1:near]]
        num.nearest[j,] = corj
      }
      for (oi in 1:dim(num.nearest)[1]) {
        for (oj in 1:dim(num.nearest)[2]) {
          i=oi;j=num.nearest[oi,oj]
          ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
          subi[i,j]=ij;subi[j,i]=ij
        }
      }
      for (i in 1:pp) {
        subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+0.1
      }
      submatrix3[[iii]] = subi
    }
  }
  A=submatrix[[1]]
  for (ss in 2:s) {
    A=bdiag(A,submatrix[[ss]])
  }
  A = as.matrix(A)
  A2=submatrix2[[1]]
  for (ss in 2:s) {
    A2=bdiag(A2,submatrix2[[ss]])
  }
  A2 = as.matrix(A2)
  A3=submatrix3[[1]]
  for (ss in 2:s) {
    A3=bdiag(A3,submatrix3[[ss]])
  }
  A3 = as.matrix(A3)
  return(list(A1=A,A2=A2,A3=A3))
}

################################################################################
# Functions for evaluating performances of proposed methods
################################################################################
Esti.error = function(data, mu_hat, Theta_hat, Mu0, Theta0, K0, L.mat, L0, prob){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: Esti.error
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Gauging performance of the proposed and alternative approaches
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: 
  ##            R functions: f.den.vec() 
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: ( (q+1)*(p+1)-1 ) * s matrix, all regression coefficients corresponding s choices of given tuning parameters.
  ## @ mu_hat: K0_hat * p matrix, the estimated mean vectors of K0_hat subgroups.
  ## @ Theta_hat: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ other input parameters: Similar to the previous definition.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## The vector including:
  ## @ K: The estimated number of subgroups.
  ## @ CE: The sub-grouping error
  ## @ CME: The mean squared error (MSE) for the mean vectors.
  ## @ PME: The mean squared error (MSE) for the precision matrices.
  ## @ TPR/FPR: The true and false positive rates for the off-diagonal elements of the precision matrices.
  ## -----------------------------------------------------------------------------------------------------------------
  
  p = dim(mu_hat)[2]
  K_hat = dim(mu_hat)[1]
  n_all = dim(data)[1]
  if(K_hat == K0){
    num = rep(0,K0)
    numk = NULL
    for (k in 1:K0) {
      errork = apply((t(mu_hat) - Mu0[k,])^2,2,sum)
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }
    mu_hat = mu_hat[num,]
    Theta_hat = Theta_hat[,,num]
    #  CME & PME
    CME = sqrt(sum((mu_hat - Mu0)^2))/K0
    PME = sqrt(sum((Theta_hat - Theta0)^2))/K0
    
    #  TPR & FPR
    TPk = (apply((Theta_hat!=0) + (Theta0!=0) == 2,3,sum) - p) / (apply(Theta0!=0,3,sum) - p)
    TPk[(is.na(TPk))] = 1
    TPR = sum(TPk) / K0
    FPk = apply((Theta_hat!=0) + (Theta0==0) == 2,3,sum) / apply(Theta0==0,3,sum)
    FPk[(is.na(FPk))] = 0
    FPR = sum(FPk[(!is.na(FPk))]) / K0
    
    #  CE
    f.mat = matrix(0,n_all,K0)
    L.mat = matrix(0,n_all,K0)
    for(k.ind in 1:K_hat) {
      f.mat[,k.ind]=f.den.vec( data, as.numeric(mu_hat[k.ind,]), Theta_hat[,,k.ind] )                  
    }
    for(k.ind in 1:K0) {
      for(i in 1:n_all) {
        L.mat[i,k.ind] = prob[k.ind] * f.mat[i,k.ind] / prob %*% f.mat[i,]
      }
    }
    member = apply(L.mat,1,which.max)
    
    aa = L0
    cap_matrix0 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = member
    cap_matrix1 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n_all*(n_all-1)/2)
  } else{
    num = rep(0,K_hat)
    for (k in 1:K_hat) {
      mu_hatk = mu_hat[k,]
      errork.mu = apply((t(Mu0) - mu_hatk)^2,2,sum)
      Theta_hatk = Theta_hat[,,k]
      errork.Theta = apply((Theta0 - rep(Theta_hatk,K0))^2,3,sum)
      errork = errork.mu + errork.Theta
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }
    Mu0.re = Mu0[num,]
    Theta0.re = Theta0[,,num]
    if(K_hat == 1){
      Mu0.re = as.matrix(t(Mu0.re))
      Theta0.re = as.array(Theta0.re)
      dim(Theta0.re) <- c(p,p,K_hat)
    }
    
    #  CME & PME
    CME = sqrt(sum((mu_hat - Mu0.re)^2))/K_hat
    PME = sqrt(sum((Theta_hat - Theta0.re)^2))/K_hat
    
    #  TPR & FPR
    TPk = (apply((Theta_hat!=0) + (Theta0.re!=0) == 2,3,sum) - p) / (apply(Theta0.re!=0,3,sum) - p)
    TPk[(is.na(TPk))] = 1
    TPR = sum(TPk) / K_hat
    FPk = apply((Theta_hat!=0) + (Theta0.re==0) == 2,3,sum) / apply(Theta0.re==0,3,sum)
    FPk[(is.na(FPk))] = 0
    FPR = sum(FPk[(!is.na(FPk))]) / K_hat
    
    #  CE
    num = rep(0,K0)
    numk = NULL
    L.hat = rep(0,n_all)
    for (k in 1:K0) {
      errork = apply((t(mu_hat) - Mu0[k,])^2,2,sum)
      numk = which(errork == min(errork))[1]
      num[k] = numk
      L.hat[which(L0 == k)] = numk
    }
    if(K_hat < K0){L.hat=L0}
    
    f.mat = matrix(0,n_all,K_hat)
    L.mat = matrix(0,n_all,K_hat)
    for(k.ind in 1:K_hat) {
      f.mat[,k.ind]=f.den.vec( data, as.numeric(mu_hat[k.ind,]), Theta_hat[,,k.ind] )                  
    }
    for(k.ind in 1:K_hat) {
      for(i in 1:n_all) {
        L.mat[i,k.ind] = prob[k.ind] * f.mat[i,k.ind] / prob %*% f.mat[i,]
      }
    }
    member = apply(L.mat,1,which.max)
    
    aa = L.hat
    cap_matrix0 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = member
    cap_matrix1 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n_all*(n_all-1)/2)
  }
  index = as.data.frame(t(c(K_hat,CE,CME,PME,TPR,FPR)))
  names(index) = c("K","CE","CME","PME","TPR","FPR")
  
  return(index)
}

