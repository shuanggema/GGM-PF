##################################################################################
# This document includes main functions of proposed methods, 
# which are used to support numerical simulation studies and real data analysis
# in the paper
# "Gaussian Graphical Model-based Heterogeneity Analysis via Penalized Fusion"
##################################################################################


############################# Functions for main algorithms ############################
FGGM = function(data, K, lambda1 = 0.5, lambda2 = 0.2, lambda3 = 2, a = 3, rho = 1, 
                eps = 5e-2, niter = 20, maxiter=10, maxiter.AMA=5, initialization=T, initialize, 
                average=F, asymmetric=T, local_appro=T){

  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: FGGM
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            The key function of Gaussian graphical model-based heterogeneity analysis via penalized fusion:
  ##            identifying K_0 and reconstructing the network structure.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: Update_Theta()    AMA_XI()   mcp_d()   S_soft()    MCP_soft()    f.den.vec()    Symmetrize()    cut_diff_ama()
  ##            R packages: NbClust
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ K: int, a selected upper bound of K_0.
  ## @ lambda1: a float value, the tuning parameter controlling the sparse of the mean parameter. 
  ## @ lambda2: a float value, the tuning parameter controlling the sparse of the precision matrix.
  ## @ lambda3: a float value, the tuning parameter controlling the number of subgroup.
  ## @ a: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ rho: a float value, the penalty parameter in ADMM algorithm of updating precision matrix Theta, the default setting is 1.
  ## @ eps: a float value, algorithm termination threshold.
  ## @ niter: int, Maximum number of cycles of the EM algorithm, the default setting is 20.
  ## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
  ## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
  ## @ initialization: the logical variable, whether to calculate the initial value, the default setting is T,
  ##                                     if initialization = F, the initial value uses initialize.
  ## @ initialize: a given initial value used if initialization = F.
  ## @ average: the logical variable, whether to use averaging when integrating parameters that are identified as identical subgroups,
  ## @          the default setting is F, which means the estimated parameters for the subgroup with the largest sample size among 
  ## @          the subgroups identified as identical subgroups is used as the final parameter for this subgroup.
  ## @ asymmetric: the logical variable, symmetry of the precision matrices or not, the default setting is T.
  ## @ local_appro: the logical variable, whether to use local approximations when updating mean parameters, the default setting is T.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list FGGM_res including:
  ## @ mu: K0_hat * p matrix, the estimated mean vectors of K0_hat subgroups (K0_hat is the number of identified subgroups).
  ## @ Theta: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ Xi: p * p * K0_hat array, the estimated dual variables of the precision matrices in ADMM.
  ## @ niter: int, the actual number of cycles of the algorithm.
  ## @ diff_Xi: a float value, the L2-norm difference between subgroups.
  ## @ prob0: K0_hat * 1 vector, the estimated mixture probabilities of subgroups.
  ## @ L.mat0: n * K0_hat matrix, the estimated probability that each sample belongs to each subgroup.
  ## @ Theta0: p * p * K array, the estimated original precision matrices of given K subgroups.
  ## @ Xi0: p * p * K array, the estimated original dual variables of the precision matrices of given K subgroups in ADMM.
  ## @ mu0: K * p matrix, the estimated original mean vectors of K subgroups.
  ## @ group: a list, the sequence number partition of original K subgroups compressed to K0_hat.
  ## @ bic: a float value, the BIC value corresponding the choice of given tuning parameters.
  ## @ fit.error: a float value, the value of the loss function (without penalty function) corresponding the choice of given tuning parameters.
  ## @ df: a float value, the penalty value for non-zero parameters corresponding the choice of given tuning parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])
  K_c <- combn(K,2)
  if (initialization) {
    out.initial = initialize_fuc.dbscan(data,K)
    prob = out.initial$prob
    mu = out.initial$Mu
    Theta = out.initial$Theta
    memb = out.initial$memb
    L.mat = matrix(0,n,K)
    for(jj in 1:n) L.mat[jj, memb[jj]]=1
  } else {
    Theta = initialize$Theta
    mu = initialize$Mu
    prob = initialize$prob
    L.mat = initialize$L.mat
    memb = apply(L.mat,1,which.max)
  }
  
  # EM algorithm 
  f.mat = matrix(0,n,K)
  L.mat.old = matrix(0,n,K)
  mu.old = matrix(0, K, p)
  Theta.old = array(0, dim = c(p, p, K))
  
  t = 0
  diff_mu = 10
  diff_theta = 10
  while( diff_theta >= eps  && t < niter )
  {  
    prob.old = prob
    mu.old = mu
    Theta.old = Theta
    L.mat.old = L.mat
    
    # calculate pdfs
    for(k.ind in 1:K) {
      f.mat[,k.ind]=f.den.vec( data, as.numeric(mu.old[k.ind,]), Theta.old[,,k.ind] )                  
    }
    
    # update L and pi
    for(k.ind in 1:K) {
      for(i in 1:n) {
        L.mat[i,k.ind] = prob.old[k.ind] * f.mat[i,k.ind] / prob.old %*% f.mat[i,]
      }
      prob[k.ind] = mean(L.mat[,k.ind])
    }
    nK = apply(L.mat,2,sum)
    
    Theta_kk_diff <- rep(0,dim(K_c)[2])
    for (l in 1:dim(K_c)[2]) {
      Theta_kk_diff[l] <- sum((Theta.old[,,K_c[1,l]] - Theta.old[,,K_c[2,l]])^2)
    }
    
    # update mean vectors                       
    for(j in 1:p){   
      for(k.ind in 1:K) {  
        tmp = t(t(data) - as.numeric(mu[k.ind,])) %*% Theta.old[,j,k.ind] + mu[k.ind,j] * Theta.old[j,j,k.ind]
        hj = t(L.mat[,k.ind])%*%tmp
        tau_k = sqrt(apply((t(mu[k.ind,] - t(mu[-k.ind,])))^2,1,sum) + Theta_kk_diff[ceiling(which(K_c == k.ind)/2)])
        v_k = sum(mcp_d(tau_k, lambda3, a) / (tau_k+0.00001) * mu[-k.ind,j])
        mcp_lambda1 = mcp_d(mu[k.ind,j], lambda1, a)
        v_k_hat = sum(mcp_d(tau_k, lambda3, a) / (tau_k+0.00001))
        if(local_appro){
          mu[k.ind,j] = (hj + n*v_k) / (nK[k.ind] * Theta.old[j,j,k.ind] + n*v_k_hat + n*mcp_lambda1/(abs(mu[k.ind,j])+0.00001))
        } else {
          if(n*mcp_lambda1 >= abs(hj + n*v_k)){
            mu[k.ind,j] = 0
          }else{
            mu[k.ind,j] = (hj + n*v_k - n*mcp_lambda1*sign(mu[k.ind,j])) / (nK[k.ind] * Theta.old[j,j,k.ind] + n*v_k_hat)
          }
        }
      }
    }
    mu[abs(mu) < 1e-3] <- 0
    
    # update precision matrices
    mu_kk_diff <- rep(0,dim(K_c)[2])
    for (l in 1:dim(K_c)[2]) {mu_kk_diff[l] <- sum((mu[K_c[1,l],] - mu[K_c[2,l],])^2)}
    # define the pseudo sample covariance matrices \tilde{S}
    S = array(0, dim = c(p, p, K))
    for (k.ind in 1:K) {
      nk = nK[k.ind]
      L_ikx = sqrt(L.mat[,k.ind])*t(t(data) - mu[k.ind,])
      S[,,k.ind] = t(L_ikx) %*% L_ikx / nK[k.ind]
    }
    
    Theta_out = Update_Theta(S,nK,lambda2,lambda3, mu_kk_diff, K_c, tol=eps, maxiter=maxiter, maxiter.AMA=maxiter.AMA)
    Theta = Theta_out$Theta
    Xi = Theta_out$Xi
    V_kk = Theta_out$V_kk
    
    t = t + 1
    diff_mu = norm(mu.old-mu,type="2")/(norm(mu,type="2")+0.001)
    diff_theta = norm(Theta.old-Theta,type="2")/(norm(Theta,type="2")+0.001)
  }
  
  group_final = cut_diff_ama(V_kk,K_c,K,cutoff=0.01)
  K_0 = length(group_final)
  mu_final = matrix(0, K_0, p)
  Theta_final = array(0, dim = c(p, p, K_0))
  Xi_final = array(0, dim = c(p, p, K_0))
  prob0 = rep(0,K_0)
  L.mat0 = matrix(0,n,K_0)
  for (l in 1:K_0) {
    gg = group_final[[l]]
    prob0[l] = sum(prob[gg])
    if(length(gg) > 1){
      L.mat0[,l] = apply(L.mat[,gg],1,sum)
      mu_final[l,] = apply(mu[gg,],2,mean)
      
      if(!average){
        Theta_final[,,l] = Theta[,,gg[which.max(nK[gg])]] 
        Xi_final[,,l] = Xi[,,gg[which.max(nK[gg])]]       
      } else {
        Theta_final[,,l] = Theta[,,gg[1]]/length(gg)
        Xi_final[,,l] = Xi[,,gg[1]]/length(gg)
        for (gi in gg[-1]) {
          Theta_final[,,l] = Theta_final[,,l] + Theta[,,gi]/length(gg)
          Xi_final[,,l] = Xi_final[,,l] + Xi[,,gi]/length(gg)
        }
      }
      
    }else{ mu_final[l,] = mu[gg,]; Theta_final[,,l] = Theta[,,gg]; Xi_final[,,l] = Xi[,,gg]; L.mat0 = matrix(1,n,K_0) }
  }
  
  if(asymmetric){
    for(k in 1:K_0) {
      Theta_final[,,k] = Symmetrize(Theta_final[,,k])
      Xi_final[,,k] = Symmetrize(Xi_final[,,k])
    }
  }
  
  Theta_final[abs(Theta_final) < 1e-3] <- 0
  member = apply(L.mat,1,which.max)
  BIC.res = BIC(data, mu_final, Xi_final, L.mat)
  
  FGGM_res <- list();FGGM_res$mu <-  mu_final;FGGM_res$Theta <- Theta_final
  FGGM_res$Xi <- Xi_final; FGGM_res$niter <- t; FGGM_res$diff_Xi <- V_kk
  FGGM_res$prob0 <- prob0; FGGM_res$L.mat0 <- L.mat0;
  FGGM_res$Theta0 <- Theta; FGGM_res$member <- member;
  FGGM_res$Xi0 <- Xi;FGGM_res$mu0 <- mu;FGGM_res$group <- group_final
  FGGM_res$bic <- BIC.res$bic;FGGM_res$fit.error <- BIC.res$fit.error;FGGM_res$df <- BIC.res$df
  return(FGGM_res)
}

Update_Theta = function(S, nK, lambda2, lambda3, mu_kk_diff, K_c, a = 3, rho=1, 
                        maxiter=10, maxiter.AMA=5, tol=1e-2, rho.increment=1){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: Update_Theta
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Updating the precision matrices using the ADMM algorithm after updating the mean vectors in the EM algorithm. 
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: AMA_XI()
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ S: p * p * K array, the estimated pseudo sample covariance matrices of given K subgroups in the iteration of the EM algorithm.
  ## @ nK: K * 1 vector, the estimated sample sizes of K subgroups in the iteration of the EM algorithm.
  ## @ lambda2: a float value, the tuning parameter controlling the sparse of the precision matrix.
  ## @ lambda3: a float value, the tuning parameter controlling the number of subgroup.
  ## @ mu_kk_diff: L2-norm differences in the mean vector between different subgroups.
  ## @ K_c: All combinations of natural numbers within K.
  ## @ a: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ rho: a float value, the penalty parameter in ADMM algorithm of updating precision matrix Theta, the default setting is 1.
  ## @ tol: a float value, algorithm termination threshold.
  ## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
  ## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list Theta_out including:
  ## @ Theta: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ Xi: p * p * K0_hat array, the estimated dual variables of the precision matrices in ADMM.
  ## @ V_kk: L2-norm differences in the precision matrices between different subgroups.
  ## @ iters: int, the actual number of cycles of the algorithm.
  ## @ diff: a float value, the L1-norm difference between Theta and Xi.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  p = dim(S)[1]
  K = dim(S)[3]
  nkk = rep(1,K)
  # initialize Theta:
  Theta = array(0, dim = c(p, p, K))
  for (k in 1:K) { Theta[,,k] =  diag(1,p)}
  # initialize Xi:
  Xi = array(0, dim = c(p, p, K))
  # initialize Phi:
  Phi = array(0, dim = c(p, p, K))
  
  iter=0
  diff_value = 10
  diff_value_Xi = 10
  while( iter<maxiter && !(diff_value < tol || diff_value_Xi < tol) )
  {
    Theta.prev = Theta
    Xi.prev = Xi
    Theta=Theta.prev
    # update Theta:
    for(k.ind in 1:K){
      Sa=S[,,k.ind] - rho*Xi[,,k.ind]/nkk[k.ind] + rho*Phi[,,k.ind]/nkk[k.ind]
      edecomp = eigen(Sa)
      D = edecomp$values
      if(is.complex(D)){Sa=Symmetrize(Sa);edecomp = eigen(Sa);D = edecomp$values}
      V = edecomp$vectors
      D2 = nkk[k.ind]/(2*rho) * ( -D + sqrt(D^2 + 4*(rho/nkk[k.ind]) ) )
      Theta[,,k.ind] = V %*% diag(D2) %*% t(V)
    }
    # update Xi:
    # define B matrices:
    B = array(0, dim = c(p, p, K))
    for(k in 1:K){ B[,,k] = Theta[,,k] + Phi[,,k] }
    
    Xi_out_list = AMA_XI(B,K_c,lambda2,lambda3,mu_kk_diff,maxiter=maxiter.AMA)
    Xi = Xi_out_list$Xi
    Xi[abs(Xi) < 1e-3] <- 0
    V = Xi_out_list$V
    V_kk = round(apply(V^2,3,sum),4)
    
    # update the dual variable Phi:
    Phi = Phi + Theta - Xi
    iter = iter+1
    diff_value = sum(abs(Theta - Theta.prev)) / sum(abs(Theta.prev))
    diff_value_Xi = sum(abs(Xi - Xi.prev)) / (sum(abs(Xi.prev))+0.001)
    # increment rho by a constant factor:
    rho = rho*rho.increment
  }
  diff = sum(abs(Theta-Xi))
  Theta_out = list(Theta=Theta,Xi=Xi,V_kk=V_kk,diff=diff,iters=iter)
  return(Theta_out)
}

AMA_XI = function(B, K_c, lambda2, lambda3, mu_kk_diff, a = 3, 
                  kappa=1, maxiter=5, tol=1e-2){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: AMA_XI
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Solving (A.11) in Supplementary Materials using S-AMA algorithm, which is the key step of updating the precision matrices. 
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: S_soft()   MCP_soft()
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ B: The output of the previous step in ADMM algorithm (see (A.11) in Supplementary Materials).
  ## @ K_c: All combinations of natural numbers within K.
  ## @ lambda2: a float value, the tuning parameter controlling the sparse of the precision matrix.
  ## @ lambda3: a float value, the tuning parameter controlling the number of subgroup.
  ## @ mu_kk_diff: L2-norm differences in the mean vector between different subgroups.
  ## @ a: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ kappa: a float value, the penalty parameter in S-AMA algorithm, the default setting is 1.
  ## @ tol: a float value, algorithm termination threshold.
  ## @ maxiter: int, Maximum number of cycles of the algorithm.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list Xi_out including:
  ## @ Xi: p * p * K0_hat array, the estimated dual variables of the precision matrices in ADMM.
  ## @ V: p * p * K0_hat array, the estimated dual variables of Xis in S-AMA.
  ## @ Delta: Lagrange multiplier in in S-AMA.
  ## @ iters: int, the actual number of cycles of the algorithm.
  ## @ diff: a float value, the convergence criteria value.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  p = dim(B)[1]
  K = dim(B)[3]
  num_Kc = dim(K_c)[2]
  # initialize Xi:
  Xi = array(0, dim = c(p, p, K))
  for (k in 1:K) { Xi[,,k] =  diag(1,p)}
  diag_index = Xi
  # initialize V:
  V = array(0, dim = c(p, p, num_Kc))
  # initialize Delta:
  Delta = array(0, dim = c(p, p, num_Kc))
  
  Z = array(0, dim = c(p, p, K))
  Omega = array(0, dim = c(p, p, num_Kc))
  e_k12 = matrix(0,K,num_Kc)
  for (i in 1:K) {e_k12[i,which(K_c[1,] == i)] = 1;e_k12[i,which(K_c[2,] == i)] = -1}
  
  # iterations
  iter=0
  diff_value = 10
  while( iter<maxiter && diff_value > tol )
  {
    V.prev = V
    Xi.prev = Xi
    Delta.prev = Delta
    for (i in 1:p) {
      for (j in 1:p) {
        Z[i,j,] = B[i,j,] + apply(t(Delta[i,j,] * t(e_k12)),1,sum)
      }
    }
    Xi = S_soft(Z,lambda2)/(1-1/a) * (abs(Z) <= a*lambda2 & diag_index==0) + Z * (abs(Z) > a*lambda2 | diag_index==1)

    for (l in 1:num_Kc) {
      Omega[,,l] = Xi[,,K_c[1,l]] - Xi[,,K_c[2,l]] - Delta[,,l]/kappa
      Omegal2 = sum(Omega[,,l]^2)
      V[,,l] = Omega[,,l] * MCP_soft(sqrt(Omegal2 + mu_kk_diff[l]),lambda3/kappa) / (sqrt(Omegal2+mu_kk_diff[l])+1e-8) 
      Delta[,,l] = Delta[,,l] + kappa * ( V[,,l] - Xi[,,K_c[1,l]] + Xi[,,K_c[2,l]] )* (sum(V[,,l]^2) > 0 )
    }
    
    iter = iter+1
    diff_value = sum(abs(Xi - Xi.prev)^2) / sum(abs(Xi.prev)^2)
  }
  Xi_out = list(Xi=Xi,V=V,Delta=Delta,diff=diff_value,iters=iter)
  return(Xi_out)
}

FGGM.refit = function(data, K, lambda1 = 0.5, lambda2 = 0.2, lambda3 = 2, a = 3, rho = 1, 
                      eps = 5e-2, niter = 20, maxiter=10, maxiter.AMA=5, initialization=T, initialize, 
                      average=F, asymmetric=T, local_appro=T){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: FGGM.refit
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Refitting when K0_hat is identified.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: FGGM()    Update_Theta()    AMA_XI()
  ##            R packages: NbClust
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ K: int, a selected upper bound of K_0.
  ## @ lambda1: a float value, the tuning parameter controlling the sparse of the mean parameter. 
  ## @ lambda2: a float value, the tuning parameter controlling the sparse of the precision matrix.
  ## @ lambda3: a float value, the tuning parameter controlling the number of subgroup.
  ## @ a: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ rho: a float value, the penalty parameter in ADMM algorithm of updating precision matrix Theta, the default setting is 1.
  ## @ eps: a float value, algorithm termination threshold.
  ## @ niter: int, Maximum number of cycles of the algorithm, the default setting is 20.
  ## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
  ## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
  ## @ initialization: the logical variable, whether to calculate the initial value, the default setting is T,
  ##                                     if initialization = F, the initial value uses initialize.
  ## @ initialize: A given initial value used if initialization = F.
  ## @ average: the logical variable, whether to use averaging when integrating parameters that are identified as identical subgroups,
  ## @          the default setting is F, which means the estimated parameters for the subgroup with the largest sample size among 
  ## @          the subgroups identified as identical subgroups is used as the final parameter for this subgroup.
  ## @ asymmetric: the logical variable, symmetry of the precision matrices or not, the default setting is T.
  ## @ local_appro: the logical variable, whether to use local approximations when updating mean parameters, the default setting is T.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list PP.refit including:
  ## @ mu: K0_hat * p matrix, the estimated mean vectors of K0_hat subgroups (K0_hat is the number of identified subgroups).
  ## @ Theta: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ Xi: p * p * K0_hat array, the estimated dual variables of the precision matrices in ADMM.
  ## @ niter: int, the actual number of cycles of the algorithm.
  ## @ diff_Xi: a float value, the L2-norm difference between subgroups.
  ## @ prob0: K0_hat * 1 vector, the estimated mixture probabilities of subgroups.
  ## @ L.mat0: n * K0_hat matrix, the estimated probability that each sample belongs to each subgroup.
  ## @ Theta0: p * p * K array, the estimated original precision matrices of given K subgroups.
  ## @ Xi0: p * p * K array, the estimated original dual variables of the precision matrices of given K subgroups in ADMM.
  ## @ mu0: K * p matrix, the estimated original mean vectors of K subgroups.
  ## @ group: a list, the sequence number partition of original K subgroups compressed to K0_hat.
  ## @ bic: a float value, the BIC value corresponding the choice of given tuning parameters.
  ## @ fit.error: a float value, the value of the loss function (without penalty function) corresponding the choice of given tuning parameters.
  ## @ df: a float value, the penalty value for non-zero parameters corresponding the choice of given tuning parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  set.seed(1)
  PP = FGGM(data, K, lambda1, lambda2, lambda3, eps=eps, maxiter=maxiter, maxiter.AMA=maxiter.AMA, 
            initialization=initialization, initialize=initialize, average=average, asymmetric=asymmetric, local_appro=local_appro)
  K_hat = length(PP$group)
  if(K_hat == 1){
    return(PP)
  } else {
    set.seed(1)
    PP.refit = FGGM(data, K_hat, lambda1, lambda2, 0, eps=eps, maxiter=maxiter, maxiter.AMA=maxiter.AMA, 
                    initialization=T, initialize=initialize, average=average, asymmetric=asymmetric, local_appro=local_appro)
  }
  return(PP.refit)
}


############################## Functions for tuning lambdas ############################
BIC = function(data, mu_hat, Theta_hat, L.mat){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: BIC
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the adaptive BIC-type criterion.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: f.den.vec()
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ mu_hat: K0_hat * p matrix, the estimated mean vectors of K0_hat subgroups.
  ## @ Theta_hat: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ L.mat: n * K0_hat matrix, the estimated probability that each sample belongs to each subgroup.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list P including:
  ## @ fit.error: a float value, the value of the loss function (without penalty function).
  ## @ df: a float value, the penalty value for non-zero parameters corresponding the choice of given tuning parameters.
  ## @ bic: a float value, the BIC value corresponding the choice of given tuning parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  n = nrow(data)
  K = nrow(mu_hat)
  
  # fitting error
  pi_vec = apply(L.mat, 2, sum)/n
  fit.error_mat = matrix(0, n, K)
  for(k in 1:K) {
    fit.error_mat[,k] = pi_vec[k] * f.den.vec( data, as.numeric(mu_hat[k,]), Theta_hat[,,k] )                  
  }
  fit0 = apply(fit.error_mat, 1, sum)
  fit.error = sum(log( fit0 + min(fit0[fit0>0]) ))
  fit.error = - 2*fit.error
  
  # degrees of freedom
  for(i in 1:K){
    Theta_hat[upper.tri(Theta_hat[, , i], diag = T)] = 0
  }
  
  df =  log(n) * length(which(mu_hat != 0)) + 2 * length(which(Theta_hat != 0))
  bic = fit.error + df
  P = list()
  P$fit.error = fit.error
  P$df = df
  P$bic = bic
  return(P)
}

genelambda.obo = function(nlambda1=10,lambda1_max=1,lambda1_min=0.05,
                          nlambda2=10,lambda2_max=1,lambda2_min=0.01,
                          nlambda3=10,lambda3_max=5,lambda3_min=0.5){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: genelambda.obo
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating a sequence of the tuning parameters (lambda1, lambda2, and lambda3).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ nlambda1, nlambda2, and nlambda3: The numbers of lambda 1 2 3.
  ## @ lambda1_min, lambda2_min, and lambda3_min: The minimum values of lambda 1 2 3.
  ## @ lambda1_max, lambda2_max, and lambda3_max: The maximum values of lambda 1 2 3.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ lambda: a sequence of the tuning parameters (lambda1, lambda2, and lambda3).   
  ## -----------------------------------------------------------------------------------------------------------------
  
  lambda1 = exp(seq(log(lambda1_max),log(lambda1_min),len= nlambda1))
  lambda2 =exp(seq(log(lambda2_max),log(lambda2_min),len= nlambda2))
  lambda3 =exp(seq(log(lambda3_max),log(lambda3_min),len= nlambda3))
  lambda = list(lambda1=lambda1,lambda2=lambda2,lambda3=lambda3)
  return(lambda)
}

tuning.lambda.FGGM = function(lambda, data, K, initial.selection="K-means", 
                              initialize, average=F, asymmetric=T, eps = 5e-2, maxiter=10, maxiter.AMA=5,
                              local_appro=T, trace = T){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: tuning.lambda.FGGM
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Searching and selecting the optional tuning parameters under the adaptive BIC-type criterion
  ##            using the proposed method.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: FGGM.refit()
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ lambda: a list, the sequences of the tuning parameters (lambda1, lambda2, and lambda3).
  ## @ data: n * p matrix, the design matrix.
  ## @ K: int, a selected upper bound of K_0.
  ## @ initial.selection: the different initial values from two clustering methods, which can be selected from c("K-means","dbscan").
  ## @ initialize: A given initial values, which should be given when initial.selection is not in c("K-means","dbscan").
  ## @ average: the logical variable, whether to use averaging when integrating parameters that are identified as identical subgroups,
  ## @          the default setting is F, which means the estimated parameters for the subgroup with the largest sample size among 
  ## @          the subgroups identified as identical subgroups is used as the final parameter for this subgroup.
  ## @ asymmetric: the logical variable, symmetry of the precision matrices or not, the default setting is T.
  ## @ eps: a float value, algorithm termination threshold.
  ## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
  ## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
  ## @ local_appro: the logical variable, whether to use local approximations when updating mean parameters, the default setting is T.
  ## @ trace: the logical variable, whether or not to output the number of identified subgroups during the search for parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list "result" including:
  ## @ Opt_lambda: the selected optional tuning parameters.
  ## @ Mu_hat.list: the estimated mean vectors of K0_hat subgroups corresponding all choices of given tuning parameters.
  ## @ Theta_hat.list: the estimated precision matrices of K0_hat subgroups corresponding all choices of given tuning parameters.
  ## @ prob.list: the estimated mixture probabilities of subgroups corresponding all choices of given tuning parameters.
  ## @ member.list: subgroup labels to which each sample belongs corresponding all choices of given tuning parameters.
  ## @ L.mat.list: the estimated probability that each sample belongs to each subgroup corresponding all choices of given tuning parameters.
  ## @ Opt_aBIC: the optional BIC value.
  ## @ Opt_num: the position of the optimal parameter for all given parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  lambda1 = lambda$lambda1
  lambda2 = lambda$lambda2
  lambda3 = lambda$lambda3
  L1 = length(lambda1)
  L2 = length(lambda2)
  L3 = length(lambda3)
  L = L1+L2+L3

  aBIC = rep(0,L)
  n_all = dim(data)[1]
  # initialize
  if(initial.selection=="K-means"){
    set.seed(1)
    out.initial = initialize_fuc(data,K)
    memb = out.initial$memb
    L.mat = matrix(0,n_all,K)
    for(jj in 1:n_all) L.mat[jj, memb[jj]]=1
    out.initial$L.mat = L.mat
  } else if(initial.selection=="dbscan"){
    out.initial = initialize_fuc.dbscan(data,K)
    memb = out.initial$memb
    L.mat = matrix(0,n_all,K)
    for(jj in 1:n_all) L.mat[jj, memb[jj]]=1
    out.initial$L.mat = L.mat
  }
  else {out.initial = initialize}
  
  if(L == 3){
    l=1
    aBIC = rep(0,l)
    Mu_hat.list = list()
    Theta_hat.list = list()
    prob.list = list()
    L.mat.list = list()
    member.list = list()
    lam1 = lambda1;lam2 = lambda2;lam3 = lambda3;
    PP = FGGM.refit(data, K, lam1, lam2, lam3, initialization=F, initialize=out.initial, average=average, asymmetric=asymmetric, local_appro=local_appro)
    mu_hat=PP$mu;Theta_hat=PP$Xi;L.mat = PP$L.mat0;group = PP$group;prob = PP$prob0;aBIC[l] = PP$bic; member = PP$member
    Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat; prob.list[[l]]=prob; L.mat.list[[l]] = L.mat; member.list[[l]]=member
  } else {
    Mu_hat.list = list()
    Theta_hat.list = list()
    prob.list = list()
    L.mat.list = list()
    member.list = list()
    # search lam3
    lam1 = median(lambda1);lam2 = median(lambda2)
    for (l in 1:L3) {
      lam3 = lambda3[l]
      PP = FGGM.refit(data, K, lam1, lam2, lam3, initialization=F, initialize=out.initial, average=average, asymmetric=asymmetric, local_appro=local_appro)
      mu_hat=PP$mu;Theta_hat=PP$Xi;L.mat = PP$L.mat0;group = PP$group;prob = PP$prob0;aBIC[l] = PP$bic; member = PP$member
      Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat; prob.list[[l]]=prob; L.mat.list[[l]] = L.mat; member.list[[l]]=member
      if(trace){
        print(c(cat(paste(l,"lam1 lam2 lam3 ="),c(round(lam1,2),round(lam2,2),round(lam3,2)),":"),paste("K =",as.numeric(dim(Theta_hat)[3]))))
      }
    }
    aBIC[is.na(aBIC)] = 10^10
    aBIC[aBIC==0] = abs(aBIC[aBIC!=0][1])*10
    n_lam3 = which(aBIC[1:L3] == min(aBIC[1:L3]))[1];lam3 = lambda3[n_lam3]
    
    # search lam2
    for (l2 in 1:L2) {
      lam2 = lambda2[l2];l = L3+l2
      PP = FGGM.refit(data, K, lam1, lam2, lam3, initialization=F, initialize=out.initial, average=average, asymmetric=asymmetric, local_appro=local_appro)
      mu_hat=PP$mu;Theta_hat=PP$Xi;L.mat = PP$L.mat0;group = PP$group;prob = PP$prob0;aBIC[l] = PP$bic; member = PP$member
      Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat; prob.list[[l]]=prob; L.mat.list[[l]] = L.mat; member.list[[l]]=member
      if(trace){
        print(c(cat(paste(l,"lam1 lam2 lam3 ="),c(round(lam1,2),round(lam2,2),round(lam3,2)),":"),paste("K =",as.numeric(dim(Theta_hat)[3]))))
      }
    }
    aBIC[is.na(aBIC)] = 10^10
    n_lam2 = which(aBIC[(L3+1):(L3+L2)] == min(aBIC[(L3+1):(L3+L2)]))[1];lam2 = lambda2[n_lam2]
    # search lam1
    for (l1 in 1:L1) {
      lam1 = lambda1[l1];l = L3+L2+l1
      PP = FGGM.refit(data, K, lam1, lam2, lam3, initialization=F, initialize=out.initial, average=average, asymmetric=asymmetric, local_appro=local_appro)
      mu_hat=PP$mu;Theta_hat=PP$Xi;L.mat = PP$L.mat0;group = PP$group;prob = PP$prob0;aBIC[l] = PP$bic; member = PP$member
      Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat; prob.list[[l]]=prob; L.mat.list[[l]] = L.mat; member.list[[l]]=member
      if(trace){
        print(c(cat(paste(l,"lam1 lam2 lam3 ="),c(round(lam1,2),round(lam2,2),round(lam3,2)),":"),paste("K =",as.numeric(dim(Theta_hat)[3]))))
      }
    }
    aBIC[is.na(aBIC)] = 10^10
    aBIC[aBIC==0] = abs(aBIC[aBIC!=0][1])*10
    n_lam1 = which(aBIC[(L3+L2+1):(L3+L2+L1)] == min(aBIC[(L3+L2+1):(L3+L2+L1)]))[1];lam1 = lambda1[n_lam1]
  }
  K.list <- rep(1,length(Theta_hat.list))
  for (l in 1:length(Theta_hat.list)) {
    K.list[l] <- as.numeric(dim(Theta_hat.list[[l]])[3])
  }
  aBIC[which(K.list == 1)] <- 10^10
  
  n_lam = which(aBIC == min(aBIC))[1]
  Opt_aBIC = min(aBIC)
  Opt_lambda = c(lam1,lam2,lam3)
  result = list(Opt_lambda=Opt_lambda,Mu_hat.list=Mu_hat.list,Theta_hat.list=Theta_hat.list,prob.list=prob.list,member.list=member.list,L.mat.list=L.mat.list,Opt_aBIC=Opt_aBIC,BIC=aBIC,Opt_num=n_lam)
  return(result)
}

############################# Some fundamental supporting functions ############################
mcp_d = function(x,lambda,a=3){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: mcp_d
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the derivative of the MCP
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ x: a float value or a vector, the independent variable in the MCP.
  ## @ lambda: a float value, the tuning parameter in the MCP.
  ## @ a: a float value, regularization parameter in the MCP, the default setting is 3.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ rho: the derivative of the MCP.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  if(lambda!=0){
    rho <- lambda*( 1 > abs(x)/( lambda*a ) )*( 1 - abs(x)/( lambda*a ))
  } else{
    rho=0
  }
  return(rho)
}

S_soft = function(z,lambda){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: S_soft
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Define the soft threshold operator.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ z: a float value or a vector, the independent variable.
  ## @ lambda: a float value, the tuning parameter.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ The result of the soft threshold operator.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  return((abs(z) - lambda)*(abs(z) - lambda > 0)*sign(z))
}

MCP_soft = function(z,lambda,a=3){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: MCP_soft
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Define the MCP threshold operator.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ z: a float value or a vector, the independent variable in the MCP.
  ## @ lambda: a float value, the tuning parameter in the MCP.
  ## @ a: a float value, regularization parameter in the MCP, the default setting is 3.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ The result of the MCP threshold operator.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  return( S_soft(z,lambda)/(1-1/a) * (abs(z) - a*lambda <= 0) + z * (abs(z) - a*lambda > 0) )
}

f.den.vec = function(data, mu, Theta){ 
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: f.den.vec
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            calculate the density function values at each sample point.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ mu1: p * 1 vector, the mean vector.
  ## @ Omega1: p * p matrix, the precision matrix.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ fdensity: The density function values at each sample point.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  p = length(mu)
  fden = as.numeric( (2*pi)^(-p/2) * (det(Theta))^(1/2) * exp(-1/2*diag(t(t(data) - as.numeric(mu)) %*% Theta %*% (t(data) - as.numeric(mu)))) )
  return(fden)
}

Symmetrize = function(X){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: Symmetrize
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Symmetrize the precision matrices using the symmetrization strategy of Cai et al. (2016).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ X: p * p matrix.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ X: a symmetrical matrix.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  p = dim(X)[1]
  for(i in 1:p){
    for(j in i:p){
      if(X[i,j] < X[j, i]){
        X[j, i] = X[i, j]
      }else{
        X[i, j] = X[j, i]
      }
    }
  }
  return(X)
}

cut_diff_ama = function(V_kk,K_c,K,cutoff=0.01){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: cut_diff_ama
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Truncation of small differences in parameter estimates between subgroups.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ V_kk: L2-norm differences in the precision matrices between different subgroups.
  ## @ K_c: All combinations of natural numbers within K.
  ## @ K: int, a selected upper bound of K_0.
  ## @ cutoff: a float value, a given cut-off value.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ K_group_final: the partition of the original K subgroups.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  V_kk_num = which(V_kk < cutoff)
  K_group_final = list()
  if(length(V_kk_num) > 0){
    K_group = list()
    for (j in 1:length(V_kk_num)) {
      K_group[[j]] = K_c[,V_kk_num[j]]
    }
    outnum = setdiff(1:K,Reduce(union,K_group))
    if(length(outnum) > 0){
      for (j in 1:length(outnum)) {
        K_group[[length(V_kk_num)+j]] = outnum[j]
      }
      K_group[[length(V_kk_num)+j+1]] = K
    } else{
      K_group[[length(V_kk_num)+1]] = K
    }
    
    kk = 1
    repeat{
      repeat{
        K_group_old=K_group
        k_del = NULL
        for (kkk in setdiff(1:length(K_group),1) ) {
          if(length(Reduce(intersect,list(K_group[[1]],K_group_old[[kkk]]))) > 0){
            K_group[[1]] = sort(unique(c(K_group[[1]],K_group_old[[kkk]])))
            k_del = c(k_del,kkk)
          }
        }
        if(length(k_del) > 0){
          for (j in sort(k_del,decreasing = T)) {
            K_group[[j]] = NULL
          }
        }
        if(length(K_group_old) == length(K_group)){break}
      }
      K_group_final[[kk]] = K_group[[1]]
      if(kk==1 && length(K_group) == 1){print("Warning: Only one cluster!");break}
      if(length(K_group) == 2){K_group_final[[kk+1]] = K_group[[2]];break}
      if(kk>1 && length(K_group) == 1){break}
      K_group[[1]] = NULL
      kk = kk+1
    }
  }else {
    for (k in 1:K) {
      K_group_final[[k]] = k
    }
  }
  
  return(K_group_final)
}


############################# Functions for generating the initial values ############################
initialize_fuc = function(data, K, n.start = 100){  
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: initialize_fuc
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the initial values using K-means.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  Mu <- matrix(0, K, p)
  kmeans.clust <- kmeans(data, K, nstart = n.start)
  memb <- kmeans.clust$cluster
  prob <- kmeans.clust$size/n
  Theta <- array(0, dim = c(p, p, K))
  S <- array(0, dim = c(p, p, K))
  for(k in 1:K)
  {
    Mu[k,] <- t(colMeans(data[memb == k, , drop = FALSE]) )
    S[,,k]  <- cov(data[memb == k, , drop = FALSE])
    Theta[,,k] <- solve(S[,,k])
  }
  
  int.res <- list()
  int.res$prob <-  prob
  int.res$Mu <-  Mu
  int.res$Theta <- Theta
  int.res$S <- S
  int.res$memb <- memb
  return(int.res)
}

initialize_fuc.dbscan = function(data, K){  
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: initialize_fuc.dbscan
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the initial values using hierarchical clustering.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: 
  ##            R packages: NbClust
  ## -----------------------------------------------------------------------------------------------------------------
  
  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])
  
  # initialization via dbscan clustering## 
  Mu <- matrix(0, K, p)
  hc <- hclust(dist(data,method = "euclidean"),method = "ward.D2")
  memb <- cutree(hc,k=K)
  
  prob <- rep(0,K)
  Theta <- array(0, dim = c(p, p, K))
  S <- array(0, dim = c(p, p, K))
  for(k in 1:K)
  {
    Mu[k,] <- t(colMeans(data[memb == k, , drop = FALSE]) )
    S[,,k]  <- cov(data[memb == k, , drop = FALSE])
    Theta[,,k] <- solve(S[,,k] + diag(1,p))
    prob[k] <- sum(memb == k)/n
  }
  
  P <- list()
  P$prob <-  prob
  P$Mu <-  Mu
  P$Theta <- Theta
  P$S <- S
  P$memb <- memb
  return(P)
}
