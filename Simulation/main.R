################################################################################################
#Codes for conducting simulation studies
################################################################################################
rm(list = ls(all = TRUE))
ls()
###############################

library(MASS)
library(NbClust)
library(JGL)
library(Matrix)
source("function.R")
source("sim-func.R")
n <- 200              # The sample size of each subgroup
p <- 20               # The dimension of the precision matrix
K0 <- 3               # The true number of subgroups
N <- rep(n,K0)        # The sample sizes of K0 subgroups
K <- 6                # The given upper bound of K0.


################ The true parameters ################
mue <- 1.5
nonnum <- 4
mu01 <- c(rep(mue,nonnum),rep(-mue,nonnum),rep(0,p-2*nonnum))
mu02 <- c(rep(mue,2*nonnum),rep(0,p-2*nonnum))
mu03 <- c(rep(-mue,2*nonnum),rep(0,p-2*nonnum))
# Power law network
set.seed(2)
A.list <- Power.law.network(p,s=5,I2=c(1),I3=c(2))
Theta01 <- A.list$A1
Theta02 <- A.list$A2
Theta03 <- A.list$A3
sigma01 <- solve(Theta01)
sigma02 <- solve(Theta02)
sigma03 <- solve(Theta03)
Mu0.list <- list(mu01,mu02,mu03)
Sigma0.list <- list(sigma01,sigma02,sigma03)
Theta0.list <- list(Theta01,Theta02,Theta03)

################ Generating simulated data ################
whole.data <- generate.data(N,Mu0.list,Theta0.list,Sigma0.list)


################ The implementation and evaluation ################
# ------------- The propose method --------------
lambda <- genelambda.obo(nlambda1=5,lambda1_max=0.5,lambda1_min=0.1,
                         nlambda2=15,lambda2_max=1.5,lambda2_min=0.1,
                         nlambda3=10,lambda3_max=3.5,lambda3_min=0.5)
t1 <- proc.time()
res <- tuning.lambda.FGGM(lambda, whole.data$data, K, initial.selection="K-means")
proc.time() - t1
Theta_hat.list <- res$Theta_hat.list
Mu_hat.list <- res$Mu_hat.list
prob.list <- res$prob.list
L.mat.list <- res$L.mat.list
opt_num <- res$Opt_num
opt_Theta_hat <- Theta_hat.list[[opt_num]]
opt_Mu_hat <- Mu_hat.list[[opt_num]]
opt_L.mat <- L.mat.list[[opt_num]]
opt_prob <- prob.list[[opt_num]]
K_hat <- dim(opt_Theta_hat)[3]
index.FGGM <- Esti.error(whole.data$data, opt_Mu_hat, opt_Theta_hat, whole.data$Mu0, whole.data$Theta0, K0, opt_L.mat, whole.data$L0, opt_prob)
index.FGGM

