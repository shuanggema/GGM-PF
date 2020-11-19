################################################################################################
#Codes for Heterogeneity analysis of regulatory T cells in non-small-cell lung cancer
################################################################################################
rm(list = ls(all = TRUE))
ls()
###############################

library(NbClust)
library(igraph)
source("function.R")

# ------------- load data ----------------
cell <- read.csv("Treg cell factor.csv")
data0 <- read.csv("Wnt pathway.csv")
data <- data0
n = dim(data)[1]
a = data[,-c(1:5)]
a.sd = apply(a,2,sd)
zero.num = apply(a,2,function(a){sum(a==0)})
a = a[,-union(as.numeric(which(a.sd < 0.25)),as.numeric(which(zero.num > 0.95*n)))]
data=a

# ------------- identifying the number of subgroups ----------------
K=10
lambda=genelambda.obo(nlambda1=5,lambda1_max=3,lambda1_min=0.1,
                      nlambda2=15,lambda2_max=2.5,lambda2_min=0.2,
                      nlambda3=10,lambda3_max=8,lambda3_min=4)
t1 <- proc.time()
set.seed(1)
res = tuning.lambda.FGGM(lambda, data, K, initial.selection="dbscan")
t10 <- proc.time() - t1
opt_num = res$Opt_num
Theta_hat.list = res$Theta_hat.list
opt_Theta_hat = Theta_hat.list[[opt_num]]
K = dim(opt_Theta_hat)[3]
# save.image("T-cell-FGGM-K10.RData")

# ------------- refitting can be conducted to improve empirical accuracy ----------------
lambda=genelambda.obo(nlambda1=8,lambda1_max=0.4,lambda1_min=0.15,
                      nlambda2=15,lambda2_max=0.4,lambda2_min=0.15,
                      nlambda3=1,lambda3_max=8,lambda3_min=4)
lambda$lambda3 = 0
t1 <- proc.time()
set.seed(1)
res = tuning.lambda.FGGM(lambda, data, K, initial.selection="dbscan")
t20 <- proc.time() - t1
opt_num = res$Opt_num
Theta_hat.list = res$Theta_hat.list
Mu_hat.list = res$Mu_hat.list
prob.list = res$prob.list
member.list = res$member
opt_Theta_hat = Theta_hat.list[[opt_num]]
opt_Mu_hat = Mu_hat.list[[opt_num]]
opt_member = member.list[[opt_num]]
K = dim(opt_Theta_hat)[3]
# save.image("T-cell-FGGM-K10_refit.RData")



#### The sample sizes of three identified subgroups
table(opt_member)
sum(table(opt_member))
#### The numbers of edges of three identified subgroups
p = dim(opt_Theta_hat)[2]
(sum(opt_Theta_hat[,,1]!=0) - p) / 2
(sum(opt_Theta_hat[,,2]!=0) - p) / 2
(sum(opt_Theta_hat[,,3]!=0) - p) / 2


# ----------------- The first part of Figure 1: identified subgrouping network structures ------------------
gene.name = names(data)
network.list <- list()
for (k in 1:K) {
  network = as.data.frame(opt_Theta_hat[,,k])
  names(network) = gene.name
  row.names(network) = gene.name
  network = as.matrix(network)
  network[which(network != 0)] <- 1
  network.list[[k]] <- network
}
network1 = network.list[[1]]
network2 = network.list[[2]]
network3 = network.list[[3]]
net.plot1 = graph_from_adjacency_matrix(network1,mode = "undirected",diag = FALSE)
net.plot2 = graph_from_adjacency_matrix(network2,mode = "undirected",diag = FALSE)
net.plot3 = graph_from_adjacency_matrix(network3,mode = "undirected",diag = FALSE)


# ------------------ finding differences of expressions of immune related genes among subgroups ----------------
member = opt_member
plot.data <- cell[6:dim(cell)[2]]
plot.data <- scale(plot.data)
plot.data <- as.data.frame(cbind(plot.data,as.factor(member)))
P.value <- rep(1,dim(plot.data)[2]-1)
for (i in 1:(dim(plot.data)[2]-1)) {
  y <- plot.data[,i]
  x <- as.factor(plot.data[,dim(plot.data)[2]])
  fit <- aov(y~x)
  fit <- summary(fit)
  P.value[i] <- fit[[1]][1,5]
}
alpha <- 0.05
sig.n <- c(which(P.value < alpha),match("IRF4",names(plot.data)))
plot.data.sig <- as.data.frame(plot.data[,c(sig.n,dim(plot.data)[2])])
names(plot.data.sig)[dim(plot.data.sig)[2]] <- "group"
n <- dim(data0)[1]
mean.list <- as.data.frame(matrix(0,K,dim(cell)[2]-5))
names(mean.list) <- names(cell)[6:dim(cell)[2]]
mu.list <- mean.list
sd.list <- mean.list
cell.list <- list()
for (k in 1:K) {
  cell.list[[k]] <- cell[match(data0[which(member == k),1],cell[,1]),6:dim(cell)[2]]
  mean.list[k,] <- paste(round(apply(cell.list[[k]],2,mean),2),round(apply(cell.list[[k]],2,sd),2),sep = "+")
  mu.list[k,] <- round(apply(cell.list[[k]],2,mean),2)
  sd.list[k,] <- round(apply(cell.list[[k]],2,sd),2)
}
mu.sig <- mu.list[,sig.n]
mu.sig = mu.sig[,order(mu.sig[2,],decreasing = T)]


#####################################################################################
# In the R package "igraph", the node positions are randomized.
# In order to fully reproduce the Figure 1 in the article,
# the following ".RData" file stores the node location parameters used.
load("FGGM_net.RData")
#####################################################################################

# ------------------ Figure 1 ----------------
pdf(file=paste0("Tregs.pdf"),width=10,height=10)
par(mfrow = c(2,2), mar = c(1, 1, 1, 1)) 
plot.igraph(net.plot1,vertex.size=2,vertex.label.cex=0.7,vertex.label.dist=0.75, layout = l, edge.width = 0.1)
plot.igraph(net.plot2,vertex.size=2,vertex.label.cex=0.7,vertex.label.dist=0.75, layout = l, edge.width = 0.1)
plot.igraph(net.plot3,vertex.size=2,vertex.label.cex=0.7,vertex.label.dist=0.75, layout = l, edge.width = 0.1)
par(las = 2, mar = c(5, 3, 1, 1))
barplot(as.matrix(mu.sig),ylab="Mean normalized gene expression",
        col=c("brown1","green2","blue1"),legend=row.names(mu.sig),beside = T,
        mgp=c(1, 0.25, 0),tck=0.01)
dev.off()


