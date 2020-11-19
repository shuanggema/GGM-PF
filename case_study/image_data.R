################################################################################################
#Codes for LUAD heterogeneity analysis using histopathological imaging data
################################################################################################
rm(list = ls(all = TRUE))
ls()
###############################

library(HDtest)
library(igraph)
library(survival)
library(survminer)
library(Hmisc)
source("function.R")

# ------------- load data ----------------
clinical1 <- read.csv("data_bcr_clinical_data_patient_match.csv")
clinical2 <- read.csv("data_bcr_clinical_data_sample_match.csv")
data0 <- read.csv("Luad_feat_surv_data_307.csv",row.names= 1)
data <- data0
data <- data[,-c(222:229)]
n = dim(data)[1]
p = dim(data)[2]
data.sd = as.numeric(apply(data,2,sd))
data.mean = as.numeric(apply(data,2,mean))
data = data[,which(data.sd > 1)]
data = t(t(data) / as.numeric(apply(data,2,sd)))

# ------------- identifying the number of subgroups ----------------
K=10
lambda=genelambda.obo(nlambda1=5,lambda1_max=4,lambda1_min=0.75,
                      nlambda2=20,lambda2_max=5,lambda2_min=0.2,
                      nlambda3=20,lambda3_max=10,lambda3_min=5)
t1 <- proc.time()
set.seed(1)
res = tuning.lambda.FGGM(lambda, data, K, initial.selection="dbscan")
t10 <- proc.time() - t1
opt_num = res$Opt_num
Theta_hat.list = res$Theta_hat.list
opt_Theta_hat = Theta_hat.list[[opt_num]]
K = dim(opt_Theta_hat)[3]
# save.image("image-FGGM-K10.RData")

# ------------- refitting can be conducted to improve empirical accuracy ----------------
lambda=genelambda.obo(nlambda1=15,lambda1_max=5,lambda1_min=0.5,
                      nlambda2=29,lambda2_max=5,lambda2_min=0.25,
                      nlambda3=1,lambda3_max=10,lambda3_min=5)
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
# save.image("image-FGGM-K10_refit.RData")


#### The sample sizes of three identified subgroups
table(opt_member)
sum(table(opt_member))
#### The numbers of edges of three identified subgroups
p = dim(opt_Theta_hat)[2]
(sum(opt_Theta_hat[,,1]!=0) - p) / 2
(sum(opt_Theta_hat[,,2]!=0) - p) / 2
#### two subgroup differ significantly ( p-value from the test of Li and Chen (2012) )
X1 = as.matrix(data[which(opt_member == 1),])
X2 = as.matrix(data[which(opt_member == 2),])
equalCovs(X1, X2, 0.05, DNAME="0")




# ----------------- Figure 2 ------------------
data = as.data.frame(data)
gene.name = names(data)
k=1
network = as.data.frame(opt_Theta_hat[,,k])
names(network) = 1:dim(network)[1]
row.names(network) = 1:dim(network)[1]
network = as.matrix(network)
network[which(network != 0)] <- 1
network1 = network

k=2
network = as.data.frame(opt_Theta_hat[,,k])
names(network) = 1:dim(network)[1]
row.names(network) = 1:dim(network)[1]
network = as.matrix(network)
network[which(network != 0)] <- 1
network2 = network

net.plot1 = graph_from_adjacency_matrix(network1,mode = "undirected",diag = FALSE)
net.plot2 = graph_from_adjacency_matrix(network2,mode = "undirected",diag = FALSE)

network.joint = network1+network2
network.joint[which(network.joint != 0)] <- 1
set.seed(12)
net.plot = graph_from_adjacency_matrix(network.joint,mode = "undirected",diag = FALSE)
l = layout.fruchterman.reingold(net.plot)

pdf(file=paste0("image_","network",".pdf"),width=10,height=5)
par(mfrow = c(1,2), mar = c(0.1, 1, 0.1, 1)) 
plot.igraph(net.plot1,vertex.size=2,vertex.label.cex=0.7,vertex.label.dist=0.75, layout = l, edge.width = 0.1)
plot.igraph(net.plot2,vertex.size=2,vertex.label.cex=0.7,vertex.label.dist=0.75, layout = l, edge.width = 0.1)
dev.off()




# ----------------- Table 3 ------------------
## find differences
na.split1 = strsplit(clinical1[,1],split='-')
na.split2 = strsplit(clinical2[,1],split='-')
for (j in 1:length(na.split1)) {
  clinical1[j,1] = paste(na.split1[[j]][1],na.split1[[j]][2],na.split1[[j]][3],sep = ".")
}
for (j in 1:length(na.split2)) {
  clinical2[j,1] = paste(na.split2[[j]][1],na.split2[[j]][2],na.split2[[j]][3],sep = ".")
}
clinical1 = clinical1[match(rownames(data),clinical1[,1]),]
clinical2 = clinical2[match(rownames(data),clinical2[,1]),]
clinical3 = data0[,c(222:226)]
row.names(clinical1) = clinical1[,1]
row.names(clinical2) = clinical2[,1]
clinical1 = clinical1[,-1]
clinical2 = clinical2[,-1]
clinical = as.data.frame(cbind(clinical1,clinical2,clinical3))
clinical = clinical[,unique(names(clinical))]
t.c.n = c()
for (j in 1:dim(clinical)[2]) {
  t.c = as.data.frame((table(clinical[,j])))
  t.c.n[j] = dim(t.c)[1]
}
clinical = clinical[,which(t.c.n > 1)]
t.c.n = t.c.n[t.c.n!=1]
names(clinical)[which(t.c.n > 10)]
a1 = clinical[,which(t.c.n > 10)]
name.cont = names(clinical)[which(t.c.n > 10)]
clinical.class = clinical[,which(t.c.n <= 10)]
clinical.cont = clinical[,which(t.c.n > 10)]


#### Tumor status
namej = "TUMOR_STATUS"
namejj = "TUMOR STATUS"
data.tab = as.data.frame(cbind(opt_member,clinical.class[,namej]))
data.tab = data.tab[(as.character(data.tab[,2]) != "[Not Available]"),]
names(data.tab) = c("memb","x")
c = t(xtabs(~memb+x,data=data.tab))
c = c[-1,]
incidence = round(c[2,] / apply(c, 2, sum),3)
tab = as.data.frame(rbind(c,incidence))
tab[-dim(tab)[1],] = apply(tab[-dim(tab)[1],],2,function(a){return(as.character(a))})
cla = tolower(row.names(tab))
rowna = c(paste0(capitalize(tolower(namejj)),"(n=",sum(c),")"),rep("",dim(tab)[1]-1))
c.sum = fisher.test(c)
P.vec = c(rep("",dim(tab)[1]-1),as.character(signif(c.sum$p.value, 3)))
sum.tab1 = cbind(rowna,cla,tab,P.vec)
sum.tab1


#### Kras indicator
namej = "KRAS_GENE_ANALYSIS_INDICATOR"
namejj = "KRAS INDICATOR"
data.tab = as.data.frame(cbind(opt_member,clinical.class[,namej]))
data.tab = data.tab[(as.character(data.tab[,2]) != "[Not Available]"),]
names(data.tab) = c("memb","x")
c = t(xtabs(~memb+x,data=data.tab))
incidence = round(c[2,] / apply(c, 2, sum),3)
tab = as.data.frame(rbind(c,incidence))
tab[-dim(tab)[1],] = apply(tab[-dim(tab)[1],],2,function(a){return(as.character(a))})
cla = c("negative","positive","positive rate")
rowna = c(paste0(capitalize(tolower(namejj)),"(n=",sum(c),")"),rep("",dim(tab)[1]-1))
c.sum = fisher.test(c)
P.vec = c(rep("",dim(tab)[1]-1),as.character(signif(c.sum$p.value, 3)))
sum.tab2 = cbind(rowna,cla,tab,P.vec)
sum.tab2


#### ICD 10
namej = "ICD_10"
namejj = "ICD 10"
data.tab = as.data.frame(cbind(opt_member,clinical.class[,namej]))
data.tab = data.tab[(as.character(data.tab[,2]) != "[Not Available]"),]
names(data.tab) = c("memb","x")
c = t(xtabs(~memb+x,data=data.tab))
others = apply(c[-c(2:4),],2,sum)
c = as.data.frame(rbind(c[2:4,],others))
rowna = c(paste0(namejj,"(n=",sum(c),")"),rep("",dim(c)[1]-1))
cla = row.names(c)
c.sum = fisher.test(c)
P.vec = c(rep("",dim(c)[1]-1),as.character(signif(c.sum$p.value, 3)))
sum.tab3 = cbind(rowna,cla,c,P.vec)
sum.tab3

#### Overall survival status
namej = c("OS_STATUS","OS_MONTHS")
namejj = "Overall Survival STATUS"
data.tab = as.data.frame(cbind(opt_member,clinical.class[,namej[1]],clinical.cont[,namej[2]]))
names(data.tab)[2:3] = c("status","time")
data.tab = data.tab[(as.character(data.tab[,2]) != "[Not Available]"),]
table(data.tab[,2])
data.tab[(data.tab[,2] == "DECEASED"),2] = "2"
data.tab[!(data.tab[,2] == "DECEASED"),2] = "1"
data.tab[,2] = as.numeric(data.tab[,2])
data.tab[,3] = as.numeric(data.tab[,3])
surv_diff <- survdiff(Surv(time, status) ~ opt_member, data = data.tab)
p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
namej = "OS_STATUS"
namejj = "Overall Survival STATUS"
data.tab = as.data.frame(cbind(opt_member,clinical.class[,namej]))
data.tab = data.tab[(as.character(data.tab[,2]) != "[Not Available]"),]
names(data.tab) = c("memb","x")
c = t(xtabs(~memb+x,data=data.tab))
tab = c
cla = tolower(row.names(tab))
P.vec = c(rep("",dim(tab)[1]-1),as.character(signif(p.value, 3)))
rowna = c(paste0(capitalize(tolower(namejj)),"(n=",sum(c),")"),rep("",dim(tab)[1]-1))
sum.tab4 = as.data.frame(cbind(rowna,cla,tab,P.vec))
sum.tab4



#### Per(pre-broncholiator)
namej = "FEV1_PERCENT_REF_PREBRONCHOLIATOR"
namejj = "PER(PREBRONCHOLIATOR)"
data.tab = as.data.frame(cbind(opt_member,clinical.cont[,namej]))
data.tab = data.tab[(as.character(data.tab[,2]) != "[Not Available]"),]
names(data.tab) = c("memb","x")
fit <- aov(x~memb,data = data.tab)
fit <- summary(fit)
signif(fit[[1]][1,5],3)
x1 = as.numeric(data.tab$x[(data.tab$memb==1)])
x2 = as.numeric(data.tab$x[(data.tab$memb==2)])
sum.tab5 = c(paste0(capitalize(tolower(namejj)),"(n=",dim(data.tab)[1],")"),"",paste0(round(mean(x1),2),"(",round(sd(x1),2),")"),paste0(round(mean(x2),2),"(",round(sd(x2),2),")"),signif(fit[[1]][1,5],3))
sum.tab5


#### Per(post-broncholiator)
namej = "FEV1_PERCENT_REF_POSTBRONCHOLIATOR"
namejj = "PER(POSTBRONCHOLIATOR)"
data.tab = as.data.frame(cbind(opt_member,clinical.cont[,namej]))
data.tab = data.tab[(as.character(data.tab[,2]) != "[Not Available]"),]
names(data.tab) = c("memb","x")
fit <- aov(x~memb,data = data.tab)
fit <- summary(fit)
signif(fit[[1]][1,5],3)
x1 = as.numeric(data.tab$x[(data.tab$memb==1)])
x2 = as.numeric(data.tab$x[(data.tab$memb==2)])
sum.tab6 = c(paste0(capitalize(tolower(namejj)),"(n=",dim(data.tab)[1],")")," ",paste0(round(mean(x1),2),"(",round(sd(x1),2),")"),paste0(round(mean(x2),2),"(",round(sd(x2),2),")"),signif(fit[[1]][1,5],3))
sum.tab6



summary.stat = as.data.frame(rbind(sum.tab1,sum.tab2,sum.tab3,sum.tab4,sum.tab5,sum.tab6))
row.names(summary.stat) = NULL
summary.stat
write.csv(summary.stat,"summary.stat.csv")
