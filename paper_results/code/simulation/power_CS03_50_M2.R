# Power plots, 9 methods

library(ggplot2)
library(reshape2)

# other methods
library(MTAR)
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/compare_methods/aMAT.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/compare_methods/metaUSAT.R")

# mtafS
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/mtaf_davies_cauchy.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/cauchy_combination.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/simulation.R")
#source("/fs/project/PAS1149/qldeng/MTAFSC/functions/mtafS.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/mtafS_ratio.R")

# load real data
#load("/fs/project/PAS1149/qldeng/MTAFSC/data/T1_FAST_ROIs/corr_t1_fast_rois.rda")
#load("~/PAS1149/MTAFSC/data/Freesurfer_volume/corr_freesurfer_volume.rda")
#load("~/PAS1149/MTAFSC/data/thickness/corr_thickness.rda")
corr_z = cs_cor(50,0.3)


# set parameters
K=k=ncol(corr_z)
B=1e3
R = corr_z

# estimate true R by sample correlation
z_samples <- mvtnorm::rmvnorm(1e5, mean = rep(0,K), sigma = R)
R_hat <- cor(z_samples)
sigma_hat <- cov(z_samples)

r_eig <- eigen(sigma_hat)

# SUM: sum_denom
a <- rep(1, K)
sum_denom <- (sqrt(as.numeric(t(a) %*% R_hat %*% (a))))

# SSU:
cr <- eigen(R_hat, only.values = TRUE)$values
## approximate the distri by alpha Chisq_d + beta:
alpha1 <- sum(cr * cr * cr) / sum(cr * cr)
beta1 <- sum(cr) - (sum(cr * cr)^2) / (sum(cr * cr * cr))
d1 <- (sum(cr * cr)^3) / (sum(cr * cr * cr)^2)
alpha1 <- as.double(alpha1)
beta1 <- as.double(beta1)
d1 <- as.double(d1)

# S_HOM
Wi <- matrix(rep(1e3,K), nrow = 1)
sumW <- sqrt(sum(Wi^2))
W <- Wi / sumW
Sigma <- ginv(R_hat)


# T=50
#seq(2,6,length.out=10) 0.3 4%
#seq(1,3,length.out=10) 0.3 20%
#seq(2,4,length.out=10) 0.7 4%
#seq(1,2,length.out=10) 0.7 20%

# T=100
#seq(4,7,length.out=10) 0.3 2%
#seq(1.5,2.5,length.out=10) 0.3 20%
#seq(2.5,4.5,length.out=10) 0.7 2%
#seq(1,1.5,length.out=10) 0.7 20%

for (sp in c(0.04,0.2)) {
  print(sp)
  
  if (sp == 0.04) {
    w = seq(2,6,length.out=10)
  } else if (sp == 0.2) {
    w = seq(1,3,length.out=10)
  } 
  
  
  power_plot <- matrix(NA,nrow = 10,ncol = 9) # 10 effect sizes, 9 methods
  l=1
  
  for (u in w) {
    
    mu <- rep(0,K)
    mu[1:(sp*k)] <- u
    
    # generate data
    Z <- mvtnorm::rmvnorm(B, mean = mu, sigma = R)
    
    
    # metaUSAT and metaMANOVA
    #p_metaUSAT <- matrix(0,nrow = B,ncol = 2)
    p_metaUSAT <- matrix(as.numeric(unlist(apply(Z, 1, metausat,R=R_hat,metamanova=TRUE))),ncol=7,byrow=TRUE)[,c(2,5)]
    colnames(p_metaUSAT) <- c("p.metamanova","p.metausat")
    colMeans(p_metaUSAT < 1e-5, na.rm = TRUE)
    
    # minP
    #p_minP <- apply(Z, 1, UminPd, CovS=sigma_hat)
    
    # SUM
    p_sum <- Sum_fast_m(Z,sum_denom)
    mean(p_sum < 1e-5)
    
    # SSU
    p_ssu <- SumSqU_fast_m(Z,beta1,alpha1,d1)
    mean(p_ssu < 1e-5)  
    
    
    # S_HOM
    x1 <- matrix(Z, ncol = K)
    T_hom <- W %*% Sigma %*% t(x1)
    T_hom <-  as.vector(T_hom ** 2) / as.numeric(W %*% Sigma %*% t(W))
    p_Shom <- pchisq(T_hom, df = 1, ncp = 0, lower.tail = F)
    mean(p_Shom < 1e-5)
    
    # Cauchy
    p_cauchy <- as.vector(cct(2*pnorm(abs(Z),lower.tail = FALSE),weights = NULL))
    mean(p_cauchy < 1e-5)
    
    # aMAT
    p_amat <- aMAT(Z,R_hat)[,5]
    mean(p_amat < 1e-5)
    
    # MTAR
    p_mtar <- length(as.numeric(Lemats(Z,R_hat,alpha=5e-8)$idA))/B
    p_mtar
    
    # mtafS
    p_mtafS <- mtafS_ratio(Z,r_eig)
    colnames(p_mtafS) <- c("r1","r2","r3","r4","r5","MTAFSO")
    colMeans(p_mtafS < 1e-5)
    
    
    
    p_results <- cbind(p_mtafS[,"MTAFSO"],p_metaUSAT,p_amat,p_cauchy,p_sum,p_ssu,p_Shom)
    p_results <- c(colMeans(p_results < 5e-8, na.rm = TRUE),p_mtar)
    names(p_results) <- c("MTAFS","metaMANOVA","metaUSAT","aMAT","Cauchy","SUM","SSU","HOM","MTAR")
    
    power_plot[l,] <- p_results
    colnames(power_plot) <- names(p_results)
    
    l=l+1
  }
  
  fname <- paste0("/fs/project/PAS1149/qldeng/MTAFSC/results/data_plots/power_EJHG/ep2/power_cs03_50_M2_",sp,".rds")
  saveRDS(power_plot, file = fname)
  
}

