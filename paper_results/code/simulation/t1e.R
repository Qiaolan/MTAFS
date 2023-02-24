# Type 1 error, 9 methods

# other methods
library(emulator)
library(MTAR)
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/compare_methods/aMAT.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/compare_methods/metaUSAT.R")

# mtafS
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/mtaf_davies_cauchy.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/cauchy_combination.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/simulation.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/mtafS_ratio.R")


# load real data
#load("~/PAS1149/MTAFSC/data/Freesurfer_volume/corr_freesurfer_volume.rda")
#load("/fs/project/PAS1149/qldeng/MTAFSC/data/T1_FAST_ROIs/sample_corr/corr_t1_fast_rois.rda")
#corr_z = cs_cor(50,0.7)
#corr_z = ar1_cor(50,0.3)
#corr_z = diag(100)

# batches
batch <- commandArgs("TRUE")
print(paste0("batch ",batch))
#set.seed(as.integer(batch))

# set parameters
K=k=ncol(corr_z)
B=1e6
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

# inverse
R_inv <- mpinv(sigma_hat)

# generate data
Z <- mvtnorm::rmvnorm(B, mean = rep(0,K), sigma = R)


# metaUSAT and metaMANOVA
print("metaUSAT")
p_metaUSAT <- matrix(as.numeric(unlist(apply(Z, 1, metausat,R=R_hat,metamanova=TRUE, AbsTol=1e-10))),ncol=7,byrow=TRUE)[,c(2,5)]
colnames(p_metaUSAT) <- c("p.metamanova","p.metausat")


# SUM
p_sum <- Sum_fast_m(Z,sum_denom)


# SSU
p_ssu <- SumSqU_fast_m(Z,beta1,alpha1,d1)



# S_HOM
x1 <- matrix(Z, ncol = K)
T_hom <- W %*% Sigma %*% t(x1)
T_hom <-  as.vector(T_hom ** 2) / as.numeric(W %*% Sigma %*% t(W))
p_Shom <- pchisq(T_hom, df = 1, ncp = 0, lower.tail = F)


# Cauchy
p_cauchy <- as.vector(cct(2*pnorm(abs(Z),lower.tail = FALSE),weights = NULL))


# aMAT
p_amat <- aMAT(Z,R_hat)[,5]


# MTAR
p_mtar <- length(as.numeric(Lemats(Z,R_hat,alpha=1e-5)$idA))/B



# MTAFS
print("MTAFS")
p_mtafS <- mtafS_ratio(Z,r_eig)
colnames(p_mtafS) <- c("first2","q1","q2","q3","all","MTAFS")


p_results <- cbind(p_metaUSAT,p_sum,p_ssu,p_Shom,p_cauchy,p_amat,p_mtar,p_mtafS)


# save results
#save_path = paste0("~/PAS1149/MTAFSC/results/t1e_5e-8/UKCOR1/","t1e_UKCOR1_",batch,".rds")
#save_path = paste0("~/PAS1149/MTAFSC/results/t1e_5e-8/UKCOR2/","t1e_UKCOR2_",batch,".rds")
#save_path = paste0("~/PAS1149/MTAFSC/results/t1e_5e-8/CS07/","t1e_CS07_",batch,".rds")
save_path = paste0("~/PAS1149/MTAFSC/results/t1e_5e-8/AR03/","t1e_AR03_",batch,".rds")

saveRDS(p_results, file = save_path)
