# EJHG revision: real data analysis - Volumne

# batch index
b = as.numeric(commandArgs("TRUE"))
print(paste0("batch ",b))

# MTAFS functions
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/mtaf_davies_cauchy.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/cauchy_combination.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/mtafS_ratio.R")

# other methods
library(emulator)
library(MTAR)
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/compare_methods/aMAT.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/compare_methods/metaUSAT.R")


# Load LDSC estimated trait correlation matrix
sigma_hat <- readRDS("~/PAS1149/MTAFSC/data/ejhg_volume_R4/sample_cov/cov_volume.rds")


# setup
# eigen-decompostion
r_eig <- eigen(sigma_hat)
# inverse
R_inv <- mpinv(sigma_hat)

# set parameters
K=58


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


# batch info
if (b == 20) {
  n = 9500001:9971805
  bs = length(n) # batch size
} else {
  bs = 5e5 # batch size
  n = (1+(b-1)*bs):(b*bs)  
}

# snp info
pos_qc_all <- readRDS("~/PAS1149/MTAFSC/data/ejhg_volume/pos_qc_all.rds")
pos_qc_all <- pos_qc_all[n,]

# Z scores for one batch
Z_volume <- readRDS("~/PAS1149/MTAFSC/data/ejhg_volume_R4/sumstats/z_volume_all.rds")
Z_volume <- Z_volume[n,]
Z_volume <- as.matrix(Z_volume)




# MTAFS
print(c("Start MTAFS"))
p_mtafS <- mtafS_ratio(Z_volume, r_eig)
colnames(p_mtafS) <- c("first2","q1","q2","q3","all","MTAFS")
p_mtafS <- cbind(pos_qc_all, p_mtafS)

# save results
fname <- paste0("~/PAS1149/MTAFSC/data/ejhg_volume_R4/mtafs_batch/mtafs_batch_",b,".rds")
saveRDS(p_mtafS, file = fname)


# metaUSAT & metaMANOVA
print(c("Start metaUSAT"))
p_metaUSAT <- matrix(as.numeric(unlist(apply(Z_volume, 1, metausat,R=R_hat, metamanova=TRUE))),ncol=7,byrow=TRUE)[,c(2,5)]
colnames(p_metaUSAT) <- c("p.metamanova","p.metausat")


# SUM
p_sum <- Sum_fast_m(Z_volume,sum_denom)


# SSU
p_ssu <- SumSqU_fast_m(Z_volume,beta1,alpha1,d1)



# S_HOM
x1 <- matrix(Z_volume, ncol = K)
T_hom <- W %*% Sigma %*% t(x1)
T_hom <-  as.vector(T_hom ** 2) / as.numeric(W %*% Sigma %*% t(W))
p_Shom <- pchisq(T_hom, df = 1, ncp = 0, lower.tail = F)


# Cauchy
p_cauchy <- as.vector(cct(2*pnorm(abs(Z_volume),lower.tail = FALSE),weights = NULL))


# aMAT
p_amat <- aMAT(Z_volume, R_hat)[,5]


# result
result <- as.data.frame(cbind(p_metaUSAT,p_sum,p_ssu,p_Shom,p_cauchy,p_amat))
colnames(result) <- c("metaMANOVA","metaUSAT","SUM","SSU","HOM","Cauchy","aMAT")
result <- cbind(pos_qc_all, result)

fname <- paste0("~/PAS1149/MTAFSC/data/ejhg_volume_R4/competing_batch/competing_batch_",b,".rds")
saveRDS(result, file = fname)



# MTAR
p_mtar <- Lemats(Z_volume, R_hat,alpha=5e-8)
result_mtar <- pos_qc_all[p_mtar$idA,]

fname <- paste0("~/PAS1149/MTAFSC/data/ejhg_volume_R4/competing_batch/mtar_batch_",b,".rds")
saveRDS(result_mtar, file = fname)
