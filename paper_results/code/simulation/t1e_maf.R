# EJHG: evaluate t1e inflation


# batch index
b = as.numeric(commandArgs("TRUE"))
print(paste0("batch ",b))


#
# mtafS
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/mtaf_davies_cauchy.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/cauchy_combination.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/simulation.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/mtafS_ratio.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/compare_methods/aMAT.R")
# simulate 
#set.seed(1)

# set parameters
m = 1e6 # number of SNPs
N = 1000 # number of individuals
k = 58 # number of traits


# save sumstats
sumstat <- matrix(NA, nrow = m, ncol = k)

# 1. generate genotypes
#maf <- runif(m, 1e-3, 1e-2)
#maf <- runif(m, 1e-2, 5e-2)
maf <- runif(m, 5e-2, 5e-1)
G <- sapply(maf, rbinom, n=N, size=2)


# 2. generate Y
corY <- readRDS("~/PAS1149/MTAFSC/data/Freesurfer_volume/cor_volume_ldsc.rds")

Y <- mvtnorm::rmvnorm(N, rep(0,k), corY)


# 3. get summary statistics
print("Association test starts")

for (i in 1:m) {
  if (i %% 1e4 ==0) {
    print(i/1e4) 
  }
  assoc <- summary(lm(Y~G[,i]))
  z <- rep(NA, k)
  for (j in 1:k) {
    z[j] <- assoc[[paste0("Response Y",j)]][["coefficients"]][6]
  }
  sumstat[i,] <- z
}

Z1 <- as.matrix(cbind(maf, sumstat))

fname <- paste0("~/PAS1149/MTAFSC/results/inflation/zscores/z_",b,".rds")
saveRDS(Z1, fname)
#rm(Z1)

Z1 <- na.omit(Z1)

maf <- Z1[,1]
Z1 <- Z1[,-1]

# 4. MTAFS
print("MTAFS is processing")

sigma_hat <- readRDS("~/PAS1149/MTAFSC/data/ejhg_volume_R4/sample_cov/cov_volume.rds")
r_eig <- eigen(sigma_hat)
p_mtafS <- mtafS_ratio(Z1,r_eig)
colnames(p_mtafS) <- c("first2","q1","q2","q3","all","MTAFS")


# 5. aMAT
p_amat <- aMAT(Z1,cov2cor(sigma_hat))[,5]



sumstat <- as.data.frame(cbind(maf, p_mtafS[,6], p_amat))
colnames(sumstat) <- c("maf","MTAFS", "aMAT")

fname <- paste0("~/PAS1149/MTAFSC/results/inflation/sumstats2/sumstats_",b,".rds")
saveRDS(sumstat, fname)

