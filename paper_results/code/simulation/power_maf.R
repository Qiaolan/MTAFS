# EJHG: evaluate power give mafs


# batch index
b = as.numeric(commandArgs("TRUE"))
print(paste0("batch ",b))


#
# mtafS
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/mtaf_davies_cauchy.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/cauchy_combination.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/simulation.R")
source("/fs/project/PAS1149/qldeng/MTAFSC/functions/mtafS_ratio.R")

# simulate 
#set.seed(1)

# set parameters
m = 1e5 # number of SNPs
N = 1000 # number of individuals
k = 58 # number of traits


# save sumstats
sumstat <- matrix(NA, nrow = m, ncol = k)

# 1. generate genotypes
maf <- runif(m, 1e-3, 0.5)
G <- sapply(maf, rbinom, n=N, size=2)


sigma_hat <- readRDS("~/PAS1149/MTAFSC/data/ejhg_volume_R4/sample_cov/cov_volume.rds")


# loop over each SNP

for (i in 1:m) {
  #print(i)
  p <- maf[i]
  x <- G[,i]

  # 1. effect with 20% heritability
  eff <- sqrt(diag(sigma_hat)/(198*p*(1-p)))
  eff[1:55] <- 0  # sparse


  # 2. generate Y
  Y <-  x %o% eff + mvtnorm::rmvnorm(N, rep(0,k), sigma_hat)
  colnames(Y) <- paste0("Y",1:58)

  # 3. get summary statistics
  assoc <- summary(lm(Y~x))
  z <- rep(NA, k)
  for (j in 1:k) {
    z[j] <- assoc[[paste0("Response Y",j)]][["coefficients"]][6]
  }
  
  sumstat[i,] <- z

}


# save sumstats
Z1 <- as.matrix(cbind(maf, sumstat))
fname <- paste0("~/PAS1149/MTAFSC/results/power_maf/ukcor1_m2_h1/sumstats/z_",b,".rds")
saveRDS(Z1, fname)
#rm(Z1)



#
Z1 <- na.omit(Z1)
maf <- Z1[,1]
Z1 <- Z1[,-1]

# 4. MTAFS
print("MTAFS is processing")

r_eig <- eigen(sigma_hat)
p_mtafS <- mtafS_ratio(Z1,r_eig)
colnames(p_mtafS) <- c("first2","q1","q2","q3","all","MTAFS")

sumstat <- as.data.frame(cbind(maf, p_mtafS[,6]))
colnames(sumstat) <- c("maf","MTAFS")

fname <- paste0("~/PAS1149/MTAFSC/results/power_maf/ukcor1_m2_h1/mtafs/mtafs_",b,".rds")
saveRDS(sumstat, fname)

