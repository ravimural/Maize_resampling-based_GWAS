library(rMVP)
library(bestNormalize)
library(plyr)
library(doparallel)

args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
n <- sub("-","",myargument)


print(n)


############MVP.Data(fileBed = "/scratch/pruneGeno", fileKin = T, filePC = T, out = "/scratch/rMVP2", priority = "speed", ncpus = 10) #, pcs.keep = 5
#############print("pc and kinship done")

MVP.Data(fileBed = "/scratch/geno", fileKin = F, filePC = F, out = "/scratch/rMVP", priority = "speed", ncpus = 10)
print("genotype matrix done")

genoList <- read.table("/scratch/rMVP.geno.ind", header = F)
colnames(genoList) <- "ID"

pheno <- read.table(paste("pheno_", n, ".txt", sep = ""), header = T)


ph <- plyr::join(genoList, pheno, by="ID")

#########################################

geno <- attach.big.matrix("/scratch/rMVP.geno.desc")
#K <- attach.big.matrix("/scratch/rMVP2.kin.desc")
#pc <- bigmemory::as.matrix(attach.big.matrix("/scratch/rMVP2.pc.desc"))
map <- data.table::fread("/scratch/rMVP.geno.map", data.table = F)

for (k in (2:ncol(ph))){
trait <- ph[,c(1, k)]
print(colnames(trait)[2])
m <- nrow(ph)*0.2
set.seed(44)
for (i in 3:103){ #3:53 and 54 to 103
  z <- sample(1:nrow(ph), m)
  trait[,i] <- trait[,2]
  trait[z,i] <- NA
}
RMIP <- c()
for(j in 3:ncol(trait)){
  imMVP <- MVP(
    phe=trait[, c(1, j)],
    geno=geno,
    map=map,
    #CV.FarmCPU=pc,
    nPC.FarmCPU=3,
    priority="speed",
    #K = K,
    ncpus=12,
    maxLoop=10,
    method.bin="FaST-LMM",
    method=c("FarmCPU"), 
    file.output = F, 
    p.threshold = 0.05/nrow(map)
  )
  farm <- cbind(imMVP$map, imMVP$farmcpu.results)
  farm <- na.omit(farm[farm[,8]<(0.05/nrow(map)),])
  print(dim(farm)[1])
  colnames(farm)[8] <- "pvalue"
  RMIP <- rbind(RMIP, farm)
  rm(farm, imMVP)
  gc()
}


########### Save RMIP ########### 
sigSNps <- table(RMIP$SNP) # count SNP
sigSNps2 <- as.data.frame(sigSNps) 
colnames(sigSNps2)[1] <- 'SNP'
RMIPsigSNP <- plyr::join(RMIP, sigSNps2, by="SNP")

RMIPsigSNP <- RMIPsigSNP[order(RMIPsigSNP$SNP, RMIPsigSNP$pvalue, decreasing=TRUE),] # FALSE = sort in ascending/increasing order to keep SNP with lowest pvalue # TRUE = sort in descending or decreasing order tio keep SNP with highest pvalue
RMIPsigSNP <- RMIPsigSNP[!duplicated(RMIPsigSNP$SNP),] # remove duplicates except the first one which should be lowest (FALSE) or highest (TRUE) pvalue as sorted in above step
RMIPsigSNP <- RMIPsigSNP[order(RMIPsigSNP$SNP, RMIPsigSNP$pvalue),] # Sort again so it looks better organized

write.csv(RMIPsigSNP, file = paste("/scratch/", colnames(ph)[2], ".csv", sep = ""), row.names = F, quote = F)
}