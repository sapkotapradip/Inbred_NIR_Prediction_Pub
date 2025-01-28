################################### Sapkota, Pradip ####################################################
################################# NIRBLUP Pradip 2025 ##################################################
################################# Agronomic traits #####################################################
################################# multi- environment ###################################################

####################### Across-Environment Prediction ##################################################
########################################################################################################

library(dplyr)
library(prospectr)
#### reading pheno data for hybrids
pheno_19CS = read.csv("blues_19CS_hybrids_spectra.csv", check.names = FALSE)
pheno_19TA = read.csv("blues_19TA_hybrids_spectra.csv", check.names = FALSE)
pheno_20CS = read.csv("blues_20CS_hybrids_spectra.csv", check.names = FALSE)
pheno_20MA = read.csv("blues_20MA_hybrids_spectra.csv", check.names = FALSE)

#### reading pheno data for inbreds
parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
parents_20MA = read.csv("blues_20MA_inbreds_spectra.csv", check.names = FALSE)

pheno_combined = rbind(pheno_19CS,pheno_19TA, pheno_20CS, pheno_20MA)
ne <- as.vector(table(pheno_combined$env)) ## counting the number of observations

#### loading marker data 
T92I012 <- read.delim(file = "T92I012.txt", header = TRUE, sep = "\t", dec = ".")
parents = read.delim(file = "parents.txt", header = TRUE, sep = "\t", dec = ".")$pedigree

geno = T92I012 %>% filter(
     X1 %in% parents)

rownames(geno) = geno$X1  #filtered available parents
geno = geno[,-1]

Z = scale(geno,center=TRUE) 
G = tcrossprod(Z)/ncol(Z) #GRM using vanraden method

# Constructing marker data
tmp<-strsplit(pheno_combined$pedigree,"/")

P1<-rep(NA,length(tmp))
P2<-rep(NA,length(tmp))

for(i in 1:length(tmp))
{
  P1[i]=tmp[[i]][1]
  P2[i]=tmp[[i]][2]
}

P1<-as.factor(P1)
P2<-as.factor(P2)
cross<-as.factor(pheno_combined$pedigree)

Site=as.factor(pheno_combined$env)

ZE=model.matrix(~Site-1)

Z1=model.matrix(~P1-1)
dim(Z1)

Z2=model.matrix(~P2-1)
dim(Z2)

Z3=model.matrix(~cross-1)
dim(Z3)


K1=G[levels(P1), levels(P1)]
dim(K1)
rownames(K1)
K2=G[levels(P2), levels(P2)]
dim(K2)
rownames(K2)

load("K3.Rdata") # loading hybrid relationshp created via kronecker product


#### modeling effects for priors
K1star=Z1%*%K1%*%t(Z1) #female
K2star=Z2%*%K2%*%t(Z2) #male
K3star=Z3%*%K3%*%t(Z3) #female * male

ZEZEt=tcrossprod(ZE) #env

K4=K1star*ZEZEt #female *env
K5=K2star*ZEZEt #male*env
K6=K3star*ZEZEt #female*male*env

## preparing phenomic data

NIR_19CS = scale(savitzkyGolay(pheno_19CS[,-c(1:7)], m=1, p=1, w=11))
NIR_19TA = scale(savitzkyGolay(pheno_19TA[,-c(1:7)], m=1, p=1, w=11))
NIR_20CS = scale(savitzkyGolay(pheno_20CS[,-c(1:7)], m=1, p=1, w=11))
NIR_20MA = scale(savitzkyGolay(pheno_20MA[,-c(1:7)], m=1, p=1, w=11))

rownames(NIR_19CS) = pheno_19CS$pedigree
rownames(NIR_19TA) = pheno_19TA$pedigree
rownames(NIR_20CS) = pheno_20CS$pedigree
rownames(NIR_20MA) = pheno_20MA$pedigree

NIR.d1 = rbind(NIR_19CS, NIR_19TA, NIR_20CS, NIR_20MA) 
ZN1 = tcrossprod(as.matrix(NIR.d1)/ncol(as.matrix(NIR.d1))) #phenomic relationship matrices
dim(ZN1) #phenomic relationship matrices from hybrids




