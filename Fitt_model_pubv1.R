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

rownames(pheno_19CS) = pheno_19CS$pedigree
rownames(pheno_19TA) = pheno_19TA$pedigree
rownames(pheno_20CS) = pheno_20CS$pedigree
rownames(pheno_20MA) = pheno_20MA$pedigree

NIR_19CS = scale(savitzkyGolay(pheno_19CS[,-c(1:7)], m=1, p=1, w=11))
NIR_19TA = scale(savitzkyGolay(pheno_19TA[,-c(1:7)], m=1, p=1, w=11))
NIR_20CS = scale(savitzkyGolay(pheno_20CS[,-c(1:7)], m=1, p=1, w=11))
NIR_20MA = scale(savitzkyGolay(pheno_20MA[,-c(1:7)], m=1, p=1, w=11))



NIR.d1 = rbind(NIR_19CS, NIR_19TA, NIR_20CS, NIR_20MA) 
ZN1 = tcrossprod(as.matrix(NIR.d1)/ncol(as.matrix(NIR.d1))) #phenomic relationship matrices
dim(ZN1) #phenomic relationship matrices from hybrids

## phenomic relationship matrices from parental NIRS
# use 20CS mid parent values for hybrids
NIR_pheno_parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
NIR_parents_20CS = NIR_pheno_parents_20CS[,-c(2:5)]
rownames(NIR_parents_20CS) = NIR_parents_20CS$pedigree
NIR_parents_20CS = NIR_parents_20CS[,-1]
F = NIR_parents_20CS[1:45,]
M = NIR_parents_20CS[46:89,]


## cakculating averages 

library(dplyr)
available.hybrid

averages = list()
# row_F =1; row_M =1
for (row_F in rownames(F)) {
  for (row_M in rownames(M)) {
    # Compute the average for the current pair of rows
    avg <- (F[row_F, ] + M[row_M, ]) / 2
    
    # Store the result in the list with a descriptive name
    averages[[paste(row_F, row_M, sep = "/")]] <- avg
  }
}

# Convert the list of averages into a data frame
averages_df <- do.call(rbind, averages)

# Ensure the row names of averages_df are accessible
averages_pedigrees <- rownames(averages_df)
pedigrees_list = pheno_combined$pedigree
filtered_df <- averages_df[rownames(averages_df) %in% pedigrees_list, ]
final_pedigree_df <- filtered_df[match(pedigrees_list, rownames(filtered_df)), ]

mid_parent_19CS = final_pedigree_df[1:364,]
mid_parent_19TA = final_pedigree_df[365:543,]
mid_parent_20CS = final_pedigree_df[544:907,]
mid_parent_20LY = final_pedigree_df[-c(1:907),]

## estimate first derivative of mid-parent NIR spectra

NIR_mid_parent_19CS = scale(savitzkyGolay(mid_parent_19CS, m=1, p=1, w=11))
NIR_mid_parent_19TA = scale(savitzkyGolay(mid_parent_19TA, m=1, p=1, w=11))
NIR_mid_parent_20CS = scale(savitzkyGolay(mid_parent_20CS, m=1, p=1, w=11))
NIR_mid_parent_20LY = scale(savitzkyGolay(mid_parent_20LY, m=1, p=1, w=11))

NIR_mid_parent.d1 = rbind(NIR_mid_parent_19CS, NIR_mid_parent_19TA, NIR_mid_parent_20CS, NIR_mid_parent_20LY) 
NIR_mid_parent.ZN1 = tcrossprod(as.matrix(NIR_mid_parent.d1)/ncol(as.matrix(NIR_mid_parent.d1))) #phenomic relationship matrices
dim(NIR_mid_parent.ZN1) #phenomic relationship matrices from hybrids

# Fix MA and LY discrepency
### Mid-parent heterosis ###
# heterosis = F1-MP/MP

combined_NIR_F1 = rbind(pheno_19CS[,-c(1:7)],  #hybrid 
                     pheno_19TA[,-c(1:7)], 
                     pheno_20CS[,-c(1:7)],
                     pheno_20MA[,-c(1:7)])

combined_NIR_MP = final_pedigree_df # mid-parent
MP_heterosis = (combined_NIR_F1- combined_NIR_MP)/combined_NIR_MP
MP_heterosis

MP_heterosis_19CS = MP_heterosis[1:364,]
MP_heterosis_19TA = MP_heterosis[365:543,]
MP_heterosis_20CS = MP_heterosis[544:907,]
MP_heterosis_20LY = MP_heterosis[-c(1:907),]


## estimate first derivative of mid-parent_heterosis NIR spectra

NIR_mid_parent_heterosis_19CS = scale(savitzkyGolay(MP_heterosis_19CS, m=1, p=1, w=11))
NIR_mid_parent_heterosis_19TA = scale(savitzkyGolay(MP_heterosis_19TA, m=1, p=1, w=11))
NIR_mid_parent_heterosis_20CS = scale(savitzkyGolay(MP_heterosis_20CS, m=1, p=1, w=11))
NIR_mid_parent_heterosis_20LY = scale(savitzkyGolay(MP_heterosis_20LY, m=1, p=1, w=11))

NIR_mid_parent_heterosis.d1 = rbind(NIR_mid_parent_heterosis_19CS, 
                                    NIR_mid_parent_heterosis_19TA, 
                                    NIR_mid_parent_heterosis_20CS, 
                                    NIR_mid_parent_heterosis_20LY) 

NIR_mid_parent_heterosis.ZN1 = tcrossprod(as.matrix(NIR_mid_parent_heterosis.d1)/ncol(as.matrix(NIR_mid_parent_heterosis.d1))) #phenomic relationship matrices
dim(NIR_mid_parent_heterosis.ZN1) #phenomic relationship matrices from NIR_mid_parent_heterosis


######### High parent heterosis 
F;M
available.hybrid
write.csv(pheno_combined$pedigree, "total.hybrids.csv")

hybrids = read.csv("total.hybrids.csv")
F$female = rownames(F)
M$male = rownames(M)

females_hybrids <- hybrids %>%
  inner_join(F, by = c("female"="female"))

males_hybrids <- hybrids %>%
  inner_join(M, by = c("male"="male"))

females_numeric = females_hybrids[,-c(1:3)]
males_numeric = males_hybrids[,-c(1:3)]
High_parent_values = pmax(males_numeric, females_numeric)
High_parent_heterosis = (combined_NIR_F1-High_parent_values)/High_parent_values

High_perent_heterosis_19CS = High_parent_heterosis[1:364,]
High_parent_heterosis_19TA = High_parent_heterosis[365:543,]
High_parent_heterosis_20CS = High_parent_heterosis[544:907,]
High_parent_heterosis_20LY = High_parent_heterosis[-c(1:907),]


# calculation of heterosis using high parent values i.e. high parent hterosis

NIR_high_parent_heterosis_19CS = scale(savitzkyGolay(High_perent_heterosis_19CS, m=1, p=1, w=11))
NIR_high_parent_heterosis_19TA = scale(savitzkyGolay(High_parent_heterosis_19TA, m=1, p=1, w=11))
NIR_high_parent_heterosis_20CS = scale(savitzkyGolay(High_parent_heterosis_20CS, m=1, p=1, w=11))
NIR_high_parent_heterosis_20LY = scale(savitzkyGolay(High_parent_heterosis_20LY, m=1, p=1, w=11))

NIR_high_parent_heterosis.d1 = rbind(NIR_high_parent_heterosis_19CS, 
                                     NIR_high_parent_heterosis_19TA, 
                                     NIR_high_parent_heterosis_20CS, 
                                     NIR_high_parent_heterosis_20LY) 

NIR_high_parent_heterosis.ZN1 = tcrossprod(as.matrix(NIR_high_parent_heterosis.d1)/ncol(as.matrix(NIR_high_parent_heterosis.d1))) #phenomic relationship matrices
dim(NIR_high_parent_heterosis.ZN1) 


# Set the predictors/priors to run the model
# mid-parent

mid_parent.ZNZE1 = NIR_mid_parent.ZN1 * ZEZEt
mid_parent_heterosis_ZNZE1 = NIR_mid_parent_heterosis.ZN1 * ZEZEt
high_parent.ZNZE1 = NIR_high_parent_heterosis.ZN1 * ZEZEt



# Set ETAs predictors for models

# Genomic model
Eta1 <- list(list(X = ZE, model = "BRR"),     # Env
             list(K = K1star, model = "RKHS"), # Female
             list(K = K2star, model = "RKHS"), # male
             list(K = K3star, model = "RKHS"), # Femalexmale
             list(K = K4, model = "RKHS"), # Female x environment
             list(K = K5, model = "RKHS"), # male x environemnt
             list(K = K6, model = "RKHS")) #female x male x environment

# Mid-parent model

##### first derivative of mid_parent NIR

Eta2<-list(list(X=ZE,model="BRR"),      #Env
           list(K=NIR_mid_parent.ZN1, model="RKHS"),    #NIR1
           list(K=mid_parent.ZNZE1, model="RKHS")) #NIR1 x Env


##### first derivative of mid_parent heterosis NIR

Eta3<-list(list(X=ZE,model="BRR"),      #Env
           list(K=NIR_mid_parent_heterosis.ZN1, model="RKHS"),    #NIR1
           list(K=mid_parent_heterosis_ZNZE1, model="RKHS"))  #NIR1 x Env

##### first derivative of high_parent heterosis NIR

Eta4<-list(list(X=ZE,model="BRR"),      #Env
           list(K=NIR_high_parent_heterosis.ZN1, model="RKHS"),    #NIR1
           list(K=high_parent.ZNZE1, model="RKHS")) #NIR1 x Env


##### first derivative of mid_parent + high-parent heterosis NIR

Eta5<-list(list(X=ZE,model="BRR"),      #Env
           list(K=NIR_mid_parent_heterosis.ZN1, model="RKHS"),    #NIR1
           list(K=mid_parent_heterosis_ZNZE1, model="RKHS"),
           list(K=NIR_high_parent_heterosis.ZN1, model="RKHS"),    #NIR1
           list(K=high_parent.ZNZE1, model="RKHS")) #NIR1 x Env


##### first derivative of mid_parent NIR + genomic
Eta6 <- list(list(X = ZE, model = "BRR"),     # Env
             list(K = K1star, model = "RKHS"), # Female
             list(K = K2star, model = "RKHS"), # male
             list(K = K3star, model = "RKHS"), #female xmale 
             list(K = K4, model = "RKHS"), #female x env
             list(K = K5, model = "RKHS"), #male x env
             list(K = K6, model = "RKHS"), #female x male x env
             list(K=NIR_mid_parent.ZN1, model="RKHS"),    #NIR1_mid_parent
             list(K=mid_parent.ZNZE1, model="RKHS")) #NIR x env






MP_heterosis = list()
# row_F =1; row_M =1
for (row_F in rownames(F)) {
  for (row_M in rownames(M)) {
    # Compute the average for the current pair of rows
    avg <- (F[row_F, ] + M[row_M, ]) / 2
    
    # Store the result in the list with a descriptive name
    averages[[paste(row_F, row_M, sep = "/")]] <- avg
  }
}


# Next-step
# extract mid-parents only for evaluated hybrids and match it with respective 
# phenotypic values 


### calculate mid-parent heterosis
# hybrid - mid-parent value (should be straightforward)
# hybrid - high parent -- for this: pick the one having higher values


