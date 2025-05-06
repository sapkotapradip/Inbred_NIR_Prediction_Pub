################################### Sapkota, Pradip ####################################################
################################# NIRBLUP Pradip 2025 ##################################################
################################# Agronomic traits #####################################################
################################# multi- environment ###################################################

####################### Across-Environment Prediction ##################################################
########################################################################################################

library(dplyr)
library(prospectr)
library(BGLR)
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
ne

#### loading marker data 
T92I012 <- read.delim(file = "T92I012.txt", header = TRUE, sep = "\t", dec = ".")
parents  = read.delim(file = "parents.txt", header = TRUE, sep = "\t", dec = ".")$pedigree

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

load("K3.Rdata") # loading hybrid relationship created via kronecker product


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

midparent_20CS = read.csv("mid_parent_NIR_20CS_all.csv", check.names = FALSE)
midparent_20MA = read.csv("mid_parent_NIR_20MA_all.csv", check.names = FALSE)

rownames(midparent_20CS) = midparent_20CS[,1]
rownames(midparent_20MA) = midparent_20MA[,1]
midparent_20CS = midparent_20CS[,-1]
midparent_20MA = midparent_20MA[,-1]

# Ensure the row names of averages_df are accessible
# averages_pedigrees <- rownames(averages_df) this one is when calculating using above loop

pedigrees_list = pheno_combined$pedigree
filtered_df_20CS <- midparent_20CS[rownames(midparent_20CS) %in% pedigrees_list, ]
final_pedigree_20CS <- filtered_df_20CS[match(pedigrees_list, rownames(filtered_df_20CS)), ]

filtered_df_20MA <- midparent_20MA[rownames(midparent_20MA) %in% pedigrees_list, ]
final_pedigree_20MA <- filtered_df_20MA[match(pedigrees_list, rownames(filtered_df_20MA)), ]

mid_parent_19CS_20CS = final_pedigree_20CS[1:364,]
mid_parent_19TA_20CS = final_pedigree_20CS[365:543,]
mid_parent_20CS_20CS = final_pedigree_20CS[544:907,]
mid_parent_20MA_20CS = final_pedigree_20CS[-c(1:907),]

mid_parent_19CS_20MA = final_pedigree_20MA[1:364,]
mid_parent_19TA_20MA = final_pedigree_20MA[365:543,]
mid_parent_20CS_20MA = final_pedigree_20MA[544:907,]
mid_parent_20MA_20MA = final_pedigree_20MA[-c(1:907),]

## estimate first derivative of mid-parent NIR spectra
# using 20CS

NIR_mid_parent_19CS_20CS = scale(savitzkyGolay(mid_parent_19CS_20CS, m=1, p=1, w=11))
NIR_mid_parent_19TA_20CS = scale(savitzkyGolay(mid_parent_19TA_20CS, m=1, p=1, w=11))
NIR_mid_parent_20CS_20CS = scale(savitzkyGolay(mid_parent_20CS_20CS, m=1, p=1, w=11))
NIR_mid_parent_20MA_20CS = scale(savitzkyGolay(mid_parent_20MA_20CS, m=1, p=1, w=11))

NIR_mid_parent.d1_20CS = rbind(NIR_mid_parent_19CS_20CS, 
                               NIR_mid_parent_19TA_20CS, 
                               NIR_mid_parent_20CS_20CS, 
                               NIR_mid_parent_20MA_20CS)

NIR_mid_parent.ZN1_20CS = tcrossprod(as.matrix(NIR_mid_parent.d1_20CS)/ncol(as.matrix(NIR_mid_parent.d1_20CS))) #phenomic relationship matrices
dim(NIR_mid_parent.ZN1_20CS) #phenomic relationship matrices from hybrids


# using 20MA

NIR_mid_parent_19CS_20MA = scale(savitzkyGolay(mid_parent_19CS_20MA, m=1, p=1, w=11))
NIR_mid_parent_19TA_20MA = scale(savitzkyGolay(mid_parent_19TA_20MA, m=1, p=1, w=11))
NIR_mid_parent_20CS_20MA = scale(savitzkyGolay(mid_parent_20CS_20MA, m=1, p=1, w=11))
NIR_mid_parent_20MA_20MA = scale(savitzkyGolay(mid_parent_20MA_20MA, m=1, p=1, w=11))

NIR_mid_parent.d1_20MA = rbind(NIR_mid_parent_19CS_20MA, 
                               NIR_mid_parent_19TA_20MA, 
                               NIR_mid_parent_20CS_20MA, 
                               NIR_mid_parent_20MA_20MA)

NIR_mid_parent.ZN1_20MA = tcrossprod(as.matrix(NIR_mid_parent.d1_20MA)/ncol(as.matrix(NIR_mid_parent.d1_20MA))) #phenomic relationship matrices
dim(NIR_mid_parent.ZN1_20MA) #phenomic relationship matrices from hybrids

### Mid-parent heterosis ###
# heterosis = F1-MP/MP

combined_NIR_F1 = rbind(pheno_19CS[,-c(1:7)],  #hybrid 
                        pheno_19TA[,-c(1:7)], 
                        pheno_20CS[,-c(1:7)],
                        pheno_20MA[,-c(1:7)])

combined_NIR_MP = final_pedigree_20CS # mid-parent
MP_heterosis_20CS = (combined_NIR_F1- combined_NIR_MP)/combined_NIR_MP
MP_heterosis_20CS

MP_heterosis_20MA = (combined_NIR_F1- final_pedigree_20MA)/final_pedigree_20MA
MP_heterosis_20MA


MP_heterosis_19CS_20CS = MP_heterosis_20CS[1:364,]
MP_heterosis_19TA_20CS = MP_heterosis_20CS[365:543,]
MP_heterosis_20CS_20CS = MP_heterosis_20CS[544:907,]
MP_heterosis_20MA_20CS = MP_heterosis_20CS[-c(1:907),]

MP_heterosis_19CS_20MA = MP_heterosis_20MA[1:364,]
MP_heterosis_19TA_20MA = MP_heterosis_20MA[365:543,]
MP_heterosis_20CS_20MA = MP_heterosis_20MA[544:907,]
MP_heterosis_20MA_20MA = MP_heterosis_20MA[-c(1:907),]


## estimate first derivative of mid-parent_heterosis NIR spectra
NIR_MP_heterosis_19CS_20CS = scale(savitzkyGolay(MP_heterosis_19CS_20CS, m=1, p=1, w=11))
NIR_MP_heterosis_19TA_20CS = scale(savitzkyGolay(MP_heterosis_19TA_20CS, m=1, p=1, w=11))
NIR_MP_heterosis_20CS_20CS = scale(savitzkyGolay(MP_heterosis_20CS_20CS, m=1, p=1, w=11))
NIR_MP_heterosis_20MA_20CS = scale(savitzkyGolay(MP_heterosis_20MA_20CS, m=1, p=1, w=11))

NIR_MP_heterosis.d1_20CS = rbind(NIR_MP_heterosis_19CS_20CS, #combine all
                            NIR_MP_heterosis_19TA_20CS, 
                            NIR_MP_heterosis_20CS_20CS, 
                            NIR_MP_heterosis_20MA_20CS) 

NIR_MP_heterosis.ZN1_20CS = tcrossprod(as.matrix(NIR_MP_heterosis.d1_20CS)/ncol(as.matrix(NIR_MP_heterosis.d1_20CS))) #phenomic relationship matrices
dim(NIR_MP_heterosis.ZN1_20CS) #phenomic relationship matrices from NIR_mid_parent_heterosis

#20MA

NIR_MP_heterosis_19CS_20MA = scale(savitzkyGolay(MP_heterosis_19CS_20MA, m=1, p=1, w=11))
NIR_MP_heterosis_19TA_20MA = scale(savitzkyGolay(MP_heterosis_19TA_20MA, m=1, p=1, w=11))
NIR_MP_heterosis_20CS_20MA = scale(savitzkyGolay(MP_heterosis_20CS_20MA, m=1, p=1, w=11))
NIR_MP_heterosis_20MA_20MA = scale(savitzkyGolay(MP_heterosis_20MA_20MA, m=1, p=1, w=11))

NIR_MP_heterosis.d1_20MA = rbind(NIR_MP_heterosis_19CS_20MA, #combine all
                                 NIR_MP_heterosis_19TA_20MA, 
                                 NIR_MP_heterosis_20CS_20MA, 
                                 NIR_MP_heterosis_20MA_20MA) 

NIR_MP_heterosis.ZN1_20MA = tcrossprod(as.matrix(NIR_MP_heterosis.d1_20MA)/ncol(as.matrix(NIR_MP_heterosis.d1_20MA))) #phenomic relationship matrices
dim(NIR_MP_heterosis.ZN1_20MA) #phenomic relationship matrices from NIR_mid_parent_heterosis


######### High parent heterosis 
# F;M
# available.hybrid
# write.csv(pheno_combined$pedigree, "total.hybrids.csv")

hybrids = read.csv("total.hybrids.csv")

NIR_pheno_parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
NIR_parents_20CS = NIR_pheno_parents_20CS[,-c(2:5)]
rownames(NIR_parents_20CS) = NIR_parents_20CS$pedigree
NIR_parents_20CS = NIR_parents_20CS[,-1]
F_20CS = NIR_parents_20CS[1:45,] #female parents
M_20CS = NIR_parents_20CS[46:89,] #male parents
F_20CS$female = rownames(F_20CS)
M_20CS$male = rownames(M_20CS)

NIR_pheno_parents_20MA = read.csv("blues_20MA_inbreds_spectra.csv", check.names = FALSE)
NIR_parents_20MA = NIR_pheno_parents_20MA[,-c(2:5)]
rownames(NIR_parents_20MA) = NIR_parents_20MA$pedigree
NIR_parents_20MA = NIR_parents_20MA[,-1]
F_20MA = NIR_parents_20MA[1:45,] #female parents
M_20MA = NIR_parents_20MA[46:89,] #male parents
F_20MA$female = rownames(F_20MA)
M_20MA$male = rownames(M_20MA)


females_hybrids_20CS <- hybrids %>%
  inner_join(F_20CS, by = c("female"="female"))

males_hybrids_20CS <- hybrids %>%
  inner_join(M_20CS, by = c("male"="male"))


females_hybrids_20MA <- hybrids %>%
  inner_join(F_20MA, by = c("female"="female"))

males_hybrids_20MA <- hybrids %>%
  inner_join(M_20MA, by = c("male"="male"))


females_numeric_20CS = females_hybrids_20CS[,-c(1:3)]
males_numeric_20CS = males_hybrids_20CS[,-c(1:3)]

females_numeric_20MA = females_hybrids_20MA[,-c(1:3)]
males_numeric_20MA = males_hybrids_20MA[,-c(1:3)]

High_parent_values_20CS = pmax(males_numeric_20CS, females_numeric_20CS)
High_parent_heterosis_20CS = (combined_NIR_F1-High_parent_values_20CS)/High_parent_values_20CS

High_parent_values_20MA = pmax(males_numeric_20MA, females_numeric_20MA)
High_parent_heterosis_20MA = (combined_NIR_F1-High_parent_values_20MA)/High_parent_values_20MA


HP_heterosis_19CS_20CS = High_parent_heterosis_20CS[1:364,]
HP_heterosis_19TA_20CS = High_parent_heterosis_20CS[365:543,]
HP_heterosis_20CS_20CS = High_parent_heterosis_20CS[544:907,]
HP_heterosis_20MA_20CS = High_parent_heterosis_20CS[-c(1:907),]

HP_heterosis_19CS_20MA = High_parent_heterosis_20MA[1:364,]
HP_heterosis_19TA_20MA = High_parent_heterosis_20MA[365:543,]
HP_heterosis_20CS_20MA = High_parent_heterosis_20MA[544:907,]
HP_heterosis_20MA_20MA = High_parent_heterosis_20MA[-c(1:907),]

# calculation of heterosis using high parent values i.e. high parent hterosis
# using 20CS
NIR_HP_heterosis_19CS_20CS = scale(savitzkyGolay(HP_heterosis_19CS_20CS, m=1, p=1, w=11))
NIR_HP_heterosis_19TA_20CS = scale(savitzkyGolay(HP_heterosis_19TA_20CS, m=1, p=1, w=11))
NIR_HP_heterosis_20CS_20CS = scale(savitzkyGolay(HP_heterosis_20CS_20CS, m=1, p=1, w=11))
NIR_HP_heterosis_20MA_20CS = scale(savitzkyGolay(HP_heterosis_20MA_20CS, m=1, p=1, w=11))

NIR_HP_heterosis.d1_20CS = rbind(NIR_HP_heterosis_19CS_20CS, 
                            NIR_HP_heterosis_19TA_20CS, 
                            NIR_HP_heterosis_20CS_20CS, 
                            NIR_HP_heterosis_20MA_20CS) 

NIR_HP_heterosis.ZN1_20CS = tcrossprod(as.matrix(NIR_HP_heterosis.d1_20CS)/ncol(as.matrix(NIR_HP_heterosis.d1_20CS))) #phenomic relationship matrices
dim(NIR_HP_heterosis.ZN1_20CS) 

# using 20MA
NIR_HP_heterosis_19CS_20MA = scale(savitzkyGolay(HP_heterosis_19CS_20MA, m=1, p=1, w=11))
NIR_HP_heterosis_19TA_20MA = scale(savitzkyGolay(HP_heterosis_19TA_20MA, m=1, p=1, w=11))
NIR_HP_heterosis_20CS_20MA = scale(savitzkyGolay(HP_heterosis_20CS_20MA, m=1, p=1, w=11))
NIR_HP_heterosis_20MA_20MA = scale(savitzkyGolay(HP_heterosis_20MA_20MA, m=1, p=1, w=11))

NIR_HP_heterosis.d1_20MA = rbind(NIR_HP_heterosis_19CS_20MA, 
                                 NIR_HP_heterosis_19TA_20MA, 
                                 NIR_HP_heterosis_20CS_20MA, 
                                 NIR_HP_heterosis_20MA_20MA) 

NIR_HP_heterosis.ZN1_20MA = tcrossprod(as.matrix(NIR_HP_heterosis.d1_20MA)/ncol(as.matrix(NIR_HP_heterosis.d1_20MA))) #phenomic relationship matrices
dim(NIR_HP_heterosis.ZN1_20MA) 

## use high parent values
HP_19CS_20CS = High_parent_values_20CS[1:364,]
HP_19TA_20CS = High_parent_values_20CS[365:543,]
HP_20CS_20CS = High_parent_values_20CS[544:907,]
HP_20MA_20CS = High_parent_values_20CS[-c(1:907),]


HP_19CS_20MA = High_parent_values_20MA[1:364,]
HP_19TA_20MA = High_parent_values_20MA[365:543,]
HP_20CS_20MA = High_parent_values_20MA[544:907,]
HP_20MA_20MA = High_parent_values_20MA[-c(1:907),]

# using 20CS
NIR_HP_19CS_20CS = scale(savitzkyGolay(HP_19CS_20CS, m=1, p=1, w=11))
NIR_HP_19TA_20CS = scale(savitzkyGolay(HP_19TA_20CS, m=1, p=1, w=11))
NIR_HP_20CS_20CS = scale(savitzkyGolay(HP_20CS_20CS, m=1, p=1, w=11))
NIR_HP_20MA_20CS = scale(savitzkyGolay(HP_20MA_20CS, m=1, p=1, w=11))

NIR_HP.d1_20CS = rbind(NIR_HP_19CS_20CS, 
                       NIR_HP_19TA_20CS, 
                       NIR_HP_20CS_20CS, 
                       NIR_HP_20MA_20CS) 

NIR_HP_ZN1_20CS = tcrossprod(as.matrix(NIR_HP.d1_20CS)/ncol(as.matrix(NIR_HP.d1_20CS))) #phenomic relationship matrices
dim(NIR_HP_ZN1_20CS) 

# using 20MA
NIR_HP_19CS_20MA = scale(savitzkyGolay(HP_19CS_20MA, m=1, p=1, w=11))
NIR_HP_19TA_20MA = scale(savitzkyGolay(HP_19TA_20MA, m=1, p=1, w=11))
NIR_HP_20CS_20MA = scale(savitzkyGolay(HP_20CS_20MA, m=1, p=1, w=11))
NIR_HP_20MA_20MA = scale(savitzkyGolay(HP_20MA_20MA, m=1, p=1, w=11))

NIR_HP.d1_20MA = rbind(NIR_HP_19CS_20MA, 
                       NIR_HP_19TA_20MA, 
                       NIR_HP_20CS_20MA, 
                       NIR_HP_20MA_20MA) 

NIR_HP_ZN1_20MA = tcrossprod(as.matrix(NIR_HP.d1_20MA)/ncol(as.matrix(NIR_HP.d1_20MA))) #phenomic relationship matrices
dim(NIR_HP_ZN1_20MA) 

# calculating NIR x E interaction for different models using Hadamard product
# 20CS
mid_parent.ZNZE1.CS           = NIR_mid_parent.ZN1_20CS * ZEZEt
mid_parent_heterosis_ZNZE1.CS = NIR_MP_heterosis.ZN1_20CS * ZEZEt
high_parent_heterosis_ZNZE1.CS          = NIR_HP_heterosis.ZN1_20CS * ZEZEt
high_parent_ZNZE1.CS          = NIR_HP_ZN1_20CS * ZEZEt


# calculating NIR x E interaction for different models using Hadamard product
#20MA

mid_parent.ZNZE1.MA           = NIR_mid_parent.ZN1_20MA * ZEZEt
mid_parent_heterosis_ZNZE1.MA = NIR_MP_heterosis.ZN1_20MA * ZEZEt
high_parent_heterosis.ZNZE1.MA          = NIR_HP_heterosis.ZN1_20MA * ZEZEt
high_parent_ZNZE1.MA          = NIR_HP_ZN1_20MA * ZEZEt

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

Eta2<-list(list(X = ZE, model = "BRR"),      #Env
           list(K= NIR_mid_parent.ZN1_20CS, model = "RKHS"),    #NIR1
           list(K = mid_parent.ZNZE1.CS, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent heterosis NIR
Eta3<-list(list(X = ZE,model="BRR"),      #Env
           list(K = NIR_MP_heterosis.ZN1_20CS, model="RKHS"),    #NIR1
           list(K = mid_parent_heterosis_ZNZE1.CS, model="RKHS"))  #NIR1 x Env

##### first derivative of high_parent heterosis NIR
Eta4<-list(list(X = ZE,model = "BRR"),      #Env
           list(K = NIR_HP_heterosis.ZN1_20CS, model = "RKHS"),    #NIR1
           list(K = high_parent_heterosis_ZNZE1.CS, model = "RKHS")) #NIR1 x Env


#### first derivative of high parent NIR
Eta5<-list(list(X = ZE,model = "BRR"),      #Env
           list(K = high_parent_ZNZE1.CS, model = "RKHS"),    #NIR1
           list(K = NIR_HP_ZN1_20CS, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent + high-parent heterosis NIR

Eta6<-list(list(X = ZE,model="BRR"),      #Env
           list(K = NIR_MP_heterosis.ZN1_20CS, model="RKHS"),    #NIR1
           list(K = mid_parent_heterosis_ZNZE1.CS, model="RKHS"),
           list(K = NIR_HP_heterosis.ZN1_20CS, model="RKHS"),    #NIR1
           list(K = high_parent_heterosis_ZNZE1.CS, model="RKHS")) #NIR1 x Env


##### first derivative of mid_parent NIR + genomic
Eta7 <- list(list(X = ZE, model = "BRR"),     # Env
             list(K = K1star, model = "RKHS"), # Female
             list(K = K2star, model = "RKHS"), # male
             list(K = K3star, model = "RKHS"), #female xmale 
             list(K = K4, model = "RKHS"), #female x env
             list(K = K5, model = "RKHS"), #male x env
             list(K = K6, model = "RKHS"), #female x male x env
             list(K=NIR_mid_parent.ZN1_20CS, model="RKHS"),    #NIR1_mid_parent
             list(K=mid_parent.ZNZE1.CS, model="RKHS")) #NIR x env


####Use 20MA

##### first derivative of mid_parent NIR

Eta8<-list(list(X = ZE, model = "BRR"),      #Env
           list(K= NIR_mid_parent.ZN1_20MA, model = "RKHS"),    #NIR1
           list(K = mid_parent.ZNZE1.MA, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent heterosis NIR

Eta9<-list(list(X = ZE,model="BRR"),      #Env
           list(K = NIR_MP_heterosis.ZN1_20MA, model="RKHS"),    #NIR1
           list(K = mid_parent_heterosis_ZNZE1.MA, model="RKHS"))  #NIR1 x Env

##### first derivative of high_parent heterosis NIR

Eta10<-list(list(X = ZE,model = "BRR"),      #Env
           list(K = NIR_HP_heterosis.ZN1_20MA, model = "RKHS"),    #NIR1
           list(K = High_parent_heterosis_20MA, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent + high-parent heterosis NIR

Eta11<-list(list(X = ZE,model="BRR"),      #Env
           list(K = NIR_MP_heterosis.ZN1_20MA, model="RKHS"),    #NIR1
           list(K = mid_parent_heterosis_ZNZE1.MA, model="RKHS"),
           list(K = NIR_HP_heterosis.ZN1_20MA, model="RKHS"),    #NIR1
           list(K = high_parent_heterosis.ZNZE1.MA, model="RKHS")) #NIR1 x Env


#### first derivative of high parent NIR
Eta12<- list(list(X = ZE,model = "BRR"),      #Env
           list(K = NIR_HP_ZN1_20MA, model = "RKHS"),    #NIR1
           list(K = high_parent_ZNZE1.MA, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent NIR + genomic 20MA

Eta13 <- list(list(X = ZE, model = "BRR"),     # Env
             list(K = K1star, model = "RKHS"), # Female
             list(K = K2star, model = "RKHS"), # male
             list(K = K3star, model = "RKHS"), #female xmale 
             list(K = K4, model = "RKHS"), #female x env
             list(K = K5, model = "RKHS"), #male x env
             list(K = K6, model = "RKHS"), #female x male x env
             list(K=NIR_mid_parent.ZN1_20MA, model="RKHS"),    #NIR1_mid_parent
             list(K=mid_parent.ZNZE1.MA, model="RKHS")) #NIR x env


## Just GCA with no SCA
Eta14 <- list(list(X = ZE, model = "BRR"),     # Env
             list(K = K1star, model = "RKHS"), # Female
             list(K = K2star, model = "RKHS"), # male
             list(K = K4, model = "RKHS"), # Female x environment
             list(K = K5, model = "RKHS")) # male x environemnt


# also add high parent values
             

###Missing for multi-trait models using DA and PH


### 6 different models 
# MODEL 1: Genomic model
# MODEL 2: Mid-parent model-20CS
# MODEL 3: mid_parent heterosis NIR-20CS
# MODEL 4: high_parent heterosis NIR-20CS
# MODEL 5: mid_parent + high-parent heterosis NIR-20CS
# MODEL 6: first derivative of mid_parent NIR + genomic-20CS
# MODEL 7: Mid-parent model-20MA
# MODEL 8: mid_parent heterosis NIR-20MA
# MODEL 9: high_parent heterosis NIR-20MA
# MODEL 10: mid-parent + high_parent heterosis NIR-20CS
# MODEL 11: first derivative of mid_parent NIR + genomic-20MA

# tr=1
Models <- list(Eta1, Eta2, Eta3, Eta4, Eta5, Eta6, Eta7, Eta8, Eta9, Eta10, Eta11, Eta12, Eta13, Eta14)
traitnames <- c("yield", "da", "ph", "starch", "protein", "fat", "fiber", "ash")
pheno_combined[1:10,1:10]
#parents = read.delim("parents.txt")$"pedigree"


# pheno_yield = pheno_combined[,1:5]
# pheno_ph = pheno_combined[,c(1:4,6)]
# pheno_da = pheno_combined[,c(1:4,7)]
# colnames(pheno_yield)[5] = "blue"
# colnames(pheno_ph)[5] = "blue"
# colnames(pheno_da)[5] = "blue"
# 
# 
# write.csv(pheno_yield, "pheno_yield.csv")
# write.csv(pheno_da, "pheno_da.csv")
# write.csv(pheno_ph, "pheno_ph.csv")


# 70:30 partition testing

# tr=3
set.seed(1234)
for (tr in 1:length(traitnames)) {
  
  if (tr == 1) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==2) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==3) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==4) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==5) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==6) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==7) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  }
  
  
  # cross-validation #
  parents = c(unique(pheno$female), unique(pheno$male))
  hybrid = as.character(unique(pheno$pedigree))
  Phenotype_data1 = pheno
  
  set.seed(123)
  cycles = 10
  CV1 = list()
  CV2 = list()
  CV3_20LY = list()
  CV3_19TA = list()
  CV4 = list()
  
  #MODEL =2; rep_num=1
  for (MODEL in 1:length(Models)) {  
    
    for (rep_num in 1:5) {
      set.seed(123)
      CVa = sample(parents[1:89], 15, replace = FALSE)
      train_geno <- setdiff(parents[1:89], CVa)
      pheno$y = pheno$blue
      
      for (a in 1:15) {
        pheno <- pheno %>% mutate(y = replace(y, female == CVa[a], NA))
      }
      
      for (b in 1:15) {
        pheno <- pheno %>% mutate(y = replace(y, male == CVa[b], NA))
      }
      
      test_pheno <- subset(pheno, is.na(y))
      test_geno = unique(test_pheno$pedigree) # 30% of hybrid were removed
      train_geno = setdiff(hybrid, test_geno)
      
      # Preparing for CV1
      CV_Data_1_2<-Phenotype_data1
      CV_Data_1_2$Y<-NA
      CV_Data_1_2$Y[CV_Data_1_2$pedigree%in%train_geno]<-CV_Data_1_2$blue[CV_Data_1_2$pedigree%in%train_geno] 
      
      y_t<-as.numeric(CV_Data_1_2$Y)
      
      fit<-BGLR(y=y_t,
                ETA=Models[[MODEL]],
                nIter=5000,
                burnIn=1000, 
                thin=10) #nIter=5000,burnIn=1000, thin =10
      
      CV_Data_1_2$yhat <- fit$yHat
      
      
      # CV1 # untested genotypes in observed environment
      df_test <- subset(CV_Data_1_2, CV_Data_1_2$pedigree %in% test_geno)
      CV1[[(rep_num)]] <- as.data.frame(df_test %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat,use = "complete.obs")))
      
      # Preparing for CV2
      # CV2 simulates sparse testing where 30% of hybrids were sampled from each environments
      
      ne <- as.vector(table(pheno$env)) ## counting the number of observations
      ne
      
      # CV_Data_1_2<-Phenotype_data1
      # CV_Data_1_2$Y = CV_Data_1_2$blue
      # 
      # index = c(sample(1:364,round(0.30*364)), sample(365:543,round(0.30*179)), 
      #           sample(544:907,round(0.30*364)), sample(908:1094,round(0.30*187)))
      # 
      # CV_Data_1_2$Y[index] <- NA
      
      test_geno = sample(hybrids, length(unique(pedigrees_list))*0.3)
      
      CV_Data_1_2<-Phenotype_data1
      CV_Data_1_2$Y<-NA
      CV_Data_1_2$Y[CV_Data_1_2$pedigree%in%train_geno]<-CV_Data_1_2$blue[CV_Data_1_2$pedigree%in%train_geno] 
      
      y_t1<-as.numeric(CV_Data_1_2$Y)
      
      fit1<-BGLR(y=y_t1,
                 ETA=Models[[MODEL]],
                 nIter=5000,
                 burnIn=1000, 
                 thin=10) #nIter=5000,burnIn=1000, thin =10
      
      CV_Data_1_2$yhat1 <- fit1$yHat
      
      df_test1 <- subset(CV_Data_1_2, CV_Data_1_2$pedigree %in% test_geno)
      CV2[[(rep_num)]] <- as.data.frame(df_test1 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat1,use = "complete.obs")))
      
      # For CV3_20LY
      # leave environments out #20LY
      CV_Data_1_2<-Phenotype_data1
      CV_Data_1_2$Y = CV_Data_1_2$blue
      
      CV_Data_1_2$Y[CV_Data_1_2$env == "20LY"] <- NA
      
      y_t2<-as.numeric(CV_Data_1_2$Y)
      
      fit2<-BGLR(y=y_t2,
                 ETA=Models[[MODEL]],
                 nIter=5000,
                 burnIn=1000, 
                 thin=10) #nIter=5000,burnIn=1000, thin =10
      
      CV_Data_1_2$yhat2 <- fit2$yHat
      
      df_test2 <- subset(CV_Data_1_2, env == "20LY")
      CV3_20LY[[(rep_num)]] <- as.data.frame(df_test2 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat2,use = "complete.obs")))
      
      # For CV3_19TA
      # leave environments out #19TA
      CV_Data_1_2<-Phenotype_data1
      CV_Data_1_2$Y = CV_Data_1_2$blue
      
      CV_Data_1_2$Y[CV_Data_1_2$env == "19TA"] <- NA
      
      y_t3<-as.numeric(CV_Data_1_2$Y)
      
      fit3<-BGLR(y=y_t3,
                 ETA=Models[[MODEL]],
                 nIter=5000,
                 burnIn=1000, thin=10) #nIter=5000,burnIn=1000, thin =10
      
      CV_Data_1_2$yhat3 <- fit3$yHat
      
      df_test3 <- subset(CV_Data_1_2, CV_Data_1_2$env == "19TA")
      CV3_19TA[[(rep_num)]] <- as.data.frame(df_test3 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat3,use = "complete.obs")))
      
      # For CV4 #remove 30% hybrids
      CV_Data_1_2<-Phenotype_data1
      CV_Data_1_2$Y = CV_Data_1_2$blue

      sample()

      y_t3<-as.numeric(CV_Data_1_2$Y)
      fit3<-BGLR(y=y_t3,ETA=Models[[MODEL]],nIter=5000,burnIn=1000, thin=10) #nIter=5000,burnIn=1000, thin =10
      CV_Data_1_2$yhat3 <- fit3$yHat

      df_test3 <- subset(CV_Data_1_2, CV_Data_1_2$env == "19TA")
      CV3_19TA[[(rep_num)]] <- as.data.frame(df_test3 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat3,use = "complete.obs")))

      
    }
    
    #rep_num =1
    if (rep_num == 5) {
      CV1out <- plyr::ldply(CV1, data.frame)
      CV2out <- plyr::ldply(CV2, data.frame)
      CV3out_20LY <- plyr::ldply(CV3_20LY, data.frame)
      CV3out_19TA <- plyr::ldply(CV3_19TA, data.frame)
      write.csv(CV1out, file = paste("ACC_", traitnames[tr],"_CV1_", MODEL, ".csv", sep=""), row.names = FALSE)
      write.csv(CV2out, file = paste("ACC_", traitnames[tr],"_CV2_", MODEL, ".csv", sep=""), row.names = FALSE)
      write.csv(CV3out_20LY, file = paste("ACC_", traitnames[tr],"_CV3_20LY_", MODEL, ".csv", sep=""), row.names = FALSE)
      write.csv(CV3out_19TA, file = paste("ACC_", traitnames[tr],"_CV3_19TA_", MODEL, ".csv", sep=""), row.names = FALSE)
    }
  }
}

# Plot results
setwd("output/results_4_29/")
###lets plot
library(plyr)
library(readr)
library(dplyr)
library(stringr)

list_csv_files <- list.files("../results_4_29/")
df2 <- readr::read_csv(list_csv_files, id = "file_name") %>% as.data.frame()
df2
write.csv(df2, "Pred.ability_11models_3cv_5rep.csv")
df2 = read.csv("Pred.ability_11models_3cv_5rep.csv")

df1 = df2 %>%
  mutate(split_file_name = str_split(file_name, "_", simplify = TRUE)) %>%
  mutate(trait = split_file_name[, 2],
         CV = split_file_name[, 3],
         model = split_file_name[, 3],
         ".csv") 

df1 =
  df1  %>%
  mutate(
    model = case_when(
      model == "1.csv" ~ "G",
      model == "2.csv" ~ "MP",
      model == "3.csv" ~ "MPH",
      model == "4.csv" ~ "HPH",
      model == "5.csv" ~ "MPH + HPH",
      model == "6.csv" ~ "G + MP",
      TRUE~model),
    trait = case_when(
      trait == "yield" ~ "Grain Yield",
      trait == "ph" ~ "Plant Height",
      trait == "da" ~ "Days to Anthesis",
      TRUE~trait),
    CV_scheme = case_when(
      CV_scheme == "CV1" ~ "single_trait",
      CV_scheme == "CV2" ~ "GY + PH",
      CV_scheme == "CV3" ~ "DA + PH"
    )
  )

df1 = df1[,c(3,4,6,7)]


df1 <- as.data.frame(df2 %>%  dplyr::group_by(traits,Model,CV) %>% 
                       dplyr::summarise(M = mean(cor, na.rm=TRUE),
                                        SD = sd(cor, na.rm=TRUE)))
#df1 = read.csv("df1.csv")
df1$traits <- factor(df1$traits, levels =  c("yield", "ph", "da"))
df1$Model <- factor(df1$Model, levels =  c("M1", "M2", "M3", "M4", "M5", "M6",
                                           "M7", "M8", "M9", "M10", "M11"))
df1$CV <- factor(df1$CV, levels =  c("CV1", "CV2", "CV3"))


library(tidyr)
library(ggplot2)



p = ggplot(df1, aes(CV, y=M, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=round(M,2)  ), hjust=3, color="white",
            position = position_dodge(0.9), angle = 90,size=3.5)+
  geom_errorbar(aes(ymin=M, ymax=M+SD), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  facet_grid(~traits)+
  scale_y_continuous("Prediction ability")+
  xlab("Cross Validation Schemes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
  
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=5, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=12, face="bold", colour = "black"),    
    axis.title.y = element_text(size=12, face="bold", colour = "black"),    
    axis.text.x = element_text(size=7, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=7, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 8, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 8, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    legend.position = "none"
    #axis.text.x.bottom = element_blank()
  )

jpeg("Pred.ability.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()

####################################################################################################







#################use parental phenotypic data:
# multi-trait GS and PS
# use parents data from 20CS
#setwd("../../")
parents_pheno_20CS = read.csv("blues_20CS_pheno_parents.csv")
parents_pheno_20MA = read.csv("blues_20MA_pheno_parents.csv")

rownames(parents_pheno_20CS) = parents_pheno_20CS$pedigree
parents_pheno_20CS = parents_pheno_20CS[,-c(1,2)]

F_20CS = parents_pheno_20CS[1:45,] #female parents
M_20CS = parents_pheno_20CS[46:89,] #male parents


library(dplyr)
# available.hybrid

averages = list()
# row_F =1; row_M =1
for (row_F in rownames(F_20CS)) {
  for (row_M in rownames(M_20CS)) {
    # Compute the average for the current pair of rows
    avg <- (F_20CS[row_F, ] + M_20CS[row_M, ]) / 2

    # Store the result in the list with a descriptive name
    averages[[paste(row_F, row_M, sep = "/")]] <- avg
  }
}

averages

# Convert the list of averages into a data frame
averages_df <- do.call(rbind, averages)
write.csv(averages_df, "mid_parent_pheno_20CS_all.csv")

######### mid_parent values

rownames(parents_pheno_20MA) = parents_pheno_20MA$pedigree
parents_pheno_20MA = parents_pheno_20MA[,-c(1,2)]

F_20MA = parents_pheno_20MA[1:45,] #female parents
M_20MA = parents_pheno_20MA[46:89,] #male parents


library(dplyr)
# available.hybrid

averages = list()
# row_F =1; row_M =1
for (row_F in rownames(F_20MA)) {
  for (row_M in rownames(M_20MA)) {
    # Compute the average for the current pair of rows
    avg <- (F_20MA[row_F, ] + M_20MA[row_M, ]) / 2
    
    # Store the result in the list with a descriptive name
    averages[[paste(row_F, row_M, sep = "/")]] <- avg
  }
}

averages

# Convert the list of averages into a data frame
averages_df <- do.call(rbind, averages)
write.csv(averages_df, "mid_parent_pheno_20MA_all.csv")

## combining with 20CS parents values
pheno_hybrids_yield = read.csv("pheno_yield.csv")
pheno_hybrids_da = read.csv("pheno_da.csv")
pheno_hybrids_ph = read.csv("pheno_ph.csv")

## using mid parent of 20CS
mid_parent_pheno_20CS = read.csv("mid_parent_pheno_20CS_all.csv")


pedigrees_list = pheno_hybrids_yield$pedigree
filtered_df <- mid_parent_pheno_20CS[mid_parent_pheno_20CS$X %in% pedigrees_list, ]
final_pedigree_df <- filtered_df[match(pedigrees_list, filtered_df$X), ]

colnames(final_pedigree_df)[2] = "gy_mid_parent"
colnames(final_pedigree_df)[3] = "ph_mid_parent"
colnames(final_pedigree_df)[4] = "dy_mid_parent"
pheno_gy_20CS = cbind(pheno_hybrids_yield,final_pedigree_df[,c(2,3,4)])
pheno_gy_20CSMP = pheno_gy_20CS[,-1]

pheno_da_20CSMP = cbind(pheno_hybrids_da, final_pedigree_df[,c(2,3,4)])
pheno_ph_20CSMP = cbind(pheno_hybrids_ph, final_pedigree_df[,c(2,3,4)])

pheno_da_20CSMP = pheno_da_20CSMP[,-1]
pheno_ph_20CSMP = pheno_ph_20CSMP[,-1]


write.csv(pheno_gy_20CSMP, "pheno_yield_20CSMP.csv")
write.csv(pheno_da_20CSMP, "pheno_da_20CSMP.csv")
write.csv(pheno_ph_20CSMP, "pheno_ph_20CSMP.csv")

## using mid parent of 20MA
mid_parent_pheno_20MA = read.csv("mid_parent_pheno_20MA_all.csv")

pedigrees_list = pheno_hybrids_yield$pedigree
filtered_df <- mid_parent_pheno_20MA[mid_parent_pheno_20MA$X %in% pedigrees_list, ]
final_pedigree_df <- filtered_df[match(pedigrees_list, filtered_df$X), ]

colnames(final_pedigree_df)[2] = "gy_mid_parent"
colnames(final_pedigree_df)[3] = "ph_mid_parent"
colnames(final_pedigree_df)[4] = "dy_mid_parent"
pheno_gy_20MA = cbind(pheno_hybrids_yield,final_pedigree_df[,c(2,3,4)])
pheno_gy_20MAMP = pheno_gy_20MA[,-1]

pheno_da_20MAMP = cbind(pheno_hybrids_da, final_pedigree_df[,c(2,3,4)])
pheno_ph_20MAMP = cbind(pheno_hybrids_ph, final_pedigree_df[,c(2,3,4)])

pheno_da_20MAMP = pheno_da_20MAMP[,-1]
pheno_ph_20MAMP = pheno_ph_20MAMP[,-1]


write.csv(pheno_gy_20MAMP, "pheno_yield_20MAMP.csv")
write.csv(pheno_da_20MAMP, "pheno_da_20MAMP.csv")
write.csv(pheno_ph_20MAMP, "pheno_ph_20MAMP.csv")

### run the model with multi-trait PS models


# tr=1
set.seed(1234)
mpnames = c("20CSMP", "20MAMP")

## adding multi_trait models also

# tr =1; mp =1
for (tr in 1:length(traitnames)) {
  for (mp in 1:length(mpnames)){
    
  if (tr == 1) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], "_", mpnames[mp], ".csv", sep = ""))
  } else if (tr ==2) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], "_", mpnames[mp], ".csv", sep = ""))
  } else {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr],"_", mpnames[mp], ".csv", sep = ""))
  }
  
  # cross-validation #
  hybrid = as.character(unique(pheno$pedigree))
  Phenotype_data1 = pheno
  
  cycles = 10
  CV3 = list()
  CV2 = list()
  CV1 = list()
  
  
  #MODEL =2; rep_num=1
  for (MODEL in 1:length(Models)) {  
    
    for (rep_num in 1:5) {
      test_geno <- sample(hybrid, round(0.3*length(hybrid)))  
      train_geno <-setdiff(hybrid, test_geno)
      
      CV_Data_1_2<-Phenotype_data1
      CV_Data_1_2$Y<-NA
      CV_Data_1_2$Y[CV_Data_1_2$pedigree%in%train_geno]<-CV_Data_1_2$blue[CV_Data_1_2$pedigree%in%train_geno] 
      
      y_t<-as.numeric(CV_Data_1_2$Y)
      fit<-BGLR(y=y_t,ETA=Models[[MODEL]],nIter=500,burnIn=100, thin=10) #nIter=5000,burnIn=1000, thin =10
      CV_Data_1_2$yhat <- fit$yHat
      
      
      # CV1
      df_test <- subset(CV_Data_1_2, CV_Data_1_2$pedigree %in% test_geno)
      CV1[[(rep_num)]] <- as.data.frame(df_test %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat,use = "complete.obs")))
      
      # multi_trait using yield, anthesis and height
      CV_Data_1_2<-Phenotype_data1
      CV_Data_1_2$Y<-NA
      CV_Data_1_2$Y[CV_Data_1_2$pedigree%in%train_geno]<-CV_Data_1_2$blue[CV_Data_1_2$pedigree%in%train_geno]
      
      y_t1 = as.matrix(CV_Data_1_2[,c(7,8,9,10)]) #all the traits
      fit1<- Multitrait(y = y_t1, ETA = Models[[MODEL]], nIter = 500, burnIn = 100, thin = 10)
      whichNA = fit1$missing_records[fit1$patterns[,4]]
      pred = fit1$ETAHat[whichNA,] #predicted values
      
      df_test1 <- subset(pheno, pheno$pedigree %in% test_geno)
      df_test1 = cbind(df_test1, pred[,4])
      colnames(df_test1)[10] = "yhat1"
      CV2[[rep_num]] <- as.data.frame(df_test1 %>% group_by(env) %>% dplyr::summarize(cor = cor(blue, yhat1, use = "complete.obs")))
      
      
      # multi_trait using anthesis and height
      CV_Data_1_2<-Phenotype_data1
      CV_Data_1_2$Y<-NA
      CV_Data_1_2$Y[CV_Data_1_2$pedigree%in%train_geno]<-CV_Data_1_2$blue[CV_Data_1_2$pedigree%in%train_geno]
      
      y_t2 = as.matrix(CV_Data_1_2[,c(8,9,10)]) # days to anthesis and plant height
      fit2<- Multitrait(y = y_t2, ETA = Models[[MODEL]], nIter = 500, burnIn = 100, thin = 10)
      whichNA = fit2$missing_records[fit2$patterns[,3]]
      pred = fit2$ETAHat[whichNA,] #predicted values
      
      df_test2 <- subset(pheno, pheno$pedigree %in% test_geno)
      df_test2 = cbind(df_test2, pred[,3])
      colnames(df_test2)[10] = "yhat2"
      CV3[[rep_num]] <- as.data.frame(df_test2 %>% group_by(env) %>% dplyr::summarize(cor = cor(blue, yhat2, use = "complete.obs")))
    }
    #rep_num =1
    if (rep_num == 5) {
      CV1out <- plyr::ldply(CV1, data.frame)
      CV2out <- plyr::ldply(CV2, data.frame)
      CV3out <- plyr::ldply(CV3, data.frame)
      write.csv(CV1out, file = paste("ACC_singltrait_", traitnames[tr], "_", mpnames[mp], "_CV1_", MODEL, "_", ".csv", sep=""), row.names = FALSE)
      write.csv(CV2out, file = paste("ACC_multitrait_", traitnames[tr], "_", mpnames[mp], "_CV2_", MODEL, "_", ".csv", sep=""), row.names = FALSE)
      write.csv(CV3out, file = paste("ACC_multitrait_", traitnames[tr], "_", mpnames[mp], "_CV3_", MODEL, "_", ".csv", sep=""), row.names = FALSE)
      
    }
  }
  }
}

########## wow; just wow


rownames(pheno_hybrids_yield_20MAyield) = NULL
pheno_hybrids_yield_20MAparents = pheno_hybrids_yield_20MAyield[,-1]

### add parent's pheno from 20CS


averages_pedigrees = rownames(averages_df)
pedigrees_list = pheno_combined$pedigree
filtered_df <- averages_df[rownames(averages_df) %in% pedigrees_list, ]
final_pedigree_df <- filtered_df[match(pedigrees_list, rownames(filtered_df)), ]



# calculate mid_parent values

#leave parents out model
# random forest and svm for GS and PS
# multi-trait GS

####### plotting the results
# quick view of result

setwd("output/11_models/")
###lets plot
library(plyr)
library(readr)
library(dplyr)
library(stringr)

list_csv_files <- list.files("../11_models/")
df2 <- readr::read_csv(list_csv_files, id = "file_name") %>% as.data.frame()
df2
write.csv(df2, "Pred.ability.70_30_5rep_11models.csv")
df2 = read.csv("Pred.ability.70_30_5rep_11models.csv")

df1 = df2 %>%
  mutate(split_file_name = str_split(file_name, "_", simplify = TRUE)) %>%
  mutate(trait = split_file_name[, 2],
         model = split_file_name[, 4], ".csv") 

df1 =
  df1  %>%
  mutate(
    model = case_when(
      model == "1.csv" ~ "G",
      model == "2.csv" ~ "MP",
      model == "3.csv" ~ "MPH",
      model == "4.csv" ~ "HPH",
      model == "5.csv" ~ "MPH + HPH",
      model == "6.csv" ~ "G + MP",
      TRUE~model),
    trait = case_when(
      trait == "yield" ~ "Grain Yield",
      trait == "ph" ~ "Plant Height",
      trait == "da" ~ "Days to Anthesis",
      TRUE~trait),
    CV_scheme = case_when(
      CV_scheme == "CV1" ~ "single_trait",
      CV_scheme == "CV2" ~ "GY + PH",
      CV_scheme == "CV3" ~ "DA + PH"
    )
  )

df1 = df1[,c(3,4,6,7)]


df1 <- as.data.frame(df1 %>%  dplyr::group_by(trait,model) %>% 
                       dplyr::summarise(M = mean(cor, na.rm=TRUE),
                                        SD = sd(cor, na.rm=TRUE)))
#df1 = read.csv("df1.csv")
df1$trait <- factor(df1$trait, levels =  c("Grain Yield", "Plant Height", "Days to Anthesis"))
df1$model <- factor(df1$model, levels =  c("G", "MP", "G + MP", "MPH", "HPH", "MPH + HPH"))

library(tidyr)
library(ggplot2)



p = ggplot(df1, aes(model, y=M, fill=trait)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=round(M,2)  ), hjust=3, color="white",
            position = position_dodge(0.9), angle = 90,size=3.5)+
  geom_errorbar(aes(ymin=M, ymax=M+SD), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  facet_grid(~trait)+
  scale_y_continuous("Prediction ability")+
  xlab("Agronomic traits") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
  
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=5, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=12, face="bold", colour = "black"),    
    axis.title.y = element_text(size=12, face="bold", colour = "black"),    
    axis.text.x = element_text(size=7, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=7, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 8, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 8, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    legend.position = "none"
    #axis.text.x.bottom = element_blank()
  )

jpeg("Pred.ability.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()

df1 = read.csv("Pred.ability.11models_3CV.csv")

df1 = df1 %>%
  mutate(split_file_name = str_split(file_name, "_", simplify = TRUE)) %>%
  mutate(trait = split_file_name[, 2],
         NIR = split_file_name[, 3],
         CV_scheme = split_file_name[, 4],
         model = split_file_name[, 5], ".csv") 


## for 11 models 
df1 =
  df1  %>%
  mutate(
    model = case_when(
      model == "1" ~ "G",
      model == "2" ~ "MP_20CS",
      model == "3" ~ "MPH_20CS",
      model == "4" ~ "HPH_20CS",
      model == "5" ~ "MPH + HPH_20CS",
      model == "6" ~ "G + MP_20CS",
      model == "7" ~ "MP_20LY",
      model == "8" ~ "MPH_20LY",
      model == "9" ~ "HPH_20LY",
      model == "10" ~ "MPH + HPH_20LY",
      model == "11" ~ "G + MP_20LY",
      
      TRUE~model),
    trait = case_when(
      trait == "yield" ~ "Grain Yield",
      trait == "ph" ~ "Plant Height",
      trait == "da" ~ "Days to Anthesis",
      TRUE~trait),
    CV_scheme = case_when(
      CV_scheme == "CV1" ~ "single_trait",
      CV_scheme == "CV2" ~ "GY + PH",
      CV_scheme == "CV3" ~ "DA + PH"
    )
  )

df2 = df1[,c(3,4,6,7,8,9)]


df1 <- as.data.frame(df2 %>%  dplyr::group_by(trait,NIR, CV_scheme, model) %>% 
                       dplyr::summarise(M = mean(cor, na.rm=TRUE),
                                        SD = sd(cor, na.rm=TRUE)))
#df1 = read.csv("df1.csv")
df1$trait <- factor(df1$trait, levels =  c("Grain Yield", "Plant Height", "Days to Anthesis"))
df1$model <- factor(df1$model, levels =  c("G", "MP_20CS", "G + MP_20CS", "MPH_20CS", "HPH_20CS", "MPH + HPH_20CS",
                                           "MP_20LY", "G + MP_20LY", "MPH_20LY", "HPH_20LY", "MPH + HPH_20LY"))

df1$CV_scheme =  factor(df1$CV_scheme, levels =  c("single_trait",
                                               "GY + PH",
                                               "DA + PH"))

library(tidyr)
library(ggplot2)



p = ggplot(df1, aes(model, y=M, fill=trait)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=round(M,2)  ), hjust=3, color="white",
            position = position_dodge(1.2), angle = 90,size=5)+ #size = 3.5, position_dodge(0.9) 
  geom_errorbar(aes(ymin=M, ymax=M+SD), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  facet_grid(~trait~NIR~CV_scheme)+
  scale_y_continuous("Prediction ability")+
  xlab("Agronomic traits") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
  
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=5, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=12, face="bold", colour = "black"),    
    axis.title.y = element_text(size=12, face="bold", colour = "black"),    
    axis.text.x = element_text(size=7, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=7, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 8, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 8, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    legend.position = "none"
    #axis.text.x.bottom = element_blank()
  )

jpeg("Pred.ability_combined.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()


jpeg("Pred.ability_combined.jpeg",width = 9,height =6,units = "in", res=600)
p
dev.off()



#############DPAC
library(adegenet)
dapc1 = dapc(pheno_combined[,-c(1:7)], pheno_combined$env)
p = scatter(dapc1)

jpeg("dapc1.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()

inbreds_combined = rbind(parents_20CS, parents_20MA)
dapc2 = dapc(inbreds_combined[,-c(1:5)], inbreds_combined$env)
scatter(dapc2)






jpeg("pred_acc.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()



#### plot the output ###
## run machine learning models 
## do dapc





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
# if this didnot improve, i might add them as covariate the way daniel did in his paper

