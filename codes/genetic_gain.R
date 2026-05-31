## Using methods from Wang et al., 2025 The Plant Phenome ####################################################
## Codes for assessing breeder's equation ####################################################################
library(dplyr)
library(prospectr)
library(BGLR)

#### reading pheno data for hybrids
pheno_19CS = read.csv("data/blues_19CS_hybrids_spectra.csv", check.names = FALSE)
pheno_19TA = read.csv("data/blues_19TA_hybrids_spectra.csv", check.names = FALSE)
pheno_20CS = read.csv("data/blues_20CS_hybrids_spectra.csv", check.names = FALSE)
pheno_20LY = read.csv("data/blues_20MA_hybrids_spectra.csv", check.names = FALSE)

#### reading pheno data for inbreds
parents_20CS = read.csv("data/blues_20CS_inbreds_spectra.csv", check.names = FALSE)
parents_20LY = read.csv("data/blues_20MA_inbreds_spectra.csv", check.names = FALSE)

pheno_combined = rbind(pheno_19CS,pheno_19TA, pheno_20CS, pheno_20LY)
ne <- as.vector(table(pheno_combined$env)) ## counting the number of observations
ne

#### loading marker data 
T92I012 <- read.table(file = "data/T92I012_clean.txt", header = TRUE)
parents  = read.delim(file = "data/parents.txt", header = TRUE, sep = "\t", dec = ".")$pedigree

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

#### added 

# K3=kronecker(K1,K2)
# dim(K3)
# 
# 
# ################################
# K3.data = as.data.frame(K3)
# 
# 
# complete.list = row.names(G)
# numF = dim(K1)[1]
# Fnames = complete.list[1:numF]
# numM= dim(K2)[1]
# Mnames = complete.list[numF+1:numM]
# 
# total.hybrids = c()
# #loop
# for (f in Fnames) {
#   for (m in Mnames) {
#     total.hybrids = c(total.hybrids, paste(f,"/",m, sep = ""))
#   }
# }
# 
# # 
# # temp = as.data.frame(matrix(0, ncol = 1980, nrow = 1980))
# # colnames(temp) <- total.hybrids
# # rownames(temp) <- total.hybrids
# 
# available.hybrid = unique(pheno_combined$pedigree)
# # 
# # for(av.hyb1 in available.hybrid){
# #   for(av.hyb2 in available.hybrid){
# #     temp[av.hyb1, av.hyb2] = 1
# #   }
# # }
# 
# dim(K3.data)
# # dim(temp)
# 
# 
# # tmpf.k3 = temp*K3.data
# 
# 
# final.k3 = as.data.frame(matrix(1, ncol = 366, nrow = 366))
# colnames(final.k3) <- available.hybrid 
# rownames(final.k3) <- available.hybrid
# 
# 
# colnames(K3.data) <- total.hybrids
# rownames(K3.data) <- total.hybrids
# 
# 
# for(av.hyb1 in available.hybrid){
#   for(av.hyb2 in available.hybrid){
#     # print( K3.data[av.hyb1, av.hyb2])
#     final.k3[av.hyb1, av.hyb2] = K3.data[av.hyb1, av.hyb2]
#   }
# }
# 
# final.k3
# K3 = as.matrix(final.k3)
# 
# save(K3, file = "K3.Rdata")


###

load("data/K3.Rdata") # loading hybrid relationship created via kronecker product


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
rownames(pheno_20LY) = pheno_20LY$pedigree

NIR_19CS = scale(savitzkyGolay(pheno_19CS[,-c(1:7)], m=1, p=1, w=11))
NIR_19TA = scale(savitzkyGolay(pheno_19TA[,-c(1:7)], m=1, p=1, w=11))
NIR_20CS = scale(savitzkyGolay(pheno_20CS[,-c(1:7)], m=1, p=1, w=11))
NIR_20LY = scale(savitzkyGolay(pheno_20LY[,-c(1:7)], m=1, p=1, w=11))

NIR.d1 = rbind(NIR_19CS, NIR_19TA, NIR_20CS, NIR_20LY) 
ZN1 = tcrossprod(as.matrix(NIR.d1)/ncol(as.matrix(NIR.d1))) #phenomic relationship matrices
dim(ZN1) #phenomic relationship matrices from hybrids

midparent_20CS = read.csv("data/mid_parent_NIR_20CS_all.csv", check.names = FALSE)
midparent_20LY = read.csv("data/mid_parent_NIR_20MA_all.csv", check.names = FALSE)

rownames(midparent_20CS) = midparent_20CS[,1]
rownames(midparent_20LY) = midparent_20LY[,1]
midparent_20CS = midparent_20CS[,-1]
midparent_20LY = midparent_20LY[,-1]

# Ensure the row names of averages_df are accessible
# averages_pedigrees <- rownames(averages_df) this one is when calculating using above loop

pedigrees_list = pheno_combined$pedigree
filtered_df_20CS <- midparent_20CS[rownames(midparent_20CS) %in% pedigrees_list, ]
final_pedigree_20CS <- filtered_df_20CS[match(pedigrees_list, rownames(filtered_df_20CS)), ]

filtered_df_20LY <- midparent_20LY[rownames(midparent_20LY) %in% pedigrees_list, ]
final_pedigree_20LY <- filtered_df_20LY[match(pedigrees_list, rownames(filtered_df_20LY)), ]

mid_parent_19CS_20CS = final_pedigree_20CS[1:364,]
mid_parent_19TA_20CS = final_pedigree_20CS[365:543,]
mid_parent_20CS_20CS = final_pedigree_20CS[544:907,]
mid_parent_20LY_20CS = final_pedigree_20CS[-c(1:907),]

mid_parent_19CS_20LY = final_pedigree_20LY[1:364,]
mid_parent_19TA_20LY = final_pedigree_20LY[365:543,]
mid_parent_20CS_20LY = final_pedigree_20LY[544:907,]
mid_parent_20LY_20LY = final_pedigree_20LY[-c(1:907),]

## estimate first derivative of mid-parent NIR spectra
# using 20CS

NIR_mid_parent_19CS_20CS = scale(savitzkyGolay(mid_parent_19CS_20CS, m=1, p=1, w=11))
NIR_mid_parent_19TA_20CS = scale(savitzkyGolay(mid_parent_19TA_20CS, m=1, p=1, w=11))
NIR_mid_parent_20CS_20CS = scale(savitzkyGolay(mid_parent_20CS_20CS, m=1, p=1, w=11))
NIR_mid_parent_20LY_20CS = scale(savitzkyGolay(mid_parent_20LY_20CS, m=1, p=1, w=11))

NIR_mid_parent.d1_20CS = rbind(NIR_mid_parent_19CS_20CS, 
                               NIR_mid_parent_19TA_20CS, 
                               NIR_mid_parent_20CS_20CS, 
                               NIR_mid_parent_20LY_20CS)

NIR_mid_parent.ZN1_20CS = tcrossprod(as.matrix(NIR_mid_parent.d1_20CS)/ncol(as.matrix(NIR_mid_parent.d1_20CS))) #phenomic relationship matrices
dim(NIR_mid_parent.ZN1_20CS) #phenomic relationship matrices from hybrids

# using 20LY
NIR_mid_parent_19CS_20LY = scale(savitzkyGolay(mid_parent_19CS_20LY, m=1, p=1, w=11))
NIR_mid_parent_19TA_20LY = scale(savitzkyGolay(mid_parent_19TA_20LY, m=1, p=1, w=11))
NIR_mid_parent_20CS_20LY = scale(savitzkyGolay(mid_parent_20CS_20LY, m=1, p=1, w=11))
NIR_mid_parent_20LY_20LY = scale(savitzkyGolay(mid_parent_20LY_20LY, m=1, p=1, w=11))

NIR_mid_parent.d1_20LY = rbind(NIR_mid_parent_19CS_20LY, 
                               NIR_mid_parent_19TA_20LY, 
                               NIR_mid_parent_20CS_20LY, 
                               NIR_mid_parent_20LY_20LY)

NIR_mid_parent.ZN1_20LY = tcrossprod(as.matrix(NIR_mid_parent.d1_20LY)/ncol(as.matrix(NIR_mid_parent.d1_20LY))) #phenomic relationship matrices
dim(NIR_mid_parent.ZN1_20LY) #phenomic relationship matrices from hybrids

### Mid-parent heterosis ###
# heterosis = F1-MP/MP

combined_NIR_F1 = rbind(pheno_19CS[,-c(1:7)],  #hybrid 
                        pheno_19TA[,-c(1:7)], 
                        pheno_20CS[,-c(1:7)],
                        pheno_20LY[,-c(1:7)])

combined_NIR_MP = final_pedigree_20CS # mid-parent
MP_heterosis_20CS = (combined_NIR_F1- combined_NIR_MP)/combined_NIR_MP
MP_heterosis_20CS

MP_heterosis_20LY = (combined_NIR_F1- final_pedigree_20LY)/final_pedigree_20LY
MP_heterosis_20LY


MP_heterosis_19CS_20CS = MP_heterosis_20CS[1:364,]
MP_heterosis_19TA_20CS = MP_heterosis_20CS[365:543,]
MP_heterosis_20CS_20CS = MP_heterosis_20CS[544:907,]
MP_heterosis_20LY_20CS = MP_heterosis_20CS[-c(1:907),]

MP_heterosis_19CS_20LY = MP_heterosis_20LY[1:364,]
MP_heterosis_19TA_20LY = MP_heterosis_20LY[365:543,]
MP_heterosis_20CS_20LY = MP_heterosis_20LY[544:907,]
MP_heterosis_20LY_20LY = MP_heterosis_20LY[-c(1:907),]


## estimate first derivative of mid-parent_heterosis NIR spectra
NIR_MP_heterosis_19CS_20CS = scale(savitzkyGolay(MP_heterosis_19CS_20CS, m=1, p=1, w=11))
NIR_MP_heterosis_19TA_20CS = scale(savitzkyGolay(MP_heterosis_19TA_20CS, m=1, p=1, w=11))
NIR_MP_heterosis_20CS_20CS = scale(savitzkyGolay(MP_heterosis_20CS_20CS, m=1, p=1, w=11))
NIR_MP_heterosis_20LY_20CS = scale(savitzkyGolay(MP_heterosis_20LY_20CS, m=1, p=1, w=11))

NIR_MP_heterosis.d1_20CS = rbind(NIR_MP_heterosis_19CS_20CS, #combine all
                                 NIR_MP_heterosis_19TA_20CS, 
                                 NIR_MP_heterosis_20CS_20CS, 
                                 NIR_MP_heterosis_20LY_20CS) 

NIR_MP_heterosis.ZN1_20CS = tcrossprod(as.matrix(NIR_MP_heterosis.d1_20CS)/ncol(as.matrix(NIR_MP_heterosis.d1_20CS))) #phenomic relationship matrices
dim(NIR_MP_heterosis.ZN1_20CS) #phenomic relationship matrices from NIR_mid_parent_heterosis

#20LY
NIR_MP_heterosis_19CS_20LY = scale(savitzkyGolay(MP_heterosis_19CS_20LY, m=1, p=1, w=11))
NIR_MP_heterosis_19TA_20LY = scale(savitzkyGolay(MP_heterosis_19TA_20LY, m=1, p=1, w=11))
NIR_MP_heterosis_20CS_20LY = scale(savitzkyGolay(MP_heterosis_20CS_20LY, m=1, p=1, w=11))
NIR_MP_heterosis_20LY_20LY = scale(savitzkyGolay(MP_heterosis_20LY_20LY, m=1, p=1, w=11))

NIR_MP_heterosis.d1_20LY = rbind(NIR_MP_heterosis_19CS_20LY, #combine all
                                 NIR_MP_heterosis_19TA_20LY, 
                                 NIR_MP_heterosis_20CS_20LY, 
                                 NIR_MP_heterosis_20LY_20LY) 

NIR_MP_heterosis.ZN1_20LY = tcrossprod(as.matrix(NIR_MP_heterosis.d1_20LY)/ncol(as.matrix(NIR_MP_heterosis.d1_20LY))) #phenomic relationship matrices
dim(NIR_MP_heterosis.ZN1_20LY) #phenomic relationship matrices from NIR_mid_parent_heterosis

# calculating NIR x E interaction for different models using Hadamard product
# 20CS
mid_parent.ZNZE1.CS             = NIR_mid_parent.ZN1_20CS * ZEZEt
mid_parent_heterosis_ZNZE1.CS   = NIR_MP_heterosis.ZN1_20CS * ZEZEt

# calculating NIR x E interaction for different models using Hadamard product
#20LY

mid_parent.ZNZE1.LY             = NIR_mid_parent.ZN1_20LY * ZEZEt
mid_parent_heterosis_ZNZE1.LY   = NIR_MP_heterosis.ZN1_20LY * ZEZEt

# Set ETAs predictors for models
# Genomic model
Eta1_CS <- list(list(X = ZE, model = "BRR"),     # Env
             list(K = K1star, model = "RKHS"), # Female
             list(K = K2star, model = "RKHS"), # male
             list(K = K3star, model = "RKHS"), # Femalexmale
             list(K = K4, model = "RKHS"), # Female x environment
             list(K = K5, model = "RKHS"), # male x environemnt
             list(K = K6, model = "RKHS")) #female x male x environment

##### first derivative of mid_parent NIR
Eta2_CS <-list(list(X = ZE, model = "BRR"),      #Env
           list(K= NIR_mid_parent.ZN1_20CS, model = "RKHS"),    #NIR1
           list(K = mid_parent.ZNZE1.CS, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent heterosis NIR
Eta3_CS <-list(list(X = ZE,model="BRR"),      #Env
           list(K = NIR_MP_heterosis.ZN1_20CS, model="RKHS"),    #NIR1
           list(K = mid_parent_heterosis_ZNZE1.CS, model="RKHS"))  #NIR1 x Env

##### first derivative of mid_parent NIR + genomic
Eta4_CS <- list(list(X = ZE, model = "BRR"),     # Env
             list(K = K1star, model = "RKHS"), # Female
             list(K = K2star, model = "RKHS"), # male
             list(K = K3star, model = "RKHS"), #female xmale 
             list(K = K4, model = "RKHS"), #female x env
             list(K = K5, model = "RKHS"), #male x env
             list(K = K6, model = "RKHS"), #female x male x env
             list(K = NIR_mid_parent.ZN1_20CS, model="RKHS"),    #NIR1_mid_parent
             list(K = mid_parent.ZNZE1.CS, model="RKHS")) #NIR x env


####Use 20LY

# Genomic model
Eta1_LY <- list(list(X = ZE, model = "BRR"),     # Env
                list(K = K1star, model = "RKHS"), # Female
                list(K = K2star, model = "RKHS"), # male
                list(K = K3star, model = "RKHS"), # Femalexmale
                list(K = K4, model = "RKHS"), # Female x environment
                list(K = K5, model = "RKHS"), # male x environemnt
                list(K = K6, model = "RKHS")) #female x male x environment


##### first derivative of mid_parent NIR
Eta2_LY <- list(list(X = ZE, model = "BRR"),      #Env
           list(K= NIR_mid_parent.ZN1_20LY, model = "RKHS"),    #NIR1
           list(K = mid_parent.ZNZE1.LY, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent heterosis NIR
Eta3_LY <- list(list(X = ZE,model="BRR"),      #Env
           list(K = NIR_MP_heterosis.ZN1_20LY, model="RKHS"),    #NIR1
           list(K = mid_parent_heterosis_ZNZE1.LY, model="RKHS"))  #NIR1 x Env

##### first derivative of mid_parent NIR + genomic 20LY
Eta4_LY <- list(list(X = ZE, model = "BRR"),     # Env
                list(K = K1star, model = "RKHS"), # Female
                list(K = K2star, model = "RKHS"), # male
                list(K = K3star, model = "RKHS"), #female xmale 
                list(K = K4, model = "RKHS"), #female x env
                list(K = K5, model = "RKHS"), #male x env
                list(K = K6, model = "RKHS"), #female x male x env
                list(K=NIR_mid_parent.ZN1_20LY, model="RKHS"),    #NIR1_mid_parent
                list(K=mid_parent.ZNZE1.LY, model="RKHS")) #NIR x env


### 4 different models 
# MODEL 1: Genomic model
# MODEL 2: Mid-parent model
# MODEL 3: mid_parent heterosis
# MODEL 4: high_parent heterosis
# MODEL 5: mid_parent + high-parent heterosis
# MODEL 6: first derivative of mid_parent NIR + genomic

Models <- list(Eta1_CS, Eta2_CS, Eta3_CS, Eta4_CS)
traitnames <- c("yield", "da", "ph", "starch", "protein", "fat", "fiber", "ash")

pheno_combined[1:10,1:10]

source("codes/Estimate_gcor_prediction.R")
# tr=1
#set.seed(1234)
for (tr in 1:length(traitnames)) {
  
  if (tr == 1) {
    pheno <- read.csv(file = paste("data/pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==2) {
    pheno <- read.csv(file = paste("data/pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==3) {
    pheno <- read.csv(file = paste("data/pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==4) {
    pheno <- read.csv(file = paste("data/pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==5) {
    pheno <- read.csv(file = paste("data/pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==6) {
    pheno <- read.csv(file = paste("data/pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==7) {
    pheno <- read.csv(file = paste("data/pheno_", traitnames[tr], ".csv", sep = ""))
  } else {
    pheno <- read.csv(file = paste("data/pheno_", traitnames[tr], ".csv", sep = ""))
  }
  
  
  # cross-validation #
  parents = c(unique(pheno$female), unique(pheno$male))
  hybrid = as.character(unique(pheno$pedigree))
  Phenotype_data1 = pheno
  
  set.seed(123)
  cycles = 10
  CV1 = list()
  CV2 = list()
  #CV3_20LY = list()
  #CV3_19TA = list()
  #CV4 = list()
  
  #MODEL =1; rep_num=1
  for (MODEL in 1:length(Models)) {  
    
    for (rep_num in 1:5) {
      #set.seed(123)
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
      CV1_adjusted = mean(CV1[[1]]$cor)
      h2 = fit$ETA[[1]]$varB / var(CV_Data_1_2$blue)
      Acc = CV1_adjusted/h2
      estimate_gcor(data = df_test, Knn = df_test, method = "MCMCglmm", normalize = F)[1]


      # Preparing for CV2
      # CV2 simulates sparse testing where 30% of hybrids were sampled from each environments
      
      test_geno = sample(unique(pedigrees_list), round(length(unique(pedigrees_list))*0.3))
      train_geno = setdiff(hybrid, test_geno)
      
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
      # CV_Data_1_2<-Phenotype_data1
      # CV_Data_1_2$Y = CV_Data_1_2$blue
      # 
      # CV_Data_1_2$Y[CV_Data_1_2$env == "20LY"] <- NA
      # 
      # y_t2<-as.numeric(CV_Data_1_2$Y)
      # 
      # fit2<-BGLR(y=y_t2,
      #            ETA=Models[[MODEL]],
      #            nIter=5000,
      #            burnIn=1000, 
      #            thin=10) #nIter=5000,burnIn=1000, thin =10
      # 
      # CV_Data_1_2$yhat2 <- fit2$yHat
      # 
      # df_test2 <- subset(CV_Data_1_2, env == "20LY")
      # CV3_20LY[[(rep_num)]] <- as.data.frame(df_test2 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat2,use = "complete.obs")))
      # 
      # # For CV3_19TA
      # # leave environments out #19TA
      # CV_Data_1_2<-Phenotype_data1
      # CV_Data_1_2$Y = CV_Data_1_2$blue
      # 
      # CV_Data_1_2$Y[CV_Data_1_2$env == "19TA"] <- NA
      # 
      # y_t3<-as.numeric(CV_Data_1_2$Y)
      # 
      # fit3<-BGLR(y=y_t3,
      #            ETA=Models[[MODEL]],
      #            nIter=5000,
      #            burnIn=1000, thin=10) #nIter=5000,burnIn=1000, thin =10
      # 
      # CV_Data_1_2$yhat3 <- fit3$yHat
      # 
      # df_test3 <- subset(CV_Data_1_2, CV_Data_1_2$env == "19TA")
      # CV3_19TA[[(rep_num)]] <- as.data.frame(df_test3 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat3,use = "complete.obs")))
    }
    
    #rep_num =1
    if (rep_num == 5) {
      CV1out <- plyr::ldply(CV1, data.frame)
      CV2out <- plyr::ldply(CV2, data.frame)
      #CV3out_20LY <- plyr::ldply(CV3_20LY, data.frame)
      #CV3out_19TA <- plyr::ldply(CV3_19TA, data.frame)
      write.csv(CV1out, file = paste("ACC_", traitnames[tr],"_CV1_20CS_", MODEL, ".csv", sep=""), row.names = FALSE)
      write.csv(CV2out, file = paste("ACC_", traitnames[tr],"_CV2_20CS_", MODEL, ".csv", sep=""), row.names = FALSE)
      #write.csv(CV3out_20LY, file = paste("ACC_", traitnames[tr],"_CV3_20CS_20LY_", MODEL, ".csv", sep=""), row.names = FALSE)
      #write.csv(CV3out_19TA, file = paste("ACC_", traitnames[tr],"_CV3_20CS_19TA_", MODEL, ".csv", sep=""), row.names = FALSE)
    }
  }
}


