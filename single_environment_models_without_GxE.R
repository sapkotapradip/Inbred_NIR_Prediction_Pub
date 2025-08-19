####################### Within-Environment Prediction ##################################################
########################################################################################################
# looping for single environment trials #####

library(dplyr)
library(prospectr)
library(BGLR)

traitnames <- c("yield", "da", "ph", "starch", "protein", "fat", "fiber", "ash")
envnames <- c("19CS", "19TA", "20CS", "20LY")
inbreds <- c("20CS", "20LY")

# tr = 1; env =1; inb = 1

for (tr in 1:length(traitnames)){
  for (env in 1:length(envnames)) {
    for(inb in 1: length(inbreds))
    tryCatch({
      cat("Processing environment:", envnames[env], "and trait:", traitnames[tr], "\n")
      
      if (tr == 1) {
        pheno <- read.csv("pheno_yield.csv")
      } else if (tr == 2) {
        pheno <- read.csv("pheno_da.csv")
      } else if (tr == 3) {
        pheno <- read.csv("pheno_ph.csv")
      } else if (tr == 4) {
        pheno <- read.csv("pheno_starch.csv")
      } else if (tr == 5) {
        pheno <- read.csv("pheno_protein.csv")
      } else if (tr == 6) {
        pheno <- read.csv("pheno_fat.csv")
      } else if (tr == 7) {
        pheno <- read.csv("pheno_fiber.csv")
      } else {
        pheno <- read.csv("pheno_ash.csv")
      }
      
      cat("Phenotype data dimensions for env:", envnames[env], "are", dim(pheno), "\n")
      
      Phenotype_data1 = pheno[pheno$env == envnames[env],]
      
      NIR_hybrids = read.csv(file = paste0("blues_",
                                          envnames[env],
                                          "_hybrids_spectra.csv"), 
                                          check.names = FALSE)
      
      NIR_inbreds =  read.csv(file = paste0("blues_",
                                            inbreds[inb],
                                            "_inbreds_spectra.csv"), 
                              check.names = FALSE)
      
      NIR_mid_parent = read.csv(file = paste0("mid_parent_NIR_",
                                                   inbreds[inb],
                                                   "_all.csv"), 
                                     check.names = FALSE)
      
      rownames(NIR_mid_parent) = NIR_mid_parent[,1]
      mid_parent_ = NIR_mid_parent[,-1]
      
      # Ensure the row names of averages_df are accessible
      # averages_pedigrees <- rownames(averages_df) this one is when calculating using above loop
      pheno_combined = read.csv("pheno_yield.csv" )
      
      pedigrees_list = pheno_combined$pedigree
      filtered_df_ <- NIR_mid_parent[rownames(NIR_mid_parent) %in% pedigrees_list, ]
      final_pedigree_ <- filtered_df_[match(pedigrees_list, rownames(filtered_df_)), ]
      
      final_pedigree_ = final_pedigree_[,-1]
      
      mid_parent_19CS = final_pedigree_[1:364,]
      mid_parent_19TA = final_pedigree_[365:543,]
      mid_parent_20CS = final_pedigree_[544:907,]
      mid_parent_20LY = final_pedigree_[-c(1:907),]
      
      ## estimate first derivative of mid-parent NIR spectra
      # using 20CS
      colnames = envnames[env]
      
      mid_parent_ [[paste(envnames[env], sep = "_")]]
      
      NIR_mid_parent_ = scale(savitzkyGolay(get(paste0("mid_parent_", envnames[env])), m=1, p=1, w=11))
  
      NIR_mid_parent.ZN1_ = tcrossprod(as.matrix(NIR_mid_parent_)/ncol(as.matrix(NIR_mid_parent_))) #phenomic relationship matrices
      dim(NIR_mid_parent.ZN1_) #phenomic relationship matrices from hybrids
      
      
      ### Mid-parent heterosis ###
      # heterosis = F1-MP/MP
      
      pheno_19CS = read.csv("blues_19CS_hybrids_spectra.csv", check.names = FALSE)
      pheno_19TA = read.csv("blues_19TA_hybrids_spectra.csv", check.names = FALSE)
      pheno_20CS = read.csv("blues_20CS_hybrids_spectra.csv", check.names = FALSE)
      pheno_20LY = read.csv("blues_20MA_hybrids_spectra.csv", check.names = FALSE)
      pheno_combined_ = rbind(pheno_19CS,pheno_19TA, pheno_20CS, pheno_20LY)
      
      combined_NIR_F1 = pheno_combined_[,-c(1:7)]
      NIR_19CS = pheno_19CS[,-c(1:7)]
      NIR_19TA = pheno_19TA[,-c(1:7)]
      NIR_20CS = pheno_20CS[,-c(1:7)]
      NIR_20LY = pheno_20LY[,-c(1:7)]
      
      
      MPH = (combined_NIR_F1- final_pedigree_)/final_pedigree_
      
      MPH_19CS = MPH[1:364,]
      MPH_19TA = MPH[365:543,]
      MPH_20CS = MPH[544:907,]
      MPH_20LY = MPH[-c(1:907),]
      
      NIR_MPH = scale(savitzkyGolay(get(paste0("MPH_", envnames[env])), m=1, p=1, w=11))
      
      NIR_MPH.ZN1_ = tcrossprod(as.matrix(NIR_MPH)/ncol(as.matrix(NIR_MPH))) #phenomic relationship matrices
      dim(NIR_MPH.ZN1_) 
      
      # High parent heterosis
      
      NIR_parents = NIR_inbreds[,-c(2:5)]
      rownames(NIR_parents) = NIR_parents$pedigree
      NIR_parents = NIR_parents[,-1]
      F_ = NIR_parents[1:45,] #female parents
      M_ = NIR_parents[46:89,] #male parents
      F_$female = rownames(F_)
      M_$male = rownames(M_)
      
      hybrids = read.csv("total.hybrids.csv")
      
      females_hybrids <- hybrids %>%
        inner_join(F_, by = c("female"="female"))
      
      males_hybrids <- hybrids %>%
        inner_join(M_, by = c("male"="male"))
      
      
      females_numeric = females_hybrids[,-c(1:3)]
      males_numeric = males_hybrids[,-c(1:3)]
      
      
      High_parent_values_ = pmax(females_numeric, males_numeric)
      High_parent_heterosis_ = ( combined_NIR_F1 - High_parent_values_)/High_parent_values_
      
      HPH_19CS = High_parent_heterosis_[1:364,]
      HPH_19TA = High_parent_heterosis_[365:543,]
      HPH_20CS = High_parent_heterosis_[544:907,]
      HPH_20MA = High_parent_heterosis_[-c(1:907),]
      
      NIR_HPH = scale(savitzkyGolay(get(paste0("HPH_", envnames[env])), m=1, p=1, w=11))
      NIR_HPH.ZN1_ = tcrossprod(as.matrix(NIR_HPH)/ncol(as.matrix(NIR_HPH))) #phenomic relationship matrices
      dim(NIR_HPH.ZN1_) 
      
      
      # Phenotype_data1
      
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
      tmp<-strsplit(Phenotype_data1$pedigree,"/")
      
      P1<-rep(NA,length(tmp))
      P2<-rep(NA,length(tmp))
      
      for(i in 1:length(tmp))
      {
        P1[i]=tmp[[i]][1]
        P2[i]=tmp[[i]][2]
      }
      
      P1<-as.factor(P1)
      P2<-as.factor(P2)
      cross<-as.factor(Phenotype_data1$pedigree)
          
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
      K3 = K3[Phenotype_data1$pedigree, Phenotype_data1$pedigree]
      
      #### modeling effects for priors
      
      K1star=Z1%*%K1%*%t(Z1) #female
      K2star=Z2%*%K2%*%t(Z2) #male
      K3star=Z3%*%K3%*%t(Z3) #female * male
      
    
      #### Predictors for models
      
      # Genomic model
      Eta1 <- list(list(K = K1star, model = "RKHS"), # Female
                   list(K = K2star, model = "RKHS"), # male
                   list(K = K3star, model = "RKHS")) # Femalexmale
              
      ##### first derivative of mid_parent NIR
      Eta2<- list(list(K = NIR_mid_parent.ZN1_, model = "RKHS"))
      
      
      ##### first derivative of mid_parent heterosis NIR
      Eta3<- list(list(K = NIR_MPH.ZN1_, model = "RKHS"))
      
      ##### first derivative of high_parent heterosis NIR
      Eta4<- list(list(K = NIR_HPH, model = "RKHS"))
      
      
      #### first derivative of high parent NIR
      Eta5<-list(list(K = NIR_MPH.ZN1_, model = "RKHS"),
                 list(K = NIR_HPH, model = "RKHS"))
      
      
      ##### first derivative of mid_parent + high-parent heterosis NIR
      Eta6 <- list(list(K = K1star, model = "RKHS"), # Female
                   list(K = K2star, model = "RKHS"), # male
                   list(K = K3star, model = "RKHS"),
                   list(K = NIR_mid_parent.ZN1_, model = "RKHS")) # Femalexmale
      
      Models <- list(Eta1, Eta2, Eta3, Eta4, Eta5, Eta6)
      
      hybrid = as.character(unique(Phenotype_data1$pedigree))
      
      # cross-validation #
      parents = c(unique(Phenotype_data1$female), unique(Phenotype_data1$male))
      hybrid = as.character(unique(Phenotype_data1$pedigree))
      
      set.seed(123)
      cycles = 10
      CV1 = list()
      CV2 = list()
      
      #MODEL =2; rep_num=1
      for (MODEL in 1:length(Models)) {  
        
        for (rep_num in 1:5) {
          set.seed(123)
          CVa = sample(parents[1:89], 15, replace = FALSE)
          train_geno <- setdiff(parents[1:89], CVa)
          Phenotype_data1$y = Phenotype_data1$blue
          
          for (a in 1:15) {
            Phenotype_data1 <- Phenotype_data1 %>% mutate(y = replace(y, female == CVa[a], NA))
          }
          
          for (b in 1:15) {
            Phenotype_data1 <- Phenotype_data1 %>% mutate(y = replace(y, male == CVa[b], NA))
          }
          
          test_pheno <- subset(Phenotype_data1, is.na(y))
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
          
          test_geno = sample(unique(hybrid), round(length(hybrid)*0.3))
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
        
      }
        if (rep_num == 5) {
          CV1out <- plyr::ldply(CV1, data.frame)
          CV2out <- plyr::ldply(CV2, data.frame)
          write.csv(CV1out, file = paste("ACC_", envnames[env], "_", traitnames[tr], "_CV1_model", MODEL, ".csv", sep = ""), row.names = F)
          write.csv(CV2out, file = paste("ACC_", envnames[env], "_", traitnames[tr], "_CV2_model", MODEL, ".csv", sep = ""), row.names = F)
        }
      }
    }, error = function(e) {
      cat("Error in processing environment:", envnames[env], "and trait:", traitnames[tr], "\n")
      cat("Error message:", e$message, "\n")
    })
  }
}

















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

eta = NIR_mid_parent.ZN1_20CS[1:364,1:364]

# Set ETAs predictors for models

# Genomic model
Eta1_CS <- list(list(X = eta, model = "RKHS"))

##### first derivative of mid_parent NIR
Eta2_CS <-list(list(X = ZE, model = "BRR"),      #Env
               list(K= NIR_mid_parent.ZN1_20CS, model = "RKHS"),    #NIR1
               list(K = mid_parent.ZNZE1.CS, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent heterosis NIR
Eta3_CS <-list(list(X = ZE,model="BRR"),      #Env
               list(K = NIR_MP_heterosis.ZN1_20CS, model="RKHS"),    #NIR1
               list(K = mid_parent_heterosis_ZNZE1.CS, model="RKHS"))  #NIR1 x Env

##### first derivative of high_parent heterosis NIR
Eta4_CS <-list(list(X = ZE,model = "BRR"),      #Env
               list(K = NIR_HP_heterosis.ZN1_20CS, model = "RKHS"),    #NIR1
               list(K = high_parent_heterosis_ZNZE1.CS, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent + high-parent heterosis NIR
Eta5_CS <-list(list(X = ZE,model="BRR"),      #Env
               list(K = NIR_MP_heterosis.ZN1_20CS, model="RKHS"),    #NIR1
               list(K = mid_parent_heterosis_ZNZE1.CS, model="RKHS"),
               list(K = NIR_HP_heterosis.ZN1_20CS, model="RKHS"),    #NIR1
               list(K = high_parent_heterosis_ZNZE1.CS, model="RKHS")) #NIR1 x Env


##### first derivative of mid_parent NIR + genomic
Eta6_CS <- list(list(X = ZE, model = "BRR"),     # Env
                list(K = K1star, model = "RKHS"), # Female
                list(K = K2star, model = "RKHS"), # male
                list(K = K3star, model = "RKHS"), #female xmale 
                list(K = K4, model = "RKHS"), #female x env
                list(K = K5, model = "RKHS"), #male x env
                list(K = K6, model = "RKHS"), #female x male x env
                list(K = NIR_mid_parent.ZN1_20CS, model="RKHS"),    #NIR1_mid_parent
                list(K = mid_parent.ZNZE1.CS, model="RKHS")) #NIR x env


####Use 20MA

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
                list(K= NIR_mid_parent.ZN1_20MA, model = "RKHS"),    #NIR1
                list(K = mid_parent.ZNZE1.MA, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent heterosis NIR
Eta3_LY <- list(list(X = ZE,model="BRR"),      #Env
                list(K = NIR_MP_heterosis.ZN1_20MA, model="RKHS"),    #NIR1
                list(K = mid_parent_heterosis_ZNZE1.MA, model="RKHS"))  #NIR1 x Env

##### first derivative of mid_parent + high-parent heterosis NIR
Eta4_LY <- list(list(X = ZE,model="BRR"),      #Env
                list(K = NIR_MP_heterosis.ZN1_20MA, model="RKHS"),    #NIR1
                list(K = mid_parent_heterosis_ZNZE1.MA, model="RKHS"),
                list(K = NIR_HP_heterosis.ZN1_20MA, model="RKHS"),    #NIR1
                list(K = high_parent_heterosis.ZNZE1.MA, model="RKHS")) #NIR1 x Env


#### first derivative of high parent NIR
Eta5_LY <- list(list(X = ZE,model = "BRR"),      #Env
                list(K = NIR_HP_ZN1_20MA, model = "RKHS"),    #NIR1
                list(K = high_parent_ZNZE1.MA, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent NIR + genomic 20MA
Eta6_LY <- list(list(X = ZE, model = "BRR"),     # Env
                list(K = K1star, model = "RKHS"), # Female
                list(K = K2star, model = "RKHS"), # male
                list(K = K3star, model = "RKHS"), #female xmale 
                list(K = K4, model = "RKHS"), #female x env
                list(K = K5, model = "RKHS"), #male x env
                list(K = K6, model = "RKHS"), #female x male x env
                list(K=NIR_mid_parent.ZN1_20MA, model="RKHS"),    #NIR1_mid_parent
                list(K=mid_parent.ZNZE1.MA, model="RKHS")) #NIR x env


### 6 different models 
# MODEL 1: Genomic model
# MODEL 2: Mid-parent model
# MODEL 3: mid_parent heterosis
# MODEL 4: high_parent heterosis
# MODEL 5: mid_parent + high-parent heterosis
# MODEL 6: first derivative of mid_parent NIR + genomic

Models <- list(Eta1_CS, Eta2_CS, Eta3_CS, Eta4_CS, Eta5_CS, Eta6_CS)
traitnames <- c("yield", "da", "ph", "starch", "protein", "fat", "fiber", "ash")
traitnames <- c("fiber", "ash")


pheno_combined[1:10,1:10]


# tr=1
#set.seed(1234)
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
  #CV3_20LY = list()
  #CV3_19TA = list()
  #CV4 = list()
  
  #MODEL =6; rep_num=1
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

# started on 2:14 PM 5/12

### 20MA
Models <- list(Eta1_LY, Eta2_LY, Eta3_LY, Eta4_LY, Eta5_LY, Eta6_LY)
traitnames <- c("yield", "da", "ph", "starch", "protein", "fat", "fiber", "ash")
pheno_combined[1:10,1:10]


# tr=1
#set.seed(1234)
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
  
  # set.seed(123)
  cycles = 10
  CV1 = list()
  CV2 = list()
  # CV3_20LY = list()
  # CV3_19TA = list()
  # CV4 = list()
  # 
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
      
      # CV_Data_1_2$Y[CV_Data_1_2$env == "20LY"] <- NA
      
      # y_t2<-as.numeric(CV_Data_1_2$Y)
      
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
      # For CV3_19TA
      # leave environments out #19TA
      #   CV_Data_1_2<-Phenotype_data1
      #   CV_Data_1_2$Y = CV_Data_1_2$blue
      #   
      #   CV_Data_1_2$Y[CV_Data_1_2$env == "19TA"] <- NA
      #   
      #   y_t3<-as.numeric(CV_Data_1_2$Y)
      #   
      #   fit3<-BGLR(y=y_t3,
      #              ETA=Models[[MODEL]],
      #              nIter=5000,
      #              burnIn=1000, thin=10) #nIter=5000,burnIn=1000, thin =10
      #   
      #   CV_Data_1_2$yhat3 <- fit3$yHat
      #   
      #   df_test3 <- subset(CV_Data_1_2, CV_Data_1_2$env == "19TA")
      #   CV3_19TA[[(rep_num)]] <- as.data.frame(df_test3 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat3,use = "complete.obs")))
    }
    
    #rep_num =1
    if (rep_num == 5) {
      CV1out <- plyr::ldply(CV1, data.frame)
      CV2out <- plyr::ldply(CV2, data.frame)
      #CV3out_20LY <- plyr::ldply(CV3_20LY, data.frame)
      #CV3out_19TA <- plyr::ldply(CV3_19TA, data.frame)
      write.csv(CV1out, file = paste("ACC_", traitnames[tr],"_CV1_20LY_", MODEL, ".csv", sep=""), row.names = FALSE)
      write.csv(CV2out, file = paste("ACC_", traitnames[tr],"_CV2_20LY_", MODEL, ".csv", sep=""), row.names = FALSE)
      #write.csv(CV3out_20LY, file = paste("ACC_", traitnames[tr],"_CV3_20LY_20LY_", MODEL, ".csv", sep=""), row.names = FALSE)
      #write.csv(CV3out_19TA, file = paste("ACC_", traitnames[tr],"_CV3_20LY_19TA_", MODEL, ".csv", sep=""), row.names = FALSE)
    }
  }
}
