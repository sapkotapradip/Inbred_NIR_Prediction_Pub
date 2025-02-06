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
F = NIR_parents_20CS[1:45,] #female parents
M = NIR_parents_20CS[46:89,] #male parents


## cakculating averages 

# library(dplyr)
# # available.hybrid
# 
# averages = list()
# # row_F =1; row_M =1
# for (row_F in rownames(F)) {
#   for (row_M in rownames(M)) {
#     # Compute the average for the current pair of rows
#     avg <- (F[row_F, ] + M[row_M, ]) / 2
#     
#     # Store the result in the list with a descriptive name
#     averages[[paste(row_F, row_M, sep = "/")]] <- avg
#   }
# }
# 
# averages
# 
# # Convert the list of averages into a data frame
# averages_df <- do.call(rbind, averages)
# write.csv(averages_df, "mid_parent_NIR.csv")

averages_df = read.csv("mid_parent_NIR_20CS_all.csv", check.names = FALSE)
rownames(averages_df) = averages_df[,1]
averages_df = averages_df[,-1]

# Ensure the row names of averages_df are accessible
# averages_pedigrees <- rownames(averages_df) this one is when calculating using above loop
averages_pedigrees = rownames(averages_df)
pedigrees_list = pheno_combined$pedigree
filtered_df <- averages_df[rownames(averages_df) %in% pedigrees_list, ]
final_pedigree_df <- filtered_df[match(pedigrees_list, rownames(filtered_df)), ]

mid_parent_19CS = final_pedigree_df[1:364,]
mid_parent_19TA = final_pedigree_df[365:543,]
mid_parent_20CS = final_pedigree_df[544:907,]
mid_parent_20MA = final_pedigree_df[-c(1:907),]

## estimate first derivative of mid-parent NIR spectra

NIR_mid_parent_19CS = scale(savitzkyGolay(mid_parent_19CS, m=1, p=1, w=11))
NIR_mid_parent_19TA = scale(savitzkyGolay(mid_parent_19TA, m=1, p=1, w=11))
NIR_mid_parent_20CS = scale(savitzkyGolay(mid_parent_20CS, m=1, p=1, w=11))
NIR_mid_parent_20MA = scale(savitzkyGolay(mid_parent_20LY, m=1, p=1, w=11))

NIR_mid_parent.d1 = rbind(NIR_mid_parent_19CS, NIR_mid_parent_19TA, NIR_mid_parent_20CS, NIR_mid_parent_20MA) 
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
MP_heterosis_20MA = MP_heterosis[-c(1:907),]


## estimate first derivative of mid-parent_heterosis NIR spectra

NIR_MP_heterosis_19CS = scale(savitzkyGolay(MP_heterosis_19CS, m=1, p=1, w=11))
NIR_MP_heterosis_19TA = scale(savitzkyGolay(MP_heterosis_19TA, m=1, p=1, w=11))
NIR_MP_heterosis_20CS = scale(savitzkyGolay(MP_heterosis_20CS, m=1, p=1, w=11))
NIR_MP_heterosis_20MA = scale(savitzkyGolay(MP_heterosis_20MA, m=1, p=1, w=11))

NIR_MP_heterosis.d1 = rbind(NIR_MP_heterosis_19CS, #combine all
                            NIR_MP_heterosis_19TA, 
                            NIR_MP_heterosis_20CS, 
                            NIR_MP_heterosis_20MA) 

NIR_MP_heterosis.ZN1 = tcrossprod(as.matrix(NIR_MP_heterosis.d1)/ncol(as.matrix(NIR_MP_heterosis.d1))) #phenomic relationship matrices
dim(NIR_MP_heterosis.ZN1) #phenomic relationship matrices from NIR_mid_parent_heterosis


######### High parent heterosis 
# F;M
# available.hybrid
# write.csv(pheno_combined$pedigree, "total.hybrids.csv")

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

HP_heterosis_19CS = High_parent_heterosis[1:364,]
HP_heterosis_19TA = High_parent_heterosis[365:543,]
HP_heterosis_20CS = High_parent_heterosis[544:907,]
HP_heterosis_20MA = High_parent_heterosis[-c(1:907),]


# calculation of heterosis using high parent values i.e. high parent hterosis

NIR_HP_heterosis_19CS = scale(savitzkyGolay(HP_heterosis_19CS, m=1, p=1, w=11))
NIR_HP_heterosis_19TA = scale(savitzkyGolay(HP_heterosis_19TA, m=1, p=1, w=11))
NIR_HP_heterosis_20CS = scale(savitzkyGolay(HP_heterosis_20CS, m=1, p=1, w=11))
NIR_HP_heterosis_20MA = scale(savitzkyGolay(HP_heterosis_20MA, m=1, p=1, w=11))

NIR_HP_heterosis.d1 = rbind(NIR_HP_heterosis_19CS, 
                            NIR_HP_heterosis_19TA, 
                            NIR_HP_heterosis_20CS, 
                            NIR_HP_heterosis_20MA) 

NIR_HP_heterosis.ZN1 = tcrossprod(as.matrix(NIR_HP_heterosis.d1)/ncol(as.matrix(NIR_HP_heterosis.d1))) #phenomic relationship matrices
dim(NIR_HP_heterosis.ZN1) 

# Set the predictors/priors to run the model
# mid-parent

# calculating NIR x E interaction for different models using Hadamard product
mid_parent.ZNZE1           = NIR_mid_parent.ZN1 * ZEZEt
mid_parent_heterosis_ZNZE1 = NIR_MP_heterosis.ZN1 * ZEZEt
high_parent.ZNZE1          = NIR_HP_heterosis.ZN1 * ZEZEt



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
           list(K= NIR_mid_parent.ZN1, model = "RKHS"),    #NIR1
           list(K = mid_parent.ZNZE1, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent heterosis NIR

Eta3<-list(list(X = ZE,model="BRR"),      #Env
           list(K = NIR_MP_heterosis.ZN1, model="RKHS"),    #NIR1
           list(K = mid_parent_heterosis_ZNZE1, model="RKHS"))  #NIR1 x Env

##### first derivative of high_parent heterosis NIR

Eta4<-list(list(X = ZE,model = "BRR"),      #Env
           list(K = NIR_HP_heterosis.ZN1, model = "RKHS"),    #NIR1
           list(K = high_parent.ZNZE1, model = "RKHS")) #NIR1 x Env


##### first derivative of mid_parent + high-parent heterosis NIR

Eta5<-list(list(X = ZE,model="BRR"),      #Env
           list(K = NIR_MP_heterosis.ZN1, model="RKHS"),    #NIR1
           list(K = mid_parent_heterosis_ZNZE1, model="RKHS"),
           list(K = NIR_HP_heterosis.ZN1, model="RKHS"),    #NIR1
           list(K = high_parent.ZNZE1, model="RKHS")) #NIR1 x Env


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

###Missing for multi-trait models using DA and PH


### 6 different models 
# MODEL 1: Genomic model
# MODEL 2: Mid-parent model
# MODEL 3: mid_parent heterosis NIR
# MODEL 4: high_parent heterosis NIR
# MODEL 5: mid_parent + high-parent heterosis NIR
# MODEL 6: first derivative of mid_parent NIR + genomic


Models <- list(Eta1, Eta2, Eta3, Eta4, Eta5, Eta6)
traitnames <- c("yield", "da", "ph")
pheno_combined[1:10,1:10]

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

# tr=1
set.seed(1234)
for (tr in 1:length(traitnames)) {
  
  if (tr == 1) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==2) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  }
  
  # cross-validation #
  hybrid = as.character(unique(pheno$pedigree))
  Phenotype_data1 = pheno
  
  cycles = 10
  CV2 = list()
  CV1 = list()
  
  #MODEL =6; rep_num=1
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
      }
    #rep_num =1
    if (rep_num == 5) {
      CV1out <- plyr::ldply(CV1, data.frame)
      write.csv(CV1out, file = paste("ACC_", traitnames[tr],"_CV1_", MODEL, ".csv", sep=""), row.names = FALSE)
      }
  }
}

####### plotting the results
# quick view of result

setwd("output/")
###lets plot
library(plyr)
library(readr)
library(dplyr)
library(stringr)

list_csv_files <- list.files("../output/")
df2 <- readr::read_csv(list_csv_files, id = "file_name") %>% as.data.frame()
df2
write.csv(df2, "Pred.ability.70:30_5rep.csv")
df2 = read.csv("Pred.ability.70:30_5rep.csv")

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
      TRUE~trait)
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


#############DPAC
#K means clustering; that needs to yield lower BIC;
# runa fter transforming the data using PCA


library(adegenet)
data("dapcIllus")
class(dapcIllus)
names(dapcIllus)
x <- dapcIllus$a
grp <- find.clusters(x, max.n.clust = 40)

dapc1 = dapc(pheno_combined[,-c(1:7)], pheno_combined$env)
scatter(dapc1)


inbreds_combined = rbind(parents_20CS, parents_20MA)
dapc2 = dapc(inbreds_combined[,-c(1:5)], inbreds_combined$env)
scatter(dapc2)


scatter(dapc1, 
       scree.da=FALSE, 
       bg="white", 
       pch=20, 
       cell=0, 
       cstar=0, 
       #col=myCol, 
       solid=.4, 
       cex=3,
       clab=0, 
       leg=TRUE, 
       txt.leg=paste("Cluster",1:6))



jpeg("pred_acc.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()

#############

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


