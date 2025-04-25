#################### Phenotypic data analysis ####################################
#########
#########Distribution of BLUES ###################################################

## plotting phenotyping data 

library(ggplot2)
library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)
library(agricolae)

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

names(pheno_combined) = tolower(names(pheno_combined))
read.csv("blues_20CS_pheno_parents.csv")
pheno_combined$env = factor(pheno_combined$env, 
                            levels = c("19CS", "19TA", 
                                       "20CS", "20LY"))

########## this is just for phenotypes of hybrids #############

a <- ggplot(pheno_combined, aes(env, gy, fill = env)) +
  geom_boxplot(colour = "blue")+
  labs(y = "Grain Yield (Kg/ha)", x= "")+
  #scale_y_continuous(limits = c(0.2,1))+
  #ggtitle("b) Phenotypic distribution for kernel physical factors")+
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=10, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=10, face="bold", colour = "black"),    
    axis.title.y = element_text(size=10, face="bold", colour = "black"),    
    axis.text.x = element_text(size=10, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=10, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 5, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 5, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.4),
    axis.line.y = element_line(color="black", size = 0.4),
    panel.border = element_rect(colour = "black", fill=NA, size=0.4),
    axis.text.x.bottom = element_blank(),
    legend.text = element_text(size = rel(1))
  )


b <- ggplot(pheno_combined, aes(env, dy, fill = env)) +
  geom_boxplot(colour = "purple")+
  labs(y = "Days to Anthesis", x= "")+
  #scale_y_continuous(limits = c(0.2,1))+
  #ggtitle("b) Phenotypic distribution for kernel physical factors")+
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=10, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=10, face="bold", colour = "black"),    
    axis.title.y = element_text(size=10, face="bold", colour = "black"),    
    axis.text.x = element_text(size=10, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=10, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 5, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 5, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.4),
    axis.line.y = element_line(color="black", size = 0.4),
    panel.border = element_rect(colour = "black", fill=NA, size=0.4),
    axis.text.x.bottom = element_blank(),
    legend.text = element_text(size = rel(1))
  )


c <- ggplot(pheno_combined, aes(env, ph, fill = env)) +
  geom_boxplot(colour = "brown")+
  labs(y = "Days to Anthesis", x= "Environment")+
  #scale_y_continuous(limits = c(0.2,1))+
  #ggtitle("b) Phenotypic distribution for kernel physical factors")+
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=10, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=10, face="bold", colour = "black"),    
    axis.title.y = element_text(size=10, face="bold", colour = "black"),    
    axis.text.x = element_text(size=10, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=10, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 5, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 5, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.4),
    axis.line.y = element_line(color="black", size = 0.4),
    panel.border = element_rect(colour = "black", fill=NA, size=0.4),
    axis.text.x.bottom = element_blank(),
    legend.text = element_text(size = rel(1))
  )

merge = ggarrange(a,b,c,
                  ncol = 1, nrow = 3,
                  common.legend = TRUE,
                  legend = c("bottom"))


merge

jpeg("Blues_hybrids.jpeg",
     width = 11,height = 7,
     units = "in", res=600)

merge
dev.off()


################# Variance component estimations for hybrids ##############
hybrids = unique(pheno_combined$pedigree)
pheno = read.csv("raw_phenotypic_data.csv")
unique(pheno$Pedigree)
names(pheno) = tolower(names(pheno))

pheno = pheno %>% filter(
  pedigree %in% hybrids)

########### Variance components estimations for grain yield days to anthesis and plant height##
library(lme4)
library(lmerTest)

pheno = read.csv("raw_phenotypic_data.csv")
head(pheno)
names(pheno) = tolower(names(pheno))


## checking for NAs
colSums(is.na(pheno))

#### converting variables in factors

pheno$pedigree <- as.factor(pheno$pedigree)
pheno$env <- as.factor(pheno$env)
pheno$rep <- as.factor(pheno$rep)
pheno$female <- as.factor(pheno$female)
pheno$male <- as.factor(pheno$male)
pheno$ra <- as.factor(pheno$ra)
pheno$ro <- as.factor(pheno$ro)


##### grain yield
model <- lmer(gy ~ (1|female) + (1|male) + (1|pedigree) + (1|rep/env) +
                (1|env) + (1|ra/env) + (1|ro/env)+ (1|female:env) + 
                (1|male:env) + (1|female:male:env), 
              pheno)

# 2 reps in 19CS
# 2 reps in 19TA
# 2 reps in 20CS
# 2 reps in 20LY

ranova(model)
summary(model)
GCAmales <- ranef(model)$male
GCAfemales <- ranef(model)$female
SCA <- ranef(model)$"pedigree"
GXE = ranef(model)$"female:male:env"
GCAmales
GCAfemales
SCA

### calculation of variance heritability and CV component
varcomp <- print(VarCorr(model), com = c("Variance", "Std.Dev."))
variances = as.data.frame(varcomp)[,c(1,4,5)]
hybridvar = variances[2,2] + variances[8,2] + variances[9,2]
hybridenvvar = variances[1,2] + variances[3,2] +variances[4,2]
errorvar = variances[14,2]

addvar = variances[8,2] + variances[9,2]

cve = sqrt(errorvar)*100/mean(na.omit(pheno$gy))

repeatability = hybridvar/(hybridvar+(hybridenvvar/4)+errorvar/8) #BSH
NSH = addvar/(hybridvar+(hybridenvvar/4)+errorvar/8) #NSH
NSHf = variances[9,2]/(hybridvar+(hybridenvvar/4)+errorvar/8)
NSHm = variances[8,2]/(hybridvar+(hybridenvvar/4)+errorvar/8)


##### plant height
model <- lmer(ph ~ (1|female) + 
                (1|male) + 
                (1|pedigree) + 
                (1|rep/env) +
                (1|env) + 
                (1|ra/env) + 
                (1|ro/env)+ 
                (1|female:env) + 
                (1|male:env) + 
                (1|female:male:env), 
              pheno)

# 2 reps in 19CS
# 2 reps in 19TA
# 2 reps in 20CS
# 2 reps in 20LY

ranova(model)
summary(model)
GCAmales <- ranef(model)$male
GCAfemales <- ranef(model)$female
SCA <- ranef(model)$"pedigree"
GXE = ranef(model)$"female:male:env"
GCAmales
GCAfemales
SCA

### calculation of variance heritability and CV component
varcomp <- print(VarCorr(model), com = c("Variance", "Std.Dev."))
variances = as.data.frame(varcomp)[,c(1,4,5)]
hybridvar = variances[2,2] + variances[8,2] + variances[9,2]
hybridenvvar = variances[1,2] + variances[3,2] +variances[4,2]
errorvar = variances[14,2]

addvar = variances[8,2] + variances[9,2]

cve = sqrt(errorvar)*100/mean(na.omit(pheno$ph))

repeatability = hybridvar/(hybridvar+(hybridenvvar/4)+errorvar/8) #BSH
NSH = addvar/(hybridvar+(hybridenvvar/4)+errorvar/8) #NSH
NSHf = variances[9,2]/(hybridvar+(hybridenvvar/4)+errorvar/8)
NSHm = variances[8,2]/(hybridvar+(hybridenvvar/4)+errorvar/8)

print(cve)
print(repeatability)
print(NSH)
print(NSHf)
print(NSHm)




##### days to anthesis
model <- lmer(dy ~ (1|female) + 
                (1|male) + 
                (1|pedigree) + 
                (1|rep/env) +
                (1|env) + 
                (1|ra/env) + 
                (1|ro/env)+ 
                (1|female:env) + 
                (1|male:env) + 
                (1|female:male:env), 
              pheno)

# 2 reps in 19CS
# 2 reps in 19TA
# 2 reps in 20CS
# 2 reps in 20LY

ranova(model)
summary(model)
GCAmales <- ranef(model)$male
GCAfemales <- ranef(model)$female
SCA <- ranef(model)$"pedigree"
GXE = ranef(model)$"female:male:env"
GCAmales
GCAfemales
SCA

### calculation of variance heritability and CV component
varcomp <- print(VarCorr(model), com = c("Variance", "Std.Dev."))
variances = as.data.frame(varcomp)[,c(1,4,5)]
hybridvar = variances[2,2] + variances[8,2] + variances[9,2]
hybridenvvar = variances[1,2] + variances[3,2] +variances[4,2]
errorvar = variances[14,2]

addvar = variances[8,2] + variances[9,2]

cve = sqrt(errorvar)*100/mean(na.omit(pheno$dy))

repeatability = hybridvar/(hybridvar+(hybridenvvar/4)+errorvar/8) #BSH
NSH = addvar/(hybridvar+(hybridenvvar/4)+errorvar/8) #NSH
NSHf = variances[9,2]/(hybridvar+(hybridenvvar/4)+errorvar/8)
NSHm = variances[8,2]/(hybridvar+(hybridenvvar/4)+errorvar/8)

print(cve)
print(repeatability)
print(NSH)
print(NSHf)
print(NSHm)


