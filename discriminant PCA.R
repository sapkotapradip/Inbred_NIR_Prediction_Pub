### From youtube session
library(readxl)
library(MASS)
library(caret)
library(ggplot2)
library(dplyr)
library(pROC)

#### reading pheno data for inbreds
parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
parents_20MA = read.csv("blues_20MA_inbreds_spectra.csv", check.names = FALSE)

# 20CS
parents_20CS[1:10,1:10]
dim(parents_20CS)

group = c(rep("Seed Parents (A)",45), rep("Pollinator Parents (R)",44))
sum(is.na(parents_20CS)) # need to omit or impute missing values

lda.model = lda(group ~., data = parents_20CS[,-c(1:5)])

# Predict LDA scores
lda_scores <- predict(lda.model)$x
lda_data <- data.frame(lda_scores, Group = group)

# Print first few rows to check
head(lda_data)

print(dim(lda_scores))  # Check the dimensions of LDA output
colnames(lda_scores)    # Print available column names

lda.model$svd


fig = ggplot(lda_data, aes(x = LD1, y = Group, color = Group)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  theme_minimal() +
  # labs(title = "",
  #      x = x ,
  #      y = y,
  #      color = "Group") +
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=10, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=17, face="bold", colour = "black"),    
    axis.title.y = element_text(size=17, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 13, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 13, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    #axis.text.x.bottom = element_blank()
  )

fig <- fig + guides(color = guide_legend(title = "Heterotic Groups"))


# Combining two heterotic groups for two environments 
parents_combined = rbind(parents_20CS,parents_20MA)
group = c(c(rep("Seed Parents - CS",45), rep("Pollinator Parents - CS",44)),c(rep("Seed Parents - LY",45), rep("Pollinator Parents - LY",44)))

lda.model = lda(group ~., data = parents_combined[,-c(1:5)])

# Predict LDA scores
lda_scores <- predict(lda.model)$x
lda_data <- data.frame(lda_scores, Group = group)

# Print first few rows to check
head(lda_data)

print(dim(lda_scores))  # Check the dimensions of LDA output
colnames(lda_scores)    # Print available column names

lda.model$svd
a = lda.model$svd^2 / sum(lda.model$svd^2)

x = paste0(round(a[1]*100, digits = 1), " %")
y = paste0(round(a[2]*100, digits = 1), " %")


fig1 = ggplot(lda_data, aes(x = LD1, y = LD2, color = Group)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  theme_minimal() +
  labs(title = "A",
       x = x ,
       y = y,
       color = "Group") +
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=10, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=17, face="bold", colour = "black"),    
    axis.title.y = element_text(size=17, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 13, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 13, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    #axis.text.x.bottom = element_blank()
  )

fig1 <- fig1 + guides(color = guide_legend(title = "Heterotic Groups"))

#### combined environments for hybrids grain samples
#### reading pheno data for hybrids
pheno_19CS = read.csv("blues_19CS_hybrids_spectra.csv", check.names = FALSE)
pheno_19TA = read.csv("blues_19TA_hybrids_spectra.csv", check.names = FALSE)
pheno_20CS = read.csv("blues_20CS_hybrids_spectra.csv", check.names = FALSE)
pheno_20MA = read.csv("blues_20MA_hybrids_spectra.csv", check.names = FALSE)

pheno_combined = rbind(pheno_19CS,pheno_19TA, pheno_20CS, pheno_20MA)

group = pheno_combined$env
lda.model = lda(group ~., data = pheno_combined[,-c(1:7)])

# Predict LDA scores
lda_scores <- predict(lda.model)$x
lda_data <- data.frame(lda_scores, Group = group)
lda.model$svd
a = lda.model$svd^2 / sum(lda.model$svd^2)

x = paste0(round(a[1]*100, digits = 1), " %")
y = paste0(round(a[2]*100, digits = 1), " %")


fig2 = ggplot(lda_data, aes(x = LD1, y = LD2, color = Group)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  theme_minimal() +
  labs(title = "B",
       x = x ,
       y = y,
       color = "Group") +
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=10, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=17, face="bold", colour = "black"),    
    axis.title.y = element_text(size=17, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 13, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 13, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    #axis.text.x.bottom = element_blank()
  )

fig2 <- fig2 + guides(color = guide_legend(title = "Environments"))

library(ggpubr)
merge1 = ggarrange(fig1,fig2,
                   ncol = 1, nrow = 2,
                   common.legend = FALSE,
                   legend = c("bottom"))



jpeg("pca_heterotic_env_NIR.jpeg",width = 8,height =10,units = "in", res=600)
merge1
dev.off()


#This goes into the publication
#################################################################################
#################################################################################







############### Rough ###########################

#############DPAC
library(adegenet)

#### reading pheno data for hybrids
pheno_19CS = read.csv("blues_19CS_hybrids_spectra.csv", check.names = FALSE)
pheno_19TA = read.csv("blues_19TA_hybrids_spectra.csv", check.names = FALSE)
pheno_20CS = read.csv("blues_20CS_hybrids_spectra.csv", check.names = FALSE)
pheno_20MA = read.csv("blues_20MA_hybrids_spectra.csv", check.names = FALSE)

#### reading pheno data for inbreds
parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
parents_20MA = read.csv("blues_20MA_inbreds_spectra.csv", check.names = FALSE)


pheno_combined = rbind(pheno_19CS,pheno_19TA, pheno_20CS, pheno_20MA)
parents = rbind(parents_20CS,parents_20MA)

dapc1 = dapc(pheno_combined[,-c(1:7)], pheno_combined$env)
p = scatter(dapc1)

myCol <-c("darkblue","purple","green","orange","red","blue")

scatter(dapc1,
        scree.da = FALSE,
        bg = "white",
        pch = 20,
        cell = 0,
        cstar=0,
        col=myCol,
        solid = .4,
        cex=3,
        clab =0,
        leg = TRUE,
        txt.leg=paste("Cluster",1:3))

myPal <-colorRampPalette(c("blue","gold","red"))

scatter(dapc1,
        col = transp(myPal(6)),
        scree.da = FALSE,
        cell = 1.5,
        cex=2,
        bg = "white",
        cstar =0)


parents_20CS$group = c(rep("A",45), rep("R",44))

dapc1 = dapc(parents_20CS[,-c(1:5,4206)], parents_20CS$group)

parents_20MA$group = c(rep("A1",45), rep("R1",44))

x = rbind(parents_20CS,parents_20MA)

dapc1 = dapc(x[,-c(1:5,4206)], x$group)

p = scatter(dapc1)

scale(savitzkyGolay(x[,-c(1:5,4206)], m =1, p =3, w=11))
# MASS # LDA

scatter(dapc1,
        col = transp(myPal(6)),
        scree.da = FALSE,
        cell = 1.5,
        cex=2,
        bg = "white",
        cstar =0)


library(MASS)
p = parents_20CS[,-(1:5)]

lda = lda(group ~ ., data = p) 
g= p$group
# Predict LDA scores
lda_scores <- predict(lda)$x
lda_data <- data.frame(lda_scores, Group = g)

# Print first few rows to check
head(lda_data)

print(dim(lda_scores))  # Check the dimensions of LDA output
colnames(lda_scores)    # Print available column names

lda_result$svd
lda_result$svd^2 / sum(lda_result$svd^2)


A = ggplot(lda_data, aes(x = LD1, y = LD2, color = Group)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  theme_minimal() +
  labs(title = "Discriminant PCA (DPCA) - Single LD",
       x = "Linear Discriminant 1 (LD1)",
       y = "Linear Discriminant 2 (LD2)") +
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=10, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=17, face="bold", colour = "black"),    
    axis.title.y = element_text(size=17, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 13, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 13, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    #axis.text.x.bottom = element_blank()
  )

A = A + guides(fill = guide_legend(title = "Environments"))






jpeg("dapc1.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()

dapc1 = dapc(parents[,-c(1:5)], parents$env)
p = scatter(dapc1)

dapc1$posterior

names(dapc1)

grp = find.clusters(dapc1, max.n.clust=35)



library(adegenet) 
data(dapcIllus) 
class(dapcIllus)



inbreds_combined = rbind(parents_20CS, parents_20MA)
dapc2 = dapc(inbreds_combined[,-c(1:5)], inbreds_combined$env)
scatter(dapc2)

inbreds_20CS = parents_20CS
parents_20CS$group = c(rep("A",45), rep("R",44))

dapc2 = dapc(parents_20CS[,-c(1:5,4206)], parents_20CS$group)
scatter(dapc2)

inbreds_20MA = parents_20MA
parents_20MA$group = c(rep("A",45), rep("R",44))

dapc2 = dapc(parents_20MA[,-c(1:5,4206)], parents_20MA$group)
scatter(dapc2)



jpeg("pred_acc.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()