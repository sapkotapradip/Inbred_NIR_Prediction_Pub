# Load necessary package
# install.packages("data.table")  # If not installed
library(data.table)

# Read the CSV file
data <- fread("blues_20CS_inbreds_spectra.csv")

# Print the dimensions (rows x columns)
print(dim(data))

# Extract the outcome (first column)
outcome <- data[[1]]

# Remove the first column and next 4 columns (i.e., columns 2 to 5)
data_filtered <- data[, -c(1:5), with = FALSE]

# Print new dimensions
print(dim(data_filtered))

# Convert outcome to two categories based on the first character
outcome <- ifelse(substr(outcome, 1, 1) == "A", "A", "R")

# Convert to categorical (factor) type
outcome <- as.factor(outcome)

# Check the summary to verify
print(summary(outcome))

# 
# # Load necessary packages
# install.packages("Boruta")  # If not installed
# install.packages("randomForest")  # Boruta requires randomForest

library(Boruta)
library(randomForest)

# Run Boruta feature selection
set.seed(123)  # For reproducibility
# boruta_result <- Boruta(outcome ~ ., data = data_filtered, doTrace = 2)

boruta_result <- Boruta(outcome ~ ., data = data_filtered, doTrace = 2, pValue = 0.5)



# Print summary of results
print(boruta_result)

# Get selected important features
selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)

# Subset dataset with selected features
data_selected <- data_filtered[, selected_features, with = FALSE]

# Print new dimensions
print(dim(data_selected))












# install.packages("signal")  # If not installed
library(signal)  # For Savitzky-Golay filter





# Function to apply Savitzky-Golay filter to each column
smooth_data <- function(df, window = 11, polyorder = 3) {
  apply(df, 2, function(x) sgolayfilt(x, p = polyorder, n = window))
}

# Apply smoothing to feature columns
data_smoothed <- as.data.frame(smooth_data(as.matrix(data_filtered)))

# Print new dimensions after smoothing
print(dim(data_smoothed))


library(Boruta)
library(randomForest)

set.seed(123)  # Ensure reproducibility

# Run Boruta with smoothed data
boruta_result <- Boruta(outcome ~ ., data = data_smoothed, doTrace = 2, 
                        maxRuns = 500, pValue = 0.05, mcAdj = FALSE)

# Handle tentative features
boruta_result <- TentativeRoughFix(boruta_result)

# Get selected important features
selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)

# Subset dataset with selected features
data_selected <- data_smoothed[, selected_features, drop = FALSE]

# Print new dimensions
print(dim(data_selected))

library(MASS)
library(ggplot2)

# Ensure outcome is categorical
outcome <- as.factor(outcome)

# Apply PCA on selected features
pca_result <- prcomp(data_selected, center = TRUE, scale. = TRUE)

# sdev of PCS
sdev = pca_result$sdev
variance_explained = sdev^2
prop_variance_explained <- variance_explained / sum(variance_explained)
percentage_PC1_PC2 <- sum(prop_variance_explained[1:2]) * 100

PC1 = prop_variance_explained[1]
PC2 = prop_variance_explained[2]

# Extract principal components
pca_scores <- data.frame(pca_result$x[, 1:2])  # First 2 PCs


## added by pradip
##
pca_scores$group = c(rep("A", 45), rep("R",44))

ggplot(pca_scores, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  labs(x = PC1, y = PC2) +
  
  
  
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

#####

# Apply Linear Discriminant Analysis (LDA)
lda_result <- lda(outcome ~ ., data = pca_scores)

# Predict LDA scores
lda_scores <- predict(lda_result)$x
lda_data <- data.frame(lda_scores, Group = outcome)

# Print first few rows to check
head(lda_data)





print(dim(lda_scores))  # Check the dimensions of LDA output
colnames(lda_scores)    # Print available column names


ggplot(lda_data, aes(x = LD1, y = Group, color = Group)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  theme_minimal() +
  labs(title = "Discriminant PCA (DPCA) - Single LD",
       x = "Linear Discriminant 1 (LD1)",
       y = "Group") +
  theme(legend.title = element_blank())



# Perform PCA before LDA
pca_result <- prcomp(data_selected, center = TRUE, scale. = TRUE)

# Take top principal components (e.g., first 5)
pca_data <- as.data.frame(pca_result$x[, 1:5])  # Use first 5 PCs

# Add the outcome variable
pca_data$Group <- outcome

# Run LDA on PCA-transformed data
lda_result <- lda(Group ~ ., data = pca_data)

# Extract LDA scores
lda_scores <- as.data.frame(predict(lda_result)$x)

# Add group labels
lda_scores$Group <- outcome

# Check dimensions
print(dim(lda_scores))  # Should have LD1 and potentially LD2






ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA Scatter Plot",
       x = "Principal Component 1 (PC1)",
       y = "Principal Component 2 (PC2)") +
  theme(legend.title = element_blank())








# install.packages("e1071")  # If not installed
library(e1071)

# Train SVM classifier
svm_model <- svm(Group ~ PC1 + PC2, data = pca_data, kernel = "linear")

# Create grid for decision boundary
grid <- expand.grid(PC1 = seq(min(pca_data$PC1), max(pca_data$PC1), length.out = 100),
                    PC2 = seq(min(pca_data$PC2), max(pca_data$PC2), length.out = 100))

# Predict on grid
grid$Group <- predict(svm_model, newdata = grid)

# Plot with decision boundary
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_contour(data = grid, aes(x = PC1, y = PC2, z = as.numeric(Group)), bins = 1, color = "black") +
  theme_minimal() +
  labs(title = "PCA Scatter Plot with Decision Boundary",
       x = "Principal Component 1 (PC1)",
       y = "Principal Component 2 (PC2)") +
  theme(legend.title = element_blank())



# eigen vectors
parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
parents_20MA = read.csv("blues_20MA_inbreds_spectra.csv", check.names = FALSE)

rownames(parents_20CS) = parents_20CS$pedigree
X = scale(parents_20CS[,-c(1:6)], center = TRUE)
G = tcrossprod(X)/ncol(X)

out_eigen = eigen(G)
round(out_eigen$values/sum(out_eigen$values)*100,2)

plot(out_eigen$vectors[,1:2],
     xlab="Eigen vector 1 (X%)", ylab="Eigen vector 2 (Y%)",pch=19,
     xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))

text(0.90*out_eigen$vectors[,1],1.2*out_eigen$vectors[,2],labels=rownames(G),cex=0.5)


# Assuming out_eigen$vectors is a data frame or matrix
df <- as.data.frame(out_eigen$vectors[, 1:2])
df$label <- rownames(G)

names = as.data.frame(cbind(df$label,c(rep("A", 45), rep("R",44))))

df$Categories <- names$V2

pc1_var <- round(out_eigen$values[1] / sum(out_eigen$values) * 100, 2)
pc2_var <- round(out_eigen$values[2] / sum(out_eigen$values) * 100, 2)

A = ggplot(df, aes(x = V1, y = V2, colour = Categories)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  labs(x = pc1_var, y = pc2_var) +
  
  
  
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

################################# 20MA ##########
# eigen vectors
parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
parents_20MA = read.csv("blues_20MA_inbreds_spectra.csv", check.names = FALSE)

rownames(parents_20MA) = parents_20CS$pedigree
X = scale(parents_20MA[,-c(1:6)], center = TRUE)
G = tcrossprod(X)/ncol(X)

out_eigen = eigen(G)
round(out_eigen$values/sum(out_eigen$values)*100,2)

plot(out_eigen$vectors[,1:2],
     xlab="Eigen vector 1 (X%)", ylab="Eigen vector 2 (Y%)",pch=19,
     xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))

text(0.90*out_eigen$vectors[,1],1.2*out_eigen$vectors[,2],labels=rownames(G),cex=0.5)


# Assuming out_eigen$vectors is a data frame or matrix
df <- as.data.frame(out_eigen$vectors[, 1:2])
df$label <- rownames(G)

names = as.data.frame(cbind(df$label,c(rep("A", 45), rep("R",44))))

df$Categories <- names$V2

pc1_var <- round(out_eigen$values[1] / sum(out_eigen$values) * 100, 2)
pc2_var <- round(out_eigen$values[2] / sum(out_eigen$values) * 100, 2)

A = ggplot(df, aes(x = V1, y = V2, colour = Categories)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  labs(x = pc1_var, y = pc2_var) +
  
  
  
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
A

### raw and scaled didnot work but let do first derivative


# first derivative

# eigen vectors
parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)


rownames(parents_20CS) = parents_20CS$pedigree

X = scale(savitzkyGolay(parents_20CS[,-c(1:6)], m =1, p =1, w=11))
G = tcrossprod(X)/ncol(X)

out_eigen = eigen(G)
round(out_eigen$values/sum(out_eigen$values)*100,2)

plot(out_eigen$vectors[,1:2],
     xlab="Eigen vector 1 (X%)", ylab="Eigen vector 2 (Y%)",pch=19,
     xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))

text(0.90*out_eigen$vectors[,1],1.2*out_eigen$vectors[,2],labels=rownames(G),cex=0.5)


# Assuming out_eigen$vectors is a data frame or matrix
df <- as.data.frame(out_eigen$vectors[, 1:2])
df$label <- rownames(G)

names = as.data.frame(cbind(df$label,c(rep("A", 45), rep("R",44))))

df$Categories <- names$V2

pc1_var <- round(out_eigen$values[1] / sum(out_eigen$values) * 100, 2)
pc2_var <- round(out_eigen$values[2] / sum(out_eigen$values) * 100, 2)

A = ggplot(df, aes(x = V1, y = V2, colour = Categories)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  labs(x = pc1_var, y = pc2_var) +
  
  
  
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

################################# 20MA ##########
# eigen vectors
parents_20MA = read.csv("blues_20MA_inbreds_spectra.csv", check.names = FALSE)
rownames(parents_20MA) = parents_20MA$pedigree

X = scale(savitzkyGolay(parents_20MA[,-c(1:6)], m =1, p =1, w=11))
G = tcrossprod(X)/ncol(X)

out_eigen = eigen(G)
round(out_eigen$values/sum(out_eigen$values)*100,2)

plot(out_eigen$vectors[,1:2],
     xlab="Eigen vector 1 (X%)", ylab="Eigen vector 2 (Y%)",pch=19,
     xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))

text(0.90*out_eigen$vectors[,1],1.2*out_eigen$vectors[,2],labels=rownames(G),cex=0.5)


# Assuming out_eigen$vectors is a data frame or matrix
df <- as.data.frame(out_eigen$vectors[, 1:2])
df$label <- rownames(G)

names = as.data.frame(cbind(df$label,c(rep("A", 45), rep("R",44))))

df$Categories <- names$V2

pc1_var <- round(out_eigen$values[1] / sum(out_eigen$values) * 100, 2)
pc2_var <- round(out_eigen$values[2] / sum(out_eigen$values) * 100, 2)

A = ggplot(df, aes(x = V1, y = V2, colour = Categories)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  labs(x = pc1_var, y = pc2_var) +
  
  
  
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
A

# Second derivative; p =2 

parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)


rownames(parents_20CS) = parents_20CS$pedigree

X = scale(savitzkyGolay(parents_20CS[,-c(1:6)], m =1, p =2, w=11))
G = tcrossprod(X)/ncol(X)

out_eigen = eigen(G)
round(out_eigen$values/sum(out_eigen$values)*100,2)

plot(out_eigen$vectors[,1:2],
     xlab="Eigen vector 1 (X%)", ylab="Eigen vector 2 (Y%)",pch=19,
     xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))

text(0.90*out_eigen$vectors[,1],1.2*out_eigen$vectors[,2],labels=rownames(G),cex=0.5)


# Assuming out_eigen$vectors is a data frame or matrix
df <- as.data.frame(out_eigen$vectors[, 1:2])
df$label <- rownames(G)

names = as.data.frame(cbind(df$label,c(rep("A", 45), rep("R",44))))

df$Categories <- names$V2

pc1_var <- round(out_eigen$values[1] / sum(out_eigen$values) * 100, 2)
pc2_var <- round(out_eigen$values[2] / sum(out_eigen$values) * 100, 2)

A = ggplot(df, aes(x = V1, y = V2, colour = Categories)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  labs(x = pc1_var, y = pc2_var) +
  
  
  
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

################################# 20MA ##########
# eigen vectors
parents_20MA = read.csv("blues_20MA_inbreds_spectra.csv", check.names = FALSE)
rownames(parents_20MA) = parents_20MA$pedigree

X = scale(savitzkyGolay(parents_20MA[,-c(1:6)], m =1, p =2, w=11))
G = tcrossprod(X)/ncol(X)

out_eigen = eigen(G)
round(out_eigen$values/sum(out_eigen$values)*100,2)

plot(out_eigen$vectors[,1:2],
     xlab="Eigen vector 1 (X%)", ylab="Eigen vector 2 (Y%)",pch=19,
     xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))

text(0.90*out_eigen$vectors[,1],1.2*out_eigen$vectors[,2],labels=rownames(G),cex=0.5)


# Assuming out_eigen$vectors is a data frame or matrix
df <- as.data.frame(out_eigen$vectors[, 1:2])
df$label <- rownames(G)

names = as.data.frame(cbind(df$label,c(rep("A", 45), rep("R",44))))

df$Categories <- names$V2

pc1_var <- round(out_eigen$values[1] / sum(out_eigen$values) * 100, 2)
pc2_var <- round(out_eigen$values[2] / sum(out_eigen$values) * 100, 2)

A = ggplot(df, aes(x = V1, y = V2, colour = Categories)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  labs(x = pc1_var, y = pc2_var) +
  
  
  
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
A

# Third derivative; p = 3 

parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)


rownames(parents_20CS) = parents_20CS$pedigree

X = scale(savitzkyGolay(parents_20CS[,-c(1:6)], m =1, p =3, w=11))
G = tcrossprod(X)/ncol(X)

out_eigen = eigen(G)
round(out_eigen$values/sum(out_eigen$values)*100,2)

plot(out_eigen$vectors[,1:2],
     xlab="Eigen vector 1 (X%)", ylab="Eigen vector 2 (Y%)",pch=19,
     xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))

text(0.90*out_eigen$vectors[,1],1.2*out_eigen$vectors[,2],labels=rownames(G),cex=0.5)


# Assuming out_eigen$vectors is a data frame or matrix
df <- as.data.frame(out_eigen$vectors[, 1:2])
df$label <- rownames(G)

names = as.data.frame(cbind(df$label,c(rep("A", 45), rep("R",44))))

df$Categories <- names$V2

pc1_var <- round(out_eigen$values[1] / sum(out_eigen$values) * 100, 2)
pc2_var <- round(out_eigen$values[2] / sum(out_eigen$values) * 100, 2)

A = ggplot(df, aes(x = V1, y = V2, colour = Categories)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  labs(x = pc1_var, y = pc2_var) +
  
  
  
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

A
################################# 20MA ##########
# eigen vectors
parents_20MA = read.csv("blues_20MA_inbreds_spectra.csv", check.names = FALSE)
rownames(parents_20MA) = parents_20MA$pedigree

X = scale(savitzkyGolay(parents_20MA[,-c(1:6)], m =1, p =3, w=11))
G = tcrossprod(X)/ncol(X)

out_eigen = eigen(G)
round(out_eigen$values/sum(out_eigen$values)*100,2)

plot(out_eigen$vectors[,1:2],
     xlab="Eigen vector 1 (X%)", ylab="Eigen vector 2 (Y%)",pch=19,
     xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))

text(0.90*out_eigen$vectors[,1],1.2*out_eigen$vectors[,2],labels=rownames(G),cex=0.5)


# Assuming out_eigen$vectors is a data frame or matrix
df <- as.data.frame(out_eigen$vectors[, 1:2])
df$label <- rownames(G)

names = as.data.frame(cbind(df$label,c(rep("A", 45), rep("R",44))))

df$Categories <- names$V2

pc1_var <- round(out_eigen$values[1] / sum(out_eigen$values) * 100, 2)
pc2_var <- round(out_eigen$values[2] / sum(out_eigen$values) * 100, 2)

A = ggplot(df, aes(x = V1, y = V2, colour = Categories)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  labs(x = pc1_var, y = pc2_var) +
  
  
  
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
A




##################### scaling, first derivative, second derivative, nothing worked #####################
# Lets perform boruta algorithm or variable selction
library(Boruta)
library(randomForest)

data <- read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
data_filtered <- data[, -c(1:5)]
data_filtered$outcome = c(rep("A", 45), rep("R",44))
data_filtered$outcome = as.factor(data_filtered$outcome)
boruta_result <- Boruta(outcome ~ ., data = data_filtered, doTrace = 2, pValue = 0.5)
print(boruta_result) # no attributes were deemed to be important

# scaled 
data <- read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
data_filtered <- data[, -c(1:5)]
data_filtered1 = scale(data_filtered,center = TRUE)
data_filtered$outcome = c(rep("A", 45), rep("R",44))
data_filtered$outcome = as.factor(data_filtered$outcome)
boruta_result <- Boruta(outcome ~ ., data = data_filtered, doTrace = 2, pValue = 0.5)
print(boruta_result) # no attributes were deemed to be important

# Get selected important features
selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)

# Subset dataset with selected features
data_selected <- data_filtered[, selected_features, drop = FALSE]

# Print new dimensions
print(dim(data_selected))

library(MASS)
library(ggplot2)

# Ensure outcome is categorical
outcome <- as.factor(outcome)

# Apply PCA on selected features
pca_result <- prcomp(data_selected, center = TRUE, scale. = TRUE)

# sdev of PCS
sdev = pca_result$sdev
variance_explained = sdev^2
prop_variance_explained <- variance_explained / sum(variance_explained)
percentage_PC1_PC2 <- sum(prop_variance_explained[1:2]) * 100

PC1 = prop_variance_explained[1]
PC2 = prop_variance_explained[2]

# Extract principal components
pca_scores <- data.frame(pca_result$x)  # First 2 PCs


## added by pradip
##
pca_scores$group = c(rep("A", 45), rep("R",44))

ggplot(pca_scores, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  labs(x = PC1, y = PC2) +
  
  
  
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

## since it picked up only one wavebands; this mightnot work

## Use first derivative
# eigen vectors
parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
rownames(parents_20CS) = parents_20CS$pedigree

data_smoothed = scale(savitzkyGolay(parents_20CS[,-c(1:5)], m =1, p =1, w=11))
print(dim(data_smoothed))

library(Boruta)
library(randomForest)

set.seed(123)  # Ensure reproducibility
outcome = as.factor(c(rep("A", 45), rep("R",44)))

# Run Boruta with smoothed data
data_smoothed = as.data.frame(data_smoothed)
boruta_result <- Boruta(outcome ~ ., data = data_smoothed, doTrace = 2, 
                        maxRuns = 500, pValue = 0.05, mcAdj = FALSE)

selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)

# Subset dataset with selected features
data_selected <- data_smoothed[, selected_features, drop = FALSE]

# Print new dimensions
print(dim(data_selected))

library(MASS)
library(ggplot2)

# Ensure outcome is categorical
# outcome <- as.factor(outcome)

# Apply PCA on selected features
pca_result <- prcomp(data_selected, center = TRUE, scale. = TRUE)

# sdev of PCS
sdev = pca_result$sdev
variance_explained = sdev^2
prop_variance_explained <- variance_explained / sum(variance_explained)
percentage_PC1_PC2 <- sum(prop_variance_explained[1:2]) * 100

PC1 = prop_variance_explained[1]
PC2 = prop_variance_explained[2]

# Extract principal components
pca_scores <- data.frame(pca_result$x[, 1:2])  # First 2 PCs


## added by pradip
##
pca_scores$group = c(rep("A", 45), rep("R",44))

ggplot(pca_scores, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  labs(x = PC1, y = PC2) +
  
  
  
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


#second derivative
parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
rownames(parents_20CS) = parents_20CS$pedigree

data_smoothed = scale(savitzkyGolay(parents_20CS[,-c(1:5)], m =1, p =2, w=11))
print(dim(data_smoothed))

library(Boruta)
library(randomForest)

set.seed(123)  # Ensure reproducibility
outcome = c(rep("A",45), rep("R",44))
outcome = as.factor(outcome)
# Run Boruta with smoothed data
data_smoothed = as.data.frame(data_smoothed)
boruta_result <- Boruta(outcome ~ ., data = data_smoothed, doTrace = 2, 
                        maxRuns = 500, pValue = 0.05, mcAdj = FALSE)


boruta_result <- TentativeRoughFix(boruta_result)

# Get selected important features
selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)

# Subset dataset with selected features
data_selected <- data_smoothed[, selected_features, drop = FALSE]

# Print new dimensions
print(dim(data_selected))

library(MASS)
library(ggplot2)

# Ensure outcome is categorical
# outcome <- as.factor(outcome)

# Apply PCA on selected features
pca_result <- prcomp(data_selected, center = TRUE, scale. = TRUE)

# sdev of PCS
sdev = pca_result$sdev
variance_explained = sdev^2
prop_variance_explained <- variance_explained / sum(variance_explained)
percentage_PC1_PC2 <- sum(prop_variance_explained[1:2]) * 100

PC1 = prop_variance_explained[1]
PC2 = prop_variance_explained[2]

# Extract principal components
pca_scores <- data.frame(pca_result$x[, 1:2])  # First 2 PCs


## added by pradip
##
pca_scores$group = c(rep("A", 45), rep("R",44))

ggplot(pca_scores, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(pch = 19) +
  #geom_text(aes(label = label), vjust = 1.2, hjust = 0.9, size = 3) +
  labs(x = PC1, y = PC2) +
  
  
  
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

################### third derivative ####################
library(prospectr)
parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
rownames(parents_20CS) = parents_20CS$pedigree

data_smoothed = scale(savitzkyGolay(parents_20CS[,-c(1:5)], m =1, p =3, w=11))
print(dim(data_smoothed))

library(Boruta)
library(randomForest)

set.seed(123)  # Ensure reproducibility
outcome = c(rep("A",45), rep("R",44))
outcome = as.factor(outcome)
# Run Boruta with smoothed data
data_smoothed = as.data.frame(data_smoothed)
boruta_result <- Boruta(outcome ~ ., data = data_smoothed, doTrace = 2, 
                        maxRuns = 500, pValue = 0.05, mcAdj = FALSE)

# Get selected important features
selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)

# Subset dataset with selected features
data_selected <- data_smoothed[, selected_features, drop = FALSE]

# Print new dimensions
print(dim(data_selected))

library(MASS)
library(ggplot2)

# Ensure outcome is categorical
# outcome <- as.factor(outcome)

# Apply PCA on selected features
pca_result <- prcomp(data_selected, center = TRUE, scale. = TRUE)

# sdev of PCS
sdev = pca_result$sdev
variance_explained = sdev^2
prop_variance_explained <- variance_explained / sum(variance_explained)
percentage_PC1_PC2 <- sum(prop_variance_explained[1:2]) * 100

PC1 = prop_variance_explained[1]
PC2 = prop_variance_explained[2]

# Extract principal components
pca_scores <- data.frame(pca_result$x[,1:2])  # First 2 PCs
### LDA linear discriminant analysis 
lda_result <- lda(outcome ~ ., data = pca_scores)

# Predict LDA scores
lda_scores <- predict(lda_result)$x
lda_data <- data.frame(lda_scores, Group = outcome)

# Print first few rows to check
head(lda_data)

print(dim(lda_scores))  # Check the dimensions of LDA output
colnames(lda_scores)    # Print available column names

lda_result$svd
lda_result$svd^2 / sum(lda_result$svd^2)


ggplot(lda_data, aes(x = LD1, y = Group, color = Group)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  theme_minimal() +
  labs(title = "Discriminant PCA (DPCA) - Single LD",
       x = "Linear Discriminant 1 (LD1)",
       y = "Group") +
  theme(legend.title = element_blank())



############## This is for publication ##
## use lda for first derivative ##########

################### first derivative ####################
library(prospectr)
parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
rownames(parents_20CS) = parents_20CS$pedigree

data_smoothed = scale(savitzkyGolay(parents_20CS[,-c(1:5)], m =1, p =1, w=11))
print(dim(data_smoothed))

library(Boruta)
library(randomForest)

set.seed(123)  # Ensure reproducibility
outcome = c(rep("A",45), rep("R",44))
outcome = as.factor(outcome)
# Run Boruta with smoothed data
data_smoothed = as.data.frame(data_smoothed)
boruta_result <- Boruta(outcome ~ ., data = data_smoothed, doTrace = 2, 
                        maxRuns = 500, pValue = 0.05, mcAdj = FALSE)

# Get selected important features
selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)

# Subset dataset with selected features
data_selected <- data_smoothed[, selected_features, drop = FALSE]

# Print new dimensions
print(dim(data_selected))

library(MASS)
library(ggplot2)

# Ensure outcome is categorical
# outcome <- as.factor(outcome)

# Apply PCA on selected features
pca_result <- prcomp(data_selected, center = TRUE, scale. = TRUE)

# sdev of PCS
sdev = pca_result$sdev
variance_explained = sdev^2
prop_variance_explained <- variance_explained / sum(variance_explained)
percentage_PC1_PC2 <- sum(prop_variance_explained[1:2]) * 100

PC1 = prop_variance_explained[1]
PC2 = prop_variance_explained[2]

# Extract principal components
pca_scores <- data.frame(pca_result$x[,1:2])  # First 2 PCs
### LDA linear discriminant analysis 
lda_result <- lda(outcome ~ ., data = pca_scores)

# Predict LDA scores
lda_scores <- predict(lda_result)$x
lda_data <- data.frame(lda_scores, Group = outcome)

# Print first few rows to check
head(lda_data)

print(dim(lda_scores))  # Check the dimensions of LDA output
colnames(lda_scores)    # Print available column names

lda_result$svd
lda_result$svd^2 / sum(lda_result$svd^2)


ggplot(lda_data, aes(x = LD1, y = Group, color = Group)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  theme_minimal() +
  labs(title = "Discriminant PCA (DPCA) - Single LD",
       x = "Linear Discriminant 1 (LD1)",
       y = "Group") +
  theme(legend.title = element_blank())

# 20MA

library(prospectr)
parents_20MA = read.csv("blues_20MA_inbreds_spectra.csv", check.names = FALSE)
rownames(parents_20MA) = parents_20MA$pedigree

data_smoothed = scale(savitzkyGolay(parents_20MA[,-c(1:5)], m =1, p =1, w=11))
print(dim(data_smoothed))

library(Boruta)
library(randomForest)

set.seed(123)  # Ensure reproducibility
outcome = c(rep("A",45), rep("R",44))
outcome = as.factor(outcome)
# Run Boruta with smoothed data
data_smoothed = as.data.frame(data_smoothed)
boruta_result <- Boruta(outcome ~ ., data = data_smoothed, doTrace = 2, 
                        maxRuns = 500, pValue = 0.05, mcAdj = FALSE)

# Get selected important features
selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)

# Subset dataset with selected features
data_selected <- data_smoothed[, selected_features, drop = FALSE]

# Print new dimensions
print(dim(data_selected))

library(MASS)
library(ggplot2)

# Ensure outcome is categorical
# outcome <- as.factor(outcome)

# Apply PCA on selected features
pca_result <- prcomp(data_selected, center = TRUE, scale. = TRUE)

# sdev of PCS
sdev = pca_result$sdev
variance_explained = sdev^2
prop_variance_explained <- variance_explained / sum(variance_explained)
percentage_PC1_PC2 <- sum(prop_variance_explained[1:2]) * 100

PC1 = prop_variance_explained[1]
PC2 = prop_variance_explained[2]

# Extract principal components
pca_scores <- data.frame(pca_result$x)  # First 2 PCs
### LDA linear discriminant analysis 
lda_result <- lda(outcome ~ ., data = pca_scores)

# Predict LDA scores
lda_scores <- predict(lda_result)$x
lda_data <- data.frame(lda_scores, Group = outcome)

# Print first few rows to check
head(lda_data)

print(dim(lda_scores))  # Check the dimensions of LDA output
colnames(lda_scores)    # Print available column names

lda_result$svd
lda_result$svd^2 / sum(lda_result$svd^2)


ggplot(lda_data, aes(x = LD1, y = Group, color = Group)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  theme_minimal() +
  labs(title = "Discriminant PCA (DPCA) - Single LD",
       x = "Linear Discriminant 1 (LD1)",
       y = "Group") +
  theme(legend.title = element_blank())
