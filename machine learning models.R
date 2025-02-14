# standard normal variate (SNV) correction for NIR
# Analysis summary for paper
# perform analysis on heterotic groups using NIR spectral data
# to see if there is more or less diversity in A line or R line
# DPAC is giving one group while perofming it in R
# calculate heterosis of spectra and plot it how does it correlate with phenotypic heterosis
#

############### random forest model ############################################

# ntree = 50, 100, 200, 500, 1000
# mtry = 50, 100, 200, 300

# do the first derivative ?

# building a decision tree
install.packages("partykit")
install.packages("vip")
install.packages("radnomForestSRC")
install.packages("xgboost")
install.packages("neuralnet")
install.packages("nnet")

library(caret) # random forest
library(randomForest)
library(partykit)
library(randomForestSRC)
library(vip)
library(e1071) #svm
library(xgboost) #xgboost
library(neuralnet)
library(nnet)
library(prospectr)
#nnet for neural network


pheno_19CS = read.csv("blues_19CS_hybrids_spectra.csv", check.names = FALSE)
pheno_19CS = pheno_19CS[, -c(1:4,6:7)] 

#use first derivative of NIR for random forest model
NIR_19CS = scale(savitzkyGolay(pheno_19CS[,-c(1:7)], m=1, p=1, w=11)) 

a = cbind(pheno_19CS$gy, as.data.frame(NIR_19CS))

colnames(a)[1] = "gy"

a = as.matrix(a)
train_index <- createDataPartition(a, p = 0.8, list = FALSE)
train_data <- pheno_19CS[train_index, ]
test_data  <- pheno_19CS[-train_index, ]

train_ctrl <- trainControl(method="cv", # type of resampling in this case Cross-Validated
                           number=5, # number of folds
                           search = "grid" # we are performing a "grid-search"
)

grid_rf <- expand.grid( mtry=c(1:20) )


model_cv_grid <- train(gy ~ .,
                       data = train_data,
                       method = "rf", # this will use the randomForest::randomForest function
                       #metric = "Accuracy", # which metric should be optimized for 
                       trControl = train_ctrl, 
                       tuneGrid = grid_rf,
                       # options to be passed to randomForest
                       ntree = 741,
                       keep.forest=TRUE,
                       importance=TRUE) 


predicted = predict(model, newdata = test_data)$predicted
cor(predicted, test_data$gy)

model_cv_grid

plot(model_cv_grid)

str(model_cv_grid,1)

model_cv_grid$bestTune$mtry 

model_cv_grid$results
library(dplyr)
model_cv_grid$results %>%
  select(mtry, Accuracy, Kappa) %>%
  gather(-mtry, key="metric", value="Value") %>%
  ggplot(aes(x=mtry, y=Value, color = metric, shape=metric ) ) + 
  geom_point() + 
  geom_line()

model_cv_grid$finalModel
varImp(model_cv_grid)

model_cv_grid$finalModel$importance
randomForest::varImpPlot(model_cv_grid$finalModel)

as_tibble(model_cv_grid$finalModel$importance, rownames='Feature') %>%
  arrange(-MeanDecreaseAccuracy)

as_tibble(model_cv_grid$finalModel$importance, rownames='Feature') %>%
  arrange(-MeanDecreaseGini)
probs <- predict(model_cv_grid , TEST,"prob")
class <- predict(model_cv_grid , TEST,"raw")



model = rfsrc(gy ~.,
              data = train_data,
              importance = TRUE)

predicted = predict(model, newdata = test_data)$predicted
cor(predicted, test_data$gy)

svm_model <- svm(gy ~ ., 
                 data = train_data, 
                 kernel = "radial", cost = 1, gamma = 0.01)

predictions <- predict(svm_model, newdata = test_data)$predicted
predictions

predictions <- predict(svm_model, test_data)
cor(data.frame(predictions), test_data$gy)



# try new parameters 
# tuning

library(caret)

pheno_19CS = read.csv("blues_19CS_hybrids_spectra.csv", check.names = FALSE)
pheno_19CS = pheno_19CS[, -c(1:4,6:7)] 

train_index <- createDataPartition(pheno_19CS$gy, 
                                   p = 0.7, 
                                   list = FALSE)
train_data <- pheno_19CS[train_index, ]
test_data  <- pheno_19CS[-train_index, ]


model = rfsrc(gy ~.,
              data = train_data,
              #ntree = 200,
              #mtry = 200,
              importance = TRUE)

predicted = predict(model, newdata = test_data)$predicted
cor(predicted, test_data$gy)




##################### roughs #
library(randomForest)
randomForest()



data(iris) 
iris.rf <- randomForest(iris[,-5], 
                        iris[,5], 
                        prox=TRUE) 

iris.p <- classCenter(iris[,-5], 
                      iris[,5], 
                      iris.rf$prox) 

plot(iris[,3], iris[,4], 
     pch=21, xlab=names(iris)[3], 
     ylab=names(iris)[4], 
     bg=c("red", "blue", "green")[as.numeric(factor(iris$Species))], 
     main="Iris Data with Prototypes") 

points(iris.p[,3], iris.p[,4], pch=21, cex=2, bg=c("red", "blue", "green"))


data("iris")


iris.rf <- randomForest(iris[,-5], 
                        iris[,5], 
                        prox=TRUE) 

iris.p <- classCenter(iris[,-5], 
                      iris[,5], 
                      iris.rf$prox) 

plot(iris[,3], 
     iris[,4], 
     pch=21, 
     xlab=names(iris)[3], 
     ylab=names(iris)[4], 
     bg=c("red", "blue", "green")[as.numeric(factor(iris$Species))], 
     main="Iris Data with Prototypes") 

points(iris.p[,3], iris.p[,4], pch=21, cex=2, bg=c("red", "blue", "green"))



data(iris) 
rf1 <- randomForest(Species ~ ., 
                    iris, 
                    ntree=50, 
                    norm.votes=FALSE) 

rf2 <- randomForest(Species ~ ., 
                    iris, 
                    ntree=50, 
                    norm.votes=FALSE) 

rf3 <- randomForest(Species ~ ., 
                    iris, 
                    ntree=50, 
                    norm.votes=FALSE) 

rf.all <- combine(rf1, rf2, rf3) 
print(rf.all)



########### Montisinos lopez's book chapter 15

load("Chap_15_Data_Toy_EYT.RData")
ctree
data = read.csv("blues_19CS_hybrids_spectra.csv", check.names = FALSE)
data = data[,-c(1:4,6:8)]

library(randomForest)
model = rfs










head(model$importance)
plot(model$importance)


plot_importances <- function(importances, how_many=10,
                             color="#248a8a") {
  importances <- importances[!is.na(importances)]
  importances <- importances[importances > 0]
  n_indices <- min(length(importances), how_many)
  indices <- order(importances, decreasing = TRUE)[1:n_indices]
  best_importances <- importances[indices]
  
  X_names <- names(best_importances)
  if (is.null(X_names)) {
    X_names <- 1:length(best_importances)
  }
  X_names <- factor(X_names, levels=X_names)
  PlotData <- data.frame(value=best_importances, variable=X_names)
  Plot <- ggplot(PlotData, aes(x=variable, y=value)) +
    geom_bar(stat="identity", color=color, fill=color) +
    coord_flip()
  Plot <- apply_white_theme(Plot)
  return(Plot)
}

Plot <- plot_importances(model$importance, 10)


####################new

# Remove all variables from our workspace
rm(list=ls(all=TRUE))
library(randomForestSRC)
library(dplyr)
library(ggplot2)
# Import the own function for plotting variable importances
source("utils.R")


# Import the data set
load("Chap_15_Data_Toy_EYT.RData", verbose=TRUE)
Pheno <- Pheno_Toy_EYT
Pheno$Env <- as.factor(Pheno$Env)
Geno <- G_Toy_EYT


# Sorting data
Pheno <- Pheno[order(Pheno$Env, Pheno$GID), ]
geno_sort_lines <- sort(rownames(Geno))
Geno <- Geno[geno_sort_lines, geno_sort_lines]

### Design matrices definition ###
ZG <- model.matrix(~0 + GID, data=Pheno)


# Compute the Choleski factorization
ZL <- chol(Geno)
ZGL <- ZG %*% ZL
ZE <- model.matrix(~0 + Env, data=Pheno)

# Bind all design matrices in a single matrix to be used as predictor
X <- cbind(ZE,ZGL)
dim(X)

# Create a data frame with the information of response variable and all
# predictors
Data <- data.frame(y=Pheno$GY, X)
head(Data[, 1:5])


# Fit the model with importance=TRUE for also computing the variable
importance
model <- rfsrc(y ~ ., data=Data, importance=TRUE)

# Get the variable importance
head(model$importance,10)

# Plot the 30 most important variables with own function
Plot <- plot(model$importance)
Plot
save_plot(Plot, "plots/1.continuous.png")



model <- rfsrc(y ~ ., data=Data, 
               ntree = 200,
               mtry = 200,
               importance=TRUE)






iris.rf <- randomForest(pheno_19CS[,-c(1:7)], 
                        pheno_19CS[,5], 
                        prox=TRUE) 


iris.p <- classCenter(iris[,-5], 
                      iris[,5], 
                      iris.rf$prox) 



iris.p <- classCenter(pheno_19CS[,-c(1:7)], 
                      pheno_19CS[,5], 
                      iris.rf$prox) 


randomForest()






# DPAC
# DAPC requires the adegenet package. Let's load this package:
library("adegenet")
data(H3N2) # load the H3N2 influenza data. Type ?H3N2 for more info.
pop(H3N2) <- H3N2$other$epid
dapc.H3N2 <- dapc(H3N2, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(H3N2) - 1)
scatter(dapc.H3N2, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)



parents_20CS
parents_20CS$group = c(rep("A",45), rep("R",44))


a = dapc(parents_20CS[,-c(1:5,4206)], 
     var.contrib = TRUE, 
     scale = FALSE, 
     n.pca = 30, 
     grp = parents_20CS$group)
     
scatter(a)


n.da = nPop(H3N2) - 1)



library(ggplot2)
pheno_19CS = read.csv("blues_19CS_hybrids_spectra.csv", check.names = FALSE)
NIR_19CS = t(pheno_19CS[,-c(1:7)][1,])
NIR.d1 = t(as.data.frame(savitzkyGolay(pheno_19CS[1,-c(1:7)], m =1, p =1, w =11)))
NIR.d1.scaled = t(as.data.frame(scale(savitzkyGolay(pheno_19CS[,-c(1:7)], m=1, p=1, w=11)))[1,])

NIR_19CS = as.data.frame(cbind(as.numeric(rownames(NIR_19CS)), NIR_19CS[,1]))
NIR.d1 = as.data.frame(cbind(as.numeric(rownames(NIR.d1)), NIR.d1[,1]))
NIR.d1.scaled = as.data.frame(cbind(as.numeric(rownames(NIR.d1.scaled)), NIR.d1.scaled[,1]))

colnames(NIR_19CS)[1] = "Wavelength"
colnames(NIR_19CS)[2] = "Reflectance"

colnames(NIR.d1)[1] = "Wavelength"
colnames(NIR.d1)[2] = "Reflectance"

colnames(NIR.d1.scaled)[1] = "Wavelength"
colnames(NIR.d1.scaled)[2] = "Reflectance"

###plotting representative of bands

p1 = ggplot(NIR_19CS, aes(Wavelength, Reflectance)) +
  geom_line(stat="identity", position=position_dodge())+
  theme_bw()+
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
    legend.position = "none"
    #axis.text.x.bottom = element_blank()
  )


p2 = ggplot(NIR.d1, aes(Wavelength, Reflectance)) +
  geom_line(stat="identity", position=position_dodge())+
  theme_bw()+
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
    legend.position = "none"
    #axis.text.x.bottom = element_blank()
  )

p3 = ggplot(NIR.d1.scaled, aes(Wavelength, Reflectance)) +
  geom_line(stat="identity", position=position_dodge())+
  theme_bw()+
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
    legend.position = "none"
    #axis.text.x.bottom = element_blank()
  )

library(ggpubr)

merge = ggarrange(p1,p2,p3,
                  ncol = 3, nrow = 1,
                  common.legend = TRUE,
                  legend = c("bottom"))

jpeg("spec_smoothing.jpeg",width = 14,height =6,units = "in", res=600)
merge
dev.off()









ggplot(NIR,aes(Wavelength, Reflectance))+
  geom_line()

ggplot(NIR1,aes(Wavelength, Reflectance))+
  geom_line()


a$wave = as.numeric(rownames(a))






############
data(NIRsoil)
opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))

# plot of the 10 first spectra
matplot(as.numeric(colnames(NIRsoil$spc)),
        t(NIRsoil$spc[1:10, ]),
        type = "l",
        xlab = "",
        ylab = "Absorbance"
)

mtext("Raw spectra")
NIRsoil$spc_sg <- savitzkyGolay(
  X = NIRsoil$spc,
  m = 1,
  p = 3,
  w = 11,
  delta.wav = 2
)

matplot(as.numeric(colnames(NIRsoil$spc_sg)),
        t(NIRsoil$spc_sg[1:10, ]),
        type = "l",
        xlab = "Wavelength /nm",
        ylab = "1st derivative"
)

mtext("1st derivative spectra")
par(opar)
ggplot(a)
table(rpois(100, 5)), type = "h", col = "red", lwd = 10,
main = "rpois(100, lambda = 5)")


library(caret)
?caret

###############
# Neural Networks Using Model Averaging
avNNet()



data(BloodBrain)
## Not run: 
modelFit <- avNNet(bbbDescr, 
                   logBBB, 
                   size = 5, 
                   linout = TRUE, 
                   trace = FALSE)
modelFit

predict(modelFit, bbbDescr)

# support vector machine model 
library(e1071)
library(kernlab)

?svm() #svm function to fit svm model

#glmnet to fit ridge, lasso, and elastic net

#random forest
library(randomForest)

#model = randomForest(target ~., data = data_train)
#model
#importance(model)             # returns the importance of the variables: most significant -cp followed by thalach and so on     
#varImplot(model)              #visualizing the importance of variables of the model
#pred_test <- predict(model, newdata -data_test, type= "class")
#confusionMatrix(table(pred_test,data_test$target)) # The prediction to compute the confusion matrix and see the accuracy score 

#radnom forest
library(randomForest)
set.seed(1)
install.packages("ISLR2")
library(ISLR2)

attach(Boston)
train <- sample(1:nrow(Boston), nrow(Boston) / 2)

bag.boston = randomForest(medv ~.,
                          data = Boston,
                          subset = train,
                          mtry =12,         #12 predictors should be considered  for each split of the tree
                          importance=TRUE)

yhat.bag <- predict(bag.boston , newdata = Boston[-train , ])
boston.test <- Boston[-train, "medv"]
plot(yhat.bag, boston.test)


# support vector machines
# classification model/approach : intended for binary classification
# supervised machine learning model
# Kernel trick/ usually for pattern detection sucha as features in image
#

# support vector classfiers
library(e1071)
?svm()
# svm() can be used when kernel = "linear"


set.seed(1)
x <- matrix(rnorm(20 * 2), ncol = 2)
y <- c(rep(-1, 10), rep(1, 10))
x[y == 1, ] <- x[y == 1, ] + 1
plot(x, col = (3 - y))


# fitting the support vector classifer

dat <- data.frame(x = x, y = as.factor(y))
library(e1071)

svmfit <- svm(y ∼., 
              data = dat, 
              kernel = "linear",
              cost = 10, 
              scale = FALSE)

# need to figureout, couldnot know


## Deep learning approaches
# neural network

# single layer neural networks
library(ISLR2)
Gitters <- na.omit(Hitters)
n <- nrow(Gitters)
set.seed (13)
ntest <- trunc(n / 3)
testid <- sample(1:n, ntest)

lfit = lm(Salary ~ ., 
           data = Gitters[-testid, ])


lpred <- predict(lfit , Gitters[testid , ])
with(Gitters[testid , ], mean(abs(lpred - Salary)))


x <- scale(model.matrix(Salary ~ . - 1, data = Gitters))
y <- Gitters$Salary


library(glmnet)

cvfit <- cv.glmnet(x[-testid , ], y[-testid],
                     type.measure = "mae")
cpred <- predict(cvfit , x[testid , ], s = "lambda.min")
mean(abs(y[testid] - cpred))


## fitting neural network
#install_tensorflow()
library(keras)
library(dplyr)

modnn <- keras_model_sequential() %>%
  layer_dense(units = 50, activation = "relu",
              input_shape = ncol(x)) %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 1)

install.packages("tensorflow")
install.packages("reticulate")

library(tensorflow)
library(reticulate)
x <- scale(model.matrix(Salary ~ . - 1, data = Gitters))

x <- model.matrix(Salary ~ . - 1, data = Gitters) %>% scale()

modnn %>% compile(loss = "mse",
                  optimizer = optimizer_rmsprop(),
                  metrics = list("mean_absolute_error")
)


?smoof
install.packages("smoof")
library(smoof) # for svm package
install.packages("mlrMBO")
library(mlrMBO)
install.packages("ParamHelpers")
library(ParamHelpers)
?Param

# for svm
library(kernlab)

#hyperparameters optimization with Bayesian model based optimization approach
install.packages("tuneRanger")  
library(tuneRanger)

attach(iris)
result <- tuneRanger(Species ~ ., 
                     data = iris, 
                     num.trees = 500, 
                     iters = 10)
print(result$recommended.pars)  # Best hyperparameter set


library(tuneRanger)

result <- tuneRanger(iris, num.trees = 500, iters = 10)
print(result$recommended.pars)




######### support vector machine
# supervised learning ; 
# two types: classification (predict categories), regression (predict values)
# svm is binary classifiers,
# example inputs and desired outputs and the goal is to learn a general rule that map inputs to outputs
# detecting car vs truck 

# svm is classficiation model classify based on features
# draw lines between classes; we want a line of seperation between two groups


library(e1071)

#### Bayesian Model based optimization in R using smoof pacakge

#install.packages(c("smoof", "mlrMBO", "mlr", "DiceKriging", "rgenoud"))
install.packages("rgenoud")
library(smoof)
library(mlrMBO)
library(mlr)
library(DiceKriging)  # For Gaussian process modeling
library(rgenoud)  # Genetic algorithms for optimization


# Define objective function
obj_function <- makeSingleObjectiveFunction(
  name = "Simple Function",
  fn = function(x) (x - 2)^2 + sin(5 * x),
  par.set = makeNumericParamSet(lower = -2, upper = 5, len = 1),  # Search space
  minimize = TRUE  # We want to minimize the function
)

# Define the Bayesian Optimization Control Settings
ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters = 10)  # Number of iterations
ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())  # Use Expected Improvement (EI)


# Set up the Surrogate Model
surrogate <- makeLearner("regr.km", 
                         predict.type = "se", 
                         covtype = "matern3_2")


# Run Bayesian Optimization
set.seed(123)
result <- mbo(
  fun = obj_function,
  learner = surrogate,
  control = ctrl
)

# Results
print(result)
plot(result)  # Visualize optimization progress


# Adapting hyperparameter tuning

library(e1071)  # For SVM

svm_objective <- makeSingleObjectiveFunction(
  name = "SVM Tuning",
  fn = function(x) {
    model <- svm(Species ~ ., 
                 data = iris, 
                 cost = x[1], 
                 gamma = x[2], 
                 cross = 5)
    return(-mean(model$accuracies))  # Negate accuracy since `smoof` minimizes
  },
  par.set = makeParamSet(
    makeNumericParam("cost", lower = 0.1, upper = 10),
    makeNumericParam("gamma", lower = 0.001, upper = 1)
  ),
  minimize = TRUE
)

# Bayesian Optimization like before

result_svm <- mbo(
  fun = svm_objective,
  learner = surrogate,
  control = ctrl
)

print(result_svm)

library(e1071)
library(caret)
library(dplyr)


pheno_19CS[1:10,1:10]
pheno_19CS[]

nir_data = pheno_19CS[,-c(1:4,6,7)]

#nir_data$Class <- sample(c("Group1", "Group2"), nrow(nir_data), replace = TRUE)

set.seed(123)
nir_data$Class <- sample(c("Group1", "Group2"), nrow(nir_data), replace = TRUE)


train_index <- createDataPartition(nir_data$Class, p = 0.8, list = FALSE)
train_data <- nir_data[train_index, ]
test_data  <- nir_data[-train_index, ]


svm_model <- svm(gy ~ ., 
                 data = train_data[, -4202], 
                 kernel = "radial", cost = 1, gamma = 0.01)

predictions <- predict(svm_model, test_data[, -4202])
predictions

confusionMatrix(predictions, test_data$Class)

a = as.data.frame(predictions)
cor(a$predictions, test_data[,1])


## random forest

rf_objective <- makeSingleObjectiveFunction(
  name = "RF Hyperparameter Optimization",
  fn = function(x) {
    rf_model <- randomForest(gy ~ ., data = train_data[,-4202], 
                             ntree = as.integer(x[1]), 
                             mtry = as.integer(x[2]), 
                             nodesize = as.integer(x[3]), 
                             importance = FALSE, 
                             keep.inbag = TRUE)
    return(mean(rf_model$mse))  # OOB error
  },
  par.set = makeParamSet(
    makeIntegerParam("ntree", lower = 50, upper = 500),
    makeIntegerParam("mtry", lower = 1, upper = 4),
    makeIntegerParam("nodesize", lower = 1, upper = 10)
  ),
  minimize = TRUE
)


svm_objective <- makeSingleObjectiveFunction(
  name = "SVM Hyperparameter Optimization",
  fn = function(x) {
    svm_model <- svm(gy ~ ., data = train_data[,-4202], 
                     cost = x[1], 
                     gamma = x[2])
    
    preds <- predict(svm_model, test_data)
    rmse <- sqrt(mean((test_data$Species - preds)^2))  # Calculate RMSE
    return(rmse)  # Minimize RMSE
  },
  par.set = makeParamSet(
    makeNumericParam("cost", lower = 0.1, upper = 10),
    makeNumericParam("gamma", lower = 0.001, upper = 1)
  ),
  minimize = TRUE
)


ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters = 6)  # Optimize for 6 iterations
ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())  # Use Expected Improvement (EI)

surrogate <- makeLearner("regr.km", predict.type = "se", covtype = "matern3_2")  # Gaussian Process

library(randomForest)
set.seed(123)
rf_result <- mbo(
  fun = rf_objective,
  learner = surrogate,
  control = ctrl
)
print(rf_result)


set.seed(123)
svm_result <- mbo(
  fun = svm_objective,
  learner = surrogate,
  control = ctrl
)
print(svm_result)


