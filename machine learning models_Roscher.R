###############################################################################
# Phenomic Selection Pipeline
#
# - Loads adjusted phenotype, SNP and NIRS data
# - Builds relationship/kernel matrices
# - Defines bayesian-optimization objectives for RF, XGB and SVM
# - Loops over traits & 200 random splits, fitting 6 models on SNP, NIRS & Combined
# - Writes out per-iteration predictions and correlations
###############################################################################


###############################################################################
# 1) Environment setup
###############################################################################


# Load packages
library(sommer)        # mixed models (mmer)
library(BGLR)          # Bayesian RKHS & BL
library(ranger)        # Random Forest
library(xgboost)       # XGBoost
library(mlrMBO)        # bayesian optimization
library(smoof)         # define objective functions
library(ParamHelpers)  # parameter spaces
library(lhs)           # LHS design
library(prospectr)     # Savitzky-Golay
library(kernlab)       # SVM (ksvm)
library(DiceKriging)   # GP for MBO
library(ggplot2)       # plotting
library(dplyr)

###############################################################################
# 2) Load & prepare adjusted phenotypic data
###############################################################################

Pheno <- read.csv("pheno_yield.csv")
pheno = subset(Pheno, env == "20CS")

###############################################################################
# 3) Load & prepare SNP marker data
###############################################################################

#### loading marker data 
T92I012 <- read.delim(file = "T92I012.txt", header = TRUE, sep = "\t", dec = ".")
parents  = read.delim(file = "parents.txt", header = TRUE, sep = "\t", dec = ".")$pedigree
hybrid_marker = read.delim("T395012.txt", header = TRUE, sep = "\t", dec = ".")
hybrid_list = unique(pheno$pedigree)

marker = T92I012 %>% filter(
  X1 %in% parents)

marker2 = hybrid_marker %>% filter(
  X1 %in% hybrid_list)

rownames(marker2) = marker2$X1
marker2 = marker2[,-1]

marker_scaled <- as.data.frame(scale(marker2))

NIRS_raw = read.csv("blues_20CS_hybrids_spectra.csv", check.names = FALSE)
sg = scale(savitzkyGolay(NIRS_raw[,-c(1:7)], m=1, p=1, w=11))

NIRS_scaled <- scale(as.matrix(sg))

###############################################################################
# 5) Relationship & Distance matrices
###############################################################################

X_snp   <- as.matrix(marker_scaled);  p_snp   <- ncol(X_snp)
Dist_snp<- as.matrix(dist(X_snp))^2 / p_snp

X_nir   <- as.matrix(NIRS_scaled);    p_nir   <- ncol(X_nir)
Dist_nir<- as.matrix(dist(X_nir))^2 / p_nir

W       <- (X_nir %*% t(X_nir)) / p_nir

marker2 = as.matrix(marker2)

sg = as.data.frame(sg)
combi   <- scale(cbind(sg, marker2))

###############################################################################
# 6) Bayesian optimization functions
# objective functions: we want to minimize the RMSE by tuning hyperparameters
###############################################################################

make_rf_obj <- function(Xmat, yvec) {
  makeSingleObjectiveFunction(
    name = "rf_cv",
    fn   = function(xx) {
      rf <- ranger(x=Xmat, y=yvec,
                   num.trees=500,
                   mtry=xx["mtry"],
                   min.node.size=xx["min.node.size"],
                   num.threads=4)
      min(rf$prediction.error)
    },
    par.set = makeParamSet(
      makeIntegerParam("mtry",          lower=100, upper=floor(ncol(Xmat)/3)),
      makeIntegerParam("min.node.size", lower=3,   upper=15)
    ),
    minimize=TRUE
  )
}



obj.xgb <- makeSingleObjectiveFunction(          #seqiemto; gvm ?
  name = "xgb_cv",
  fn   = function(xx) {
    cv <- xgb.cv(
      params            = list(booster="gbtree", #?gradient boosting
                               eta=xx["eta"],
                               max_depth=xx["max_depth"],
                               min_child_weight=xx["min_child_weight"],
                               gamma=xx["gamma"],
                               subsample=xx["subsample"],
                               colsample_bytree=xx["colsample_bytree"],
                               objective="reg:squarederror"),
      data              = Z.train,
      label             = pheno.train,
      nfold             = 5,
      nrounds           = 7000,
      early_stopping_rounds = 25, # check this; maybe early ?
      verbose           = FALSE
    )
    min(cv$evaluation_log$test_rmse_mean)
  },
  par.set = makeParamSet(
    makeNumericParam("eta",               lower=0.01, upper=0.9),
    makeNumericParam("gamma",             lower=0.01, upper=5),
    makeIntegerParam("max_depth",         lower=2,    upper=20),
    makeIntegerParam("min_child_weight",  lower=1,    upper=300),
    makeNumericParam("subsample",         lower=0.1,  upper=1),
    makeNumericParam("colsample_bytree",  lower=0.1,  upper=1)
  ),
  minimize=TRUE
)

obj.svm <- makeSingleObjectiveFunction(
  name = "svm_cv",
  fn   = function(xx) {
    sv <- ksvm(x=as.matrix(train_svm), y=y.train_svm,
               kernel="rbfdot", C=xx["C"], epsilon=xx["epsilon"], cross=5)
    min(sv@cross)
  },
  par.set = makeParamSet(
    makeNumericParam("C",       lower=1e-3, upper=2^10),
    makeNumericParam("epsilon", lower=0,    upper=0.5)
  ),
  minimize=TRUE
)



###############################################################################
# 7) Bayesian optimization functions
# MBO Wrapper
###############################################################################

do_bayes <- function(n_init, n_iter, obj_fun) {
  design <- generateDesign(n_init, par.set=getParamSet(obj_fun), fun=lhs::randomLHS)
  ctrl   <- makeMBOControl(); ctrl <- setMBOControlTermination(ctrl, iters=n_iter)
  res    <- mbo(fun=obj_fun, design=design,
                learner=makeLearner("regr.km", predict.type="se", covtype="matern3_2"),
                control=ctrl, show.info=TRUE)
  df     <- as.data.frame(res$opt.path)
  df$Round<- seq_len(nrow(df))
  df$Phase<- ifelse(df$Round<=n_init, "init", "MBO")
  plot   <- ggplot(df, aes(Round, y, color=Phase)) + geom_point() +
    labs(title="MBO convergence", y="objective")
  list(result=res, plot=plot)
}


###############################################################################
# 8) Prediction Loop
###############################################################################
pheno
traits  <- colnames(pheno)[6]
set.seed(123)  # for reproducible split-ordering across traits if desired

    n         <- nrow(pheno)
    train_n   <- floor(0.8 * n)
    train_idx <- sample(n, train_n)
    test_idx  <- setdiff(seq_len(n), train_idx)
    
    pheno$y = pheno$blue
    pheno[-train_idx, ]$y <- NA
    
    # --- SNP-based predictions ---

    # (4) RF
    #marker_scaled 
    train_RF   <- marker_scaled[train_idx, ]  #using only marker
    yRF        <- pheno[train_idx,]$y
    rf_obj     <- make_rf_obj(train_RF, yRF)
    rf_run     <- do_bayes(12, 6, rf_obj)  #double check if these numbers could be altered
    bp_rf      <- rf_run$result$x
    
    rf_mod     <- ranger(y=yRF, 
                         x=train_RF,
                         num.trees=500,
                         mtry = bp_rf$mtry,
                         min.node.size = bp_rf$min.node.size)
    
    pred_rf    <- predict(rf_mod, 
                          marker2[test_idx, ])$predictions
    
    
    cor(pheno[-train_idx,]$blue, pred_rf)
    
    # marker + phenomic data 
    train_RF   <- combi[train_idx, ]  #using only marker
    yRF        <- pheno[train_idx,]$y
    rf_obj     <- make_rf_obj(train_RF, yRF)
    rf_run     <- do_bayes(12, 6, rf_obj)  #double check if these numbers could be altered
    bp_rf      <- rf_run$result$x
    
    rf_mod     <- ranger(y=yRF, 
                         x=train_RF,
                         num.trees=500,
                         mtry = bp_rf$mtry,
                         min.node.size = bp_rf$min.node.size)
    
    pred_rf    <- predict(rf_mod, 
                          marker2[test_idx, ])$predictions
    
    
    cor(pheno[-train_idx,]$blue, pred_rf)
    
    # phenomic data only 
    
    
    
    
    marker_hybrid_scaled = scale(marker2)
    # (5) XGB
    Zt         <- as.matrix(marker_hybrid_scaled[train_idx, ])
    pheno.train = pheno[train_idx, "y"]
    yXgb       <- pheno[train_idx,]$blue
    xgb_run    <- do_bayes(36, 18, obj.xgb) # error in this line
    bp_xgb     <- xgb_run$result$x
    
    xgb_mod    <- xgboost(params=bp_xgb, 
                          data=Zt, 
                          label=yXgb,
                          nrounds=7000, 
                          objective="reg:squarederror",
                          early_stopping_rounds=25, 
                          verbose=0)
    
    pred_xgb   <- predict(xgb_mod, as.matrix(marker_mother_scaled[test_idx, ]))
    
    # (6) SVM
    train_svm  <- marker_hybrid_scaled[train_idx, ]
    y.train.svm       <- Pheno[train_idx,"blue"]
    
    svm_run    <- do_bayes(12, 6, obj.svm)
    bp_svm     <- svm_run$result$x
    
    svm_mod    <- ksvm(x=as.matrix(train_svm), 
                       y=ySvm,
                       kernel="rbfdot",
                       C=bp_svm["C"], epsilon=bp_svm["epsilon"])
    
    pred_svm   <- predict(svm_mod, as.matrix(marker_mother_scaled[test_idx, ]))
    
    
    # (10) RF
    train_RF2 <- NIRS_mother_scaled[train_idx, ]
    yRF2      <- Pheno[train_idx,"trait2"]
    rf_obj2   <- make_rf_obj(train_RF2, yRF2)
    rf_run2   <- do_bayes(12,6,rf_obj2)
    bp_rf2    <- rf_run2$result$x
    rf_mod2   <- ranger(y=yRF2, x=train_RF2,
                        num.trees=500,
                        mtry=bp_rf2["mtry"],
                        min.node.size=bp_rf2["min.node.size"])
    pred_rf2  <- predict(rf_mod2, NIRS_mother_scaled[test_idx, ])$predictions
    
    # (11) XGB
    Zt3       <- as.matrix(NIRS_mother_scaled[train_idx, ])
    yXgb3     <- Pheno[train_idx, tr]
    xgb_run2  <- do_bayes(36,18,obj.xgb)
    bp_xgb2   <- xgb_run2$result$x
    xgb_mod2  <- xgboost(params=bp_xgb2, data=Zt3, label=yXgb3,
                         nrounds=7000, objective="reg:squarederror",
                         early_stopping_rounds=25, verbose=0)
    pred_xgb2 <- predict(xgb_mod2, as.matrix(NIRS_mother_scaled[test_idx, ]))
    
    # (12) SVM
    train_svm2<- NIRS_mother_scaled[train_idx, ]
    ySvm2     <- Pheno[train_idx,"trait2"]
    svm_run2  <- do_bayes(12,6,obj.svm)
    bp_svm2   <- svm_run2$result$x
    svm_mod2  <- ksvm(x=as.matrix(train_svm2), y=ySvm2,
                      kernel="rbfdot",C=bp_svm2["C"],
                      epsilon=bp_svm2["epsilon"])
    pred_svm2 <- predict(svm_mod2, as.matrix(NIRS_mother_scaled[test_idx, ]))
    
    # Assemble NIRS predictions & correlations
    Pred_NIRS <- data.frame(
      ID    = Pheno[test_idx,"ID"],
      y_obs = Pheno[test_idx,tr],
      NBLUP = pred_nblup[test_idx],
      RKHS  = pred_rkhs_nir[test_idx],
      BL    = pred_bl_nir[test_idx],
      RF    = pred_rf2,
      XGB   = pred_xgb2,
      SVM   = pred_svm2,
      seed  = seed,
      row.names = NULL
    )
    
    Corr_NIRS <- sapply(Pred_NIRS[, c("NBLUP","RKHS","BL","RF","XGB","SVM")],
                        function(v) cor(Pred_NIRS$y_obs, v, use="complete.obs"))
    
    Corr_NIRS <- data.frame(as.list(Corr_NIRS), seed = seed)
    
    write.table(Pred_NIRS,
                file = paste0("pred_NIRS/", tr, "_", seed, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(Corr_NIRS,
                file = paste0("cor_NIRS/", tr, "_", seed, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    # --- COMBINED SNP+NIRS predictions ---
    # (13) GBLUP+NBLUP
    
    combi_fit  <- mmer(trait2~1+Mother,
                       random=~vs(id,Gu=A)+vs(id2,Gu=W),
                       rcov=~units,data=Pheno,
                       getPEV=FALSE, verbose=FALSE)
    
    cg1        <- model.matrix(~id-1,Pheno) %*% combi_fit$U$`u:id`$trait2 +
      combi_fit$Beta[1,"Estimate"] +
      Pheno$Mother*combi_fit$Beta[2,"Estimate"]
    
    cg2        <- model.matrix(~id2-1,Pheno) %*% combi_fit$U$`u:id2`$trait2
    pred_com1  <- cg1 + cg2
    
    # (14) RKHS combined
    ETAc       <- c(
      list(list(K=exp(-0.5*Dist_snp),model="RKHS"),
           list(K=exp(-1  *Dist_snp),model="RKHS"),
           list(K=exp(-2.5*Dist_snp),model="RKHS")),
      list(list(K=exp(-0.5*Dist_nir),model="RKHS"),
           list(K=exp(-1  *Dist_nir),model="RKHS"),
           list(K=exp(-2.5*Dist_nir),model="RKHS"))
    )
    rkhs_com   <- BGLR(y=Pheno$trait2, ETA=ETAc,
                       nIter=5000, burnIn=1000, verbose=FALSE)
    pred_com2  <- rkhs_com$yHat
    # (15) BL combined
    bl_com     <- BGLR(y=Pheno$trait2,
                       ETA=list(
                         list(X=Pheno$Mother, model="FIXED"),
                         list(X=combi,       model="BL")
                       ),
                       nIter=5000, burnIn=1000, verbose=FALSE)
    pred_com3  <- bl_com$yHat
    # (16) RF combined
    train_RF3  <- combi[train_idx, ]
    yRF3       <- Pheno[train_idx,"trait2"]
    rf_obj3    <- make_rf_obj(train_RF3, yRF3)
    rf_run3    <- do_bayes(12,6,rf_obj3)
    bp_rf3     <- rf_run3$result$x
    rf_mod3    <- ranger(y=yRF3, x=train_RF3,
                         num.trees=500,
                         mtry=bp_rf3["mtry"],
                         min.node.size=bp_rf3["min.node.size"])
    pred_com4  <- predict(rf_mod3, combi[test_idx, ])$predictions
    # (17) XGB combined
    Zt4        <- as.matrix(combi[train_idx, ])
    ph_t4      <- Pheno[train_idx, tr]
    xgb_run3   <- do_bayes(36,18,obj.xgb)
    bp_xgb3    <- xgb_run3$result$x
    xgb_mod3   <- xgboost(params=bp_xgb3, data=Zt4, label=ph_t4,
                          nrounds=7000, objective="reg:squarederror",
                          early_stopping_rounds=25, verbose=0)
    pred_com5  <- predict(xgb_mod3, as.matrix(combi[test_idx, ]))
    # (18) SVM combined
    train_svm3 <- combi[train_idx, ]
    ySvm3      <- Pheno[train_idx,"trait2"]
    svm_run3   <- do_bayes(12,6,obj.svm)
    bp_svm3    <- svm_run3$result$x
    svm_mod3   <- ksvm(x=as.matrix(train_svm3), y=ySvm3,
                       kernel="rbfdot",C=bp_svm3["C"],
                       epsilon=bp_svm3["epsilon"])
    pred_com6  <- predict(svm_mod3, as.matrix(combi[test_idx, ]))
    
    # Assemble COMBI predictions & correlations
    Pred_COMBI <- data.frame(
      ID      = Pheno[test_idx,"ID"],
      y_obs   = Pheno[test_idx,tr],
      GBLUP   = pred_com1[test_idx],
      RKHS    = pred_com2[test_idx],
      BL      = pred_com3[test_idx],
      RF      = pred_com4,
      XGB     = pred_com5,
      SVM     = pred_com6,
      seed    = seed,
      row.names=NULL
    )
    Corr_COMBI <- sapply(Pred_COMBI[, c("GBLUP","RKHS","BL","RF","XGB","SVM")],
                         function(v) cor(Pred_COMBI$y_obs, v, use="complete.obs"))
    Corr_COMBI <- data.frame(as.list(Corr_COMBI), seed=seed)
    
    write.table(Pred_COMBI,
                file = paste0("pred_COMBI/", tr, "_", seed, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(Corr_COMBI,
                file = paste0("cor_COMBI/", tr, "_", seed, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}


######### End of the script #######