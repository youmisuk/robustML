# :::::::::::::::::::::::
# ::::: DR functions :::: ####
# :::::::::::::::::::::::

erLearnPR <- function(Z, X, ID, library=c("glm", "deeplearning"), int=T, rate=0.01, ...) {
  
  if(length(library) != sum(library %in% c("glm", "deeplearning", "gbm", "randomForests"))) {
    stop("library is not in the appropriate format. Choose methods among glm, gbm, randomForests, and deeplearning.")
  }
  
  # variables 
  IDfactor = as.factor(ID)
  IDnum = as.numeric(IDfactor)
  Zfactor = as.factor(Z)
  Z = as.numeric(as.character(Z))
  
  lvl2var <- colSums(apply(as.matrix(X), 2, function(x) tapply(x, IDnum, function(x) length(unique(x)))) == 1) == length(levels(IDfactor))
  X_lvl1 <- as.matrix(X)[, !lvl2var]
  W_lvl2 <- as.matrix(X)[, lvl2var] 
  X_lvl1_num <- as.matrix(as.matrix(X_lvl1)[, apply(as.matrix(X_lvl1), 2, function(x) sum(!unique(x) %in% c(0, 1))) > 0])
  X_lvl1_mean <- apply(as.matrix(X)[, !lvl2var], 2, function(y) ave(y, IDnum))
  
  X_lvl1_num_sq <- X_lvl1_num^2
  colnames(X_lvl1_num_sq) <- paste0(colnames(X_lvl1_num_sq), "_sq")
  
  # datasets
  h2oDataE= as.h2o(data.frame(Zfactor, X_lvl1, W_lvl2, IDfactor, X_lvl1_num_sq))   
  
  fitER_glm <- fitER_dl <- fitER_gbm <- fitER_rf <- Z.hat_glm <- Z.hat_dl <- Z.hat_gbm <- Z.hat_rf <-  NULL
  
 if (int == F) {
  if ("glm" %in% library) {
    # ensemble glm for propensity score training
    if(length(2:(ncol(h2oDataE) -ncol(X_lvl1_num) -1))>=2) {
      fitER_glm = h2o.glm(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataE,alpha=1,lambda=0,family="binomial",interactions=2:(ncol(h2oDataE) -ncol(X_lvl1_num)-1))
    } else {
      fitER_glm = h2o.glm(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataE,alpha=1,lambda=0,family="binomial")
    }
    
    Z.hat_glm = as.numeric(as.data.frame(predict(fitER_glm,h2oDataE)[,"p1"])$p1)
    
   }
  } else {
    
    X_lvl1_ID_int <- model.matrix( ~ X_lvl1_num:IDfactor-1)
    
    h2oDataE1= as.h2o(data.frame(Zfactor, X_lvl1, W_lvl2, IDfactor, X_lvl1_num_sq, X_lvl1_ID_int)) 
    
    if ("glm" %in% library) {
      # ensemble glm for propensity score training
      if(length(2:(ncol(h2oDataE1) -ncol(X_lvl1_num) -ncol(X_lvl1_ID_int) -1))>=2) {
        fitER_glm = h2o.glm(x=2:ncol(h2oDataE1),y=1,training_frame=h2oDataE1,alpha=1,lambda=0,family="binomial",interactions=2:(ncol(h2oDataE1) -ncol(X_lvl1_num)-ncol(X_lvl1_ID_int)-1))
      } else {
        fitER_glm = h2o.glm(x=2:ncol(h2oDataE1),y=1,training_frame=h2oDataE1,alpha=1,lambda=0,family="binomial")
      }
      Z.hat_glm = as.numeric(as.data.frame(predict(fitER_glm,h2oDataE1)[,"p1"])$p1)
    }
  }
  
  
  if ("deeplearning" %in% library) {
    fitER_dl = h2o.deeplearning(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataE,hidden= rep((ncol(h2oDataE) -2) + nlevels(IDfactor)+150,2), epochs=500,adaptive_rate=FALSE,
                                l1=1e-8, l2=1e-8, rate=rate, rate_annealing=1e-6, max_w2=10, ...)
    Z.hat_dl = as.numeric(as.data.frame(predict(fitER_dl,h2oDataE)[,"p1"])$p1)
  }
  
  if ("gbm" %in% library) {
    fitER_gbm = h2o.gbm(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataE,ntrees=8,max_depth=5,min_rows=min(5,min(table(IDfactor))/2))
    Z.hat_gbm = as.numeric(as.data.frame(predict(fitER_gbm,h2oDataE)$pred)$predict) 
  }
  
  if ("randomForests" %in% library) {
    fitER_rf = h2o.randomForest(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataE,ntrees = 500, max_depth = 50, stopping_rounds = 2, stopping_tolerance = 1e-2, score_each_iteration = T)
    Z.hat_rf = as.numeric(as.data.frame(predict(fitER_rf,h2oDataE)$pred)$predict)
  }
  
  
  # all predictions
  Z.hats <- cbind(Z.hat_glm, Z.hat_dl, Z.hat_gbm, Z.hat_rf)
  
  # nnls
  fit.nnlsER <- nnls(A=Z.hats, b=Z)
  initCoefER <- coef(fit.nnlsER)
  initCoefER[is.na(initCoefER)] <- 0.0
  
  # normalize so sum(coef) = 1 if possible
  if (sum(initCoefER) > 0) {
    coefER <- initCoefER / sum(initCoefER)
  } else {
    warning("All algorithms have zero weight", call. = FALSE)
    coefER <- initCoefER
  }
  
  
  Z.comb <- as.numeric(crossprod(t(Z.hats[, coefER != 0, drop = FALSE]), coefER[coefER != 0]))
  
  return(list(coef=coefER, 
              Z.hat = Z.comb,
              Z.hats = data.frame(Z.hats)))
}


orLearnPR <- function(Y, Z, X, ID, library=c("glm", "deeplearning"), int=T, rate=0.01, ...) {
  
  if(length(library) != sum(library %in% c("glm", "deeplearning", "gbm", "randomForests"))) {
    stop("library is not in the appropriate format. Choose methods among glm, gbm, randomForests, and deeplearning.")
  }
  
  # variables
  IDfactor = as.factor(ID)
  IDnum = as.numeric(IDfactor)
  
  lvl2var <- colSums(apply(as.matrix(X), 2, function(x) tapply(x, IDnum, function(x) length(unique(x)))) == 1) == length(levels(IDfactor))
  X_lvl1 <- as.matrix(X)[, !lvl2var]
  W_lvl2 <- as.matrix(X)[, lvl2var] 
  X_lvl1_num <- as.matrix(as.matrix(X_lvl1)[, apply(as.matrix(X_lvl1), 2, function(x) sum(!unique(x) %in% c(0, 1))) > 0])
  X_lvl1_mean <- apply(as.matrix(X)[, !lvl2var], 2, function(y) ave(y, IDnum))
  Z_mean <- ave(Z, IDnum)
  
  X_lvl1_num_sq <- X_lvl1_num^2
  colnames(X_lvl1_num_sq) <- paste0(colnames(X_lvl1_num_sq), "_sq")
  
  # datasets
  h2oAll = as.h2o(data.frame(Y,Z,X_lvl1, W_lvl2 ,IDfactor,X_lvl1_num_sq))  
  h2oAll_Z1 = h2oAll_Z0 = h2oAll; h2oAll_Z1$Z = 1; h2oAll_Z0$Z = 0
  
  fitOR_glm <- fitOR_dl <- fitOR_gbm <- fitOR_rf <- Y.hat_glm <- Y.hat_dl <- Y.hat_gbm <- Y.hat_rf <- 
    Y1.hat_glm <- Y1.hat_dl <- Y1.hat_gbm <- Y1.hat_rf <-
    Y0.hat_glm <- Y0.hat_dl <- Y0.hat_gbm <- Y0.hat_rf <- NULL
  
 if (int ==F) {
  if ("glm" %in% library) {
    # ensemble glm for outcome reg training
    if(length(2:(ncol(h2oAll) -ncol(X_lvl1_num)-1))>=2) {
      fitOR_glm = h2o.glm(x=2:ncol(h2oAll),y=1,training_frame=h2oAll,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oAll) -ncol(X_lvl1_num)-1))  
    } else {
      fitOR_glm = h2o.glm(x=2:ncol(h2oAll),y=1,training_frame=h2oAll,alpha=1,lambda=0,family="gaussian")
    }
    
    Y.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll)$pred)$predict)
    Y1.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll_Z1)$pred)$predict)
    Y0.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll_Z0)$pred)$predict)
   } 
  
  } else {
    
    X_lvl1_ID_int <- model.matrix( ~ X_lvl1_num:IDfactor-1)
    
    h2oAll1 = as.h2o(data.frame(Y,Z,X_lvl1, W_lvl2 ,IDfactor,X_lvl1_num_sq, X_lvl1_ID_int))
    h2oAll1_Z1 = h2oAll1_Z0 = h2oAll1; h2oAll1_Z1$Z = 1; h2oAll1_Z0$Z = 0
    
    if ("glm" %in% library) {
      # ensemble glm for outcome reg training
      if(length(2:(ncol(h2oAll1) -ncol(X_lvl1_num)-ncol(X_lvl1_ID_int)-1))>=2) {
        fitOR_glm = h2o.glm(x=2:ncol(h2oAll1),y=1,training_frame=h2oAll1,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oAll1) -ncol(X_lvl1_num)-ncol(X_lvl1_ID_int)-1))  
      } else {
        fitOR_glm = h2o.glm(x=2:ncol(h2oAll1),y=1,training_frame=h2oAll1,alpha=1,lambda=0,family="gaussian")
      }
      
      Y.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll1)$pred)$predict)
      Y1.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll1_Z1)$pred)$predict)
      Y0.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll1_Z0)$pred)$predict)
    }
   }

  
  if ("deeplearning" %in% library) {
    fitOR_dl = h2o.deeplearning(x=2:ncol(h2oAll),y=1,training_frame=h2oAll,hidden=  rep((ncol(h2oAll) -2) + nlevels(IDfactor)+150,2), epochs=500,adaptive_rate=FALSE,
                                l1=1e-8, l2=1e-8, rate=rate, rate_annealing=1e-6, max_w2=10, ...)  
    Y.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAll)$pred)$predict)
    Y1.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAll_Z1)$pred)$predict)
    Y0.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAll_Z0)$pred)$predict)
  }
  
  if ("gbm" %in% library) {
    fitOR_gbm = h2o.gbm(x=2:ncol(h2oAll),y=1,training_frame=h2oAll, ntrees=8,max_depth=5,min_rows=min(5,min(table(IDfactor))/2))
    Y.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAll)$pred)$predict)
    Y1.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAll_Z1)$pred)$predict)
    Y0.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAll_Z0)$pred)$predict)
  }
  
  if ("randomForests" %in% library) {
    fitOR_rf = h2o.randomForest(x=2:ncol(h2oAll),y=1,training_frame=h2oAll,ntrees = 500, max_depth = 50, stopping_rounds = 2, stopping_tolerance = 1e-2, score_each_iteration = T)
    Y.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAll)$pred)$predict)
    Y1.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAll_Z1)$pred)$predict)
    Y0.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAll_Z0)$pred)$predict)  
  }
  
  # all predictions
  Y.hats <- cbind(Y.hat_glm, Y.hat_dl, Y.hat_gbm, Y.hat_rf)
  Y1.hats <- cbind(Y1.hat_glm, Y1.hat_dl, Y1.hat_gbm, Y1.hat_rf)
  Y0.hats <- cbind(Y0.hat_glm, Y0.hat_dl, Y0.hat_gbm, Y0.hat_rf)
  
  # nnls
  fit.nnlsOR <- nnls(A=Y.hats, b=Y)
  initCoefOR <- coef(fit.nnlsOR)
  initCoefOR[is.na(initCoefOR)] <- 0.0
  
  # normalize so sum(coef) = 1 if possible
  if (sum(initCoefOR) > 0) {
    coefOR <- initCoefOR / sum(initCoefOR)
  } else {
    warning("All algorithms have zero weight", call. = FALSE)
    coefOR <- initCoefOR
  }
  
  Y1.comb <- as.numeric(crossprod(t(Y1.hats[, coefOR != 0, drop = FALSE]), coefOR[coefOR != 0]))
  Y0.comb <- as.numeric(crossprod(t(Y0.hats[, coefOR != 0, drop = FALSE]), coefOR[coefOR != 0]))
  
  return(list(coef=coefOR, Y1.hat = Y1.comb, Y0.hat = Y0.comb,
              Y1.hats = data.frame(Y1.hats), Y0.hats = data.frame(Y0.hats)))
  
}

DR <- function(Y, Z, interZ=formula(~1), Z.hat, Y1.hat, Y0.hat, data) {
  
  DRest.ind = (Z * (Y- Y1.hat) / Z.hat + Y1.hat) - ((1 -Z) * (Y - Y0.hat) / (1 - Z.hat) + Y0.hat)
  
  interZ.mat <- model.matrix(interZ, data=data)
  
  ATE_interZ_model = lm(DRest.ind ~ interZ.mat -1)
  SE_interZ = sqrt(diag(vcov(ATE_interZ_model)))
  ATE_interZ =coef(ATE_interZ_model)
  tau.est <- cbind(Estimate=ATE_interZ,SE=SE_interZ)
  
  colnames(tau.est) <-  c("Estimate", "Std. Error")
  rownames(tau.est) <- colnames(interZ.mat)
  
  return(tau.est)
}


# :::::::::::::::::::::::
# ::::: DD functions :::: ####
# :::::::::::::::::::::::


erLearnDD <- function(Z, X, ID, library=c("glm", "deeplearning"), int=T, rate=0.01, ...) {
  
  if(length(library) != sum(library %in% c("glm", "deeplearning", "gbm", "randomForests"))) {
    stop("library is not in the appropriate format. Choose methods among glm, gbm, randomForests, and deeplearning.")
  }
  
  # variables
  IDfactor = as.factor(ID)
  IDnum = as.numeric(IDfactor)
  
  Z_demeaned = Z - ave(Z,IDnum )
  Z_mean = ave(Z,IDnum )
  X_demeaned = apply(as.matrix(X), 2, function(y) y - ave(y, IDnum))
  
  lvl2var <-  colSums(apply(as.matrix(X), 2, function(x) tapply(x, IDnum, function(x) length(unique(x)))) == 1) == length(levels(IDfactor))
  X_lvl1 <- as.matrix(X)[, !lvl2var]
  X_demeaned_lvl1 <-  as.matrix(X_demeaned)[, !lvl2var]
  X_demeaned_lvl1_num <- as.matrix(as.matrix(X_demeaned_lvl1)[, apply(as.matrix(X_lvl1), 2, function(x) sum(!unique(x) %in% c(0, 1))) > 0])
  
  X_lvl1_mean <- apply(as.matrix(X)[, !lvl2var], 2, function(y) ave(y, IDnum))
  W_lvl2 <- as.matrix(X)[, lvl2var]
  
  # datasets
  h2oDataE = as.h2o(data.frame(Z_demeaned,X_demeaned_lvl1, X_lvl1_mean, W_lvl2, X_demeaned_lvl1_num^2))
  
  fitER_glm <- fitER_dl <- fitER_gbm <- fitER_rf <- Z.hat_glm <- Z.hat_dl <- Z.hat_gbm <- Z.hat_rf <-  NULL
  
  if (int == F) {
  if ("glm" %in% library) {
    fitER_glm = h2o.glm(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataE,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oDataE) - ncol(X_demeaned_lvl1_num)))
    Z.hat_glm = as.numeric(as.data.frame(predict(fitER_glm,h2oDataE)$pred)$predict)
   }
  } else {
    
    X_lvl1_num <- as.matrix(as.matrix(X_lvl1)[, apply(as.matrix(X_lvl1), 2, function(x) sum(!unique(x) %in% c(0, 1))) > 0])
    X_lvl1_ID_int <- model.matrix( ~ X_lvl1_num:IDfactor-1)
    X_lvl1_ID_int_dm <- apply(X_lvl1_ID_int, 2, function(y) y - ave(y, IDnum))
    
    h2oDataE1 = as.h2o(data.frame(Z_demeaned,X_demeaned_lvl1, X_lvl1_mean, W_lvl2, X_demeaned_lvl1_num^2, X_lvl1_ID_int_dm))
    if ("glm" %in% library) {
      fitER_glm = h2o.glm(x=2:ncol(h2oDataE1),y=1,training_frame=h2oDataE1,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oDataE1) - ncol(X_demeaned_lvl1_num) - ncol(X_lvl1_ID_int_dm)))
      Z.hat_glm = as.numeric(as.data.frame(predict(fitER_glm,h2oDataE1)$pred)$predict)
    }
  }
  
  if ("deeplearning" %in% library) {
    fitER_dl = h2o.deeplearning(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataE,hidden= rep((ncol(h2oDataE) -2) + nlevels(IDfactor)+150,2),epochs=500,adaptive_rate=FALSE,
                                l1=0, l2=1e-8, rate=rate, rate_annealing=1e-6, max_w2=10, ...)
    Z.hat_dl = as.numeric(as.data.frame(predict(fitER_dl,h2oDataE)$pred)$predict) 
  }
  
  if ("gbm" %in% library) {
    fitER_gbm = h2o.gbm(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataE,ntrees=8,max_depth=5,min_rows=min(5,min(table(IDfactor))/2))
    Z.hat_gbm = as.numeric(as.data.frame(predict(fitER_gbm,h2oDataE)$pred)$predict) 
  }
  
  if ("randomForests" %in% library) {
    fitER_rf = h2o.randomForest(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataE,ntrees = 500, max_depth = 50, stopping_rounds = 2, stopping_tolerance = 1e-2, score_each_iteration = T)
    Z.hat_rf = as.numeric(as.data.frame(predict(fitER_rf,h2oDataE)$pred)$predict)
  }
  
  # all predictions
  Z.hats <- cbind(Z.hat_glm, Z.hat_dl, Z.hat_gbm, Z.hat_rf)
  
  # nnls
  fit.nnlsER <- nnls(A=Z.hats, b=Z_demeaned)
  initCoefER <- coef(fit.nnlsER)
  initCoefER[is.na(initCoefER)] <- 0.0
  
  # normalize so sum(coef) = 1 if possible
  if (sum(initCoefER) > 0) {
    coefER <- initCoefER / sum(initCoefER)
  } else {
    warning("All algorithms have zero weight", call. = FALSE)
    coefER <- initCoefER
  }
  
  Z.comb <- as.numeric(crossprod(t(Z.hats[, coefER != 0, drop = FALSE]), coefER[coefER != 0]))
  
  return(list(coef = coefER, Z.hat = Z.comb, Z.hats = data.frame(Z.hats)))
}


orLearnDD <- function(Y, Z, X, ID, library=c("glm", "deeplearning"), int=T, rate=0.01, ...) {
  
  if(length(library) != sum(library %in% c("glm", "deeplearning", "gbm", "randomForests"))) {
    stop("library is not in the appropriate format. Choose methods among glm, gbm, randomForests, and deeplearning.")
  }
  
  IDfactor = as.factor(ID)
  IDnum = as.numeric(IDfactor)
  
  # group deviation & group mean
  Y_demeaned = Y - ave(Y,IDnum )
  Z_demeaned = Z - ave(Z,IDnum )
  Z_mean = ave(Z,IDnum )
  X_demeaned = apply(as.matrix(X), 2, function(y) y - ave(y, IDnum))
  
  # Data structures for h2o
  lvl2var <-  colSums(apply(as.matrix(X), 2, function(x) tapply(x, IDnum, function(x) length(unique(x)))) == 1) == length(levels(IDfactor))
  X_lvl1 <- as.matrix(X)[, !lvl2var]
  X_demeaned_lvl1 <-  as.matrix(X_demeaned)[, !lvl2var]
  X_demeaned_lvl1_num <- as.matrix(as.matrix(X_demeaned_lvl1)[, apply(as.matrix(X_lvl1), 2, function(x) sum(!unique(x) %in% c(0, 1))) > 0])
  
  X_lvl1_mean <- apply(as.matrix(X)[, !lvl2var], 2, function(y) ave(y, IDnum))
  W_lvl2 <- as.matrix(X)[, lvl2var]
 
  h2oAll = as.h2o(data.frame(Y_demeaned,Z_demeaned,X_demeaned_lvl1, X_lvl1_mean, W_lvl2, X_demeaned_lvl1_num^2))
  h2oAll_Z1 = h2oAll_Z0 = h2oAll 
  h2oAll_Z1$Z_demeaned = as.h2o(rep(1,length(Y_demeaned)) - ave(Z,IDnum)); h2oAll_Z0$Z_demeaned = as.h2o(rep(0,length(Y_demeaned)) - ave(Z,IDnum))
  
  fitOR_glm <- fitOR_dl <- fitOR_gbm <- fitOR_rf <- Y.hat_glm <- Y.hat_dl <- Y.hat_gbm <- Y.hat_rf <- 
    Y1.hat_glm <- Y1.hat_dl <- Y1.hat_gbm <- Y1.hat_rf <- Y0.hat_glm <- Y0.hat_dl <- Y0.hat_gbm <- Y0.hat_rf <- NULL
  
  if (int == F) {
  
  if ("glm" %in% library) {
    fitOR_glm = h2o.glm(x=2:ncol(h2oAll),y=1,training_frame=h2oAll,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oAll) -ncol(X_demeaned_lvl1_num)))
    Y.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll)$pred)$predict)
    Y1.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll_Z1)$pred)$predict)
    Y0.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll_Z0)$pred)$predict)
  }
  
  } else {
    
    X_lvl1_num <- as.matrix(as.matrix(X_lvl1)[, apply(as.matrix(X_lvl1), 2, function(x) sum(!unique(x) %in% c(0, 1))) > 0])
    
    X_lvl1_ID_int <- model.matrix( ~ X_lvl1_num:IDfactor-1)
    X_lvl1_ID_int_dm <- apply(X_lvl1_ID_int, 2, function(y) y - ave(y, IDnum))
    
    h2oAll1 = as.h2o(data.frame(Y_demeaned,Z_demeaned,X_demeaned_lvl1, X_lvl1_mean, W_lvl2, X_demeaned_lvl1_num^2, X_lvl1_ID_int_dm))
    h2oAll1_Z1 = h2oAll1_Z0 = h2oAll1
    h2oAll1_Z1$Z_demeaned = as.h2o(rep(1,length(Y_demeaned)) - ave(Z,IDnum)); h2oAll1_Z0$Z_demeaned = as.h2o(rep(0,length(Y_demeaned)) - ave(Z,IDnum))
    
    if ("glm" %in% library) {
      fitOR_glm = h2o.glm(x=2:ncol(h2oAll1),y=1,training_frame=h2oAll1,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oAll1) -ncol(X_demeaned_lvl1_num) - ncol(X_lvl1_ID_int_dm)))
      Y.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll1)$pred)$predict)
      Y1.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll1_Z1)$pred)$predict)
      Y0.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAll1_Z0)$pred)$predict)
    }
  }
  
  if ("deeplearning" %in% library) {
    fitOR_dl = h2o.deeplearning(x=2:ncol(h2oAll),y=1,training_frame=h2oAll,hidden= rep((ncol(h2oAll) -2) + nlevels(IDfactor)+150,2),epochs=500,adaptive_rate=FALSE,
                                l1=0, l2=1e-8, rate=rate, rate_annealing=1e-6, max_w2=10, ...)
    Y.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAll)$pred)$predict)
    Y1.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAll_Z1)$pred)$predict)
    Y0.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAll_Z0)$pred)$predict)
  }
  
  if ("gbm" %in% library) {
    fitOR_gbm = h2o.gbm(x=2:ncol(h2oAll),y=1,training_frame=h2oAll, ntrees=8,max_depth=5,min_rows=min(5,min(table(IDfactor))/2))
    Y.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAll)$pred)$predict)
    Y1.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAll_Z1)$pred)$predict)
    Y0.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAll_Z0)$pred)$predict)
  }
  
  if ("randomForests" %in% library) {
    fitOR_rf = h2o.randomForest(x=2:ncol(h2oAll),y=1,training_frame=h2oAll,ntrees = 500, max_depth = 50, stopping_rounds = 2, stopping_tolerance = 1e-2, score_each_iteration = T)
    Y.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAll)$pred)$predict)
    Y1.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAll_Z1)$pred)$predict)
    Y0.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAll_Z0)$pred)$predict)  
  }
  
  # all predictions
  Y.hats <- cbind(Y.hat_glm, Y.hat_dl, Y.hat_gbm, Y.hat_rf)
  Y1.hats <- cbind(Y1.hat_glm, Y1.hat_dl, Y1.hat_gbm, Y1.hat_rf)
  Y0.hats <- cbind(Y0.hat_glm, Y0.hat_dl, Y0.hat_gbm, Y0.hat_rf)
  
  # nnls
  fit.nnlsOR <- nnls(A=Y.hats, b=Y_demeaned)
  initCoefOR <- coef(fit.nnlsOR)
  initCoefOR[is.na(initCoefOR)] <- 0.0
  
  # normalize so sum(coef) = 1 if possible
  if (sum(initCoefOR) > 0) {
    coefOR <- initCoefOR / sum(initCoefOR)
  } else {
    warning("All algorithms have zero weight", call. = FALSE)
    coefOR <- initCoefOR
  }
  
  Y1.comb <- as.numeric(crossprod(t(Y1.hats[, coefOR != 0, drop = FALSE]), coefOR[coefOR != 0]))
  Y0.comb <- as.numeric(crossprod(t(Y0.hats[, coefOR != 0, drop = FALSE]), coefOR[coefOR != 0]))
  
  return(list(coef=coefOR, Y1.hat = Y1.comb, Y0.hat = Y0.comb,
              Y1.hats = data.frame(Y1.hats), Y0.hats = data.frame(Y0.hats)))
}



DD <- function(Y, Z, interZ=formula(~1),ID, Z.hat, Y0.hat, data) { 
  
  Ygrpcent = Y - ave(Y, ID)
  Zgrpcent = Z - ave(Z, ID)
  
  Y.resid.D = Ygrpcent - Y0.hat
  Z.resid.D = Zgrpcent - Z.hat
  
  Y.resid.DD <- Y.resid.D - ave(Y.resid.D, ID)
  Z.resid.DD <- Z.resid.D - ave(Z.resid.D, ID)
  
  Zmat <- model.matrix(interZ, data=data)
  Zmat.rr <- apply(Zmat, 2, '*', Z.resid.DD)
  Zmat.Z <- apply(Zmat, 2, '*', Z)
  Zmat.Z.grpcent <- apply(Zmat.Z, 2, function(x) x -ave(x, ID))
  
  
  ivm = summary(ivreg(Y.resid.DD ~ Zmat.Z.grpcent -1 | Zmat.rr - 1))
  tau.est = matrix(0,ncol(Zmat),2)
  tau.est[,1] = ivm$coefficients[,1]
  tau.est[,2] = ivm$coefficients[,2]
  
  rownames(tau.est) = colnames(Zmat)
  colnames(tau.est) = c("Estimate","Std. Error")
  return(tau.est)
}


# :::::::::::::::::::::::::::::::
# ::::: Additional functions :::: ####
# :::::::::::::::::::::::::::::::

overlap <- function(X, Z, bin = 20, ...)
{
  # plot a histogram of a covariate by group
  # x   ... numeric vector (covariate)
  # z   ... treatment indicator (dummy)
  # lab ... label for title and x-axis
  # bin ... number of bins for histogram
  
  r1 <- range(X)
  if (!is.numeric(Z)) Z <- as.numeric(Z) - 1
  c.dat <- hist(X[Z == 0], seq(r1[1], r1[2], length = bin), plot = F)  # histogram data for control group
  t.dat <- hist(X[Z == 1], seq(r1[1], r1[2], length = bin), plot = F)  # histogram data for treatm. group
  c.dat$counts <- -c.dat$counts
  plot(t.dat, axes = F, ylim = c(min(c.dat$counts), max(t.dat$counts)), ...)
  plot(c.dat, add = T, density = 30)
  axis(1)
  ax.loc <- axis(2, labels = F)
  axis(2, at = ax.loc, labels = abs(ax.loc))
  y <- par('usr')[3:4]
  text(rep(max(X), 2), c(y[2] - diff(y)*.05, y[1] + diff(y)*.05),
       c('Treatment', 'Control'), adj = 1, xpd = T)
}


