###############################################################################################################################
# Robust Machine Learning for Treatment Effects in Multilevel Observational Studies Under Cluster-level Unmeasured Confounding
# Youmi Suk and Hyunseung Kang 

# :::: Simulation codes ::::

# :: load packages and DGP codes
library(nnls)
library(AER)
library(h2o)
library(bit64)
library(drgee)
library(bartCause)
library(tmle)

source("DGP_twolevel_crossclassified.R") # DGP code
source("CURobustML_sepFun.R") # functions for the propensity score and outcome predictions based on our proposed methods

localH2O = h2o.init()
h2o.removeAll()

# :::: Designs 1 and 3 :::: ####
# Design 1: Two-level Data, No Cross-level Interaction, and Large Cluster Size
# Design 3: Two-level Data, No Cross-level Interaction, and Small Cluster Size

iter <- 500 # the number of replications (we use HTCondor to run simulations)

MLtestCATE_beta1.rlst <- MLtestCATE_beta2.rlst <- MLtestATE.est.rlst <- matrix(NA, nrow=iter, ncol=12)
DRtestCATE_beta1.rlst <- DRtestCATE_beta2.rlst <- DRtestATE.est.rlst <- matrix(NA, nrow=iter, ncol=16)

interZ <- formula(~ X3)
interZ_c <- formula(~ X3_c)

for (i in 1:iter) {

  dat <- twolevel.pop(Smpl.size = "50.100", crosslevel.int = F) # Design 1's the sample size (50, 100)
  # You can change the sample size conditions for Design 3: (100, 5), (100, 10), (100, 20), and (100, 30)
  
  dat$X3_c <- dat$X3 - mean(dat$X3)
  dat$term1 <- dat$X1 * dat$X2^2  
  
  ID=dat$id ###

  # ::: ML Methods :::
  # :: our proposed methods
  er_PR <- erLearnPR(Z=dat$Z, X=dat[, c("X1", "X2", "X3", "W1")], ID=ID, rate=0.01)
  or_PR <- orLearnPR(Y=dat$Y, Z=dat$Z, X=dat[, c("X1", "X2", "X3", "W1")], ID=ID, rate=0.01)
  er_DD <- erLearnDD(Z=dat$Z, X=dat[, c("X1", "X2", "X3", "W1")], ID=ID, rate=0.01)
  or_DD <- orLearnDD(Y=dat$Y, Z=dat$Z, X=dat[, c("X1", "X2", "X3", "W1")], ID=ID, rate=0.01)

  # : summarize CATE
  # - glm
  PR.CATEest_glm <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_glm, interZ =formula(~X3), Y1.hat=or_PR$Y1.hats$Y1.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)
  DD.CATEest_glm <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_DD$Z.hats$Z.hat_glm,Y0.hat=or_DD$Y0.hats$Y0.hat_glm, data=dat)
  DDPR.CATEest_glm <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_PR$Z.hats$Z.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)
  
  # - dl
  PR.CATEest_dl <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_dl, interZ =formula(~X3), Y1.hat=or_PR$Y1.hats$Y1.hat_dl, Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)
  DD.CATEest_dl <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_DD$Z.hats$Z.hat_dl,Y0.hat=or_DD$Y0.hats$Y0.hat_dl, data=dat)
  DDPR.CATEest_dl <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_PR$Z.hats$Z.hat_dl,Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)
  
  # - el
  PR.CATEest_el <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hat, interZ =formula(~X3), Y1.hat=or_PR$Y1.hat, Y0.hat=or_PR$Y0.hat, data=dat)
  DD.CATEest_el <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_DD$Z.hat,Y0.hat=or_DD$Y0.hat, data=dat)
  DDPR.CATEest_el <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_PR$Z.hat,Y0.hat=or_PR$Y0.hat, data=dat)
  
  # : summarize ATE
  PR.ATEest_glm <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_glm, interZ =formula(~1), Y1.hat=or_PR$Y1.hats$Y1.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)
  DD.ATEest_glm <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_DD$Z.hats$Z.hat_glm,Y0.hat=or_DD$Y0.hats$Y0.hat_glm, data=dat)[1,]
  DDPR.ATEest_glm <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_PR$Z.hats$Z.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)[1,]

  # - dl
  PR.ATEest_dl <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_dl, interZ =formula(~X3_c), Y1.hat=or_PR$Y1.hats$Y1.hat_dl, Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)[1,]
  DD.ATEest_dl <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_DD$Z.hats$Z.hat_dl,Y0.hat=or_DD$Y0.hats$Y0.hat_dl, data=dat)[1,]
  DDPR.ATEest_dl <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_PR$Z.hats$Z.hat_dl,Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)[1,]
  
  # - el
  PR.ATEest_el <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hat, interZ =formula(~X3_c), Y1.hat=or_PR$Y1.hat, Y0.hat=or_PR$Y0.hat, data=dat)[1,]
  DD.ATEest_el <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_DD$Z.hat,Y0.hat=or_DD$Y0.hat, data=dat)[1,]
  DDPR.ATEest_el <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_PR$Z.hat,Y0.hat=or_PR$Y0.hat, data=dat)[1,]
  
  # :: other comparisons: DRGEE, BART, and TMLE
  # : DRCGEE
  strong.man <- drgee(oformula=formula(Y~ X1 + X2 + X3 + I(X3^2)), eformula=(Z~ X1 + X2 + X3 + I(X3^2)), 
                      iaformula = interZ, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  
  strong.manATE <- drgee(oformula=formula(Y~ X1 + X2 + X3  + I(X3^2)), eformula=(Z~ X1 + X2 + X3  + I(X3^2)), 
                         iaformula = interZ_c, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  
  # : BART
  bart.fit <- bartc(response=dat$Y, treatment=dat$Z, confounders=dat[, c("X1", "X2", "X3", "W1")], method.rsp="bart", method.trt="bart", p.scoreAsCovariate=TRUE, use.rbrt=FALSE, keepTrees=TRUE)
  est.bart.ite <- fitted(bart.fit, type="icate", sample="all") 
  BART_est <- summary(lm(est.bart.ite ~ dat$X3))$coef[, 1]
  
  # : TLME
  tmle.fit <- tmle(Y=dat$Y, A=dat$Z, W=dat[, c("X1", "X2", "X3", "W1")])
  est.tmle.ite <- tmle.fit$Qstar[,2]-tmle.fit$Qstar[,1]
  TMLE_est <- summary(lm(est.tmle.ite ~ dat$X3))$coef[, 1]

  # save CATE
  CATE.rlst <- rbind(PR.CATEest_glm[,1], DD.CATEest_glm[,1], DDPR.CATEest_glm[,1], 
                     PR.CATEest_dl[,1], DD.CATEest_dl[,1], DDPR.CATEest_dl[,1],
                     PR.CATEest_el[,1], DD.CATEest_el[,1], DDPR.CATEest_el[,1],
                     summary(strong.man)$coefficients[, 1], BART_est, TMLE_est)

  MLtestCATE_beta1.rlst[i,] <- CATE.rlst[,1]
  MLtestCATE_beta2.rlst[i,] <- CATE.rlst[,2]

  # save ATE
  MLtestATE.est.rlst[i,] <- as.numeric(c(PR.ATEest_glm[1], DD.ATEest_glm[1], DDPR.ATEest_glm[1],
                                          PR.ATEest_dl[1], DD.ATEest_dl[1], DDPR.ATEest_dl[1],
                                          PR.ATEest_el[1], DD.ATEest_el[1], DDPR.ATEest_el[1],
                                          summary(strong.manATE)$coefficients[1, 1], mean(est.bart.ite), mean(est.tmle.ite)))

  # ::: glm-based methods under different specifications :::
  # -- PR.Estimator
  fit.Z.PR <- glm(Z ~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) + id, family=binomial(), data=dat)
  fit.Z.PR.mis <- glm(Z ~ X1 + X2 + X3 + id, family=binomial(), data=dat)
  fit.Y.PR <- lm(Y ~ Z + Z:X3 + X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) + id, data=dat)
  fit.Y.PR.mis <- lm(Y ~ Z + Z:X3 + X1 + X2 + X3 + id, data=dat)
  
  # -- DD Estimator
  dat$Y_dm <- dat$Y - ave(dat$Y, ID); dat$Z_dm <- dat$Z - ave(dat$Z, ID)
  dat$ZX3_dm <- dat$Z*dat$X3 - ave(dat$Z*dat$X3, ID)
  dat$X1_dm <- dat$X1 - ave(dat$X1, ID); dat$X2_dm <- dat$X2 - ave(dat$X2, ID);  dat$X3_dm <- dat$X3 - ave(dat$X3, ID) 
  dat$X3sq_dm <- dat$X3^2 - ave(dat$X3^2, ID) 
  dat$X2ind_dm <- dat$X2*I( 0<dat$X2 & dat$X2 < 1.5 ) - ave(dat$X2*I( 0<dat$X2 & dat$X2 < 1.5 ), ID) 
  dat$X1X3ind_dm <- dat$X1*I( dat$X3 < 0.2 ) - ave(dat$X1*I( dat$X3 < 0.2 ), ID) 
  dat$term1_dm <- dat$term1 - ave(dat$term1, ID)
  
  fit.Z.DD <- lm(Z_dm ~ X1_dm + X2_dm + X3_dm + term1_dm + X2ind_dm + X3sq_dm + X1X3ind_dm , data=dat)
  fit.Z.DD.mis <- glm(Z_dm ~ X1_dm + X2_dm + X3_dm, data=dat)
  fit.Y.DD <- lm(Y_dm ~ Z_dm + ZX3_dm + X1_dm + X2_dm + X3_dm + term1_dm + X2ind_dm + X3sq_dm + X1X3ind_dm , data=dat)
  fit.Y.DD.mis <- lm(Y_dm ~ Z_dm + ZX3_dm + X1_dm + X2_dm + X3_dm, data=dat)
  
  # -- drgee estimator
  # - CATE
  both.cor <- drgee(oformula=formula(Y~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), eformula=(Z~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), 
                    iaformula = interZ, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  mis.out <- drgee(oformula=formula(Y~ X1 + X2 + X3), eformula=(Z~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), 
                   iaformula = interZ, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  mis.sel <- drgee(oformula=formula(Y~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2))), eformula=(Z~ X1 + X2 + X3), 
                   iaformula = interZ, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  both.mis <- drgee(oformula=formula(Y~ X1 + X2 + X3), eformula=(Z~ X1 + X2 + X3), iaformula = interZ, olink="identity", elink="identity", 
                    data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  
  # - ATE
  both.corATE <- drgee(oformula=formula(Y~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), eformula=(Z~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), 
                       iaformula = interZ_c, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  mis.outATE <- drgee(oformula=formula(Y~ X1 + X2 + X3), eformula=(Z~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), 
                      iaformula = interZ_c, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  mis.selATE <- drgee(oformula=formula(Y~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2))), eformula=(Z~ X1 + X2 + X3), 
                      iaformula = interZ_c, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  both.misATE <- drgee(oformula=formula(Y~ X1 + X2 + X3), eformula=(Z~ X1 + X2 + X3), iaformula = interZ_c, olink="identity", elink="identity", 
                       data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  
  Z.hat.PR <- predict(fit.Z.PR, type="response")
  Z.hat.PR.mis <- predict(fit.Z.PR.mis, type="response")
  
  dat_Z1 <- dat_Z0 <- dat
  dat_Z1$Z <- 1; dat_Z0$Z <- 0
  
  Y0.hat.PR <- predict(fit.Y.PR, dat_Z0)
  Y0.hat.PR.mis <- predict(fit.Y.PR.mis, dat_Z0)
  
  Y1.hat.PR <- predict(fit.Y.PR, dat_Z1)
  Y1.hat.PR.mis <- predict(fit.Y.PR.mis, dat_Z1)
  
  Z.hat.DD <- predict(fit.Z.DD)
  Z.hat.DD.mis <- predict(fit.Z.DD.mis)
  
  dat_Z0 <- dat
  dat_Z0$Z_dm <- 0
  dat_Z0$ZX3_dm <- 0
  
  Y0.hat.DD <- predict(fit.Y.DD, dat_Z0)
  Y0.hat.DD.mis <- predict(fit.Y.DD.mis, dat_Z0)
  
  # summarize CATE
  DRtestCATE.est <- rbind(DR(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, Z.hat=Z.hat.PR, Y1.hat=Y1.hat.PR, Y0.hat=Y0.hat.PR)[,1], # both.cor - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.DD, Y0.hat=Y0.hat.DD)[,1],  # both.cor - DD
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR)[,1],  # both.cor - DDPR
                          summary(both.cor)$coefficients[, 1],  # both.cor - DRCGEE
                          
                          DR(Y=dat$Y,Z=dat$Z, interZ =interZ, data=dat, Z.hat=Z.hat.PR,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[,1], # out.mis - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ =interZ, data=dat, ID=ID, Z.hat=Z.hat.DD, Y0.hat=Y0.hat.DD.mis)[,1], # out.mis - DD
                          DD(Y=dat$Y,Z=dat$Z, interZ =interZ, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR.mis)[,1], # out.mis - DDPR
                          summary(mis.out)$coefficients[, 1], # out.mis - DRCGEE
                          
                          DR(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR,Y0.hat=Y0.hat.PR)[,1], # trt.mis - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.DD.mis, Y0.hat=Y0.hat.DD)[,1], # trt.mis - DD
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR)[,1], # trt.mis - DDPR
                          summary(mis.sel)$coefficients[, 1], # trt.mis - DRCGEE
                          
                          DR(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[,1], # both.mis - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.DD.mis, Y0.hat=Y0.hat.DD.mis)[,1], # both.mis - DD
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR.mis)[,1],  # both.mis - DDPR
                          summary(both.mis)$coefficients[, 1]) # both.mis - DRCGEE
  
  DRtestCATE_beta1.rlst[i, ] <- DRtestCATE.est[, 1]
  DRtestCATE_beta2.rlst[i, ] <- DRtestCATE.est[, 2]
  
  # summarize ATE
  DRtestATE.est.rlst[i,] <- c(DR(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, Z.hat=Z.hat.PR, Y1.hat=Y1.hat.PR, Y0.hat=Y0.hat.PR)[1,1], # both.cor - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.DD, Y0.hat=Y0.hat.DD)[1,1], # both.cor - DD
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR)[1,1], # both.cor - DDPR
                              summary(both.corATE)$coefficients[1, 1], # both.cor - DRCGEE
                              
                              DR(Y=dat$Y,Z=dat$Z, interZ =interZ_c, data=dat, Z.hat=Z.hat.PR,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[1,1],  # out.mis - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ =interZ_c, data=dat, ID=ID, Z.hat=Z.hat.DD, Y0.hat=Y0.hat.DD.mis)[1,1], # out.mis - DD
                              DD(Y=dat$Y,Z=dat$Z, interZ =interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR.mis)[1,1], # out.mis - DDPR
                              summary(mis.outATE)$coefficients[1, 1], # out.mis - DRCGEE
                              
                              DR(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR,Y0.hat=Y0.hat.PR)[1,1], # out.mis - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.DD.mis, Y0.hat=Y0.hat.DD)[1,1], # out.mis - DD
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR)[1,1], # trt.mis - DDPR
                              summary(mis.selATE)$coefficients[1, 1],  # trt.mis - DRCGEE
                              
                              DR(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[1,1], # both.mis - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.DD.mis, Y0.hat=Y0.hat.DD.mis)[1,1], # both.mis - DD
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR.mis)[1,1],  # both.mis - DDPR
                              summary(both.misATE)$coefficients[1, 1])  # both.mis - DRCGEE
  
}

MLtest_rlst <- cbind(MLtestCATE_beta1.rlst, MLtestCATE_beta2.rlst, MLtestATE.est.rlst) # 12*3
DRtest_rlst <- cbind(DRtestCATE_beta1.rlst, DRtestCATE_beta2.rlst, DRtestATE.est.rlst) # 16*3




# :::: Design 2: Two-level Data, Cross-level Interaction, and Large Cluster Size :::: ####

iter <- 500 # the number of replications (we use HTCondor to run simulations)

MLtestCATE_beta1.rlst <- MLtestCATE_beta2.rlst <- MLtestATE.est.rlst <- matrix(NA, nrow=iter, ncol=9)
DRtestCATE_beta1.rlst <- DRtestCATE_beta2.rlst <- DRtestATE.est.rlst <- matrix(NA, nrow=iter, ncol=8)

interZ <- formula(~ X3)
interZ_c <- formula(~ X3_c)

for (i in 1:iter) {
  
  dat <- twolevel.pop(Smpl.size = "50.100", crosslevel.int = T) 
  dat$X3_c <- dat$X3 - mean(dat$X3)
  dat$term1 <- dat$X1 * dat$X2^2  
  
  ID=dat$id ###
  
  # ::: ML Methods :::
  # :: our proposed methods
  er_PR <- erLearnPR(Z=dat$Z, X=dat[, c("X1", "X2", "X3", "W1")], ID=ID, rate=0.01)
  or_PR <- orLearnPR(Y=dat$Y, Z=dat$Z, X=dat[, c("X1", "X2", "X3", "W1")], ID=ID, rate=0.01)
  
  # : summarize CATE
  # - glm
  PR.CATEest_glm <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_glm, interZ =formula(~X3), Y1.hat=or_PR$Y1.hats$Y1.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)
  DDPR.CATEest_glm <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_PR$Z.hats$Z.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)
  
  # - dl
  PR.CATEest_dl <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_dl, interZ =formula(~X3), Y1.hat=or_PR$Y1.hats$Y1.hat_dl, Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)
  DDPR.CATEest_dl <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_PR$Z.hats$Z.hat_dl,Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)
  
  # - el
  PR.CATEest_el <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hat, interZ =formula(~X3), Y1.hat=or_PR$Y1.hat, Y0.hat=or_PR$Y0.hat, data=dat)
  DDPR.CATEest_el <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_PR$Z.hat,Y0.hat=or_PR$Y0.hat, data=dat)
  
  # : summarize ATE
  PR.ATEest_glm <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_glm, interZ =formula(~1), Y1.hat=or_PR$Y1.hats$Y1.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)
  DDPR.ATEest_glm <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_PR$Z.hats$Z.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)[1,]
  
  # - dl
  PR.ATEest_dl <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_dl, interZ =formula(~X3_c), Y1.hat=or_PR$Y1.hats$Y1.hat_dl, Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)[1,]
  DDPR.ATEest_dl <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_PR$Z.hats$Z.hat_dl,Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)[1,]
  
  # - el
  PR.ATEest_el <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hat, interZ =formula(~X3_c), Y1.hat=or_PR$Y1.hat, Y0.hat=or_PR$Y0.hat, data=dat)[1,]
  DDPR.ATEest_el <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_PR$Z.hat,Y0.hat=or_PR$Y0.hat, data=dat)[1,]
  
  # :: other comparisons: DRGEE, BART, and TMLE
  # : DRCGEE
  strong.man <- drgee(oformula=formula(Y~ X1 + X2 + X3 + I(X3^2)), eformula=(Z~ X1 + X2 + X3 + I(X3^2)), 
                      iaformula = interZ, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  
  strong.manATE <- drgee(oformula=formula(Y~ X1 + X2 + X3  + I(X3^2)), eformula=(Z~ X1 + X2 + X3  + I(X3^2)), 
                         iaformula = interZ_c, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  
  # : BART
  bart.fit <- bartc(response=dat$Y, treatment=dat$Z, confounders=dat[, c("X1", "X2", "X3", "W1")], method.rsp="bart", method.trt="bart", p.scoreAsCovariate=TRUE, use.rbrt=FALSE, keepTrees=TRUE)
  est.bart.ite <- fitted(bart.fit, type="icate", sample="all") 
  BART_est <- summary(lm(est.bart.ite ~ dat$X3))$coef[, 1]
  
  # : TLME
  tmle.fit <- tmle(Y=dat$Y, A=dat$Z, W=dat[, c("X1", "X2", "X3", "W1")])
  est.tmle.ite <- tmle.fit$Qstar[,2]-tmle.fit$Qstar[,1]
  TMLE_est <- summary(lm(est.tmle.ite ~ dat$X3))$coef[, 1]
  
  # save CATE
  CATE.rlst <- rbind(PR.CATEest_glm[,1], DDPR.CATEest_glm[,1], 
                     PR.CATEest_dl[,1], DDPR.CATEest_dl[,1],
                     PR.CATEest_el[,1], DDPR.CATEest_el[,1],
                     summary(strong.man)$coefficients[, 1], BART_est, TMLE_est)

  MLtestCATE_beta1.rlst[i,] <- CATE.rlst[,1]
  MLtestCATE_beta2.rlst[i,] <- CATE.rlst[,2]
  
  # save ATE
  MLtestATE.est.rlst[i,] <- as.numeric(c(PR.ATEest_glm[1], DDPR.ATEest_glm[1],
                                         PR.ATEest_dl[1], DDPR.ATEest_dl[1],
                                         PR.ATEest_el[1], DDPR.ATEest_el[1],
                                         summary(strong.manATE)$coefficients[1, 1], mean(est.bart.ite), mean(est.tmle.ite)))
  
  # ::: glm-based methods under different specifications :::
  # -- PR.Estimator
  
  intIDX3 <- model.matrix( ~ ID:X3 - 1, dat)[,-1]
  
  fit.Z.PR <- glm(Z ~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) + id + intIDX3, family=binomial(), data=dat)
  fit.Z.PR.mis <- glm(Z ~ X1 + X2 + X3 + id + intIDX3, family=binomial(), data=dat)
  fit.Y.PR <- lm(Y ~ Z + Z:X3 + X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) + id + intIDX3, data=dat)
  fit.Y.PR.mis <- lm(Y ~ Z + Z:X3 + X1 + X2 + X3 + id + intIDX3, data=dat)
  
  Z.hat.PR <- predict(fit.Z.PR, type="response")
  Z.hat.PR.mis <- predict(fit.Z.PR.mis, type="response")
  
  dat_Z1 <- dat_Z0 <- dat
  dat_Z1$Z <- 1; dat_Z0$Z <- 0
  
  Y0.hat.PR <- predict(fit.Y.PR, dat_Z0)
  Y0.hat.PR.mis <- predict(fit.Y.PR.mis, dat_Z0)
  
  Y1.hat.PR <- predict(fit.Y.PR, dat_Z1)
  Y1.hat.PR.mis <- predict(fit.Y.PR.mis, dat_Z1)
  
  # summarize CATE
  DRtestCATE.est <- rbind(DR(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, Z.hat=Z.hat.PR, Y1.hat=Y1.hat.PR, Y0.hat=Y0.hat.PR)[,1], # both.cor - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR)[,1],  # both.cor - DDPR
                          
                          DR(Y=dat$Y,Z=dat$Z, interZ =interZ, data=dat, Z.hat=Z.hat.PR,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[,1], # out.mis - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ =interZ, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR.mis)[,1], # out.mis - DDPR

                          DR(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR,Y0.hat=Y0.hat.PR)[,1], # trt.mis - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR)[,1], # trt.mis - DDPR

                          DR(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[,1], # both.mis - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR.mis)[,1])  # both.mis - DDPR

  DRtestCATE_beta1.rlst[i, ] <- DRtestCATE.est[, 1]
  DRtestCATE_beta2.rlst[i, ] <- DRtestCATE.est[, 2]
  
  # summarize ATE
  DRtestATE.est.rlst[i,] <- c(DR(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, Z.hat=Z.hat.PR, Y1.hat=Y1.hat.PR, Y0.hat=Y0.hat.PR)[1,1], # both.cor - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR)[1,1], # both.cor - DDPR

                              DR(Y=dat$Y,Z=dat$Z, interZ =interZ_c, data=dat, Z.hat=Z.hat.PR,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[1,1],  # out.mis - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ =interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR.mis)[1,1], # out.mis - DDPR

                              DR(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR,Y0.hat=Y0.hat.PR)[1,1], # out.mis - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR)[1,1], # trt.mis - DDPR

                              DR(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[1,1], # both.mis - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR.mis)[1,1])  # both.mis - DDPR
}

MLtest_rlst <- cbind(MLtestCATE_beta1.rlst, MLtestCATE_beta2.rlst, MLtestATE.est.rlst) # 9*3
DRtest_rlst <- cbind(DRtestCATE_beta1.rlst, DRtestCATE_beta2.rlst, DRtestATE.est.rlst) # 8*3


# :::: Design 4: Cross-classified Data, No Cross-level Interaction, and Large Cluster Size :::: ####

iter <- 500 # the number of replications (we use HTCondor to run simulations)

MLtestCATE_beta1.rlst <- MLtestCATE_beta2.rlst <- MLtestATE.est.rlst <- matrix(NA, nrow=iter, ncol=12)
DRtestCATE_beta1.rlst <- DRtestCATE_beta2.rlst <- DRtestATE.est.rlst <- matrix(NA, nrow=iter, ncol=16)

interZ <- formula(~ X3)
interZ_c <- formula(~ X3_c)

for (i in 1:iter) {
  
  dat <- ccrem.pop(Smpl.size = "25.25.300")
  
  dat$X3_c <- dat$X3 - mean(dat$X3)
  dat$term1 <- dat$X1 * dat$X2^2  
  
  ID=dat$f12id ###
  
  # ::: ML Methods :::
  # :: our proposed methods
  er_PR <- erLearnPR(Z=dat$Z, X=dat[, c("X1", "X2", "X3", "W1", "Q1")], ID=ID, rate=0.01)
  or_PR <- orLearnPR(Y=dat$Y, Z=dat$Z, X=dat[, c("X1", "X2", "X3", "W1", "Q1")], ID=ID, rate=0.01)
  er_DD <- erLearnDD(Z=dat$Z, X=dat[, c("X1", "X2", "X3", "W1", "Q1")], ID=ID, rate=0.01)
  or_DD <- orLearnDD(Y=dat$Y, Z=dat$Z, X=dat[, c("X1", "X2", "X3", "W1", "Q1")], ID=ID, rate=0.01)
  
  # : summarize CATE
  # - glm
  PR.CATEest_glm <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_glm, interZ =formula(~X3), Y1.hat=or_PR$Y1.hats$Y1.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)
  DD.CATEest_glm <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_DD$Z.hats$Z.hat_glm,Y0.hat=or_DD$Y0.hats$Y0.hat_glm, data=dat)
  DDPR.CATEest_glm <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_PR$Z.hats$Z.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)
  
  # - dl
  PR.CATEest_dl <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_dl, interZ =formula(~X3), Y1.hat=or_PR$Y1.hats$Y1.hat_dl, Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)
  DD.CATEest_dl <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_DD$Z.hats$Z.hat_dl,Y0.hat=or_DD$Y0.hats$Y0.hat_dl, data=dat)
  DDPR.CATEest_dl <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_PR$Z.hats$Z.hat_dl,Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)
  
  # - el
  PR.CATEest_el <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hat, interZ =formula(~X3), Y1.hat=or_PR$Y1.hat, Y0.hat=or_PR$Y0.hat, data=dat)
  DD.CATEest_el <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_DD$Z.hat,Y0.hat=or_DD$Y0.hat, data=dat)
  DDPR.CATEest_el <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3), Z.hat=er_PR$Z.hat,Y0.hat=or_PR$Y0.hat, data=dat)
  
  # : summarize ATE
  PR.ATEest_glm <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_glm, interZ =formula(~1), Y1.hat=or_PR$Y1.hats$Y1.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)
  DD.ATEest_glm <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_DD$Z.hats$Z.hat_glm,Y0.hat=or_DD$Y0.hats$Y0.hat_glm, data=dat)[1,]
  DDPR.ATEest_glm <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_PR$Z.hats$Z.hat_glm, Y0.hat=or_PR$Y0.hats$Y0.hat_glm, data=dat)[1,]
  
  # - dl
  PR.ATEest_dl <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hats$Z.hat_dl, interZ =formula(~X3_c), Y1.hat=or_PR$Y1.hats$Y1.hat_dl, Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)[1,]
  DD.ATEest_dl <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_DD$Z.hats$Z.hat_dl,Y0.hat=or_DD$Y0.hats$Y0.hat_dl, data=dat)[1,]
  DDPR.ATEest_dl <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_PR$Z.hats$Z.hat_dl,Y0.hat=or_PR$Y0.hats$Y0.hat_dl, data=dat)[1,]
  
  # - el
  PR.ATEest_el <- DR(Y=dat$Y, Z=dat$Z, Z.hat=er_PR$Z.hat, interZ =formula(~X3_c), Y1.hat=or_PR$Y1.hat, Y0.hat=or_PR$Y0.hat, data=dat)[1,]
  DD.ATEest_el <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_DD$Z.hat,Y0.hat=or_DD$Y0.hat, data=dat)[1,]
  DDPR.ATEest_el <- DD(Y=dat$Y, Z=dat$Z, ID=ID, interZ =formula(~X3_c), Z.hat=er_PR$Z.hat,Y0.hat=or_PR$Y0.hat, data=dat)[1,]
  
  # :: other comparisons: DRGEE, BART, and TMLE
  # : DRCGEE
  strong.man <- drgee(oformula=formula(Y~ X1 + X2 + X3 + I(X3^2)), eformula=(Z~ X1 + X2 + X3 + I(X3^2)), 
                      iaformula = interZ, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  
  strong.manATE <- drgee(oformula=formula(Y~ X1 + X2 + X3  + I(X3^2)), eformula=(Z~ X1 + X2 + X3  + I(X3^2)), 
                         iaformula = interZ_c, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  
  # : BART
  bart.fit <- bartc(response=dat$Y, treatment=dat$Z, confounders=dat[, c("X1", "X2", "X3", "W1", "Q1")], method.rsp="bart", method.trt="bart", p.scoreAsCovariate=TRUE, use.rbrt=FALSE, keepTrees=TRUE)
  est.bart.ite <- fitted(bart.fit, type="icate", sample="all") 
  BART_est <- summary(lm(est.bart.ite ~ dat$X3))$coef[, 1]
  
  # : TLME
  tmle.fit <- tmle(Y=dat$Y, A=dat$Z, W=dat[, c("X1", "X2", "X3", "W1", "Q1")])
  est.tmle.ite <- tmle.fit$Qstar[,2]-tmle.fit$Qstar[,1]
  TMLE_est <- summary(lm(est.tmle.ite ~ dat$X3))$coef[, 1]
  
  # save CATE
  CATE.rlst <- rbind(PR.CATEest_glm[,1], DD.CATEest_glm[,1], DDPR.CATEest_glm[,1], 
                     PR.CATEest_dl[,1], DD.CATEest_dl[,1], DDPR.CATEest_dl[,1],
                     PR.CATEest_el[,1], DD.CATEest_el[,1], DDPR.CATEest_el[,1],
                     summary(strong.man)$coefficients[, 1], BART_est, TMLE_est)
  
  MLtestCATE_beta1.rlst[i,] <- CATE.rlst[,1]
  MLtestCATE_beta2.rlst[i,] <- CATE.rlst[,2]
  
  # save ATE
  MLtestATE.est.rlst[i,] <- as.numeric(c(PR.ATEest_glm[1], DD.ATEest_glm[1], DDPR.ATEest_glm[1],
                                         PR.ATEest_dl[1], DD.ATEest_dl[1], DDPR.ATEest_dl[1],
                                         PR.ATEest_el[1], DD.ATEest_el[1], DDPR.ATEest_el[1],
                                         summary(strong.manATE)$coefficients[1, 1], mean(est.bart.ite), mean(est.tmle.ite)))
  
  # ::: glm-based methods under different specifications :::
  # -- PR.Estimator
  fit.Z.PR <- glm(Z ~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) + f12id, family=binomial(), data=dat)
  fit.Z.PR.mis <- glm(Z ~ X1 + X2 + X3 + f12id, family=binomial(), data=dat)
  fit.Y.PR <- lm(Y ~ Z + Z:X3 + X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) + f12id, data=dat)
  fit.Y.PR.mis <- lm(Y ~ Z + Z:X3 + X1 + X2 + X3 + f12id, data=dat)
  
  # -- DD Estimator
  dat$Y_dm <- dat$Y - ave(dat$Y, ID); dat$Z_dm <- dat$Z - ave(dat$Z, ID)
  dat$ZX3_dm <- dat$Z*dat$X3 - ave(dat$Z*dat$X3, ID)
  dat$X1_dm <- dat$X1 - ave(dat$X1, ID); dat$X2_dm <- dat$X2 - ave(dat$X2, ID);  dat$X3_dm <- dat$X3 - ave(dat$X3, ID) 
  dat$X3sq_dm <- dat$X3^2 - ave(dat$X3^2, ID) 
  dat$X2ind_dm <- dat$X2*I( 0<dat$X2 & dat$X2 < 1.5 ) - ave(dat$X2*I( 0<dat$X2 & dat$X2 < 1.5 ), ID) 
  dat$X1X3ind_dm <- dat$X1*I( dat$X3 < 0.2 ) - ave(dat$X1*I( dat$X3 < 0.2 ), ID) 
  dat$term1_dm <- dat$term1 - ave(dat$term1, ID)
  
  fit.Z.DD <- lm(Z_dm ~ X1_dm + X2_dm + X3_dm + term1_dm + X2ind_dm + X3sq_dm + X1X3ind_dm , data=dat)
  fit.Z.DD.mis <- glm(Z_dm ~ X1_dm + X2_dm + X3_dm, data=dat)
  fit.Y.DD <- lm(Y_dm ~ Z_dm + ZX3_dm + X1_dm + X2_dm + X3_dm + term1_dm + X2ind_dm + X3sq_dm + X1X3ind_dm , data=dat)
  fit.Y.DD.mis <- lm(Y_dm ~ Z_dm + ZX3_dm + X1_dm + X2_dm + X3_dm, data=dat)
  
  # -- drgee estimator
  # - CATE
  both.cor <- drgee(oformula=formula(Y~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), eformula=(Z~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), 
                    iaformula = interZ, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  mis.out <- drgee(oformula=formula(Y~ X1 + X2 + X3), eformula=(Z~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), 
                   iaformula = interZ, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  mis.sel <- drgee(oformula=formula(Y~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2))), eformula=(Z~ X1 + X2 + X3), 
                   iaformula = interZ, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  both.mis <- drgee(oformula=formula(Y~ X1 + X2 + X3), eformula=(Z~ X1 + X2 + X3), iaformula = interZ, olink="identity", elink="identity", 
                    data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  
  # - ATE
  both.corATE <- drgee(oformula=formula(Y~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), eformula=(Z~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), 
                       iaformula = interZ_c, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  mis.outATE <- drgee(oformula=formula(Y~ X1 + X2 + X3), eformula=(Z~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2)) ), 
                      iaformula = interZ_c, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  mis.selATE <- drgee(oformula=formula(Y~ X1 + X2 + X3 + term1 + I(X2*( 0<X2 & X2 < 1.5 )) + I(X3^2) + I(X1*(X3 < 0.2))), eformula=(Z~ X1 + X2 + X3), 
                      iaformula = interZ_c, olink="identity", elink="identity", data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  both.misATE <- drgee(oformula=formula(Y~ X1 + X2 + X3), eformula=(Z~ X1 + X2 + X3), iaformula = interZ_c, olink="identity", elink="identity", 
                       data=dat, estimation.method= "dr", clusterid = ID, cond=T)
  
  Z.hat.PR <- predict(fit.Z.PR, type="response")
  Z.hat.PR.mis <- predict(fit.Z.PR.mis, type="response")
  
  dat_Z1 <- dat_Z0 <- dat
  dat_Z1$Z <- 1; dat_Z0$Z <- 0
  
  Y0.hat.PR <- predict(fit.Y.PR, dat_Z0)
  Y0.hat.PR.mis <- predict(fit.Y.PR.mis, dat_Z0)
  
  Y1.hat.PR <- predict(fit.Y.PR, dat_Z1)
  Y1.hat.PR.mis <- predict(fit.Y.PR.mis, dat_Z1)
  
  Z.hat.DD <- predict(fit.Z.DD)
  Z.hat.DD.mis <- predict(fit.Z.DD.mis)
  
  dat_Z0 <- dat
  dat_Z0$Z_dm <- 0
  dat_Z0$ZX3_dm <- 0
  
  Y0.hat.DD <- predict(fit.Y.DD, dat_Z0)
  Y0.hat.DD.mis <- predict(fit.Y.DD.mis, dat_Z0)
  
  # summarize CATE
  DRtestCATE.est <- rbind(DR(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, Z.hat=Z.hat.PR, Y1.hat=Y1.hat.PR, Y0.hat=Y0.hat.PR)[,1], # both.cor - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.DD, Y0.hat=Y0.hat.DD)[,1],  # both.cor - DD
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR)[,1],  # both.cor - DDPR
                          summary(both.cor)$coefficients[, 1],  # both.cor - DRCGEE
                          
                          DR(Y=dat$Y,Z=dat$Z, interZ =interZ, data=dat, Z.hat=Z.hat.PR,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[,1], # out.mis - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ =interZ, data=dat, ID=ID, Z.hat=Z.hat.DD, Y0.hat=Y0.hat.DD.mis)[,1], # out.mis - DD
                          DD(Y=dat$Y,Z=dat$Z, interZ =interZ, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR.mis)[,1], # out.mis - DDPR
                          summary(mis.out)$coefficients[, 1], # out.mis - DRCGEE
                          
                          DR(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR,Y0.hat=Y0.hat.PR)[,1], # trt.mis - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.DD.mis, Y0.hat=Y0.hat.DD)[,1], # trt.mis - DD
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR)[,1], # trt.mis - DDPR
                          summary(mis.sel)$coefficients[, 1], # trt.mis - DRCGEE
                          
                          DR(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[,1], # both.mis - PR
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.DD.mis, Y0.hat=Y0.hat.DD.mis)[,1], # both.mis - DD
                          DD(Y=dat$Y,Z=dat$Z, interZ = interZ, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR.mis)[,1],  # both.mis - DDPR
                          summary(both.mis)$coefficients[, 1]) # both.mis - DRCGEE
  
  DRtestCATE_beta1.rlst[i, ] <- DRtestCATE.est[, 1]
  DRtestCATE_beta2.rlst[i, ] <- DRtestCATE.est[, 2]
  
  # summarize ATE
  DRtestATE.est.rlst[i,] <- c(DR(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, Z.hat=Z.hat.PR, Y1.hat=Y1.hat.PR, Y0.hat=Y0.hat.PR)[1,1], # both.cor - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.DD, Y0.hat=Y0.hat.DD)[1,1], # both.cor - DD
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR)[1,1], # both.cor - DDPR
                              summary(both.corATE)$coefficients[1, 1], # both.cor - DRCGEE
                              
                              DR(Y=dat$Y,Z=dat$Z, interZ =interZ_c, data=dat, Z.hat=Z.hat.PR,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[1,1],  # out.mis - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ =interZ_c, data=dat, ID=ID, Z.hat=Z.hat.DD, Y0.hat=Y0.hat.DD.mis)[1,1], # out.mis - DD
                              DD(Y=dat$Y,Z=dat$Z, interZ =interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR, Y0.hat=Y0.hat.PR.mis)[1,1], # out.mis - DDPR
                              summary(mis.outATE)$coefficients[1, 1], # out.mis - DRCGEE
                              
                              DR(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR,Y0.hat=Y0.hat.PR)[1,1], # out.mis - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.DD.mis, Y0.hat=Y0.hat.DD)[1,1], # out.mis - DD
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR)[1,1], # trt.mis - DDPR
                              summary(mis.selATE)$coefficients[1, 1],  # trt.mis - DRCGEE
                              
                              DR(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, Z.hat=Z.hat.PR.mis,Y1.hat=Y1.hat.PR.mis,Y0.hat=Y0.hat.PR.mis)[1,1], # both.mis - PR
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.DD.mis, Y0.hat=Y0.hat.DD.mis)[1,1], # both.mis - DD
                              DD(Y=dat$Y,Z=dat$Z, interZ = interZ_c, data=dat, ID=ID, Z.hat=Z.hat.PR.mis, Y0.hat=Y0.hat.PR.mis)[1,1],  # both.mis - DDPR
                              summary(both.misATE)$coefficients[1, 1])  # both.mis - DRCGEE
  
}

MLtest_rlst <- cbind(MLtestCATE_beta1.rlst, MLtestCATE_beta2.rlst, MLtestATE.est.rlst) # 12*3
DRtest_rlst <- cbind(DRtestCATE_beta1.rlst, DRtestCATE_beta2.rlst, DRtestATE.est.rlst) # 16*3



