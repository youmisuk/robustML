##############################################################################################################################
# Robust Machine Learning for Treatment Effects in Multilevel Observational Studies Under Cluster-level Unmeasured Confounding
# : The ECLS-K data analysis 
# by Youmi Suk & Hyunseung Kang

# The variables in the dataset include:

# :: ID
# S7_ID    : school ID

# :: treatment
# C7DESMTH : whether students took the algebra or a higher-level math course (= 1) or not (= 0)

# :: outcome
# C7R4MSCL : math IRT scale score

# :: student-level covariates
# C6R4MSCL : prior math IRT scale score
# P6EXPECT : expected degree of child
# GENDER : male or female
# RACE : race/ethnicity
# WKSESL : socio-economic status measure
# WKPOV_R : poverty level
# WKMOMED : mother's education level
# P7HFAMIL : family type
     
# :: school-level covariates
# S7PUPRI : public or private school
# R7REGION : census region - WEST, MIDWEST, NORTHEAST, SOUTH
# R7URBAN : school location - urban, suburb, small town
##############################################################################################################################

# load packages/sources
source("CURobustML_sepFun.R") # functions for the propensity score and outcome predictions based on our proposed methods
library(nnls)
library(AER)
library(h2o)
library(bit64)
library(drgee)
library(tmle)
library(bartCause)

localH2O = h2o.init()
h2o.removeAll()

# :: load data
dat <- read.csv("ECLSK_Algebra_complete.csv") 

dat$WHITE <- ifelse(dat$RACE=="WHITE", 1, 0)
dat$HISPANIC <- ifelse(dat$RACE=="HISPANIC", 1, 0)

dat$WKMOMED <- factor(dat$WKMOMED, levels=c("less than high school", "high school", "some college or voca program", "Bachelor's degree or higher"))
dat$R7REGION <- factor(dat$R7REGION, levels=c("WEST", "MIDWEST" , "NORTHEAST", "SOUTH"))

dat$C6R4MSCL_c <- scale(dat$C6R4MSCL, center=T, scale=F)
dat$WKSESL_c <- scale(dat$WKSESL, center=T, scale=F)

covs.dum_c <- model.matrix(~ C6R4MSCL_c + WKSESL_c + P6EXPECT + GENDER + WHITE + HISPANIC + WKPOV_R + WKMOMED + P7HFAMIL+ S7PUPRI + R7REGION + R7URBAN + 0,  data = dat)

# :: unadjusted, prima facie effect
primaFacie <- lm(C7R4MSCL ~ C7DESMTH, data=dat)
primaFacie_ATE <- summary(primaFacie)$coef[2, 1:2]

# ::  our proposed estimators
# fit erLearners and orLearners
er1_PR <- erLearnPR(Z=dat$C7DESMTH,X=covs.dum_c,ID=dat$S7_ID, rate=10^-5, int=F)
or1_PR <- orLearnPR(Y=dat$C7R4MSCL,Z=dat$C7DESMTH,X=covs.dum_c, ID=dat$S7_ID, rate=10^-5, int=F)
er1_DD <- erLearnDD(Z=dat$C7DESMTH,X=covs.dum_c, ID=dat$S7_ID, rate=10^-5, int=F)
or1_DD <- orLearnDD(Y=dat$C7R4MSCL,Z=dat$C7DESMTH,X=covs.dum_c, ID=dat$S7_ID, rate=10^-5, int=F)

# estimate ATE
WHITE_c <- dat$WHITE - mean(dat$WHITE)

PR_ATE=DR(Y=dat$C7R4MSCL, Z=dat$C7DESMTH,  Z.hat=er1_PR$Z.hat,Y1.hat=or1_PR$Y1.hat,Y0.hat=or1_PR$Y0.hat, data=dat)
DD_ATE=DD(Y=dat$C7R4MSCL, Z=dat$C7DESMTH, formula(~ WHITE_c), ID=dat$S7_ID, Z.hat=er1_DD$Z.hat,Y0.hat=or1_DD$Y0.hat, data=dat)[1,]
DDPR_ATE=DD(Y=dat$C7R4MSCL, Z=dat$C7DESMTH, formula(~ WHITE_c),  ID=dat$S7_ID,  Z.hat=er1_PR$Z.hat,Y0.hat=or1_PR$Y0.hat, data=dat)[1,]

DRDD_ATE <- rbind(PR=as.numeric(PR_ATE), DD=DD_ATE, DDPR=DDPR_ATE)

# estimate CATE
PR_WHITE=DR(Y=dat$C7R4MSCL, Z=dat$C7DESMTH, interZ=formula(~ as.factor(WHITE) -1), Z.hat=er1_PR$Z.hat,Y1.hat=or1_PR$Y1.hat,Y0.hat=or1_PR$Y0.hat, data=dat)
DD_WHITE=DD(Y=dat$C7R4MSCL, Z=dat$C7DESMTH, interZ=formula(~ as.factor(WHITE) -1), ID=dat$S7_ID, Z.hat=er1_DD$Z.hat,Y0.hat=or1_DD$Y0.hat, data=dat)
DDPR_WHITE=DD(Y=dat$C7R4MSCL, Z=dat$C7DESMTH, interZ=formula(~ as.factor(WHITE) -1), ID=dat$S7_ID,  Z.hat=er1_PR$Z.hat,Y0.hat=or1_PR$Y0.hat, data=dat)

DRDD_WHITE <- rbind(PR=c(PR_WHITE[1,], PR_WHITE[2,]), DD=c(DD_WHITE[1,], DD_WHITE[2,]), DDPR=c(DDPR_WHITE[1,], DDPR_WHITE[2,]))

# :: DRCGEE or BART/TMLE
# DRCGEE
drgee_ATE <- summary(drgee(oformula=formula(C7R4MSCL ~ C6R4MSCL_c + P6EXPECT + GENDER +  WHITE + HISPANIC + WKSESL_c + WKPOV_R + WKMOMED + P7HFAMIL), eformula=formula(C7DESMTH ~ C6R4MSCL_c + P6EXPECT + GENDER +  WHITE  + HISPANIC + WKSESL_c + WKPOV_R + WKMOMED + P7HFAMIL), iaformula = formula(~ WHITE_c), olink="identity", elink="identity", data=dat, estimation.method="dr", clusterid=dat$S7_ID, cond=T))$coef[1, 1:2]
drgee_WHITE <- summary(drgee(oformula=formula(C7R4MSCL ~ C6R4MSCL_c + P6EXPECT + GENDER +  WHITE + HISPANIC + WKSESL_c + WKPOV_R + WKMOMED + P7HFAMIL), eformula=formula(C7DESMTH ~ C6R4MSCL_c + P6EXPECT + GENDER +  WHITE  + HISPANIC + WKSESL_c + WKPOV_R + WKMOMED + P7HFAMIL), iaformula = formula(~ as.factor(WHITE) - 1 ), olink="identity", elink="identity", data=dat, estimation.method="dr", clusterid=dat$S7_ID, cond=T))$coef[, 1:2]

# BART
covs.dum_c_bart <- covs.dum_c
colnames(covs.dum_c_bart) <- paste0("X_", 1:ncol(covs.dum_c))
bart.fit_WHITE <- bartc(response=dat$C7R4MSCL, treatment=dat$C7DESMTH, confounders=covs.dum_c_bart, method.rsp="bart", method.trt="bart", p.scoreAsCovariate=TRUE, use.rbrt=FALSE, keepTrees=TRUE, group.by = dat$WHITE,group.effects = TRUE)
BART_ALL <- summary(bart.fit_WHITE)$estimate # with subgroups

# TMLE
tmle.fit <- tmle(Y=dat$C7R4MSCL, A=dat$C7DESMTH, W=covs.dum_c)
est.tmle.ite <- tmle.fit$Qstar[,2]-tmle.fit$Qstar[,1]
TMLE_ate <- summary(tmle.fit)$estimates$ATE
TMLE_WHITE <- summary(lm(est.tmle.ite ~ as.factor(dat$WHITE) -1))$coef[, 1:2] # we use cluster bootstrap SE
TMLE_ATE <- c(TMLE_ate$psi, sqrt(TMLE_ate$var.psi))

# summmary
# note that you'll see slightly different estimates for PR, DD, DDPR, BART, and TMLE, except for DRCGEE estimator.
ATE_sum <- rbind(primaFacie=primaFacie_ATE, DRDD_ATE, DRCGEE=drgee_ATE, BART=as.numeric(BART_ALL[3, 1:2]), TMLE=TMLE_ATE)
CATE_sum <- rbind(DRDD_WHITE, DRCGEE=as.numeric(c(drgee_WHITE[1,], drgee_WHITE[2,])), BART=as.numeric(c(BART_ALL[1, 1:2], BART_ALL[2, 1:2])), TMLE=as.numeric(c(TMLE_WHITE[1, 1:2], TMLE_WHITE[2, 1:2])))

round(ATE_sum, 5)
round(CATE_sum, 5) # we use cluster bootstrap SE for TMLE


