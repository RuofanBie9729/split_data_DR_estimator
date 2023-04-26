simulation <- function(seed){

#####one trial and one target without missing data
set.seed(seed)
library(mvtnorm)
library(ranger)
source("beta0_solver.R")

p=5 # number of covariates
rho=0.5 # covariance between covariates
N=3000 #total sample size
marg=0.5
B=100

###coefficients for trial participation deciding model
beta_vec <- c(0.7, -0.7, 0.2, 0.7, -0.7, 0.2, 0.2) ### x1+x2+x3+x4+x5+x1:x2+x2^2

###coefficients for outcome generating model
theta1 <- c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.2) ### x1+x2+x3+x4+x5+x2^2
theta0 <- c(1, -0.5, -0.5, -0.5, -0.5, -0.5, -0.2)

###coefficients for trial assignment
tau_coef <- c(log(1.3), log(1.3), log(1.3), log(1.3), log(1.3))
tau_coef <- rbind(tau_coef, c(log(0.8), log(0.8), log(0.8), log(0.8), log(0.8)))
tau_coef <- rbind(tau_coef, rep(0, 5)) 

###coefficients for missing data R=0
xi1_R_arm1 <- c(0.6, 0.3, 0.3, 0.3, 0.3)
xi1_R_arm0 <- c(0.6, 0.3, 0.3, 0.3, 0.3)

xi2_R_arm1 <- c(0.6, 0.4, 0.4, 0.4, 0.4)
xi2_R_arm0 <- c(0.6, 0.4, 0.4, 0.4, 0.4)

xi3_R_arm1 <- c(0.6, 0.5, 0.5, 0.5, 0.5)
xi3_R_arm0 <- c(0.6, 0.5, 0.5, 0.5, 0.5)

Sigma <- matrix(rho, 3, 3)
diag(Sigma) <- 1

expit <- function(x){return(exp(x)/(1+exp(x)))}
paticipat_assign <- function(x){
  n <- dim(x)[1]
  x <- cbind(rep(1,n), x)
  expR <- sapply(x%*%beta, expit)
  random_num <- runif(n)
  return(random_num<expR)
}
Y_gen <- function(Y){
  n <- length(Y)
  expY <- sapply(Y, expit)
  random_num <- runif(n)
  return(abs(random_num<expY))
}

trial_assign <- function(x){
  expTau <- apply(x%*%t(tau_coef), 1, exp)
  prob <- apply(expTau, 2, function(x){x/sum(x)})
  return(prob)
}

Subjects <- rmvnorm(N, rep(0, 3), Sigma) 
#+ rmvnorm(N, rep(3, 3), Sigma)
Subjects <- cbind(Subjects, matrix(rbinom(N*2, 1, 0.5), ncol=2))
external_dat <- cbind(rep(1, N), Subjects, (Subjects[,1]*Subjects[,2]), (Subjects[,2])^2)
Intcpt <- determine_intercept(beta_vec=beta_vec, 
  marg=marg, lower_bound=-20, upper_bound=20, 
  external_dataset=as.data.frame(external_dat))
beta <- c(as.numeric(Intcpt[1]), beta_vec)

betaX <- cbind(Subjects, (Subjects[,1]*Subjects[,2]), (Subjects[,2])^2)
paticipat_prob <- paticipat_assign(betaX)
TrialsData <- Subjects[which(paticipat_prob==1), ]
Target <- Subjects[which(paticipat_prob==0), ]

### trial assignment
trial_prob <- trial_assign(TrialsData)
trial_onehot <- apply(trial_prob, 2, rmultinom, n=1, size=1)
trial <- apply(trial_onehot, 2, function(x)(which(x==1)))
TrialsData <- cbind(TrialsData, trial)
TrialsData <- TrialsData[order(TrialsData[,p+1]),]  ##should be changed to 10

### treatment assignment
Trt_assign <- sapply(table(trial), function(x){rbinom(x, 1, 0.5)}, simplify=F)
TrialsData <- cbind(TrialsData, unlist(Trt_assign))
Target <- cbind(Target, rbinom(nrow(Target), 1, 0.5))

X_Trl <- cbind(rep(1, nrow(TrialsData)), TrialsData[,1:p], (TrialsData[,2])^2)
Trial_Y1 <- Y_gen(X_Trl%*%theta1)
Trial_Y0 <- Y_gen(X_Trl%*%theta0)
Trial_Y <- unlist(Trt_assign)*Trial_Y1+(1-unlist(Trt_assign))*Trial_Y0
TrialsData <- cbind(TrialsData, Trial_Y)

X_Tgt <- cbind(rep(1, nrow(Target)), Target[,1:p], (Target[,2])^2)
Target_Y1 <- Y_gen(X_Tgt%*%theta1)
Target_Y0 <- Y_gen(X_Tgt%*%theta0)
Target_Y <- Target[,p+1]*Target_Y1+(1-Target[,p+1])*Target_Y0
Target <- cbind(Target, Target_Y)

TrialData <- as.data.frame(TrialsData)
colnames(TrialData) <- c("V1", "V2", "V3", "V4", "V5", "trial", "trt", "Y")

Target <- as.data.frame(Target)
colnames(Target) <- c("V1", "V2", "V3", "V4", "V5", "trt", "Y")

TrialData_list <- split(TrialData, TrialData$trial)

for(j in 1:length(TrialData_list)){
  if(j==1){
    X_trl_A1 <- as.matrix(TrialData_list[[j]][which(TrialData_list[[j]]$trt==1),c("V2", "V3", "V4", "V5")])
    X_trl_A0 <- as.matrix(TrialData_list[[j]][which(TrialData_list[[j]]$trt==0),c("V2", "V3", "V4", "V5")])
    X_trl_A1 <- cbind(rep(1, nrow(X_trl_A1)), X_trl_A1)
    X_trl_A0 <- cbind(rep(1, nrow(X_trl_A0)), X_trl_A0)
    Trial_R1 <- X_trl_A1 %*% xi1_R_arm1
    Trial_R0 <- X_trl_A0 %*% xi1_R_arm0
    TrialData_list[[j]][which(TrialData_list[[j]]$trt==1),"R"] <- Y_gen(Trial_R1)
    TrialData_list[[j]][which(TrialData_list[[j]]$trt==0),"R"] <- Y_gen(Trial_R0)
    TrialData_list[[j]][which(TrialData_list[[j]]$R==0),"V1"] <- NA
  }else if(j==2){
    X_trl_A1 <- as.matrix(TrialData_list[[j]][which(TrialData_list[[j]]$trt==1),c("V1", "V3", "V4", "V5")])
    X_trl_A0 <- as.matrix(TrialData_list[[j]][which(TrialData_list[[j]]$trt==0),c("V1", "V3", "V4", "V5")])
    X_trl_A1 <- cbind(rep(1, nrow(X_trl_A1)), X_trl_A1)
    X_trl_A0 <- cbind(rep(1, nrow(X_trl_A0)), X_trl_A0)
    Trial_R1 <- X_trl_A1 %*% xi2_R_arm1
    Trial_R0 <- X_trl_A0 %*% xi2_R_arm0   
    TrialData_list[[j]][which(TrialData_list[[j]]$trt==1),"R"] <- Y_gen(Trial_R1)
    TrialData_list[[j]][which(TrialData_list[[j]]$trt==0),"R"] <- Y_gen(Trial_R0)
    TrialData_list[[j]][which(TrialData_list[[j]]$R==0),"V2"] <- NA
  }else if(j==3){
    X_trl_A1 <- as.matrix(TrialData_list[[j]][which(TrialData_list[[j]]$trt==1),c("V1", "V2", "V4", "V5")])
    X_trl_A0 <- as.matrix(TrialData_list[[j]][which(TrialData_list[[j]]$trt==0),c("V1", "V2", "V4", "V5")])
    X_trl_A1 <- cbind(rep(1, nrow(X_trl_A1)), X_trl_A1)
    X_trl_A0 <- cbind(rep(1, nrow(X_trl_A0)), X_trl_A0)
    Trial_R1 <- X_trl_A1 %*% xi3_R_arm1
    Trial_R0 <- X_trl_A0 %*% xi3_R_arm0
    TrialData_list[[j]][which(TrialData_list[[j]]$trt==1),"R"] <- Y_gen(Trial_R1)
    TrialData_list[[j]][which(TrialData_list[[j]]$trt==0),"R"] <- Y_gen(Trial_R0)
    TrialData_list[[j]][which(TrialData_list[[j]]$R==0),"V3"] <- NA
  }
}
Target$R <- rep(1, nrow(Target))



Phi_hat <- function(a, Trials, Target){
  Trials[,"Y"] <- as.factor(Trials[,"Y"]) 
  Trials[,"trt"] <- as.factor(Trials[,"trt"])
  Trials[,"R"] <- as.factor(Trials[,"R"])
  Target[,"R"] <- as.factor(Target[,"R"])
  variable_name <- names(Trials)[1:p]
  cpt_variable_name <- variable_name[!is.na(apply(Trials[,1:p], 2, mean))] ### used for complete-variable estimators 
                                                                        ### when between trial missingness exists
  covariates <-  paste(variable_name, collapse="+")
  covariates_outcome <- paste0(covariates, "+I(V2^2)")
  cpt_covariates <- paste(cpt_variable_name, collapse="+")
  Trial <- Trials
  Trials <- Trials[Trials$R==1,]
  Whole <- rbind(Trials[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y", "R")], Target)
  In_Trial <- c(rep(1, nrow(Trials)), rep(0, nrow(Target)))
  Whole <- cbind(Whole, In_Trial)
  Trial_whole <- rbind(Trial[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y", "R")], Target)
  Whole[,"In_Trial"] <- as.factor(Whole[,"In_Trial"])

  #### compute hat_p
  hat_p_fit <- glm(as.formula(paste0("In_Trial~", covariates)), data=Whole, family="binomial")
  hat_p <- predict(hat_p_fit, newdata=Whole, type="response")
  hat_p_crt_fit <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial_whole, family="binomial") 
  hat_p_crt <- predict(hat_p_crt_fit, newdata=Whole, type="response")
  
  #### ML prediction
  w <- 1/table(Whole$In_Trial)
  w <- w/sum(w)
  weights <- rep(0, nrow(Whole))
  weights[Whole$In_Trial==0] <- w['0']
  weights[Whole$In_Trial==1] <- w['1']
  hat_p_fit_rf <- ranger(as.formula(paste0("In_Trial~", covariates)), data=Whole, case.weights = weights, probability=T)
  hat_p_rf <- predict(hat_p_fit_rf, data=Whole, type="response")$predictions[,2]
  w <- 1/table(Trial_whole$R)
  w <- w/sum(w)
  weights <- rep(0, nrow(Trial_whole))
  weights[Trial_whole$R==0] <- w['0']
  weights[Trial_whole$R==1] <- w['1']
  hat_p_crt_fit_rf <- ranger(as.formula(paste0("R~", cpt_covariates)), data=Trial_whole, case.weights = weights, probability=T) 
  hat_p_crt_rf <- predict(hat_p_crt_fit_rf, data=Whole, type="response")$predictions[,2]

  ### compute ga
  Trials_arm <- Trials[which(Trials$trt==a),]
  ga_fit <- glm(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm, family="binomial")
  ga <- predict(ga_fit, newdata=Whole, type="response")

  ### ML ga prediction
  w <- 1/table(Trials_arm$Y)
  w <- w/sum(w)
  weights <- rep(0, nrow(Trials_arm))
  weights[Trials_arm$Y==0] <- w['0']
  weights[Trials_arm$Y==1] <- w['1']
  Trials_arm[,"V22"] <- Trials_arm[,"V2"]^2
  Whole[,"V22"] <- Whole[,"V2"]^2
  ga_fit_rf <- ranger(Y ~ V1+V2+V3+V4+V5+V22, data=Trials_arm, case.weights = weights, probability=T)
  ga_rf <- predict(ga_fit_rf, data=Whole, type="response")$predictions[,2]

  ### compute ea1
  w <- 1/table(Trials$trt)
  w <- w/sum(w)
  weights <- rep(0, nrow(Trials))
  weights[Trials$trt==0] <- w['0']
  weights[Trials$trt==1] <- w['1']
  ea1_fit <- glm(as.formula(paste0("trt~", covariates)), data=Trials, family="binomial")
  ea1_fit_rf <- ranger(as.formula(paste0("trt~", covariates)), data=Trials, case.weights = weights, probability=T)
  ea1 <- predict(ea1_fit, Whole, type="response") 
  ea1_rf <- predict(ea1_fit_rf, data=Whole, type="response")$predictions[,2]
  ea0 <- 1-ea1
  ea0_rf <- 1-ea1_rf
  ea_crt_fit <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial, family="binomial")
  ea_crt <- predict(ea_crt_fit, newdata=Whole, type="response")
  w <- 1/table(Trial$R)
  w <- w/sum(w)
  weights <- rep(0, nrow(Trial))
  weights[Trial$R==0] <- w['0']
  weights[Trial$R==1] <- w['1']
  ea_crt_fit_rf <- ranger(as.formula(paste0("R~", cpt_covariates)), data=Trial, case.weights = weights, probability=T)
  ea_crt_rf <- predict(ea_crt_fit_rf, data=Whole, type="response")$predictions[,2]

  ### compute pr
  Trial_arm <- Trial[which(Trial$trt==a),]
  w <- 1/table(Trial_arm$R)
  w <- w/sum(w)
  weights <- rep(0, nrow(Trial_arm))
  weights[Trial_arm$R==0] <- w['0']
  weights[Trial_arm$R==1] <- w['1']
  pr_fit <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial_arm, family="binomial")
  pr <- predict(pr_fit, Whole, type="response")
  pr_fit_rf <- ranger(as.formula(paste0("R~", cpt_covariates)), data=Trial_arm, case.weights = weights, probability=T)
  pr_rf <- predict(pr_fit_rf, Whole, type="response")$predictions[,2]

  if(a==1){
    denominator <- (hat_p/hat_p_crt)*(ea1/ea_crt)*pr
    denominator_rf <- (hat_p_rf/hat_p_crt_rf)*(ea1_rf/ea_crt_rf)*pr_rf
    wa <- (1-hat_p/hat_p_crt)/denominator
    wa_rf <- (1-hat_p_rf/hat_p_crt_rf)/denominator_rf
  }else{
    denominator <- (hat_p/hat_p_crt)*(ea0/ea_crt)*pr
    denominator_rf <- (hat_p_rf/hat_p_crt_rf)*(ea0_rf/ea_crt_rf)*pr_rf
    wa <- (1-hat_p/hat_p_crt)/denominator
    wa_rf <- (1-hat_p_rf/hat_p_crt_rf)/denominator_rf
  }
  
  Whole[,"trt"] <- as.numeric(Whole[,"trt"])-1
  Whole[,"Y"] <- as.numeric(Whole[,"Y"])-1
   ### outcome estimator
  out <- sum(ga*(In_Trial==0))/nrow(Target)
  out_rf <- sum(ga_rf*(In_Trial==0))/nrow(Target)

  ### IPW estimator
  IPW <- sum((In_Trial==1)*(Whole$trt==a)*wa*Whole$Y)/nrow(Target)
  IPW_rf <- sum((In_Trial==1)*(Whole$trt==a)*wa_rf*Whole$Y)/nrow(Target)

  ### normalized DR estimator
  DR <- out + sum((In_Trial==1)*(Whole$trt==a)*wa*(Whole$Y-ga))/nrow(Target)
  DR_rf <- out_rf + sum((In_Trial==1)*(Whole$trt==a)*wa_rf*(Whole$Y-ga_rf))/nrow(Target)

 return(c(out, IPW, DR, out_rf, IPW_rf, DR_rf))
}

Phi <- function(a){
  ### True Y
  if(a==1){
    TrueY <- mean(Target_Y1)
  }else{
    TrueY <- mean(Target_Y0)
  }
  return(TrueY)
}

DR_split <- function(a, Trials, Target){
  variable_name <- names(Trials)[1:p]
  cpt_variable_name <- variable_name[!is.na(apply(Trials[,1:p], 2, mean))] ### used for complete-variable estimators 
                                                                        ### when between trial missingness exists
  covariates <-  paste(variable_name, collapse="+")
  covariates_outcome <- paste0(covariates, "+I(V2^2)")
  cpt_covariates <- paste(cpt_variable_name, collapse="+")

  ### split trials and target
  n_trl <- nrow(Trials)
  n_tgt <- nrow(Target)
  S1_trl_id <- sample(n_trl, round(n_trl/2))
  S1_tgt_id <- sample(n_tgt, round(n_tgt/2))
  Trials_S1 <- Trials[S1_trl_id,]
  Target_S1 <- Target[S1_tgt_id,]
  Trials_S2 <- Trials[-S1_trl_id,]
  Target_S2 <- Target[-S1_tgt_id,]
  Trial_S1 <- Trials_S1
  Trial_S2 <- Trials_S2
  Trials_S1 <- Trials_S1[Trials_S1$R==1,]
  Whole_S1 <- rbind(Trials_S1[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y", "R")], Target_S1)
  In_Trial_S1 <- c(rep(1, nrow(Trials_S1)), rep(0, nrow(Target_S1)))
  Whole_S1 <- cbind(Whole_S1, In_Trial_S1)
  Trial_whole_S1 <- rbind(Trial_S1[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y", "R")], Target_S1)

  Trials_S2 <- Trials_S2[Trials_S2$R==1,]
  Whole_S2 <- rbind(Trials_S2[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y", "R")], Target_S2)
  In_Trial_S2 <- c(rep(1, nrow(Trials_S2)), rep(0, nrow(Target_S2)))
  Whole_S2 <- cbind(Whole_S2, In_Trial_S2)
  Trial_whole_S2 <- rbind(Trial_S2[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y", "R")], Target_S2)

  #### compute hat_p
  hat_p_fit_S1 <- glm(as.formula(paste0("Y~", covariates)), data=Whole_S1, family="binomial")
  hat_p_S2 <- predict(hat_p_fit_S1, newdata=Whole_S2, type="response")
  hat_p_crt_fit_S1 <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial_whole_S1, family="binomial") 
  hat_p_crt_S2 <- predict(hat_p_crt_fit_S1, newdata=Whole_S2, type="response")

  hat_p_fit_S2 <- glm(as.formula(paste0("Y~", covariates)), data=Whole_S2, family="binomial")
  hat_p_S1 <- predict(hat_p_fit_S2, newdata=Whole_S1, type="response")
  hat_p_crt_fit_S2 <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial_whole_S2, family="binomial") 
  hat_p_crt_S1 <- predict(hat_p_crt_fit_S2, newdata=Whole_S1, type="response")
  
  #### ML prediction
  hat_p_fit_rf_S1 <- randomForest(as.formula(paste0("Y~", covariates_outcome)), data=Whole_S1, mtry=5)
  hat_p_rf_S2 <- predict(hat_p_fit_rf_S1, newdata=Whole_S2, type="response")
  hat_p_crt_fit_rf_S1 <- randomForest(as.formula(paste0("R~", cpt_covariates)), data=Trial_whole_S1, mtry=4) 
  hat_p_crt_rf_S2 <- predict(hat_p_crt_fit_rf_S1, newdata=Whole_S2, type="response")

  hat_p_fit_rf_S2 <- randomForest(as.formula(paste0("Y~", covariates_outcome)), data=Whole_S2, mtry=5)
  hat_p_rf_S1 <- predict(hat_p_fit_rf_S2, newdata=Whole_S1, type="response")
  hat_p_crt_fit_rf_S2 <- randomForest(as.formula(paste0("R~", cpt_covariates)), data=Trial_whole_S2, mtry=4) 
  hat_p_crt_rf_S1 <- predict(hat_p_crt_fit_rf_S2, newdata=Whole_S1, type="response")

  ### compute ga
  Trials_arm_S1 <- Trials_S1[which(Trials_S1$trt==a),]
  Trials_arm_S2 <- Trials_S2[which(Trials_S2$trt==a),]
  ga_fit_S1 <- glm(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm_S1, family="binomial")
  ga_S2 <- predict(ga_fit_S1, newdata=Whole_S2, type="response")
  ga_fit_S2 <- glm(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm_S2, family="binomial")
  ga_S1 <- predict(ga_fit_S2, newdata=Whole_S1, type="response")

  #### ML ga prediction
  ga_fit_rf_S1 <- randomForest(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm_S1, mtry=5)
  ga_rf_S2 <- predict(ga_fit_rf_S1, newdata=Whole_S2, type="response")

  ga_fit_rf_S2 <- randomForest(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm_S2, mtry=5)
  ga_rf_S1 <- predict(ga_fit_rf_S2, newdata=Whole_S1, type="response")

  ### compute ea1
  ea1_fit_S1 <- glm(as.formula(paste0("trt~", covariates)), data=Trials_S1, family="binomial")
  ea1_fit_rf_S1 <- randomForest(as.formula(paste0("trt~", covariates)), data=Trials_S1, mtry=5)
  ea1_S2 <- predict(ea1_fit_S1, Whole_S2, type="response")
  ea1_rf_S2 <- predict(ea1_fit_rf_S1, newdata=Whole_S2, type="response")
  ea0_S2 <- 1-ea1_S2
  ea0_rf_S2 <- 1-ea1_rf_S2

  ea_crt_fit_S1 <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial_S1, family="binomial")
  ea_crt_S2 <- predict(ea_crt_fit_S1, newdata=Whole_S2, type="response")
  ea_crt_fit_rf_S1 <- randomForest(as.formula(paste0("R~", cpt_covariates)), data=Trial_S1, mtry=4)
  ea_crt_rf_S2 <- predict(ea_crt_fit_rf_S1, newdata=Whole_S2, type="response")

  ea1_fit_S2 <- glm(as.formula(paste0("trt~", covariates)), data=Trials_S2, family="binomial")
  ea1_fit_rf_S2 <- randomForest(as.formula(paste0("trt~", covariates)), data=Trials_S2, mtry=5)
  ea1_S1 <- predict(ea1_fit_S2, Whole_S1, type="response")
  ea1_rf_S1 <- predict(ea1_fit_rf_S2, newdata=Whole_S1, type="response")
  ea0_S1 <- 1-ea1_S1
  ea0_rf_S1 <- 1-ea1_rf_S1

  ea_crt_fit_S2 <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial_S2, family="binomial")
  ea_crt_S1 <- predict(ea_crt_fit_S2, newdata=Whole_S1, type="response") 
  ea_crt_fit_rf_S2 <- randomForest(as.formula(paste0("R~", cpt_covariates)), data=Trial_S2, mtry=4)
  ea_crt_rf_S1 <- predict(ea_crt_fit_rf_S2, newdata=Whole_S1, type="response")

  ### compute pr
  Trial_arm_S1 <- Trial_S1[which(Trial_S1$trt==a),]  
  pr_fit_S1 <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial_arm_S1, family="binomial")
  pr_S2 <- predict(pr_fit_S1, Whole_S2, type="response")
  pr_fit_rf_S1 <- randomForest(as.formula(paste0("R~", cpt_covariates)), data=Trial_arm_S1, mtry=4)
  pr_rf_S2 <- predict(pr_fit_rf_S1, Whole_S2, type="response")

  Trial_arm_S2 <- Trial_S2[which(Trial_S2$trt==a),]  
  pr_fit_S2 <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial_arm_S2, family="binomial")
  pr_S1 <- predict(pr_fit_S2, Whole_S1, type="response")
  pr_fit_rf_S2 <- randomForest(as.formula(paste0("R~", cpt_covariates)), data=Trial_arm_S2, mtry=4)
  pr_rf_S1 <- predict(pr_fit_rf_S2, Whole_S1, type="response")

  if(a==1){
    denominator_S1 <- (hat_p_S1/hat_p_crt_S1)*(ea1_S1/ea_crt_S1)*pr_S1
    denominator_rf_S1 <- (hat_p_rf_S1/hat_p_crt_rf_S1)*(ea1_rf_S1/ea_crt_rf_S1)*pr_rf_S1
    wa_S1 <- (1-hat_p_S1/hat_p_crt_S1)/denominator_S1
    wa_rf_S1 <- (1-hat_p_rf_S1/hat_p_crt_rf_S1)/denominator_rf_S1
    denominator_S2 <- (hat_p_S2/hat_p_crt_S2)*(ea1_S2/ea_crt_S2)*pr_S2
    denominator_rf_S2 <- (hat_p_rf_S2/hat_p_crt_rf_S2)*(ea1_rf_S2/ea_crt_rf_S2)*pr_rf_S2
    wa_S2 <- (1-hat_p_S2/hat_p_crt_S2)/denominator_S2
    wa_rf_S2 <- (1-hat_p_rf_S2/hat_p_crt_rf_S2)/denominator_rf_S2
  }else{
    denominator_S1 <- (hat_p_S1/hat_p_crt_S1)*(ea0_S1/ea_crt_S1)*pr_S1
    denominator_rf_S1 <- (hat_p_rf_S1/hat_p_crt_rf_S1)*(ea0_rf_S1/ea_crt_rf_S1)*pr_rf_S1
    wa_S1 <- (1-hat_p_S1/hat_p_crt_S1)/denominator_S1
    wa_rf_S1 <- (1-hat_p_rf_S1/hat_p_crt_rf_S1)/denominator_rf_S1
    denominator_S2 <- (hat_p_S2/hat_p_crt_S2)*(ea0_S2/ea_crt_S2)*pr_S2
    denominator_rf_S2 <- (hat_p_rf_S2/hat_p_crt_rf_S2)*(ea0_rf_S2/ea_crt_rf_S2)*pr_rf_S2
    wa_S2 <- (1-hat_p_S2/hat_p_crt_S2)/denominator_S2
    wa_rf_S2 <- (1-hat_p_rf_S2/hat_p_crt_rf_S2)/denominator_rf_S2
  }
  
   ### outcome estimator
  out_S1 <- sum(ga_S1*(In_Trial_S1==0))/nrow(Target_S1)
  out_S2 <- sum(ga_S2*(In_Trial_S2==0))/nrow(Target_S2)

  out_rf_S1 <- sum(ga_rf_S1*(In_Trial_S1==0))/nrow(Target_S1)
  out_rf_S2 <- sum(ga_rf_S2*(In_Trial_S2==0))/nrow(Target_S2)

  ### normalized DR estimator
  DR_S1 <- out_S1 + sum((In_Trial_S1==1)*(Whole_S1$trt==a)*wa_S1*(Whole_S1$Y-ga_S1))/nrow(Target_S1)
  DR_S2 <- out_S2 + sum((In_Trial_S2==1)*(Whole_S2$trt==a)*wa_S2*(Whole_S2$Y-ga_S2))/nrow(Target_S2)
  DR_rf_S1 <- out_rf_S1 + sum((In_Trial_S1==1)*(Whole_S1$trt==a)*wa_rf_S1*(Whole_S1$Y-ga_rf_S1))/nrow(Target_S1)
  DR_rf_S2 <- out_rf_S2 + sum((In_Trial_S2==1)*(Whole_S2$trt==a)*wa_rf_S2*(Whole_S2$Y-ga_rf_S2))/nrow(Target_S2)
  return(c((DR_S1+DR_S2)/2, (DR_rf_S1+DR_rf_S2)/2))
}
ak <- sapply(TrialData_list, nrow, simplify=T)/nrow(TrialData)

boot_trl <- function(Trial){
  boot_Trial <- Trial[sample(nrow(Trial), nrow(Trial), replace=T),]
  return(boot_Trial)
}

Bootstrap <- function(b){
  boot_TrialData <- sapply(TrialData_list, boot_trl, simplify=F)
  boot_Target <- Target[sample(nrow(Target), nrow(Target), replace=T),]
  boot_est1 <- sapply(boot_TrialData, Phi_hat, a=1, Target=boot_Target, simplify=T)
  boot_est0 <- sapply(boot_TrialData, Phi_hat, a=0, Target=boot_Target, simplify=T)
  boot_split1 <- sapply(boot_TrialData, DR_split, Target=boot_Target, a=1, simplify=T)
  boot_split0 <- sapply(boot_TrialData, DR_split, Target=boot_Target, a=0, simplify=T)
  return(c(boot_est1 %*% ak - boot_est0 %*% ak, boot_split1 %*% ak - boot_split0 %*% ak))
}


est1 <- sapply(TrialData_list, Phi_hat, a=1, Target=Target, simplify=T)
est0 <- sapply(TrialData_list, Phi_hat, a=0, Target=Target, simplify=T)

true <- Phi(1)-Phi(0)
split_est1 <- sapply(TrialData_list, DR_split, Target=Target, a=1, simplify=T)
split_est0 <- sapply(TrialData_list, DR_split, Target=Target, a=0, simplify=T)

#boot_res <- sapply(1:B, Bootstrap)

Res <- c(true, est1%*%ak-est0%*%ak, split_est1%*%ak-split_est0%*%ak)
Res[2:8]-Res[1]
return(Res)
#return(list(Res, boot_res))
}

results <- pbmcapply::pbmclapply(
	1:10,
	function(seed){
		simulation(seed)
	},
	mc.cores = 10
)
#saveRDS(results, "hat_p_misspecified_n1500.rds")

bias <- sapply(results, function(x) x)
bias <- apply(bias, 1, mean)
bias[2:8]-bias[1]
