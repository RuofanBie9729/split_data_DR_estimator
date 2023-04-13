#####one trial and one target without missing data
set.seed(123456)
library(mvtnorm)
library(randomForest)
source("beta0_solver.R")

p=5 # number of covariates
rho=0.5 # covariance between covariates
N=3000 #total sample size
marg=0.5

###coefficients for trial participation deciding model
beta_vec <- c(0.7, -0.7, 0.2, 0.7, -0.7, 0.2, 0.2) ### x1+x2+x3+x4+x5+x1:x2+x2^2

###coefficients for outcome generating model
theta1 <- c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.2) ### x1+x2+x3+x4+x5+x2^2
theta0 <- c(1, -0.5, -0.5, -0.5, -0.5, -0.5, -0.2)

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
  return(random_num<expY)
}

Res <- NULL
for(rep in 1:100){

Subjects <- rmvnorm(N, rep(0, 3), Sigma)
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

### treatment assignment
Trt_assign <- rbinom(nrow(TrialsData), 1, 0.5)
TrialsData <- cbind(TrialsData, Trt_assign)
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
colnames(TrialData) <- c("V1", "V2", "V3", "V4", "V5", "trt", "Y")

Target <- as.data.frame(Target)
colnames(Target) <- c("V1", "V2", "V3", "V4", "V5", "trt", "Y")

Phi_hat <- function(a, Trials){
  variable_name <- names(Trials)[1:p]
  cpt_variable_name <- variable_name[!is.na(apply(Trials[,1:p], 2, mean))] ### used for complete-variable estimators 
                                                                        ### when between trial missingness exists
  covariates <-  paste(variable_name, collapse="+")
  covariates_outcome <- paste0(covariates, "+I(V2^2)")
  cpt_covariates <- paste(cpt_variable_name, collapse="+")
  Whole <- rbind(Trials[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y")], Target)
  In_Trial <- c(rep(1, nrow(Trials)), rep(0, nrow(Target)))
  Whole <- cbind(Whole, In_Trial)
  #### compute hat_p
  hat_p_fit <- glm(as.formula(paste0("Y~", covariates_outcome)), data=Whole, family="binomial")
  hat_p <- sapply(predict(hat_p_fit, newdata=Whole), expit)
  
  #### ML prediction
  hat_p_fit_rf <- randomForest(as.formula(paste0("Y~", covariates_outcome)), data=Whole, mtry=5)
  hat_p_rf <- predict(hat_p_fit_rf, newdata=Whole, type="response")

  ### compute ga
  Trials_arm <- Trials[which(Trials$trt==a),]
  ga_fit <- glm(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm, family="binomial")
  ga <- sapply(predict(ga_fit, newdata=Whole), expit)

  ### ML ga prediction
  ga_fit_rf <- randomForest(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm, mtry=5)
  ga_rf <- predict(ga_fit_rf, newdata=Whole, type="response")

  ### compute ea1
  ea1_fit <- glm(as.formula(paste0("trt~", covariates)), data=Trials, family="binomial")
  ea1_fit_rf <- randomForest(as.formula(paste0("trt~", covariates)), data=Trials, mtry=5)
  ea1 <- sapply(predict(ea1_fit, Whole), expit)  
  ea1_rf <- predict(ea1_fit_rf, newdata=Whole, type="response")
  ea0 <- 1-ea1
  ea0_rf <- 1-ea1_rf

  if(a==1){
    wa <- (1-hat_p)/hat_p*ea1
    wa_rf <- (1-hat_p_rf)/hat_p_rf*ea1_rf
  }else{
    wa <- (1-hat_p)/hat_p*ea0
    wa_rf <- (1-hat_p_rf)/hat_p_rf*ea0_rf
  }
  
   ### outcome estimator
  out <- sum(ga*(In_Trial==0))/nrow(Target)
  out_rf <- sum(ga_rf*(In_Trial==0))/nrow(Target)

  ### IPW estimator
  IPW <- sum((In_Trial==1)*(Whole$trt==a)*wa*Whole$Y)/nrow(Target)
  IPW_rf <- sum((In_Trial==1)*(Whole$trt==a)*wa_rf*Whole$Y)/nrow(Target)

  ### normalized DR estimator
  DR <- out + sum((In_Trial==1)*(Whole$trt==a)*wa*(Whole$Y-ga))/nrow(Target)
  DR_rf <- out_rf + sum((In_Trial==1)*(Whole$trt==a)*wa_rf*(Whole$Y-ga_rf))/nrow(Target)

 return(c(out, IPW, DR, IPW_rf, DR_rf))
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

DR_split <- function(a, Trials){
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
  Whole_S1 <- rbind(Trials_S1[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y")], Target_S1)
  Whole_S2 <- rbind(Trials_S2[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y")], Target_S2)
  In_Trial_S1 <- c(rep(1, nrow(Trials_S1)), rep(0, nrow(Target_S1)))
  In_Trial_S2 <- c(rep(1, nrow(Trials_S2)), rep(0, nrow(Target_S2)))
  Whole_S1 <- cbind(Whole_S1, In_Trial_S1)
  Whole_S2 <- cbind(Whole_S2, In_Trial_S2)
  #### compute hat_p
  hat_p_fit_S1 <- glm(as.formula(paste0("Y~", covariates_outcome)), data=Whole_S1, family="binomial")
  hat_p_S2 <- sapply(predict(hat_p_fit_S1, newdata=Whole_S2), expit)

  hat_p_fit_S2 <- glm(as.formula(paste0("Y~", covariates_outcome)), data=Whole_S2, family="binomial")
  hat_p_S1 <- sapply(predict(hat_p_fit_S2, newdata=Whole_S1), expit)
  
  #### ML prediction
  hat_p_fit_rf_S1 <- randomForest(as.formula(paste0("Y~", covariates_outcome)), data=Whole_S1, mtry=5)
  hat_p_rf_S2 <- predict(hat_p_fit_rf_S1, newdata=Whole_S2, type="response")

  hat_p_fit_rf_S2 <- randomForest(as.formula(paste0("Y~", covariates_outcome)), data=Whole_S2, mtry=5)
  hat_p_rf_S1 <- predict(hat_p_fit_rf_S2, newdata=Whole_S1, type="response")

  ### compute ga
  Trials_arm_S1 <- Trials_S1[which(Trials_S1$trt==a),]
  Trials_arm_S2 <- Trials_S2[which(Trials_S2$trt==a),]
  ga_fit_S1 <- glm(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm_S1, family="binomial")
  ga_S2 <- sapply(predict(ga_fit_S1, newdata=Whole_S2), expit)
  ga_fit_S2 <- glm(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm_S2, family="binomial")
  ga_S1 <- sapply(predict(ga_fit_S2, newdata=Whole_S1), expit)

  #### ML ga prediction
  ga_fit_rf_S1 <- randomForest(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm_S1, mtry=5)
  ga_rf_S2 <- predict(ga_fit_rf_S1, newdata=Whole_S2, type="response")

  ga_fit_rf_S2 <- randomForest(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm_S2, mtry=5)
  ga_rf_S1 <- predict(ga_fit_rf_S2, newdata=Whole_S1, type="response")

  ### compute ea1
  ea1_fit_S1 <- glm(as.formula(paste0("trt~", covariates)), data=Trials_S1, family="binomial")
  ea1_fit_rf_S1 <- randomForest(as.formula(paste0("trt~", covariates)), data=Trials_S1, mtry=5)
  ea1_S2 <- sapply(predict(ea1_fit_S1, Whole_S2), expit)  
  ea1_rf_S2 <- predict(ea1_fit_rf_S1, newdata=Whole_S2, type="response")
  ea0_S2 <- 1-ea1_S2
  ea0_rf_S2 <- 1-ea1_rf_S2

  ea1_fit_S2 <- glm(as.formula(paste0("trt~", covariates)), data=Trials_S2, family="binomial")
  ea1_fit_rf_S2 <- randomForest(as.formula(paste0("trt~", covariates)), data=Trials_S2, mtry=5)
  ea1_S1 <- sapply(predict(ea1_fit_S2, Whole_S1), expit)  
  ea1_rf_S1 <- predict(ea1_fit_rf_S2, newdata=Whole_S1, type="response")
  ea0_S1 <- 1-ea1_S1
  ea0_rf_S1 <- 1-ea1_rf_S1

  if(a==1){
    wa_S1 <- (1-hat_p_S1)/hat_p_S1*ea1_S1
    wa_rf_S1 <- (1-hat_p_rf_S1)/hat_p_rf_S1*ea1_rf_S1
    wa_S2 <- (1-hat_p_S2)/hat_p_S2*ea1_S2
    wa_rf_S2 <- (1-hat_p_rf_S2)/hat_p_rf_S2*ea1_rf_S2
  }else{
    wa_S1 <- (1-hat_p_S1)/hat_p_S1*ea0_S1
    wa_rf_S1 <- (1-hat_p_rf_S1)/hat_p_rf_S1*ea0_rf_S1
    wa_S2 <- (1-hat_p_S2)/hat_p_S2*ea0_S2
    wa_rf_S2 <- (1-hat_p_rf_S2)/hat_p_rf_S2*ea0_rf_S2
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

est <- Phi_hat(1, TrialData)-Phi_hat(0, TrialData)
true <- Phi(1)-Phi(0)
split_est <- DR_split(1, TrialData)-DR_split(0, TrialData)
Res <- rbind(Res, c(true, est, split_est))
}

write.csv(Res, "toyexample_res.R")
