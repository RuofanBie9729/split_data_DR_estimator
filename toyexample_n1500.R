simulation <- function(seed){
#####one trial and one target without missing data
set.seed(seed)
library(mvtnorm)
library(randomForest)
source("beta0_solver.R")

p=5 # number of covariates
rho=0.5 # covariance between covariates
N=3000 #total sample size
marg=0.5
B=100

###coefficients for trial participation deciding model
beta_vec <- c(1.2, 1.2, -1.2, 1.2, -1.2, 1.2, -1.2) ### x1+x2+x3+x4+x5+x1:x2+x2^2

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


Subjects <- rmvnorm(N, rep(0, 3), Sigma)
Subjects <- cbind(Subjects, matrix(rbinom(N*2, 1, 0.5), ncol=2))
external_dat <- cbind(rep(1, N), sin(Subjects[,1]), cos(Subjects[,2]), sin(Subjects[,3]), 
    Subjects[,4]*sin(Subjects[,1]), Subjects[,5]*cos(Subjects[,2]), sin(Subjects[,1])^2, cos(Subjects[,3])^2)
Intcpt <- determine_intercept(beta_vec=beta_vec, 
  marg=marg, lower_bound=-20, upper_bound=20, 
  external_dataset=as.data.frame(external_dat))
beta <- c(as.numeric(Intcpt[1]), beta_vec)

betaX <- cbind(sin(Subjects[,1]), cos(Subjects[,2]), sin(Subjects[,3]), 
    Subjects[,4]*sin(Subjects[,1]), Subjects[,5]*cos(Subjects[,2]), sin(Subjects[,1])^2, cos(Subjects[,3])^2)
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

case_weight <- function(data, Y){
  w <- 1/table(data[,"Y"])
  w <- w/sum(w)
  weights <- rep(0, nrow(data))
  weights[data[,"Y"]==0] = w["0"]
  weights[data[,"Y"]==1] = w["1"]
  return(weights)
}

get_prob <- function(data, Y, k=50){
  featureMatrix <- data[,1:5]
  KNNInds <- FNN::get.knn(featureMatrix, k = k)$nn.index
  est_prob <- sapply(
    1:nrow(featureMatrix),
    function(i){
      nnObs <- data[KNNInds[i, ],]
      mean(nnObs[,Y])
    })
  return(est_prob)
}

#TrialData$trt_rf <- get_prob(TrialData, "trt")
#TrialData_arm1 <- TrialData[TrialData$trt==1,]
#TrialData_arm0 <- TrialData[TrialData$trt==0,]
#Y_rf1 <- get_prob(TrialData_arm1, "Y")
#Y_rf0 <- get_prob(TrialData_arm0, "Y")
#Y_rf <- rep(0, nrow(TrialData))
#Y_rf[TrialData$trt==1] <- Y_rf1
#Y_rf[TrialData$trt==0] <- Y_rf0
#TrialData$Y_rf <- Y_rf

#TrialData$trt_rf <- rep(0.5, nrow(TrialData))
#Y_rf1 <- sapply(X_Trl%*%theta1, expit)
#Y_rf0 <- sapply(X_Trl%*%theta0, expit)
#Y_rf <- Trt_assign*Y_rf1+(1-Trt_assign)*Y_rf0
#TrialData$Y_rf <- Y_rf

#Whole <- rbind(TrialData[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y")], Target)
#Whole_X <- as.matrix(cbind(rep(1, nrow(Whole)), Whole[,1:5], Whole[,2]^2, Whole[,1]*Whole[,2]))
#In_Trial_rf <- sapply(Whole_X %*% beta, expit)

#Whole$In_Trial <- c(rep(1, nrow(TrialData)), rep(0, nrow(Target)))
#In_Trial_rf <- get_prob(Whole, "In_Trial")
#TrialData$In_Trial_rf <- In_Trial_rf[1:nrow(TrialData)]
#Target$In_Trial_rf <- In_Trial_rf[(1+nrow(TrialData)):(nrow(Whole))]


Phi_hat <- function(a, Trials, Target){
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
  hat_p_fit <- glm(as.formula(paste0("In_Trial~", covariates)), data=Whole, family="binomial")
  hat_p <- predict(hat_p_fit, newdata=Whole, type="response")

  #### ML prediction
  #In_Trial_rf <- get_prob(Whole, "In_Trial")
  hat_p_fit_rf <- randomForest(as.factor(In_Trial)~V1+V2+V3+V4+V5, data=Whole)
  hat_p_rf <- predict(hat_p_fit_rf, newdata=Whole, type="prob")[,2]
  hat_p_rf[hat_p_rf<=0.005] <- 0.005

  ### compute ga
  Trials_arm <- Trials[which(Trials$trt==a),]
  ga_fit <- glm(as.formula(paste0("Y~", covariates)), data=Trials_arm, family="binomial")
  ga <- predict(ga_fit, newdata=Whole, type="response")

  ### ML ga prediction
  #Y_rf <- get_prob(Trials_arm, "Y")
  #Trials_arm$V22 <- Trials_arm[,"V2"]^2
  #Whole$V22 <- Whole[,"V2"]^2
  ga_fit_rf <- randomForest(as.factor(Y)~V1+V2+V3+V4+V5, data=Trials_arm)
  ga_rf <- predict(ga_fit_rf, newdata=Whole, type="prob")[,2]
  ga_rf[ga_rf<=0.005] <- 0.005

  ### compute ea1
  #trt_rf <- get_prob(Trials, "trt")
  ea1_fit <- glm(as.formula(paste0("trt~", covariates)), data=Trials, family="binomial")
  ea1_fit_rf <- randomForest(as.factor(trt)~V1+V2+V3+V4+V5, data=Trials)  
  ea1 <- predict(ea1_fit, Whole, type="response") 
  ea1_rf <- predict(ea1_fit_rf, newdata=Whole, type="prob")[,2]
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
  Whole_S1 <- rbind(Trials_S1[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y")], Target_S1)
  Whole_S2 <- rbind(Trials_S2[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y")], Target_S2)
  In_Trial_S1 <- c(rep(1, nrow(Trials_S1)), rep(0, nrow(Target_S1)))
  In_Trial_S2 <- c(rep(1, nrow(Trials_S2)), rep(0, nrow(Target_S2)))
  Whole_S1 <- cbind(Whole_S1, In_Trial_S1)
  Whole_S2 <- cbind(Whole_S2, In_Trial_S2)
  #### compute hat_p
  hat_p_fit_S1 <- glm(as.formula(paste0("In_Trial_S1~", covariates)), data=Whole_S1, family="binomial")
  hat_p_S2 <- predict(hat_p_fit_S1, newdata=Whole_S2, type="response")

  hat_p_fit_S2 <- glm(as.formula(paste0("In_Trial_S2~", covariates)), data=Whole_S2, family="binomial")
  hat_p_S1 <- predict(hat_p_fit_S2, newdata=Whole_S1, type="response")
  
  #### ML prediction
  #In_Trial_rf_S1 <- get_prob(Whole_S1, "In_Trial_S1")
  hat_p_fit_rf_S1 <- randomForest(as.factor(In_Trial_S1)~V1+V2+V3+V4+V5, data=Whole_S1)
  hat_p_rf_S2 <- predict(hat_p_fit_rf_S1, newdata=Whole_S2, type="prob")[,2]
  hat_p_rf_S2[hat_p_rf_S2 <= 0.005] <- 0.005

  #In_Trial_rf_S2 <- get_prob(Whole_S2, "In_Trial_S2")
  hat_p_fit_rf_S2 <- randomForest(as.factor(In_Trial_S2)~V1+V2+V3+V4+V5, data=Whole_S2)
  hat_p_rf_S1 <- predict(hat_p_fit_rf_S2, newdata=Whole_S1, type="prob")[,2]
  hat_p_rf_S1[hat_p_rf_S1 <= 0.005] <- 0.005

  ### compute ga
  Trials_arm_S1 <- Trials_S1[which(Trials_S1$trt==a),]
  Trials_arm_S2 <- Trials_S2[which(Trials_S2$trt==a),]
  ga_fit_S1 <- glm(as.formula(paste0("Y~", covariates)), data=Trials_arm_S1, family="binomial")
  ga_S2 <- predict(ga_fit_S1, newdata=Whole_S2, type="response")
  ga_fit_S2 <- glm(as.formula(paste0("Y~", covariates)), data=Trials_arm_S2, family="binomial")
  ga_S1 <- predict(ga_fit_S2, newdata=Whole_S1, type="response")

  #### ML ga prediction
  #Y_rf_S1 <- get_prob(Trials_arm_S1, "Y")
  #Trials_arm_S1$V22 <- Trials_arm_S1[,"V2"]^2
  #Whole_S2$V22 <- Whole_S2[,"V2"]^2 
  ga_fit_rf_S1 <- randomForest(as.factor(Y)~V1+V2+V3+V4+V5, data=Trials_arm_S1)
  ga_rf_S2 <- predict(ga_fit_rf_S1, newdata=Whole_S2, type="prob")[,2]
  ga_rf_S2[ga_rf_S2 <= 0.005] <- 0.005

  #Y_rf_S2 <- get_prob(Trials_arm_S2, "Y")
  #Trials_arm_S2$V22 <- Trials_arm_S2[,"V2"]^2
  #Whole_S1$V22 <- Whole_S1[,"V2"]^2 
  ga_fit_rf_S2 <- randomForest(as.factor(Y)~V1+V2+V3+V4+V5, data=Trials_arm_S2)
  ga_rf_S1 <- predict(ga_fit_rf_S2, newdata=Whole_S1, type="prob")[,2]
  ga_rf_S1[ga_rf_S1<=0.005] <- 0.005 

  ### compute ea1
  #trt_rf_S1 <- get_prob(Trials_S1, "trt")
  ea1_fit_S1 <- glm(as.formula(paste0("trt~", covariates)), data=Trials_S1, family="binomial")
  ea1_fit_rf_S1 <- randomForest(as.factor(trt)~V1+V2+V3+V4+V5, data=Trials_S1)
  ea1_S2 <- predict(ea1_fit_S1, Whole_S2, type="response") 
  ea1_rf_S2 <- predict(ea1_fit_rf_S1, newdata=Whole_S2, type="prob")[,2]
  ea0_S2 <- 1-ea1_S2
  ea0_rf_S2 <- 1-ea1_rf_S2

  #trt_rf_S2 <- get_prob(Trials_S2, "trt")
  ea1_fit_S2 <- glm(as.formula(paste0("trt~", covariates)), data=Trials_S2, family="binomial")
  ea1_fit_rf_S2 <- randomForest(as.factor(trt)~V1+V2+V3+V4+V5, data=Trials_S2)
  ea1_S1 <- predict(ea1_fit_S2, Whole_S1, type="response") 
  ea1_rf_S1 <- predict(ea1_fit_rf_S2, newdata=Whole_S1, type="prob")[,2]
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

Bootstrap <- function(b, Trial, Target){
  set.seed(b)
  boot_Trial <- Trial[sample(nrow(Trial), nrow(Trial), replace=T),]
  boot_Target <- Target[sample(nrow(Target), nrow(Target), replace=T),]
  boot_est <- Phi_hat(a=1, boot_Trial, boot_Target) - Phi_hat(a=0, boot_Trial, boot_Target)
  boot_split <- DR_split(a=1, boot_Trial, boot_Target) - DR_split(a=0, boot_Trial, boot_Target)
  return(c(boot_est, boot_split))
}

est <- Phi_hat(1, TrialData, Target)-Phi_hat(0, TrialData, Target)
true <- Phi(1)-Phi(0)
split_est <- DR_split(1, TrialData, Target)-DR_split(0, TrialData, Target)
Res <- c(true, est, split_est)
#boot_res <- sapply(1:B, Bootstrap, Trial=TrialData, Target=Target)
return(Res)
}

results <- pbmcapply::pbmclapply(
	1:100,
	function(seed){
		simulation(seed)
	},
	mc.cores = 10
)
saveRDS(results, "toyexample_n1500.rds")

