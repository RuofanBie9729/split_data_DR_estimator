set.seed(123456)
#library(dplyr, warn.conflicts = FALSE)
#options(dplyr.summarise.inform = FALSE)
#library(mice)
library(mvtnorm)
library(randomForest)
source("beta0_solver.R")

K=3 # number of trials
p=5 # number of covariates
rho=0.5 # covariance between covariates
N=3000 #total sample size
P <- rep(1,10) #trial assignment probability when X1>0
P <- P/sum(P)
marg=0.5

###coefficients for trial participation deciding model

beta_vec <- c(0.7, -0.7, 0.2, 0.7, -0.7, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2) ### x1+x2+x3+x4+x5+x1:x2+x1:x3+x2:x3+x1^2+x2^2+x3^2

###coefficients for outcome generating model
theta1 <- c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.2) ### x1+x2+x3+x4+x5+x2^2
theta0 <- c(1, -0.5, -0.5, -0.5, -0.5, -0.5, -0.2)

###coefficients for trial assignment
tau_coef <- c(log(1.3), log(1.3), log(1.3), log(1.3), log(1.3))
tau_coef <- rbind(tau_coef, c(log(0.8), log(0.8), log(0.8), log(0.8), log(0.8)))
tau_coef <- rbind(tau_coef, rep(0, 5)) 

###coefficients for missing data R=0
#xi1_R_arm1 <- c(0.6, 0.3, 0.3, 0.3, 0.3)
#xi1_R_arm0 <- c(0.6, 0.3, 0.3, 0.3, 0.3)

#xi2_R_arm1 <- c(0.6, 0.4, 0.4, 0.4, 0.4)
#xi2_R_arm0 <- c(0.6, 0.4, 0.4, 0.4, 0.4)

#xi3_R_arm1 <- c(0.6, 0.5, 0.5, 0.5, 0.5)
#xi3_R_arm0 <- c(0.6, 0.5, 0.5, 0.5, 0.5)

Sigma <- matrix(rho, 3, 3)
diag(Sigma) <- 1


#Res <- NULL
#for(i in 1: 1000){
Subjects <- rmvnorm(N, rep(0, 3), Sigma)
Subjects <- cbind(Subjects, matrix(rbinom(N*2, 1, 0.5), ncol=2))
external_dat <- cbind(rep(1, N), Subjects, (Subjects[,1]*Subjects[,2]), 
   (Subjects[,1]*Subjects[,3]), (Subjects[,2]*Subjects[,3]), 
   (Subjects[,1])^2, (Subjects[,2])^2, (Subjects[,3])^2)
Intcpt <- determine_intercept(beta_vec=beta_vec, 
  marg=marg, lower_bound=-20, upper_bound=20, 
  external_dataset=as.data.frame(external_dat))
beta <- c(as.numeric(Intcpt[1]), beta_vec)

expit <- function(x){
  if(exp(x)==Inf){
    return(1)
  }else if(exp(x)==-Inf){
    return(-1)
  }else{
  return(exp(x)/(1+exp(x)))
  }
}

paticipat_assign <- function(x){
  n <- dim(x)[1]
  x <- cbind(rep(1,n), x)
  expR <- sapply(x%*%beta, expit)
  random_num <- runif(n)
  return(random_num<expR)
}

betaX <- cbind(Subjects, (Subjects[,1]*Subjects[,2]), 
   (Subjects[,1]*Subjects[,3]), (Subjects[,2]*Subjects[,3]), 
   (Subjects[,1])^2, (Subjects[,2])^2, (Subjects[,3])^2)
paticipat_prob <- paticipat_assign(betaX)
TrialsData <- Subjects[which(paticipat_prob==1), ]
Target <- Subjects[which(paticipat_prob==0), ]

X <- cbind(rep(1, nrow(TrialsData)), TrialsData)

trial_assign <- function(x){
  expTau <- apply(x%*%t(tau_coef), 1, exp)
  prob <- apply(expTau, 2, function(x){x/sum(x)})
  return(prob)
}

trial_prob <- trial_assign(TrialsData)
trial_onehot <- apply(trial_prob, 2, rmultinom, n=1, size=1)
trial <- apply(trial_onehot, 2, function(x)(which(x==1)))
#trial <- trial_assign(P, nrow(TrialsData))
TrialsData <- cbind(TrialsData, trial)
#library(dplyr)
TrialsData <- TrialsData[order(TrialsData[,p+1]),]  ##should be changed to 10


### treatment assignment
Trt_assign <- sapply(table(trial), function(x){rbinom(x, 1, 0.5)}, simplify=F)
TrialsData <- cbind(TrialsData, unlist(Trt_assign))
Target <- cbind(Target, rbinom(nrow(Target), 1, 0.5))

X_Trl <- cbind(rep(1, nrow(TrialsData)), TrialsData[,1:p], (TrialsData[,2])^2)
Trial_Y1 <- X_Trl%*%theta1
Trial_Y0 <- X_Trl%*%theta0
Trial_Y <- unlist(Trt_assign)*Trial_Y1+(1-unlist(Trt_assign))*Trial_Y0
Y_gen <- function(Y){
  n <- length(Y)
  expY <- sapply(Y, expit)
  random_num <- runif(n)
  return(random_num<expY)
}
TrialsData <- cbind(TrialsData, Y_gen(Trial_Y))

X_Tgt <- cbind(rep(1, nrow(Target)), Target[,1:p], (Target[,2])^2)
Target_Y1 <- X_Tgt%*%theta1
Target_Y0 <- X_Tgt%*%theta0
Target_Y <- Target[,p+1]*Target_Y1+(1-Target[,p+1])*Target_Y0
Target <- cbind(Target, Y_gen(Target_Y))

TrialData <- as.data.frame(TrialsData)
colnames(TrialData) <- c("V1", "V2", "V3", "V4", "V5", "trial", "trt", "Y")

Target <- as.data.frame(Target)
colnames(Target) <- c("V1", "V2", "V3", "V4", "V5", "trt", "Y")

TrialMiss <- split(TrialData, TrialData$trial)
 
#for(j in c(1:MissTrial)){
#TrialMiss[[j]][ , VarList[[MissVar]]] <- NA
#}

for(j in 1:length(TrialMiss)){
  if(j==1){
    X_trl_A1 <- as.matrix(TrialMiss[[j]][which(TrialMiss[[j]]$trt==1),c("V2", "V3", "V4", "V5")])
    X_trl_A0 <- as.matrix(TrialMiss[[j]][which(TrialMiss[[j]]$trt==0),c("V2", "V3", "V4", "V5")])
    X_trl_A1 <- cbind(rep(1, nrow(X_trl_A1)), X_trl_A1)
    X_trl_A0 <- cbind(rep(1, nrow(X_trl_A0)), X_trl_A0)
    Trial_R1 <- X_trl_A1 %*% xi1_R_arm1
    Trial_R0 <- X_trl_A0 %*% xi1_R_arm0
    TrialMiss[[j]][which(TrialMiss[[j]]$trt==1),"R"] <- Y_gen(Trial_R1)
    TrialMiss[[j]][which(TrialMiss[[j]]$trt==0),"R"] <- Y_gen(Trial_R0)
    TrialMiss[[j]][which(TrialMiss[[j]]$R==0),"V1"] <- NA
  }else if(j==2){
    X_trl_A1 <- as.matrix(TrialMiss[[j]][which(TrialMiss[[j]]$trt==1),c("V1", "V3", "V4", "V5")])
    X_trl_A0 <- as.matrix(TrialMiss[[j]][which(TrialMiss[[j]]$trt==0),c("V1", "V3", "V4", "V5")])
    X_trl_A1 <- cbind(rep(1, nrow(X_trl_A1)), X_trl_A1)
    X_trl_A0 <- cbind(rep(1, nrow(X_trl_A0)), X_trl_A0)
    Trial_R1 <- X_trl_A1 %*% xi2_R_arm1
    Trial_R0 <- X_trl_A0 %*% xi2_R_arm0   
    TrialMiss[[j]][which(TrialMiss[[j]]$trt==1),"R"] <- Y_gen(Trial_R1)
    TrialMiss[[j]][which(TrialMiss[[j]]$trt==0),"R"] <- Y_gen(Trial_R0)
    TrialMiss[[j]][which(TrialMiss[[j]]$R==0),"V2"] <- NA
  }else if(j==3){
    X_trl_A1 <- as.matrix(TrialMiss[[j]][which(TrialMiss[[j]]$trt==1),c("V1", "V2", "V4", "V5")])
    X_trl_A0 <- as.matrix(TrialMiss[[j]][which(TrialMiss[[j]]$trt==0),c("V1", "V2", "V4", "V5")])
    X_trl_A1 <- cbind(rep(1, nrow(X_trl_A1)), X_trl_A1)
    X_trl_A0 <- cbind(rep(1, nrow(X_trl_A0)), X_trl_A0)
    Trial_R1 <- X_trl_A1 %*% xi3_R_arm1
    Trial_R0 <- X_trl_A0 %*% xi3_R_arm0
    TrialMiss[[j]][which(TrialMiss[[j]]$trt==1),"R"] <- Y_gen(Trial_R1)
    TrialMiss[[j]][which(TrialMiss[[j]]$trt==0),"R"] <- Y_gen(Trial_R0)
    TrialMiss[[j]][which(TrialMiss[[j]]$R==0),"V3"] <- NA
  }
}
#Target$R <- rep(1, nrow(Target))

Phi_hat <- function(a, Trials){
  variable_name <- names(Trials)[1:p]
  cpt_variable_name <- variable_name[!is.na(apply(Trials[,1:p], 2, mean))] ### used for complete-variable estimators 
                                                                        ### when between trial missingness exists
  covariates <-  paste(variable_name, collapse="+")
  covariates_outcome <- paste0(covariates, "+I(V2^2)")
  cpt_covariates <- paste(cpt_variable_name, collapse="+")
  Trial <- Trials
  #Trials <- Trials[Trials$R==1,]
  Whole <- rbind(Trials[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y")], Target)
  In_Trial <- c(rep(1, nrow(Trials)), rep(0, nrow(Target)))
  Whole <- cbind(Whole, In_Trial)
  Trial_whole <- rbind(Trial[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y")], Target)
  #### compute hat_p
  hat_p_fit <- glm(as.formula(paste0("In_Trial~", covariates_outcome)), data=Whole, family="binomial")
  hat_p <- sapply(predict(hat_p_fit, newdata=Whole), expit)
  
  #### ML prediction
  hat_p_fit_rf <- randomForest(In_Trial~., data=Whole, mtry=5)
  hat_p_rf <- predict(hat_p_fit_rf, newdata=Whole, type="response")
  
  #hat_p_crt_fit <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial_whole, family="binomial") 
  #hat_p_crt <- sapply(predict(hat_p_crt_fit, newdata=Whole), expit)

  ### compute ga
  Trials_arm <- Trials[which(Trials$trt==a),]
  ga_fit <- glm(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm, family="binomial")
  ga <- sapply(predict(ga_fit, newdata=Whole), expit)


  ### compute ea1
  ea1_fit <- glm(as.formula(paste0("trt~", covariates)), data=Trials, family="binomial")
  ea1 <- sapply(predict(ea1_fit, Whole), expit)  
  ea0 <- 1-ea1
  #ea_crt_fit <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial, family="binomial")
  #ea_crt <- sapply(predict(ea_crt_fit, newdata=Whole), expit) 
  ### compute pr

  #Trial_arm <- Trial[which(Trial$trt==a),]
  #pr_fit <- glm(as.formula(paste0("R~", cpt_covariates)), data=Trial_arm, family="binomial")
  #pr <- sapply(predict(pr_fit, Whole), expit) 

  if(a==1){
    #denominator <- (hat_p/hat_p_crt)*(ea1/ea_crt)*pr
    #wa <- (1-hat_p/hat_p_crt)/denominator 
    #wa_c <- (1-hat_p)/hat_p*ea1
    wa <- (1-hat_p)/hat_p*ea1
    wa_rf <- (1-hat_p_rf)/hat_p_rf*ea1
  }else{
    #denominator <- (hat_p/hat_p_crt)*(ea0/ea_crt)*pr
    #wa <- (1-hat_p/hat_p_crt)/denominator
    #wa_c <- (1-hat_p)/hat_p*ea0
    wa <- (1-hat_p)/hat_p*ea0
    wa_rf <- (1-hat_p_rf)/hat_p_rf*ea0
  }
  
   ### outcome estimator
  out <- sum(ga*(In_Trial==0))/nrow(Target)

  ### IPW estimator
  #IPW <- sum((In_Trial==1)*(Whole$trt==a)*wa*Whole$Y)/nrow(Target)
  #IPW_c <- sum((In_Trial==1)*(Whole$trt==a)*wa_c*Whole$Y)/nrow(Target)
  IPW <- sum((In_Trial==1)*(Whole$trt==a)*wa*Whole$Y)/nrow(Target)
  IPW_rf <- sum((In_Trial==1)*(Whole$trt==a)*wa_rf*Whole$Y)/nrow(Target)

  ### normalized DR estimator
  #DR <- out + sum((In_Trial==1)*(Whole$trt==a)*wa*(Whole$Y-ga))/nrow(Target)
  #DR_c <- out + sum((In_Trial==1)*(Whole$trt==a)*wa_c*(Whole$Y-ga))/nrow(Target)
  DR <- out + sum((In_Trial==1)*(Whole$trt==a)*wa*(Whole$Y-ga))/nrow(Target)
  DR_rf <- out + sum((In_Trial==1)*(Whole$trt==a)*wa_rf*(Whole$Y-ga))/nrow(Target)

 #return(c(out, IPW, DR, IPW_c, DR_c))
 return(c(out, IPW, DR, IPW_rf, DR_rf))
}

Phi <- function(a){
  ### True Y
  if(a==1){
    TrueY <- mean(Y_gen(Target_Y1))
  }else{
    TrueY <- mean(Y_gen(Target_Y0))
  }
  return(TrueY)
}

pool <- function(a){
  if(a==1){
    poolest <- mean(TrialData[which(TrialData$trt==1),"Y"])
  }else{
    poolest <- mean(TrialData[which(TrialData$trt==0),"Y"])
  }
  return(poolest)
}

ak <- sapply(TrialMiss, nrow, simplify=T)/nrow(TrialData)
EstRes1 <- sapply(TrialMiss, Phi_hat, a=1, simplify=T)
EstRes0 <- sapply(TrialMiss, Phi_hat, a=0, simplify=T)


mi_est <- function(a, TrialsData){
  incomplete_data <- TrialsData[,c("V1", "V2", "V3", "V4", "V5")]
  imp <- mice(incomplete_data, m=5, maxit = 50, method = 'pmm', seed = 500, print=F)
  complete_data <- mice::complete(imp, action="long", include=F)
  complete_data_list <- split(complete_data, complete_data$.imp)
  Trials_list <- lapply(complete_data_list, function(x){
    cbind(x[,c("V1", "V2", "V3", "V4", "V5")],
          TrialsData[,c("trial", "trt", "Y", "R")])})
  mi_res <- sapply(Trials_list, Phi_mi, a=a, simplify=T)
  return(apply(mi_res, 1, mean))
}


Phi_mi <- function(a, Trials){
  variable_name <- names(Trials)[1:p]
  #variable_name <- variable_name[!is.na(apply(Trials[,1:p], 2, mean))] ### used for complete-variable estimators 
                                                                        ### when between trial missingness exists
  covariates <-  paste(variable_name, collapse="+")
  covariates_outcome <- paste0(covariates, "+I(V2^2)")
  Whole <- rbind(Trials[, c("V1", "V2", "V3", "V4", "V5", "trt", "Y", "R")], Target)
  In_Trial <- c(rep(1, nrow(Trials)), rep(0, nrow(Target)))
  Whole <- cbind(Whole, In_Trial)
  
  #### compute hat_p
  hat_p_fit <- glm(as.formula(paste0("In_Trial~", covariates_outcome)), data=Whole, family="binomial")
  hat_p <- sapply(predict(hat_p_fit, newdata=Whole), expit)

  ### compute ga
  Trials_arm <- Trials[which(Trials$trt==a),]
  ga_fit <- glm(as.formula(paste0("Y~", covariates_outcome)), data=Trials_arm, family="binomial")
  ga <- sapply(predict(ga_fit, newdata=Whole), expit)


  ### compute ea1
  ea1_fit <- glm(as.formula(paste0("trt~", covariates)), data=Trials, family="binomial")
  ea1 <- sapply(predict(ea1_fit, Whole), expit)  
  ea0 <- 1-ea1

  if(a==1){
    denominator <- hat_p*ea1
    #denominator <- ifelse(denominator<=0.01, 0.01, denominator)
    wa <- (1-hat_p)/denominator 
  }else{
    denominator <- hat_p*ea0
    #denominator <- ifelse(denominator<=0.01, 0.01, denominator)
    wa <- (1-hat_p)/denominator
  }
  
   ### outcome estimator
  out <- sum(ga*(In_Trial==0))/nrow(Target)

  ### IPW estimator
  IPW <- sum((In_Trial==1)*(Whole$trt==a)*wa*Whole$Y)/nrow(Target)

  ### normalized DR estimator
  DR <- out + sum((In_Trial==1)*(Whole$trt==a)*wa*(Whole$Y-ga))/nrow(Target)

  return(c(out, IPW, DR))
}

EstMIRes1 <- sapply(TrialMiss, mi_est, a=1, simplify=T)
EstMIRes0 <- sapply(TrialMiss, mi_est, a=0, simplify=T)

Res <- rbind(Res, c(pool(1)-pool(0), EstRes1%*%ak-EstRes0%*%ak, EstMIRes1%*%ak-EstMIRes0%*%ak, Phi(1)-Phi(0)))
if(i%%100==0){
  print(i)
}
}

a <- apply(Res, 2, mean)
 a[1:9]-a[10]
write.csv(Res, "benchmark.csv")
#Res2 <- read.csv("test.csv")[,-1]
# b <- apply(Res2, 2, mean)
# b[1:3]-b[4]

