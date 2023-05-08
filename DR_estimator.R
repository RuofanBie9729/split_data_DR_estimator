library(SuperLearner)
library(keras)

DR_Est <- function(Trial, Target, p = 5, Missing = FALSE, split = FALSE
			 Trial_S = NULL, Target_S = NULL){
  ##############################################################################
  # Computes single-data doubly robust estimator of 
  # average treatment effect transported from trial to target
  # Args:
  #   p: numeric number indicating the number of covariates in X
  #   Trial: the trial dataframe containing (X, A, Y) or (X, A, Y, R)
  #   Target: the target dataframe containing (X)
  #   Missing: If true, between-trial missingness exists and 
  #	         complete-case corrections is needed. Default is FALSE.
  #   split: If true, produce split data doubly robust estimators and Trial_S
  #          and Target_S must not be NULL. Default is FALSE.
  #   Trial_S: the test trial dataframe when "split = TRUE".
  #   Target_S: the test target dataframe when "split = TRUE"
  #
  # Returns:
  #   A vector of single-data/split-data DR estimators, using logistic model, 
  #   random forest, 1-layer neural network, 2-layer neural neteork,
  #   3-layer neural network, SVM with linear kernal, 
  #	SVM with non-linear kernal, gam, GBDT, Super learner
  ##############################################################################
  # 
  ### check for NAs ###
  anyNA <- apply(Trial[ ,1:p], 2, function(x){any(is.na(x))}) 
  allNA <- apply(Trial[ ,1:p], 2, function(x){all(is.na(x))})
  # get variable names
  varname <- names(Trial)[1:p]
  trtname <- names(Trial)[p+1]
  outname <- names(Trial)[p+2]
  if (Missing == TRUE){
     # If between-trial missingness exists, get observed variables
     varname <- varname[!allNA]
  }
  if ("R" %in% colnames(Trial)){
     # If within-trial missingness exsits, get complete variables 
     # for complete-case corrections
     cmp_varname <- varname[!anyNA]
  }

  ### combine Trial and Target together as Whole dataset ###
  Target[, trtname] <- rep(0, nrow(Target))
  Target[, outname] <- rep(0, nrow(Target))
  if ("R" %in% colnames(Trial)){
     # If within-trial missingness exists, get complete cases for model fitting
     cmp_Trial <- Trial[Trial$R==1, ]
     Whole <- rbind(cmp_Trial[, c(varname, trtname, outname)], 
		     	  Target[, c(varname, trtname, outname)])
     In_Trial <- c(rep(1, nrow(cmp_Trial)), rep(0, nrow(Target)))
     Whole <- cbind(Whole, In_Trial)

     # Get complete Whole dataset for complete-case corrections
     cmp_Whole <- rbind(Trial[, c(varname, trtname, outname)], 
		     		Target[, c(varname, trtname, outname)])
  }else{
     Whole <- rbind(Trial[, c(varname, trtname, outname)], 
		        Target[, c(varname, trtname, outname)])
     In_Trial <- c(rep(1, nrow(Trial)), rep(0, nrow(Target)))
     Whole <- cbind(Whole, In_Trial)
  }
  
  ### compute hat_p ###
  # logistic model
  hat_p_fit <- SuperLearner(Y = In_Trial, X = Whole[, varname], 
                            family = binomial(),
				    SL.library = "SL.glmnet")
  hat_p <- as.vector(predict(hat_p_fit, Whole[, varname])$library.predict)
  hat_p[hat_p < 0.001] <- 0.001

  # random forest
  hat_p_fit_rf <- SuperLearner(Y = In_Trial, X = Whole[, varname], 
                               family = binomial(),
					 SL.library = "SL.ranger")
  hat_p_rf <- as.vector(predict(hat_p_fit_rf, Whole[,varname])$library.predict)
  hat_p_rf[hat_p_rf < 0.001] <- 0.001

  # one-layer neural network
  hat_p_fit_nnL1 <- keras_model_sequential() 
  hat_p_fit_nnL1 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(Whole[ ,varname])) %>%
  			   layer_dense(units = 1, activation = "sigmoid")
  hat_p_fit_nnL1 %>% compile(optimizer = "rmsprop",
                             loss = "binary_crossentropy",
				     metrics = c("binary_crossentropy"))
  hat_p_fit_nnL1 %>% fit(as.matrix(Whole[ ,varname]), In_Trial, epochs = 100,
				 batch_size = 256, validation_split = 0.2,
				 verbose = 0)
  hat_p_nnL1 <- hat_p_fit_nnL1 %>% predict(as.matrix(Whole[, varname]))
  hat_p_nnL1 <- as.vector(hat_p_nnL1)
  hat_p_nnL1[hat_p_nnL1 < 0.001] <- 0.001

  # two-layer neural network
  hat_p_fit_nnL2 <- keras_model_sequential() 
  hat_p_fit_nnL2 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(Whole[ ,varname])) %>%
                     layer_dense(units = 64, activation = "relu") %>%
  			   layer_dense(units = 1, activation = "sigmoid")
  hat_p_fit_nnL2 %>% compile(optimizer = "rmsprop",
                             loss = "binary_crossentropy",
				     metrics = c("binary_crossentropy"))
  hat_p_fit_nnL2 %>% fit(as.matrix(Whole[ ,varname]), In_Trial, epochs = 100,
				 batch_size = 256, validation_split = 0.2,
				 verbose = 0)
  hat_p_nnL2 <- hat_p_fit_nnL2 %>% predict(as.matrix(Whole[, varname]))
  hat_p_nnL2 <- as.vector(hat_p_nnL2)
  hat_p_nnL2[hat_p_nnL2 < 0.001] <- 0.001  

  # three-layer neural network
  hat_p_fit_nnL3 <- keras_model_sequential() 
  hat_p_fit_nnL3 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(Whole[ ,varname])) %>%
                     layer_dense(units = 64, activation = "relu") %>%
			   layer_dense(units = 64, activation = "relu") %>%
  			   layer_dense(units = 1, activation = "sigmoid")
  hat_p_fit_nnL3 %>% compile(optimizer = "rmsprop",
                             loss = "binary_crossentropy",
				     metrics = c("binary_crossentropy"))
  hat_p_fit_nnL3 %>% fit(as.matrix(Whole[ ,varname]), In_Trial, epochs = 100,
				 batch_size = 256, validation_split = 0.2,
				 verbose = 0)
  hat_p_nnL3 <- hat_p_fit_nnL3 %>% predict(as.matrix(Whole[, varname]))
  hat_p_nnL3 <- as.vector(hat_p_nnL3)
  hat_p_nnL3[hat_p_nnL3 < 0.001] <- 0.001  

  # svm linear kernel
  SL.svm.Linear = function(...) {
    SL.svm(..., kernel = "linear")
  }
  hat_p_fit_svmL <- SuperLearner(Y = In_Trial, X = Whole[, varname], 
                                 family = binomial(),
					   SL.library = "SL.svm.Linear")
  hat_p_svmL <- predict(hat_p_fit_svmL, Whole[,varname])$library.predict
  hat_p_svmL <- as.vector(hat_p_svmL)
  hat_p_svmL[hat_p_svmL < 0.001] <- 0.001

  # svm radial kernel
  hat_p_fit_svm <- SuperLearner(Y = In_Trial, X = Whole[, varname], 
                                family = binomial(),
					  SL.library = "SL.svm")
  hat_p_svm <- predict(hat_p_fit_svm, Whole[,varname])$library.predict
  hat_p_svm <- as.vector(hat_p_svm)
  hat_p_svm[hat_p_svm < 0.001] <- 0.001

  # gam
  hat_p_fit_gam <- SuperLearner(Y = In_Trial, X = Whole[, varname], 
                                family = binomial(),
					  SL.library = "SL.gam")
  hat_p_gam <- as.vector(predict(hat_p_fit_gam, Whole[,varname])$library.predict)
  hat_p_gam[hat_p_gam < 0.001] <- 0.001

  # GBDT
  hat_p_fit_gbdt <- SuperLearner(Y = In_Trial, X = Whole[, varname], 
                                 family = binomial(),
					   SL.library = "SL.xgboost")
  hat_p_gbdt <- as.vector(predict(hat_p_fit_gbdt, Whole[,varname])$library.predict)
  hat_p_gbdt[hat_p_gbdt < 0.001] <- 0.001

  # SuperLearner
  hat_p_fit_SL <- SuperLearner(Y = In_Trial, X = Whole[, varname], 
                               family = binomial(),
					 SL.library = c("SL.glmnet", "SL.ranger", 
                                              "SL.gam", "SL.svm", 
                                              "SL.nnet", "SL.xgboost"))
  hat_p_SL <- predict(hat_p_fit_SL, Whole[,varname], onlySL = TRUE)$pred
  hat_p_SL <- as.vector(hat_p_SL)
  hat_p_SL[hat_p_SL < 0.001] <- 0.001.

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
  ga_rf[ga_rf<0.001] <- 0.001

  ga_fit_nn <- nnet(as.factor(Y)~V1+V2+V3+V4+V5, data=Trials_arm, size=8)
  ga_nn <- as.vector(predict(ga_fit_nn, newdata=Whole, type="raw"))
  ga_nn[ga_nn<0.001] <- 0.001

  ### compute ea1
  #trt_rf <- get_prob(Trials, "trt")
  ea1_fit <- glm(as.formula(paste0("trt~", covariates)), data=Trials, family="binomial")
  ea1_fit_rf <- randomForest(as.factor(trt)~V1+V2+V3+V4+V5, data=Trials)  
  ea1_fit_nn <- nnet(as.factor(trt)~V1+V2+V3+V4+V5, data=Trials, size=8)  

  ea1 <- predict(ea1_fit, Whole, type="response") 
  ea1_rf <- predict(ea1_fit_rf, newdata=Whole, type="prob")[,2]
  ea1_nn <- as.vector(predict(ea1_fit_nn, newdata=Whole, type="raw"))
  ea0 <- 1-ea1
  ea0_rf <- 1-ea1_rf
  ea0_nn <- 1-ea1_nn

  if(a==1){
    wa <- (1-hat_p)/hat_p*ea1
    wa_rf <- (1-hat_p_rf)/hat_p_rf*ea1_rf
    wa_nn <- (1-hat_p_nn)/hat_p_nn*ea1_nn
  }else{
    wa <- (1-hat_p)/hat_p*ea0
    wa_rf <- (1-hat_p_rf)/hat_p_rf*ea0_rf
    wa_nn <- (1-hat_p_nn)/hat_p_nn*ea0_nn
  }
  
   ### outcome estimator
  out <- sum(ga*(In_Trial==0))/nrow(Target)
  out_rf <- sum(ga_rf*(In_Trial==0))/nrow(Target)
  out_nn <- sum(ga_nn*(In_Trial==0))/nrow(Target)

  ### IPW estimator
  IPW <- sum((In_Trial==1)*(Whole$trt==a)*wa*Whole$Y)/nrow(Target)
  IPW_rf <- sum((In_Trial==1)*(Whole$trt==a)*wa_rf*Whole$Y)/nrow(Target)
  IPW_nn <- sum((In_Trial==1)*(Whole$trt==a)*wa_nn*Whole$Y)/nrow(Target)

  ### normalized DR estimator
  DR <- out + sum((In_Trial==1)*(Whole$trt==a)*wa*(Whole$Y-ga))/nrow(Target)
  DR_rf <- out_rf + sum((In_Trial==1)*(Whole$trt==a)*wa_rf*(Whole$Y-ga_rf))/nrow(Target)
  DR_nn <- out_nn + sum((In_Trial==1)*(Whole$trt==a)*wa_nn*(Whole$Y-ga_nn))/nrow(Target)

 return(c(DR, DR_rf, DR_nn))
}
