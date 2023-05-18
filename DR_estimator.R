library(SuperLearner)
library(keras)
#install_keras()

DR_Est <- function(Trial, Target, nu = 0.1, p = 5, Missing = FALSE, 
                   split = FALSE, Trial_S = NULL, Target_S = NULL){
  ##############################################################################
  # Computes single-data doubly robust estimator of 
  # average treatment effect transported from trial to target
  # Args:
  #   p: numeric number indicating the number of covariates in X
  #   nu: numeric number used in SVM method settings
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
  if (split == TRUE){
     Target_S[, trtname] <- rep(0, nrow(Target_S))
     Target_S[, outname] <- rep(0, nrow(Target_S))
  }
  if ("R" %in% colnames(Trial)){
     # If within-trial missingness exists, get complete cases for model fitting
     cmp_Trial <- Trial[Trial[, "R"]==1, ]
     Whole <- rbind(cmp_Trial[, c(varname, trtname, outname)], 
		     	  Target[, c(varname, trtname, outname)])
     In_Trial <- c(rep(1, nrow(cmp_Trial)), rep(0, nrow(Target)))
     Whole <- cbind(Whole, In_Trial)

     # Get complete Whole dataset for complete-case corrections
     cmp_Whole <- rbind(Trial[, c(varname, trtname, outname)], 
		     		Target[, c(varname, trtname, outname)])

     # Get whole and complete whole datasets for split-data estimators
     if (split == TRUE){
         cmp_Trial_S <- Trial_S[Trial_S[, "R"]==1, ]
         Whole_S <- rbind(cmp_Trial_S[, c(varname, trtname, outname)], 
		     	  	  Target_S[, c(varname, trtname, outname)])
         In_Trial_S <- c(rep(1, nrow(cmp_Trial_S)), rep(0, nrow(Target_S)))
         cmp_Whole_S <- rbind(Trial_S[, c(varname, trtname, outname)], 
		     			Target_S[, c(varname, trtname, outname)])
     }
  }else{
     Whole <- rbind(Trial[, c(varname, trtname, outname)], 
		        Target[, c(varname, trtname, outname)])
     In_Trial <- c(rep(1, nrow(Trial)), rep(0, nrow(Target)))
     Whole <- cbind(Whole, In_Trial)
     
     if (split == TRUE){
        Whole_S <- rbind(Trial_S[, c(varname, trtname, outname)], 
		        	 Target_S[, c(varname, trtname, outname)])
	  In_Trial_S <- c(rep(1, nrow(Trial_S)), rep(0, nrow(Target_S)))
     }
  }
  if (split == TRUE){
     testX <- Whole_S[, varname]
  }else{
     testX <- Whole[, varname]
  }

  ### define svm with hyperparameters ###
  SL.ksvm.Linear = function(...) {
    SL.ksvm(..., kernel = "vanilladot", nu = nu)
  }

  SL.ksvm.smNu = function(...) {
    SL.ksvm(..., nu = nu)
  }

  ### compute hat_p ###
  p_trainX <- Whole[, varname]
  p_trainY <- In_Trial
  # SuperLearner
  hat_p_fit_SL <- SuperLearner(Y = p_trainY, X = p_trainX, 
                               family = binomial(),
                               SL.library = c("SL.glmnet", "SL.ranger", 
                                              "SL.ksvm.Linear", "SL.ksvm.smNu", 
                                              "SL.gam", "SL.nnet", "SL.xgboost"))
  hat_p_ML <- predict(hat_p_fit_SL, testX)
  hat_p_SL <- as.vector(hat_p_ML$pred)
  hat_p_SL[hat_p_SL < 0.001] <- 0.001

  # logistic model
  hat_p <- as.vector(hat_p_ML$library.predict[,1])
  hat_p[hat_p < 0.001] <- 0.001

  # random forest
  hat_p_rf <- as.vector(hat_p_ML$library.predict[,2])
  hat_p_rf[hat_p_rf < 0.001] <- 0.001

  # one-layer neural network
  hat_p_fit_nnL1 <- keras_model_sequential() 
  hat_p_fit_nnL1 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(p_trainX)) %>%
  			   layer_dense(units = 1, activation = "sigmoid")
  hat_p_fit_nnL1 %>% compile(optimizer = "rmsprop",
                             loss = "binary_crossentropy",
				     metrics = c("binary_crossentropy"))
  hat_p_fit_nnL1 %>% fit(as.matrix(p_trainX), p_trainY, epochs = 100,
				 batch_size = 512, validation_split = 0.2,
				 verbose = 0)
  hat_p_nnL1 <- hat_p_fit_nnL1 %>% predict(as.matrix(testX))
  hat_p_nnL1 <- as.vector(hat_p_nnL1)
  hat_p_nnL1[hat_p_nnL1 < 0.001] <- 0.001

  # two-layer neural network
  hat_p_fit_nnL2 <- keras_model_sequential() 
  hat_p_fit_nnL2 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(p_trainX)) %>%
                     layer_dense(units = 64, activation = "relu") %>%
  			   layer_dense(units = 1, activation = "sigmoid")
  hat_p_fit_nnL2 %>% compile(optimizer = "rmsprop",
                             loss = "binary_crossentropy",
				     metrics = c("binary_crossentropy"))
  hat_p_fit_nnL2 %>% fit(as.matrix(p_trainX), p_trainY, epochs = 100,
				 batch_size = 512, validation_split = 0.2,
				 verbose = 0)
  hat_p_nnL2 <- hat_p_fit_nnL2 %>% predict(as.matrix(testX))
  hat_p_nnL2 <- as.vector(hat_p_nnL2)
  hat_p_nnL2[hat_p_nnL2 < 0.001] <- 0.001  

  # three-layer neural network
  hat_p_fit_nnL3 <- keras_model_sequential() 
  hat_p_fit_nnL3 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(p_trainX)) %>%
                     layer_dense(units = 64, activation = "relu") %>%
			   layer_dense(units = 64, activation = "relu") %>%
  			   layer_dense(units = 1, activation = "sigmoid")
  hat_p_fit_nnL3 %>% compile(optimizer = "rmsprop",
                             loss = "binary_crossentropy",
				     metrics = c("binary_crossentropy"))
  hat_p_fit_nnL3 %>% fit(as.matrix(p_trainX), p_trainY, epochs = 100,
				 batch_size = 512, validation_split = 0.2,
				 verbose = 0)
  hat_p_nnL3 <- hat_p_fit_nnL3 %>% predict(as.matrix(testX))
  hat_p_nnL3 <- as.vector(hat_p_nnL3)
  hat_p_nnL3[hat_p_nnL3 < 0.001] <- 0.001  

  # svm linear kernel
  hat_p_svmL <- as.vector(hat_p_ML$library.predict[,3])
  hat_p_svmL[hat_p_svmL < 0.001] <- 0.001

  # svm radial kernel
  hat_p_svm <- as.vector(hat_p_ML$library.predict[,4])
  hat_p_svm[hat_p_svm < 0.001] <- 0.001

  # gam
  hat_p_gam <- as.vector(hat_p_ML$library.predict[,5])
  hat_p_gam[hat_p_gam < 0.001] <- 0.001

  # GBDT
  hat_p_gbdt <- as.vector(hat_p_ML$library.predict[,7])
  hat_p_gbdt[hat_p_gbdt < 0.001] <- 0.001

  ### compute ga when a = 1 ###
  ga1_trainX <- Trial[which(Trial[,trtname] == 1), varname]
  ga1_trainY <- Trial[which(Trial[,trtname] == 1), outname]
  
  # SuperLearner
  ga1_fit_SL <- SuperLearner(Y = ga1_trainY, X = ga1_trainX, 
                             family = binomial(),
				             SL.library = c("SL.glmnet", "SL.ranger", 
                                            "SL.ksvm.Linear", "SL.ksvm.smNu", 
                                            "SL.gam", "SL.nnet", "SL.xgboost"))
  ga1_ML <- predict(ga1_fit_SL, testX)
  ga1_SL <- as.vector(ga1_ML$pred)
  ga1_SL[ga1_SL < 0.001] <- 0.001
  
  # logistic model
  ga1 <- as.vector(ga1_ML$library.predict[,1])
  ga1[ga1 < 0.001] <- 0.001

  # random forest
  ga1_rf <- as.vector(ga1_ML$library.predict[,2])
  ga1_rf[ga1_rf < 0.001] <- 0.001

  # one-layer neural network
  ga1_fit_nnL1 <- keras_model_sequential() 
  ga1_fit_nnL1 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(ga1_trainX)) %>%
  			 layer_dense(units = 1, activation = "sigmoid")
  ga1_fit_nnL1 %>% compile(optimizer = "rmsprop",
                           loss = "binary_crossentropy",
				   metrics = c("binary_crossentropy"))
  ga1_fit_nnL1 %>% fit(as.matrix(ga1_trainX), ga1_trainY, epochs = 100,
			     batch_size = 256, validation_split = 0.2,
			     verbose = 0)
  ga1_nnL1 <- ga1_fit_nnL1 %>% predict(as.matrix(testX))
  ga1_nnL1 <- as.vector(ga1_nnL1)
  ga1_nnL1[ga1_nnL1 < 0.001] <- 0.001

  # two-layer neural network
  ga1_fit_nnL2 <- keras_model_sequential() 
  ga1_fit_nnL2 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(ga1_trainX)) %>%
			 layer_dense(units = 64, activation = "relu") %>%
  			 layer_dense(units = 1, activation = "sigmoid")
  ga1_fit_nnL2 %>% compile(optimizer = "rmsprop",
                           loss = "binary_crossentropy",
				   metrics = c("binary_crossentropy"))
  ga1_fit_nnL2 %>% fit(as.matrix(ga1_trainX), ga1_trainY, epochs = 100,
			     batch_size = 256, validation_split = 0.2,
			     verbose = 0)
  ga1_nnL2 <- ga1_fit_nnL2 %>% predict(as.matrix(testX))
  ga1_nnL2 <- as.vector(ga1_nnL2)
  ga1_nnL2[ga1_nnL2 < 0.001] <- 0.001

  # three-layer neural network
  ga1_fit_nnL3 <- keras_model_sequential() 
  ga1_fit_nnL3 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(ga1_trainX)) %>%
			 layer_dense(units = 64, activation = "relu") %>%
                   layer_dense(units = 64, activation = "relu") %>%
  			 layer_dense(units = 1, activation = "sigmoid")
  ga1_fit_nnL3 %>% compile(optimizer = "rmsprop",
                           loss = "binary_crossentropy",
				   metrics = c("binary_crossentropy"))
  ga1_fit_nnL3 %>% fit(as.matrix(ga1_trainX), ga1_trainY, epochs = 100,
			     batch_size = 256, validation_split = 0.2,
			     verbose = 0)
  ga1_nnL3 <- ga1_fit_nnL3 %>% predict(as.matrix(testX))
  ga1_nnL3 <- as.vector(ga1_nnL3)
  ga1_nnL3[ga1_nnL3 < 0.001] <- 0.001

  # svm linear function  
  ga1_svmL <- as.vector(ga1_ML$library.predict[,3])
  ga1_svmL[ga1_svmL < 0.001] <- 0.001

  # svm radial kernel
  ga1_svm <- as.vector(ga1_ML$library.predict[,4])
  ga1_svm[ga1_svm < 0.001] <- 0.001

  # gam
  ga1_gam <- as.vector(ga1_ML$library.predict[,5])
  ga1_gam[ga1_gam < 0.001] <- 0.001

  # GBDT
  ga1_gbdt <- as.vector(ga1_ML$library.predict[,7])
  ga1_gbdt[ga1_gbdt < 0.001] <- 0.001
  

  ### compute ga when a = 0 ###
  ga0_trainX <- Trial[which(Trial[,trtname] == 0), varname]
  ga0_trainY <- Trial[which(Trial[,trtname] == 0), outname]
  
  # SuperLearner
  ga0_fit_SL <- SuperLearner(Y = ga0_trainY, X = ga0_trainX, 
                             family = binomial(),
				             SL.library = c("SL.glmnet", "SL.ranger", 
                                            "SL.ksvm.Linear", "SL.ksvm.smNu", 
                                            "SL.gam", "SL.nnet", "SL.xgboost"))
  ga0_ML <- predict(ga0_fit_SL, testX)
  ga0_SL <- as.vector(ga0_ML$pred)
  ga0_SL[ga0_SL < 0.001] <- 0.001
  
  # logistic model
  ga0 <- as.vector(ga0_ML$library.predict[,1])
  ga0[ga0 < 0.001] <- 0.001

  # random forest
  ga0_rf <- as.vector(ga0_ML$library.predict[,2])
  ga0_rf[ga0_rf < 0.001] <- 0.001

  # one-layer neural network
  ga0_fit_nnL1 <- keras_model_sequential() 
  ga0_fit_nnL1 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(ga0_trainX)) %>%
  			 layer_dense(units = 1, activation = "sigmoid")
  ga0_fit_nnL1 %>% compile(optimizer = "rmsprop",
                           loss = "binary_crossentropy",
				   metrics = c("binary_crossentropy"))
  ga0_fit_nnL1 %>% fit(as.matrix(ga0_trainX), ga0_trainY, epochs = 100,
			     batch_size = 256, validation_split = 0.2,
			     verbose = 0)
  ga0_nnL1 <- ga0_fit_nnL1 %>% predict(as.matrix(testX))
  ga0_nnL1 <- as.vector(ga0_nnL1)
  ga0_nnL1[ga0_nnL1 < 0.001] <- 0.001

  # two-layer neural network
  ga0_fit_nnL2 <- keras_model_sequential() 
  ga0_fit_nnL2 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(ga0_trainX)) %>%
			 layer_dense(units = 64, activation = "relu") %>%
  			 layer_dense(units = 1, activation = "sigmoid")
  ga0_fit_nnL2 %>% compile(optimizer = "rmsprop",
                           loss = "binary_crossentropy",
				   metrics = c("binary_crossentropy"))
  ga0_fit_nnL2 %>% fit(as.matrix(ga0_trainX), ga0_trainY, epochs = 100,
			     batch_size = 256, validation_split = 0.2,
			     verbose = 0)
  ga0_nnL2 <- ga0_fit_nnL2 %>% predict(as.matrix(testX))
  ga0_nnL2 <- as.vector(ga0_nnL2)
  ga0_nnL2[ga0_nnL2 < 0.001] <- 0.001

  # three-layer neural network
  ga0_fit_nnL3 <- keras_model_sequential() 
  ga0_fit_nnL3 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(ga0_trainX)) %>%
			 layer_dense(units = 64, activation = "relu") %>%
                   layer_dense(units = 64, activation = "relu") %>%
  			 layer_dense(units = 1, activation = "sigmoid")
  ga0_fit_nnL3 %>% compile(optimizer = "rmsprop",
                           loss = "binary_crossentropy",
				   metrics = c("binary_crossentropy"))
  ga0_fit_nnL3 %>% fit(as.matrix(ga0_trainX), ga0_trainY, epochs = 100,
			     batch_size = 256, validation_split = 0.2,
			     verbose = 0)
  ga0_nnL3 <- ga0_fit_nnL3 %>% predict(as.matrix(testX))
  ga0_nnL3 <- as.vector(ga0_nnL3)
  ga0_nnL3[ga0_nnL3 < 0.001] <- 0.001

  # svm linear function
  ga0_svmL <- as.vector(ga0_ML$library.predict[,3])
  ga0_svmL[ga0_svmL < 0.001] <- 0.001

  # svm radial kernel
  ga0_svm <- as.vector(ga0_ML$library.predict[,4])
  ga0_svm[ga0_svm < 0.001] <- 0.001

  # gam
  ga0_gam <- as.vector(ga0_ML$library.predict[,5])
  ga0_gam[ga0_gam < 0.001] <- 0.001

  # GBDT
  ga0_gbdt <- as.vector(ga0_ML$library.predict[,7])
  ga0_gbdt[ga0_gbdt < 0.001] <- 0.001
  

  ### comput ea ###
  ea1_trainX <- Trial[, varname]
  ea1_trainY <- Trial[, trtname]

  # SuperLearner
  ea1_fit_SL <- SuperLearner(Y = ea1_trainY, X = ea1_trainX, 
                             family = binomial(),
				             SL.library = c("SL.glmnet", "SL.ranger", 
                                            "SL.ksvm.Linear", "SL.ksvm.smNu", 
                                            "SL.gam", "SL.nnet", "SL.xgboost"))
  ea1_ML <- predict(ea1_fit_SL, testX)
  ea1_SL <- as.vector(ea1_ML$pred)
  ea0_SL <- 1 - ea1_SL
  ea1_SL[ea1_SL < 0.001] <- 0.001
  ea0_SL[ea0_SL < 0.001] <- 0.001

  # logistic model
  ea1 <- as.vector(ea1_ML$library.predict[,1])
  ea0 <- 1 - ea1
  ea1[ea1 < 0.001] <- 0.001
  ea0[ea0 < 0.001] <- 0.001

  # random forest
  ea1_rf <- as.vector(ea1_ML$library.predict[,2])
  ea0_rf <- 1 - ea1_rf
  ea1_rf[ea1_rf < 0.001] <- 0.001
  ea0_rf[ea0_rf < 0.001] <- 0.001

  # one-layer neural network
  ea1_fit_nnL1 <- keras_model_sequential() 
  ea1_fit_nnL1 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(ea1_trainX)) %>%
  			 layer_dense(units = 1, activation = "sigmoid")
  ea1_fit_nnL1 %>% compile(optimizer = "rmsprop",
                           loss = "binary_crossentropy",
				   metrics = c("binary_crossentropy"))
  ea1_fit_nnL1 %>% fit(as.matrix(ea1_trainX), ea1_trainY, epochs = 100,
			     batch_size = 512, validation_split = 0.2,
			     verbose = 0)
  ea1_nnL1 <- ea1_fit_nnL1 %>% predict(as.matrix(testX))
  ea1_nnL1 <- as.vector(ea1_nnL1)
  ea0_nnL1 <- 1 - ea1_nnL1
  ea1_nnL1[ea1_nnL1 < 0.001] <- 0.001
  ea0_nnL1[ea0_nnL1 < 0.001] <- 0.001

  # two-layer neural network
  ea1_fit_nnL2 <- keras_model_sequential() 
  ea1_fit_nnL2 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(ea1_trainX)) %>%
			 layer_dense(units = 64, activation = "relu") %>%
  			 layer_dense(units = 1, activation = "sigmoid")
  ea1_fit_nnL2 %>% compile(optimizer = "rmsprop",
                           loss = "binary_crossentropy",
				   metrics = c("binary_crossentropy"))
  ea1_fit_nnL2 %>% fit(as.matrix(ea1_trainX), ea1_trainY, epochs = 100,
			     batch_size = 512, validation_split = 0.2,
			     verbose = 0)
  ea1_nnL2 <- ea1_fit_nnL2 %>% predict(as.matrix(testX))
  ea1_nnL2 <- as.vector(ea1_nnL2)
  ea0_nnL2 <- 1- ea1_nnL2
  ea1_nnL2[ea1_nnL2 < 0.001] <- 0.001
  ea0_nnL2[ea0_nnL2 < 0.001] <- 0.001

  # three-layer neural network
  ea1_fit_nnL3 <- keras_model_sequential() 
  ea1_fit_nnL3 %>% layer_dense(units = 64, activation = "relu", 
					   input_shape = ncol(ea1_trainX)) %>%
			 layer_dense(units = 64, activation = "relu") %>%
                   layer_dense(units = 64, activation = "relu") %>%
  			 layer_dense(units = 1, activation = "sigmoid")
  ea1_fit_nnL3 %>% compile(optimizer = "rmsprop",
                           loss = "binary_crossentropy",
				   metrics = c("binary_crossentropy"))
  ea1_fit_nnL3 %>% fit(as.matrix(ea1_trainX), ea1_trainY, epochs = 100,
			     batch_size = 512, validation_split = 0.2,
			     verbose = 0)
  ea1_nnL3 <- ea1_fit_nnL3 %>% predict(as.matrix(testX))
  ea1_nnL3 <- as.vector(ea1_nnL3)
  ea0_nnL3 <- 1 - ea1_nnL3
  ea1_nnL3[ea1_nnL3 < 0.001] <- 0.001
  ea0_nnL3[ea0_nnL3 < 0.001] <- 0.001

  # svm linear function
  ea1_svmL <- as.vector(ea1_ML$library.predict[,3])
  ea0_svmL <- 1 - ea1_svmL
  ea1_svmL[ea1_svmL < 0.001] <- 0.001
  ea0_svmL[ea0_svmL < 0.001] <- 0.001

  # svm radial kernel
  ea1_svm <- as.vector(ea1_ML$library.predict[,4])
  ea0_svm <- 1 - ea1_svm
  ea1_svm[ea1_svm < 0.001] <- 0.001
  ea0_svm[ea0_svm < 0.001] <- 0.001

  # gam
  ea1_gam <- as.vector(ea1_ML$library.predict[,5])
  ea0_gam <- 1 - ea1_gam
  ea1_gam[ea1_gam < 0.001] <- 0.001
  ea0_gam[ea0_gam < 0.001] <- 0.001

  # GBDT
  ea1_gbdt <- as.vector(ea1_ML$library.predict[,7])
  ea0_gbdt <- 1 - ea1_gbdt
  ea1_gbdt[ea1_gbdt < 0.001] <- 0.001
  ea0_gbdt[ea0_gbdt < 0.001] <- 0.001
  

  ### comput weight ###
  wa1 <- (1 - hat_p)/hat_p * ea1
  wa1_rf <- (1 - hat_p_rf)/hat_p_rf * ea1_rf
  wa1_nnL1 <- (1 - hat_p_nnL1)/hat_p_nnL1 * ea1_nnL1
  wa1_nnL2 <- (1 - hat_p_nnL2)/hat_p_nnL2 * ea1_nnL2
  wa1_nnL3 <- (1 - hat_p_nnL3)/hat_p_nnL3 * ea1_nnL3
  wa1_svmL <- (1 - hat_p_svmL)/hat_p_svmL * ea1_svmL
  wa1_svm <- (1 - hat_p_svm)/hat_p_svm * ea1_svm
  wa1_gam <- (1 - hat_p_gam)/hat_p_gam * ea1_gam
  wa1_gbdt <- (1 - hat_p_gbdt)/hat_p_gbdt * ea1_gbdt
  wa1_SL <- (1 - hat_p_SL)/hat_p_SL * ea1_SL

  wa0 <- (1 - hat_p)/hat_p * ea0
  wa0_rf <- (1 - hat_p_rf)/hat_p_rf * ea0_rf
  wa0_nnL1 <- (1 - hat_p_nnL1)/hat_p_nnL1 * ea0_nnL1
  wa0_nnL2 <- (1 - hat_p_nnL2)/hat_p_nnL2 * ea0_nnL2
  wa0_nnL3 <- (1 - hat_p_nnL3)/hat_p_nnL3 * ea0_nnL3
  wa0_svmL <- (1 - hat_p_svmL)/hat_p_svmL * ea0_svmL
  wa0_svm <- (1 - hat_p_svm)/hat_p_svm * ea0_svm
  wa0_gam <- (1 - hat_p_gam)/hat_p_gam * ea0_gam
  wa0_gbdt <- (1 - hat_p_gbdt)/hat_p_gbdt * ea0_gbdt
  wa0_SL <- (1 - hat_p_SL)/hat_p_SL * ea0_SL

  ### outcome estimator ###
  if (split == TRUE){
     out1 <- sum(ga1*(In_Trial_S==0))/nrow(Target_S)
     out1_rf <- sum(ga1_rf*(In_Trial_S==0))/nrow(Target_S)
     out1_nnL1 <- sum(ga1_nnL1*(In_Trial_S==0))/nrow(Target_S)
     out1_nnL2 <- sum(ga1_nnL2*(In_Trial_S==0))/nrow(Target_S)
     out1_nnL3 <- sum(ga1_nnL3*(In_Trial_S==0))/nrow(Target_S)
     out1_svmL <- sum(ga1_svmL*(In_Trial_S==0))/nrow(Target_S)
     out1_svm <- sum(ga1_svm*(In_Trial_S==0))/nrow(Target_S)
     out1_gam <- sum(ga1_gam*(In_Trial_S==0))/nrow(Target_S)
     out1_gbdt <- sum(ga1_gbdt*(In_Trial_S==0))/nrow(Target_S)
     out1_SL <- sum(ga1_SL*(In_Trial_S==0))/nrow(Target_S)

     out0 <- sum(ga0*(In_Trial_S==0))/nrow(Target_S)
     out0_rf <- sum(ga0_rf*(In_Trial_S==0))/nrow(Target_S)
     out0_nnL1 <- sum(ga0_nnL1*(In_Trial_S==0))/nrow(Target_S)
     out0_nnL2 <- sum(ga0_nnL2*(In_Trial_S==0))/nrow(Target_S)
     out0_nnL3 <- sum(ga0_nnL3*(In_Trial_S==0))/nrow(Target_S)
     out0_svmL <- sum(ga0_svmL*(In_Trial_S==0))/nrow(Target_S)
     out0_svm <- sum(ga0_svm*(In_Trial_S==0))/nrow(Target_S)
     out0_gam <- sum(ga0_gam*(In_Trial_S==0))/nrow(Target_S)
     out0_gbdt <- sum(ga0_gbdt*(In_Trial_S==0))/nrow(Target_S)
     out0_SL <- sum(ga0_SL*(In_Trial_S==0))/nrow(Target_S)
  }else{
     out1 <- sum(ga1*(In_Trial==0))/nrow(Target)
     out1_rf <- sum(ga1_rf*(In_Trial==0))/nrow(Target)
     out1_nnL1 <- sum(ga1_nnL1*(In_Trial==0))/nrow(Target)
     out1_nnL2 <- sum(ga1_nnL2*(In_Trial==0))/nrow(Target)
     out1_nnL3 <- sum(ga1_nnL3*(In_Trial==0))/nrow(Target)
     out1_svmL <- sum(ga1_svmL*(In_Trial==0))/nrow(Target)
     out1_svm <- sum(ga1_svm*(In_Trial==0))/nrow(Target)
     out1_gam <- sum(ga1_gam*(In_Trial==0))/nrow(Target)
     out1_gbdt <- sum(ga1_gbdt*(In_Trial==0))/nrow(Target)
     out1_SL <- sum(ga1_SL*(In_Trial==0))/nrow(Target)

     out0 <- sum(ga0*(In_Trial==0))/nrow(Target)
     out0_rf <- sum(ga0_rf*(In_Trial==0))/nrow(Target)
     out0_nnL1 <- sum(ga0_nnL1*(In_Trial==0))/nrow(Target)
     out0_nnL2 <- sum(ga0_nnL2*(In_Trial==0))/nrow(Target)
     out0_nnL3 <- sum(ga0_nnL3*(In_Trial==0))/nrow(Target)
     out0_svmL <- sum(ga0_svmL*(In_Trial==0))/nrow(Target)
     out0_svm <- sum(ga0_svm*(In_Trial==0))/nrow(Target)
     out0_gam <- sum(ga0_gam*(In_Trial==0))/nrow(Target)
     out0_gbdt <- sum(ga0_gbdt*(In_Trial==0))/nrow(Target)
     out0_SL <- sum(ga0_SL*(In_Trial==0))/nrow(Target)
  }

  ### normalized DR estimator
  if (split == TRUE){
     DR1 <- out1 + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==1) * wa1 * 
                       (Whole_S[ ,outname] - ga1))/nrow(Target_S)
     DR1_rf <- out1_rf + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==1) * wa1_rf * 
                       (Whole_S[ ,outname] - ga1_rf))/nrow(Target_S)
     DR1_nnL1 <- out1_nnL1 + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==1) * wa1_nnL1 * 
                       (Whole_S[ ,outname] - ga1_nnL1))/nrow(Target_S)
     DR1_nnL2 <- out1_nnL2 + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==1) * wa1_nnL2 * 
                       (Whole_S[ ,outname] - ga1_nnL2))/nrow(Target_S)
     DR1_nnL3 <- out1_nnL3 + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==1) * wa1_nnL3 * 
                       (Whole_S[ ,outname] - ga1_nnL3))/nrow(Target_S)
     DR1_svmL <- out1_svmL + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==1) * wa1_svmL * 
                       (Whole_S[ ,outname] - ga1_svmL))/nrow(Target_S)
     DR1_svm <- out1_svm + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==1) * wa1_svm * 
                       (Whole_S[ ,outname] - ga1_svm))/nrow(Target_S)
     DR1_gam <- out1_gam + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==1) * wa1_gam * 
                       (Whole_S[ ,outname] - ga1_gam))/nrow(Target_S)
     DR1_gbdt <- out1_gbdt + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==1) * wa1_gbdt * 
                       (Whole_S[ ,outname] - ga1_gbdt))/nrow(Target_S)
     DR1_SL <- out1_SL + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==1) * wa1_SL * 
                       (Whole_S[ ,outname] - ga1_SL))/nrow(Target_S)

     DR0 <- out0 + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==0) * wa0 * 
                       (Whole_S[ ,outname] - ga0))/nrow(Target_S)
     DR0_rf <- out0_rf + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==0) * wa0_rf * 
                       (Whole_S[ ,outname] - ga0_rf))/nrow(Target_S)
     DR0_nnL1 <- out0_nnL1 + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==0) * wa0_nnL1 * 
                       (Whole_S[ ,outname] - ga0_nnL1))/nrow(Target_S)
     DR0_nnL2 <- out0_nnL2 + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==0) * wa0_nnL2 * 
                       (Whole_S[ ,outname] - ga0_nnL2))/nrow(Target_S)
     DR0_nnL3 <- out0_nnL3 + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==0) * wa0_nnL3 * 
                       (Whole_S[ ,outname] - ga0_nnL3))/nrow(Target_S)
     DR0_svmL <- out0_svmL + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==0) * wa0_svmL * 
                       (Whole_S[ ,outname] - ga0_svmL))/nrow(Target_S)
     DR0_svm <- out0_svm + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==0) * wa0_svm * 
                       (Whole_S[ ,outname] - ga0_svm))/nrow(Target_S)
     DR0_gam <- out0_gam + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==0) * wa0_gam * 
                       (Whole_S[ ,outname] - ga0_gam))/nrow(Target_S)
     DR0_gbdt <- out0_gbdt + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==0) * wa0_gbdt * 
                       (Whole_S[ ,outname] - ga0_gbdt))/nrow(Target_S)
     DR0_SL <- out0_SL + sum((In_Trial_S==1) * (Whole_S[ ,trtname]==0) * wa0_SL * 
                       (Whole_S[ ,outname] - ga0_SL))/nrow(Target_S)
  }else{
     DR1 <- out1 + sum((In_Trial==1) * (Whole[ ,trtname]==1) * wa1 * 
                       (Whole[ ,outname] - ga1))/nrow(Target)
     DR1_rf <- out1_rf + sum((In_Trial==1) * (Whole[ ,trtname]==1) * wa1_rf * 
                       (Whole[ ,outname] - ga1_rf))/nrow(Target)
     DR1_nnL1 <- out1_nnL1 + sum((In_Trial==1) * (Whole[ ,trtname]==1) * wa1_nnL1 * 
                       (Whole[ ,outname] - ga1_nnL1))/nrow(Target)
     DR1_nnL2 <- out1_nnL2 + sum((In_Trial==1) * (Whole[ ,trtname]==1) * wa1_nnL2 * 
                       (Whole[ ,outname] - ga1_nnL2))/nrow(Target)
     DR1_nnL3 <- out1_nnL3 + sum((In_Trial==1) * (Whole[ ,trtname]==1) * wa1_nnL3 * 
                       (Whole[ ,outname] - ga1_nnL3))/nrow(Target)
     DR1_svmL <- out1_svmL + sum((In_Trial==1) * (Whole[ ,trtname]==1) * wa1_svmL * 
                       (Whole[ ,outname] - ga1_svmL))/nrow(Target)
     DR1_svm <- out1_svm + sum((In_Trial==1) * (Whole[ ,trtname]==1) * wa1_svm * 
                       (Whole[ ,outname] - ga1_svm))/nrow(Target)
     DR1_gam <- out1_gam + sum((In_Trial==1) * (Whole[ ,trtname]==1) * wa1_gam * 
                       (Whole[ ,outname] - ga1_gam))/nrow(Target)
     DR1_gbdt <- out1_gbdt + sum((In_Trial==1) * (Whole[ ,trtname]==1) * wa1_gbdt * 
                       (Whole[ ,outname] - ga1_gbdt))/nrow(Target)
     DR1_SL <- out1_SL + sum((In_Trial==1) * (Whole[ ,trtname]==1) * wa1_SL * 
                       (Whole[ ,outname] - ga1_SL))/nrow(Target)

     DR0 <- out0 + sum((In_Trial==1) * (Whole[ ,trtname]==0) * wa0 * 
                       (Whole[ ,outname] - ga0))/nrow(Target)
     DR0_rf <- out0_rf + sum((In_Trial==1) * (Whole[ ,trtname]==0) * wa0_rf * 
                       (Whole[ ,outname] - ga0_rf))/nrow(Target)
     DR0_nnL1 <- out0_nnL1 + sum((In_Trial==1) * (Whole[ ,trtname]==0) * wa0_nnL1 * 
                       (Whole[ ,outname] - ga0_nnL1))/nrow(Target)
     DR0_nnL2 <- out0_nnL2 + sum((In_Trial==1) * (Whole[ ,trtname]==0) * wa0_nnL2 * 
                       (Whole[ ,outname] - ga0_nnL2))/nrow(Target)
     DR0_nnL3 <- out0_nnL3 + sum((In_Trial==1) * (Whole[ ,trtname]==0) * wa0_nnL3 * 
                       (Whole[ ,outname] - ga0_nnL3))/nrow(Target)
     DR0_svmL <- out0_svmL + sum((In_Trial==1) * (Whole[ ,trtname]==0) * wa0_svmL * 
                       (Whole[ ,outname] - ga0_svmL))/nrow(Target)
     DR0_svm <- out0_svm + sum((In_Trial==1) * (Whole[ ,trtname]==0) * wa0_svm * 
                       (Whole[ ,outname] - ga0_svm))/nrow(Target)
     DR0_gam <- out0_gam + sum((In_Trial==1) * (Whole[ ,trtname]==0) * wa0_gam * 
                       (Whole[ ,outname] - ga0_gam))/nrow(Target)
     DR0_gbdt <- out0_gbdt + sum((In_Trial==1) * (Whole[ ,trtname]==0) * wa0_gbdt * 
                       (Whole[ ,outname] - ga0_gbdt))/nrow(Target)
     DR0_SL <- out0_SL + sum((In_Trial==1) * (Whole[ ,trtname]==0) * wa0_SL * 
                       (Whole[ ,outname] - ga0_SL))/nrow(Target)
  }

 return(c(DR1 - DR0, DR1_rf - DR0_rf, DR1_nnL1 - DR0_nnL1, DR1_nnL2 - DR0_nnL2, 
          DR1_nnL3 - DR0_nnL3, DR1_svmL - DR0_svmL, DR1_svm - DR0_svm, 
          DR1_gam - DR0_gam, DR1_gbdt - DR0_gbdt, DR1_SL - DR0_SL))
}
