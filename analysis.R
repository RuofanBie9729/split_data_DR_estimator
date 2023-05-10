N <- c(500, 800, 1000, 1500, 2000, 2500)
DR_res <- NULL
for(n in N){
	filename <- paste0("hat_p_misspecified_n", n, ".rds")
	test <- readRDS(filename)
	Res <- sapply(test, function(x){x[[1]]})
      #Res <- na.omit(t(Res))
	ATE <- apply(Res, 2, mean)
      DR_bias <- ATE[c(4, 6:8)]-ATE[1]
	DR_res <- rbind(DR_res, cbind(rep(n,4), DR_bias, c("DR", "DR_rf", "DRSP", "DRSP_rf")))
}

DR_res <- as.data.frame(DR_res)
colnames(DR_res) <- c("N", "bias", "estimator")
library(ggplot2)
ggplot(DR_res, aes(x=as.numeric(N), y=as.numeric(bias), color=estimator, group=estimator)) +
	geom_line()+
	geom_point()+
	geom_hline(yintercept=0, linetype="dashed", color="red")

test <- readRDS("toyexample_n1500_nn.rds")
length(test)
Res <- sapply(test, function(x){x[[1]]})
boot_res <- sapply(test, function(x){x[[2]]})
Est <- apply(Res, 1, mean)
Est[2:7]-Est[1]
boot_res <- rbind(Res[1,], boot_res)
boot_ci <- function(x){
  true = x[1]
  x = x[2:length(x)]
  est_boot <- function(i){
    boot <- x[((i-1)*100+1):(i*100)]
    boot_lower <- quantile(sort(boot), 0.025)
    boot_upper <- quantile(sort(boot), 0.975)
    return((boot_lower <= true) && (boot_upper >= true))
  }
  coverage <- sapply(1:6, est_boot)
  return(coverage)
}
coverage_res <- apply(boot_res, 2, boot_ci)
apply(coverage_res, 1, mean)


ResTab <- NULL
for(i in 1:10){
  filename <- paste0("Sim", i, ".rds")
  rawRes <- readRDS(filename)
  Res <- sapply(rawRes, function(x) x)
  Est <- apply(Res, 1, mean)
  ResTab <- rbind(ResTab, c(Est[1], Est[2:21]-Est[1]))
}
colnames(ResTab) <- c("ATE", "glm", "rf", "nnL1", "nnL2", "nnL3", "svmL", "svm", 
                      "gam", "gbdt", "SL", "glm_sp", "rf_sp", "nnL1_sp", 
                      "nnL2_sp", "nnL3_sp", "svmL_sp", "svm_sp", "gam_sp", 
                      "gbdt_sp", "SL_sp")
write.csv(ResTab, "SimRes.csv")