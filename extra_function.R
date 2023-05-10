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