#' @title Adaptive shrinkage of correlation matrix from data
#'
#' @description This function performs adaptive shrinkage of sample correlation matrix obtained from data
#'
#' @param data The samples by features data matrix
#' @return The adaptive shrinkage correlation and covariance matrices, the shrink intensity and the Ledoit-Wolf adjusted adaptive shrunk versions of correlation and covraince matrices.
#'
#' @keywords adaptive shrinkage, correlation
#' @export


ash_cov_master <- function(data){
  nsamples <- dim(data)[1];
  nfeat <- dim(data)[2];
  colmean <- colMeans(data);
  covmat <- cov(data);
  cormat <- cov2cor(covmat);
  w <- array(0, c(nsamples, nfeat, nfeat));
  for(n in 1:nsamples){
    w[n,,] <- (data[n,] - colmean)%*% t(data[n,] - colmean);
  }
  var_hat_s <- (nsamples/(nsamples-1)^2) * apply(w, c(2,3), var);
  sum_var_hat_s <- sum(var_hat_s[row(var_hat_s)!=col(var_hat_s)])
  square_cor <- covmat^2;
  sum_s_square <- sum(square_cor[row(square_cor)!=col(square_cor)]);
  shrink_intensity <- sum_var_hat_s/sum_s_square;
  ash_cormat <- ash_cor(cormat, nsamples);

  if(shrink_intensity < 0){
    ash_cor_ledoit_wolf <- cormat;
  }
  if(shrink_intensity > 1){
    ash_cor_ledoit_wolf <- ash_cormat
  }
  else{
    ash_cor_ledoit_wolf <- shrink_intensity*ash_cormat + (1- shrink_intensity)*cormat;

  }
  ash_cov_ledoit_wolf <- diag(sqrt(diag(covmat)))%*%ash_cor_ledoit_wolf%*%diag(sqrt(diag(covmat)))
  ll <- list("ash_cov_ledoit_wolf"=ash_cov_ledoit_wolf,
             "ash_cor_ledoit_wolf"=ash_cor_ledoit_wolf,
             "shrink_intensity"=shrink_intensity,
             "ash_cormat"=ash_cormat,
             "sample_cormat"=cormat)
  return(ll)
}

