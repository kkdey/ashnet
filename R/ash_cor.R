#' @title Adaptive shrinkage of a correlation matrix
#'
#' @description This function performs adaptive shrinkage of a correlation matrix
#'
#' @param cormat The estimated sample correlation matrix
#' @param nsamples The number of samples over which the correlation matrix is estimated.
#' @param image if TRUE, plots an image of the shrunk and non-shrunk correlation matrices
#'
#' @return Returns a shrunk version of the sample correlation matrix (the output is also a correlation matrix)
#'
#' @keywords adaptive shrinkage, correlation
#' @export


ash_cor <- function(cormat, nsamples, image=FALSE)
{
  cor_table <- melt(cormat);
  cor_table_non_diag <- cor_table[which(cor_table[,1] !=cor_table[,2]),];

  cor_table_non_diag.val <- cor_table_non_diag[,3];
  cor_table_non_diag.val[which(cor_table_non_diag.val==1)]=(1- 1e-7);
  #cor_table_non_diag.val[which(cor_table_non_diag.val==0)]=(1e-7);

  cor_transform_mean_vec=0.5*log((1+cor_table_non_diag.val)/(1-cor_table_non_diag.val))
  cor_transform_sd_vec=rep(sqrt(1/(nsamples-3)), dim(cor_table_non_diag)[1]);
  options(warn=-1)
  fit=ashr::ash(cor_transform_mean_vec,cor_transform_sd_vec,mixcompdist="normal");
  ash_cor_vec=(exp(2*fit$PosteriorMean)-1)/(exp(2*fit$PosteriorMean)+1);

  newdata.table <- cor_table_non_diag;
  newdata.table[,3] <- ash_cor_vec;
  new_mat <- reshape2::dcast(newdata.table, X1~X2, value.var = "value")[,-1];
  new_mat[is.na(new_mat)]=1;
  pd_completion <- Matrix::nearPD(as.matrix(new_mat), conv.tol=1e-06);
  new_mat <- sweep(pd_completion$mat,diag(pd_completion$mat), MARGIN=1,"/")
  if(image) {
    image(cormat)
    image(new_mat)
  }
  return(new_mat)
}

