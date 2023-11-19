#' @useDynLib HDCD rescale_variance_R


#' @title Re-scales each row of matrix by its MAD estimate
#' @description R wrapper for C function computing the (rescaled) median absolute difference in differences for each row of the input matrix. The rescaling factor is set to 1.05 (corresponding to the Normal distribution). Each row of the input matrix then re-scaled by the corresponding noise estimate. 
#' @param X A \eqn{p \times n} matrix
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @return A list containing 
#'   \item{X}{the input matrix, variance re-scaled and flattened}
#'   \item{scales}{vector of MAD estimates of the noise level of each row of the input matrix}
#' @examples 
#' library(HDCD)
#' n = 200
#' p = 500
#' set.seed(101)
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' 
#' ret = rescale_variance(X)
#' ret$X #rescaled matrix
#' ret$scales #estimated noise level for each time series (each row)
#' 
#' # Note that the rescaled matrix is in (column wise) vector form. To transform it back to a matrix,
#' # do the following:
#' rescaled_X = matrix(ret$X, nrow = p, ncol=n)
#' @export
rescale_variance = function(X,debug=FALSE){
  p = dim(X)[1]
  n = dim(X)[2]

  res = .Call(rescale_variance_R, as.numeric(X[,]),as.integer(n), as.integer(p), 
              as.integer(debug))
  
  
  return(res)
}
