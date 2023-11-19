#' @useDynLib HDCD cInspect_single



#' @title Inspect for single change-point estimation
#' @description R wrapper for C function for single change-point estimation using Inspect \insertCite{wang_high_2018}{HDCD}. Note that the algorithm is only implemented for \eqn{\mathcal{S} = \mathcal{S}_2}, in the notation of \insertCite{wang_high_2018;textual}{HDCD}.
#' @param X Matrix of observations, where each row contains a time series
#' @param lambda Manually specified value of \eqn{\lambda} (can be \code{NULL}, in which case \eqn{\lambda \gets \sqrt{\log(p\log n)/2}})
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}
#' @param eps Threshold for declaring numerical convergence of the power method
#' @param maxiter Maximum number of iterations for the power method
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @return A list containing 
#'   \item{pos}{estimated change-point location}
#'   \item{CUSUMval}{projected CUSUM value at the estimated change-point position}
#' @examples
#' library(HDCD)
#' n = 500
#' p = 500
#' set.seed(101)
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' # Adding a single sparse change-point:
#' X[1:5, 201:500] = X[1:5, 201:500] +1
#' 
#' res = single_Inspect(X,rescale_variance=TRUE)
#' res$pos
#' 
#' # Manually setting the value of \lambda:
#' res = single_Inspect(X, lambda = 2*sqrt(log(p*log(n))/2))
#' res$pos
#' @references 
#' \insertAllCited{}
#' @export
single_Inspect = function(X, lambda = sqrt(log(p*log(n))/2), eps=1e-10, rescale_variance=FALSE,
                  maxiter=10000,debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  
  if(rescale_variance){
    rr = rescale_variance(X)
    X = rr$X
  }

  res = .Call(cInspect_single, X,as.integer(n), as.integer(p), 0,
              eps, lambda, as.integer(maxiter), as.integer(debug))
  
  
  
  res$pos = res$pos+1
  
  return(res)
}