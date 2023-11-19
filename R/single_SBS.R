#' @useDynLib HDCD cSBS_single
#' @useDynLib HDCD cSBS_single_calibrate





#' @title Sparsified Binary Segmentation for single change-point estimation
#' @description R wrapper for C function for single change-point estimation using Sparsified Binary Segmentation \insertCite{cho_multiple-change-point_2015;textual}{HDCD}.
#' @param X Matrix of observations, where each row contains a time series
#' @param threshold Manually specified value of the threshold \eqn{\pi_T}
#' @param empirical If \code{TRUE}, the threshold is based on Monte Carlo simulation
#' @param tol If \code{empirical=TRUE}, \code{tol} is the false error probability tolerance
#' @param N If \code{empirical=TRUE}, \code{N} is the number of Monte Carlo samples used
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @return A list containing 
#'   \item{pos}{estimated change-point location}
#'   \item{maxval}{maximum thresholded and aggregated CUSUM at the estimated change-point position}
#' @examples 
#' # Single SBS
#' library(HDCD)
#' n = 50
#' p = 50
#' set.seed(101)
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' # Adding a single sparse change-point:
#' X[1:5, 26:n] = X[1:5, 26:n] +1
#' 
#' res = single_SBS(X,threshold=7,rescale_variance=TRUE)
#' res$pos
#' 
#' # Choose threhsold by Monte Carlo:
#' res = single_SBS(X,empirical=TRUE,rescale_variance=TRUE)
#' res$pos
#' @references 
#' \insertAllCited{}
#' @export
single_SBS = function(X, threshold=NULL, rescale_variance = TRUE,empirical =FALSE,N=100,tol=1/100,debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  
  if(n<2){
    return(NULL)
  }
  
  if(empirical){
    ttt = single_SBS_calibrate(n=n,p=p,N=N,tol=tol,rescale_variance = rescale_variance,debug=debug)
    threshold = sqrt(ttt)
  }else if(is.null(threshold)){
    return(NULL)
  }
  
  
  res = .Call(cSBS_single, X,as.integer(n), as.integer(p), as.numeric(threshold),
              as.integer(rescale_variance),as.integer(debug))
  
  res$pos = res$pos+1

  
  return(res)
}


#' @title Generates threshold \eqn{\pi_T} for Sparsified Binary Segmentation for single change-point detection
#' @description R wrapper for function choosing empirical threshold \eqn{\pi_T} using Monte Carlo simulation for single change-point Sparsified Binary Segmentation. More specifically, the function returns the empirical upper tol quantile of CUSUMs over \eqn{p} time series, each of length \eqn{n}, based on \eqn{N} number of runs.
#' @param n Number of observations
#' @param p Number time series
#' @param tol False positive probability tolerance
#' @param N Number of Monte Carlo samples used
#' @param rescale_variance If TRUE, each row of the data is rescaled by a MAD estimate
#' @param debug If TRUE, diagnostic prints are provided during execution
#' @returns Threshold
#' @examples 
#' library(HDCD)
#' n = 50
#' p = 50
#' set.seed(101)
#' 
#' # Simulate threshold
#' pi_T_squared = single_SBS_calibrate(n=n,p=p,N=100, tol=1/100, rescale_variance = TRUE)
#' pi_T_squared
#' 
#' 
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' # Adding a single sparse change-point:
#' X[1:5, 26:n] = X[1:5, 26:n] +1
#' 
#' # Run SBS
#' res = single_SBS(X,threshold=sqrt(pi_T_squared),rescale_variance=TRUE)
#' res$pos
#' @export
single_SBS_calibrate = function(n,p,N=100, tol=1/100, rescale_variance = TRUE,debug =FALSE){

  
  if(n<2){
    return(NULL)
  }
  toln = max(round(N*tol),1)
  
  res = .Call(cSBS_single_calibrate, as.integer(n), as.integer(p), as.integer(N),as.integer(toln),
          as.integer(rescale_variance),as.integer(debug))
  
  
  
  return(res)
}