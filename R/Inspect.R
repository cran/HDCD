#' @useDynLib HDCD cInspect
#' @useDynLib HDCD cInspect_calibrate
#' @useDynLib HDCD cInspect_single
#' @useDynLib HDCD cInspect_test_calibrate





#' @title Informative sparse projection for estimating change-points (Inspect)
#' @description R wrapper for C function implementing a Narrowest-Over-Threshold variant of Inspect \insertCite{wang_high_2018;textual}{HDCD} as specified in Appendix C in \insertCite{moen2023efficient;textual}{HDCD}. Note that the algorithm is only implemented for \eqn{\mathcal{S} = \mathcal{S}_2}, in the notation of \insertCite{moen2023efficient;textual}{HDCD}.
#' @param X Matrix of observations, where each row contains a time series
#' @param alpha Parameter for generating seeded intervals
#' @param K Parameter for generating seeded intervals
#' @param empirical If \code{TRUE}, the detection threshold \eqn{xi} is based on Monte Carlo simulation using \code{\link{Inspect_calibrate}}
#' @param lambda Manually specified value of \eqn{\lambda} (can be \code{NULL}, in which case \eqn{\lambda \gets \sqrt{\log(p\log n)/2}})
#' @param xi Manually specified value of \eqn{\xi} (can be NULL, in which case \eqn{\xi \gets 4\sqrt{\log(np)}})
#' @param tol If \code{empirical=TRUE}, \code{tol} is the false error probability tolerance
#' @param N If \code{empirical=TRUE}, \code{N} is the number of Monte Carlo samples used
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}
#' @param eps Threshold for declaring numerical convergence of the power method
#' @param maxiter Maximum number of iterations for the power method
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @return A list containing 
#'   \item{changepoints}{vector of estimated change-points}
#'   \item{changepointnumber}{number of changepoints}
#'   \item{CUSUMval}{vector with the sparse projected CUSUMs corresponding to \code{changepoints}}
#'   \item{coordinates}{a matrix of zeros and ones indicating which time series are affected by a change in mean, with each row corresponding to the change-point in \code{changepoints}}
#'   \item{scales}{vector of estimated noise level for each series}
#' @examples
#' library(HDCD)
#' n = 50
#' p = 50
#' set.seed(100)
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' # Adding a single sparse change-point:
#' X[1:5, 26:n] = X[1:5, 26:n] +1
#' 
#' # Vanilla Inspect:
#' res = Inspect(X)
#' res$changepoints
#' 
#' # Manually setting leading constants for \lambda(t) and \gamma(t)
#' res = Inspect(X, 
#'               lambda = sqrt(log(p*log(n))/2),
#'               xi = 4*sqrt(log(n*p))
#' )
#' res$changepoints #estimated change-point locations
#' 
#' # Empirical choice of thresholds:
#' res = Inspect(X, empirical=TRUE, N = 100, tol = 1/100)
#' res$changepoints
#' 
#' # Manual empirical choice of thresholds (equivalent to the above)
#' thresholds_emp = Inspect_calibrate(n,p, N=100, tol=1/100)
#' res = Inspect(X, xi = thresholds_emp$max_value)
#' res$changepoints
#' @references 
#' \insertAllCited{}
#' @export
Inspect = function(X, lambda = NULL , xi = NULL , alpha = 1.5, K = 5,eps=1e-10,
                     empirical=FALSE,maxiter=10000, N = 100, tol = 1/100, rescale_variance = TRUE, debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  if(is.null(lambda)){
    lambda = sqrt(log(p*log(n))/2)
  }
  if(is.null(xi)){
    xi = 4*sqrt(log(n*p))
  }





  lens = c(1)
  last = 1
  tmp = last
  while(alpha*last<n){
    tmp = last
    last = floor(alpha*last)
    if(last==tmp){
      last = last+1
    }
    lens= c(lens, last)
  }
  if(empirical){
    ttt = Inspect_calibrate(n=n, p=p, N=N, tol=tol,lambda = lambda , alpha = alpha, K = K,eps=eps,
                                       maxiter=maxiter,rescale_variance = rescale_variance,debug =debug)
    xi = ttt$max_value
  }
  res = .Call(cInspect, X[,],as.integer(n), as.integer(p), xi,
              as.integer(lens),as.integer(length(lens)), as.integer(K), eps, lambda, as.integer(maxiter),
              as.integer(rescale_variance),as.integer(debug))

  
  res$changepoints = as.integer(res$changepoints[1:res$changepointnumber]+1)
  res$CUSUMval = res$CUSUMval[1:res$changepointnumber]
  res$depth = as.integer(res$depth[1:res$changepointnumber])
  res$coordinate = matrix(res$coordinate,nrow = p, ncol=n)
  srt_indices = as.integer(sort(res$changepoints, decreasing =FALSE, index.return=TRUE)$ix)
  res$changepoints = as.integer(res$changepoints[srt_indices])
  res$CUSUMval = res$CUSUMval[srt_indices]
  res$depth = as.integer(res$depth[srt_indices])
  res$coordinate = res$coordinate[,srt_indices]
  
  if(res$changepointnum==0){
    res$changepoints = NA
    res$CUSUMval = NA
    res$depth = NA
    res$coordinate = NA
  }

  return(res)
}


#' @title Generates empirical detection threshold \eqn{\xi} using Monte Carlo simulation
#' @description R wrapper for C function choosing empirical detection threshold \eqn{\xi} for the Narrowest-Over-Threshold variant of Inspect \insertCite{@as specified in section 4.2 in @moen2023efficient}{HDCD} using Monte Carlo simulation.
#' @param n Number of observations
#' @param p Number time series
#' @param alpha Parameter for generating seeded intervals
#' @param K Parameter for generating seeded intervals
#' @param tol False positive probability tolerance
#' @param N Number of Monte Carlo samples used
#' @param lambda Manually specified value of \eqn{\lambda} (can be \code{NULL}, in which case \eqn{\lambda \gets \sqrt{\log(p\log n)/2}})
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}
#' @param eps Threshold for declaring numerical convergence of the power method
#' @param maxiter Maximum number of iterations for the power method
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @returns A list containing 
#'    \item{max_value}{the empirical threshold}
#' @examples 
#' library(HDCD)
#' n = 50
#' p = 50
#' 
#' set.seed(100)
#' thresholds_emp = Inspect_calibrate(n,p, N=100, tol=1/100)
#' thresholds_emp$max_value # xi
#' 
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' # Adding a single sparse change-point:
#' X[1:5, 26:n] = X[1:5, 26:n] +2
#' 
#' res = Inspect(X, xi = thresholds_emp$max_value)
#' res$changepoints
#' @references 
#' \insertAllCited{}
#' @export
Inspect_calibrate = function(n, p, N=100, tol=1/100,lambda = NULL , alpha = 1.5, K = 5,eps=1e-10,
                   maxiter=10000,rescale_variance = TRUE,debug =FALSE){

  if(is.null(lambda)){
    lambda = sqrt(log(p*log(n))/2)
  }
  
  lens = c(1)
  last = 1
  tmp = last
  while(alpha*last<n){
    tmp = last
    last = floor(alpha*last)
    if(last==tmp){
      last = last+1
    }
    lens= c(lens, last)
  }

  toln = max(round(N*tol),1)

  res = .Call(cInspect_calibrate, as.integer(n), as.integer(p), as.integer(N),
              as.integer(toln),as.integer(lens), as.integer(length(lens)), as.integer(K),
              as.numeric(lambda), as.numeric(eps), as.integer(maxiter),
              as.integer(rescale_variance),as.integer(debug))
  
  return(res)
}


#' @title Inspect single change-point test
#' @description R wrapper for C function testing for a single change-point using Inspect \insertCite{wang_high_2018;textual}{HDCD}.
#' @param X Matrix of observations, where each row contains a time series
#' @param empirical If \code{TRUE}, the detection threshold \eqn{xi} is based on Monte Carlo simulation using \code{\link{Inspect_test_calibrate}}
#' @param lambda Manually specified value of \eqn{\lambda} (can be \code{NULL}, in which case \eqn{\lambda \gets \sqrt{\log(p\log n)/2}})
#' @param xi Manually specified value of \eqn{\xi} (can be NULL, in which case \eqn{\xi \gets 4\sqrt{\log(np)}})
#' @param tol If \code{empirical=TRUE}, \code{tol} is the false error probability tolerance
#' @param N If \code{empirical=TRUE}, \code{N} is the number of Monte Carlo samples used
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}
#' @param eps Threshold for declaring numerical convergence of the power method
#' @param maxiter Maximum number of iterations for the power method
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @returns 1 if a change-point is detected, 0 otherwise
#' @examples 
#' library(HDCD)
#' n = 50
#' p = 50
#' 
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' Y = matrix(rnorm(n*p), ncol = n, nrow=p)
#' 
#' # Adding a single sparse change-point to X (and not Y):
#' X[1:5, 26:n] = X[1:5, 26:n] +1
#' 
#' # Vanilla Inspect:
#' resX = Inspect_test(X)
#' resX
#' resY = Inspect_test(Y)
#' resY
#' 
#' # Manually setting \lambda and \xi:
#' resX = Inspect_test(X, 
#'                     lambda = sqrt(log(p*log(n))/2),
#'                     xi = 4*sqrt(log(n*p))
#' )
#' resX 
#' resY = Inspect_test(Y, 
#'                     lambda = sqrt(log(p*log(n))/2),
#'                     xi = 4*sqrt(log(n*p))
#' )
#' resY
#' 
#' # Empirical choice of thresholds:
#' resX = Inspect_test(X, empirical = TRUE, N = 100, tol = 1/100)
#' resX
#' resY = Inspect_test(Y, empirical = TRUE, N = 100, tol = 1/100)
#' resY
#' 
#' # Manual empirical choice of thresholds (equivalent to the above)
#' thresholds_test_emp = Inspect_test_calibrate(n,p, N=100, tol=1/100)
#' resX = Inspect_test(X, xi = thresholds_test_emp$max_value)
#' resX
#' resY = Inspect_test(Y, xi = thresholds_test_emp$max_value)
#' resY
#' @references 
#' \insertAllCited{}
#' @export
Inspect_test = function(X, lambda = NULL , xi = NULL, eps=1e-10, empirical=FALSE, N = 100,tol = 1/100,
                          maxiter=10000,rescale_variance = TRUE,debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  
  if(is.null(lambda)){
    lambda = sqrt(log(p*log(n))/2)
  }
  if(is.null(xi)){
    xi = 4*sqrt(log(n*p))
  }
  
  if(empirical){
    ttt = Inspect_test_calibrate = function(n=n, p=p, N=N, tol=tol,lambda = lambda, eps=eps,
                                            maxiter=maxiter,rescale_variance =rescale_variance,debug =debug)
    xi = ttt$max_value
  }
  if(rescale_variance){
    X = rescale_variance(X)$X
  }
  res = .Call(cInspect_single, X,as.integer(n), as.integer(p), 0,
              eps, lambda, as.integer(maxiter), as.integer(debug))
  
  
  
  if(res$cusumval > xi){
    return(1)
  }else{return(0)}
  
}


#' @title Generates empirical detection threshold \eqn{\xi} for single change-point testing using Monte Carlo simulation
#' @description R wrapper for C function choosing the empirical detection threshold \eqn{\xi} for Inspect \insertCite{wang_high_2018;textual}{HDCD} for single change-point testing using Monte Carlo simulation.
#' @param n Number of observations
#' @param p Number time series
#' @param tol False positive probability tolerance
#' @param N Number of Monte Carlo samples used
#' @param lambda Manually specified value of \eqn{\lambda} (can be \code{NULL}, in which case \eqn{\lambda \gets \sqrt{\log(p\log n)/2}})
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}
#' @param eps Threshold for declaring numerical convergence of the power method
#' @param maxiter Maximum number of iterations for the power method
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @returns A list containing 
#'    \item{max_value}{the empirical threshold}
#' @examples 
#' library(HDCD)
#' n = 50
#' p = 50
#' 
#' set.seed(100)
#' thresholds_emp = Inspect_test_calibrate(n,p,N=100, tol=1/100)
#' thresholds_emp
#' 
#' 
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' Y = matrix(rnorm(n*p), ncol = n, nrow=p)
#' 
#' # Adding a single sparse change-point to X (and not Y):
#' X[1:5, 26:n] = X[1:5, 26:n] +2
#' resX = Inspect_test(X, xi = thresholds_emp$max_value)
#' resX
#' resY = Inspect_test(Y,  xi = thresholds_emp$max_value)
#' resY
#' @references 
#' \insertAllCited{}
#' @export
Inspect_test_calibrate = function(n, p, N=100, tol=1/100,lambda = NULL, eps=1e-10,
                             maxiter=10000,rescale_variance = TRUE,debug =FALSE){
  
  if(is.null(lambda)){
    lambda = sqrt(log(p*log(n))/2)
  }
  
  toln = max(round(N*tol),1)
  
  res = .Call(cInspect_test_calibrate, as.integer(n), as.integer(p), as.integer(N),
              as.integer(toln), as.numeric(lambda), as.numeric(eps), as.integer(maxiter),
              as.integer(rescale_variance),as.integer(debug))

  return(res)
}

