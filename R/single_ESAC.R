#' @useDynLib HDCD cESAC_single



#' @title Efficient Sparsity Adaptive Change-point estimator for a single change-point
#' @description R wrapper for C function implementing ESAC for single change-point estimation, as described in section 3.1 in \insertCite{moen2023efficient;textual}{HDCD}
#' @param X Matrix of observations, where each row contains a time series
#' @param threshold_d Leading constant for \eqn{\lambda(t) \propto r(t)} for \eqn{t= p}
#' @param threshold_s Leading constant for \eqn{\lambda(t) \propto r(t)} for \eqn{t\leq \sqrt{p\log n}}. 
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}
#' @return A list containing 
#'   \item{pos}{estimated change-point location}
#'   \item{s}{the value of \eqn{t \in \mathcal{T}} at which the sparsity specific score is maximized}
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
#' res = single_ESAC(X,rescale_variance=TRUE)
#' res$pos
#' 
#' # Manually setting the leading constants for \lambda(t):
#' # here \lambda(t) = 2 (\sqrt{p \log(n^4)}  + \log (n^4)) for t=p
#' # and             = 2 (t \log (ep\log n^4 / t^2) + \log(n^4))
#' res = single_ESAC(X, threshold_d = 2, threshold_s = 2)
#' res$pos
#' @references
#' \insertAllCited{}
#' @export
single_ESAC = function(X, threshold_d=1.5, threshold_s=1, rescale_variance = FALSE,debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  
  if(n<2){
    return(NULL)
  }
  if(rescale_variance){
    X = rescale_variance(X)$X
  }
  max_s = min(p, sqrt(p*log(n)))
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*4*p*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log.p=TRUE))
  thresholds = nu_as[]
  
  thresholds[2:length(ss)] = threshold_s*(ss[2:length(ss)]*log(exp(1)*p*log(n^4)/ss[2:length(ss)]^2)+ log(n^4))
  
  thresholds[1] = threshold_d * (sqrt(p*log(n^4)) + log(n^4))
  
  res = .Call(cESAC_single, X,as.integer(n), as.integer(p), thresholds,
             as,nu_as, as.integer(length(as)), as.integer(debug))
  
  res$pos = res$pos+1
  res$s = ss[res$apos+1]
  
  
  return(res)
}




