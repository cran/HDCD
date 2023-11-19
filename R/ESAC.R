#' @useDynLib HDCD cESAC
#' @useDynLib HDCD cESAC_calibrate
#' @useDynLib HDCD cESAC_test
#' @useDynLib HDCD cESAC_test_calibrate
#' @importFrom Rdpack reprompt
#' @import stats


#' @title Efficient Sparsity Adaptive Change-point estimator
#' @description R wrapper for C function implementing the full ESAC algorithm \insertCite{@see @moen2023efficient}{HDCD}. 
#' @param X Matrix of observations, where each row contains a time series
#' @param alpha Parameter for generating seeded intervals
#' @param K Parameter for generating seeded intervals
#' @param empirical If \code{TRUE}, detection thresholds are based on Monte Carlo simulation using \code{\link{ESAC_calibrate}}
#' @param threshold_d Leading constant for \eqn{\lambda(t) \propto r(t)} for \eqn{t= p}. Only relevant when \code{thresholds=NULL}
#' @param threshold_s Leading constant for \eqn{\lambda(t) \propto r(t)} for \eqn{t\leq \sqrt{p\log n}}.  Only relevant when \code{thresholds=NULL}
#' @param threshold_d_test Leading constant for \eqn{\gamma(t) \propto r(t)} for \eqn{t=p}.  Only relevant when \code{empirical=FALSE} and \code{thresholds_test=NULL}
#' @param threshold_s_test Leading constant for \eqn{\gamma(t) \propto r(t)} for \eqn{t\leq \sqrt{p\log n}}. Only relevant when \code{empirical=FALSE} and \code{thresholds_test=NULL}
#' @param thresholds Vector of manually chosen values of \eqn{\lambda(t)} for \eqn{t \in \mathcal{T}}, decreasing in \eqn{t}
#' @param thresholds_test Vector of manually chosen values of \eqn{\gamma(t)} for \eqn{t \in \mathcal{T}}, decreasing in \eqn{t}
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @param tol If \code{empirical=TRUE}, \code{tol} is the false error probability tolerance
#' @param N If \code{empirical=TRUE}, \code{N} is the number of Monte Carlo samples used
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}
#' @param fast If \code{TRUE}, ESAC only tests for a change-point at the midpoint of each seeded interval
#' @param trim If \code{TRUE}, interval trimming is performed
#' @param NOT If \code{TRUE}, ESAC uses Narrowest-Over-Threshold selection of change-points
#' @param midpoint If \code{TRUE}, change-point positions are estimated by the mid-point of the seeded interval in which the penalized score is the largest
#' @return A list containing 
#'   \item{changepoints}{vector of estimated change-points}
#'   \item{changepointnumber}{number of changepoints}
#'   \item{CUSUMval}{the penalized score at the corresponding change-point in \code{changepoints}}
#'   \item{coordinates}{a matrix of zeros and ones indicating which time series are affected by a change in mean, with each row corresponding to the change-point in \code{changepoints}}
#'   \item{scales}{vector of estimated noise level for each series}
#'   \item{startpoints}{start point of the seeded interval detecting the corresponding change-point in \code{changepoints}}
#'   \item{endpoints}{end point of the seeded interval detecting the corresponding change-point in \code{changepoints}}
#'   \item{thresholds}{vector of values of \eqn{\lambda(t)} for \eqn{t \in \mathcal{T}} in decreasing order}
#'   \item{thresholds_test}{vector of values of \eqn{\gamma(t)} for \eqn{t \in \mathcal{T}} in decreasing order}
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
#' # Vanilla ESAC:
#' res = ESAC(X)
#' res$changepoints
#' 
#' # Manually setting leading constants for \lambda(t) and \gamma(t)
#' res = ESAC(X, 
#'            threshold_d = 2, threshold_s = 2, #leading constants for \lambda(t)
#'            threshold_d_test = 2, threshold_s_test = 2 #leading constants for \gamma(t)
#' )
#' res$changepoints #estimated change-point locations
#' 
#' # Empirical choice of thresholds:
#' res = ESAC(X, empirical = TRUE, N = 100, tol = 1/100)
#' res$changepoints
#' 
#' # Manual empirical choice of thresholds (equivalent to the above)
#' thresholds_emp = ESAC_calibrate(n,p, N=100, tol=1/100)
#' res = ESAC(X, thresholds_test = thresholds_emp[[1]])
#' res$changepoints
#' @references
#' \insertAllCited{}
#' @export
ESAC = function(X, threshold_d=1.5, threshold_s=1.0, alpha = 1.5, K = 5, debug =FALSE,
                empirical = FALSE, tol = 0.001, N = 1000, thresholds = NULL, thresholds_test = NULL,
                threshold_d_test = threshold_d, threshold_s_test = threshold_s, fast = FALSE,
                rescale_variance = TRUE,trim = FALSE,NOT=TRUE, midpoint = FALSE){
  p = dim(X)[1]
  n = dim(X)[2]

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

  
  max_s = min(sqrt(p*log(n)), p)
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*p*4*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log.p=TRUE))
  if(debug){
    print(ss)
  }
  
  ts = ss
 

  if(is.null(thresholds_test)){
    if(empirical){
      ttt = ESAC_calibrate(n=n,p=p, alpha=alpha, K=K, N=N, tol=tol , bonferroni = TRUE,
                           fast=fast, rescale_variance = rescale_variance, debug = debug)
      
      thresholds_test = ttt[[1]]
      
    }
    else{
      if(is.null(threshold_s_test)){
        threshold_s_test = threshold_s
      }
      if(is.null(threshold_d_test)){
        threshold_d_test = threshold_d
      }
      thresholds_test = nu_as[] 
      thresholds_test[2:length(ss)] = threshold_s_test*(ss[2:length(ss)]*log(exp(1)*p*log(n^4)/ss[2:length(ss)]^2)+ log(n^4))
      thresholds_test[1] = threshold_d_test* (sqrt(p*log(n^4)) + log(n^4))
    }
    
  }
  
  
  if(is.null(thresholds)){
    
    thresholds = nu_as[]
    
    thresholds[2:length(ss)] = threshold_s*(ss[2:length(ss)]*log(exp(1)*4*p*log(n)/ss[2:length(ss)]^2)+ log(n^4))
    thresholds[1] = threshold_d * (sqrt(p*log(n^4)) + log(n^4))
  }


  if(debug){
    print("thresholds:")
    print(thresholds)
    print("thresholds_test:")
    print(thresholds_test)
    print("as:")
    print(as)
  }
  
  res = .Call(cESAC, X[,],as.integer(n), as.integer(p), thresholds,thresholds_test,
              as.integer(lens),as.integer(length(lens)), as.integer(K), as,nu_as, 
              as.integer(length(as)), as.integer(0), as.integer(ts), as.integer(fast),as.integer(rescale_variance),
              as.integer(trim),as.integer(NOT),as.integer(midpoint),as.integer(debug))
  
  
  
  
  res$changepoints = as.integer(res$changepoints[1:res$changepointnum]+1)
  res$CUSUMval = res$CUSUMval[1:res$changepointnum]
  res$depth = as.integer(res$depth[1:res$changepointnum])
  res$coordinate = matrix(res$coordinate,nrow = p, ncol=n)
  srt_indices = as.integer(sort(res$changepoints, decreasing =FALSE, index.return=TRUE)$ix)
  res$changepoints = as.integer(res$changepoints[srt_indices])
  res$CUSUMval = res$CUSUMval[srt_indices]
  res$depth = as.integer(res$depth[srt_indices])
  res$coordinate = (res$coordinate[,srt_indices])
  res$startpoints = as.integer(res$startpoints[srt_indices]+1)
  res$endpoints = as.integer(res$endpoints[srt_indices]+1)
  res$maxaposes = as.integer(res$maxaposes[srt_indices])
  res$s = as.integer(ss[res$maxaposes+1])
  res$thresholds = thresholds
  res$thresholds_test = thresholds_test
  res$ts = ts
  
  if(res$changepointnum==0){
    res$changepoints = NA
    res$CUSUMval = NA
    res$depth = NA
    res$coordinate = NA
    res$startpoints = NA
    res$endpoints = NA
    res$maxaposes = NA
    res$s = NA
  }
  

  return(res)
}



#' @title ESAC single change-point test
#' @description R wrapper for C function testing for a single change-point using ESAC \insertCite{@see @moen2023efficient}{HDCD}.
#' @param X Matrix of observations, where each row contains a time series
#' @param empirical If \code{TRUE}, detection thresholds are based on Monte Carlo simulation using \code{\link{ESAC_test_calibrate}}
#' @param threshold_d Leading constant for \eqn{\gamma(t) \propto r(t)} for \eqn{t=p}.  Only relevant when \code{empirical=FALSE} and \code{thresholds=NULL}
#' @param threshold_s Leading constant for \eqn{\gamma(t) \propto r(t)} for \eqn{t\leq \sqrt{p\log n}}. Only relevant when \code{empirical=FALSE} and \code{thresholds=NULL}
#' @param thresholds Vector of manually chosen values of \eqn{\gamma(t)} for \eqn{t \in \mathcal{T}}, decreasing in \eqn{t}
#' @param fast If \code{TRUE}, ESAC only tests for a change-point at the midpoint of each seeded interval
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @param tol If \code{empirical=TRUE}, \code{tol} is the false error probability tolerance
#' @param N If \code{empirical=TRUE}, \code{N} is the number of Monte Carlo samples used
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}
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
#' # Vanilla ESAC:
#' resX = ESAC_test(X)
#' resX
#' resY = ESAC_test(Y)
#' resY
#' 
#' # Manually setting leading constants for \lambda(t) and \gamma(t)
#' resX = ESAC_test(X, 
#'                  threshold_d = 2, threshold_s = 2, #leading constants for \gamma(t)
#' )
#' resX 
#' resY = ESAC_test(Y, 
#'                  threshold_d = 2, threshold_s = 2, #leading constants for \gamma(t)
#' )
#' resY
#' 
#' # Empirical choice of thresholds:
#' resX = ESAC_test(X, empirical = TRUE, N = 100, tol = 1/100)
#' resX
#' resY = ESAC_test(Y, empirical = TRUE, N = 100, tol = 1/100)
#' resY
#' 
#' # Manual empirical choice of thresholds (equivalent to the above)
#' thresholds_test_emp = ESAC_test_calibrate(n,p, N=100, tol=1/100,bonferroni=TRUE)
#' resX = ESAC_test(X, thresholds = thresholds_test_emp[[1]])
#' resX
#' resY = ESAC_test(Y, thresholds = thresholds_test_emp[[1]])
#' resY
#' @references 
#' \insertAllCited{}
#' @export
ESAC_test = function(X, threshold_d=1.5, threshold_s=1.0,debug =FALSE,
                empirical = FALSE, thresholds = NULL, fast = FALSE,
                tol = 0.001, N = 1000, rescale_variance = TRUE){
  
  p = dim(X)[1]
  n = dim(X)[2]
  

  
  max_s = min(sqrt(p*log(n)), p)
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*p*4*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log.p=TRUE))
  if(debug){
    print(ss)
  }
  
  ts = ss
 

  
  
  if(is.null(thresholds)){
    if(empirical){
      ttt = ESAC_test_calibrate(n=n, p=p, N=N, tol=tol, fast = fast, rescale_variance = rescale_variance,
                                        debug =debug)
     
      thresholds = ttt[[1]]

        
    }
    else{
      thresholds = nu_as[] 
      thresholds[2:length(ss)] = threshold_s*(ss[2:length(ss)]*log(exp(1)*p*log(n^4)/ss[2:length(ss)]^2)+ log(n^4))
      thresholds[1] = threshold_d* (sqrt(p*log(n^4)) + log(n^4))
      
    }
    
  }
  

  if(debug){
    print("thresholds:")
    print(thresholds)
    print("as:")
    print(as)
  }

  
  
  res = .Call(cESAC_test, X[,], as.integer(n), as.integer(p), thresholds,
              as,nu_as, as.integer(length(as)),
              as.integer(0), as.integer(ts), as.integer(rescale_variance), 
              as.integer(fast), as.integer(debug))
              
            
  
  return(res)
}

#' @title Generates empirical penalty function \eqn{\gamma(t)} for single change-point testing using Monte Carlo simulation
#' @description R wrapper for C function choosing the penalty function \eqn{\gamma(t)} by Monte Carlo simulation, as described in Appendix B in \insertCite{moen2023efficient;textual}{HDCD}, for testing for a single change-point.
#' @param n Number of observations
#' @param p Number time series
#' @param bonferroni If \code{TRUE}, a Bonferroni correction applied and the empirical penalty function \eqn{\gamma(t)} is chosen by simulating leading constants of \eqn{r(t)} through Monte Carlo simulation.
#' @param tol False positive probability tolerance
#' @param N Number of Monte Carlo samples used
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}
#' @param fast If \code{TRUE}, ESAC only tests for a change-point at the midpoint of the interval \eqn{(0,\ldots,n]}
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @returns A list containing a vector of values of \eqn{\gamma(t)} for \eqn{t \in \mathcal{T}} decreasing (element #1), a vector of corresponding values of the threshold \eqn{a(t)} (element # 3), a vector of corresponding values of \eqn{\nu_{a(t)}}
#' @return A list containing 
#'   \item{without_partial}{ a vector of values of \eqn{\gamma(t)} for \eqn{t \in \mathcal{T}} decreasing in \eqn{t}}
#'   \item{with_partial}{same as \code{without_partial}}
#'   \item{as}{vector of threshold values \eqn{a(t)} for \eqn{t \in \mathcal{T}} decreasing in \eqn{t}}
#'   \item{nu_as}{vector of conditional expectations \eqn{\nu_{a(t)}} of a thresholded Gaussian, for \eqn{t \in \mathcal{T}} decreasing in \eqn{t}}
#' @examples
#' library(HDCD)
#' n = 50
#' p = 50
#' 
#' set.seed(100)
#' thresholds_emp = ESAC_test_calibrate(n,p, bonferroni=TRUE,N=100, tol=1/100)
#' set.seed(100)
#' thresholds_emp_without_bonferroni = ESAC_test_calibrate(n,p, bonferroni=FALSE,N=100, tol=1/100)
#' thresholds_emp[[1]] # vector of \gamma(t) for t = p,...,1
#' thresholds_emp_without_bonferroni[[1]] # vector of \gamma(t) for t = p,...,1
#' 
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' Y = matrix(rnorm(n*p), ncol = n, nrow=p)
#' 
#' # Adding a single sparse change-point to X (and not Y):
#' X[1:5, 26:n] = X[1:5, 26:n] +2
#' resX = ESAC_test(X, thresholds = thresholds_emp[[1]])
#' resX
#' resY = ESAC_test(Y,  thresholds = thresholds_emp[[1]])
#' resY
#' @references
#' \insertAllCited{}
#' @export
ESAC_test_calibrate = function(n, p, bonferroni=TRUE, N=1000, tol=1/1000, fast = FALSE, rescale_variance = TRUE,
                               debug =FALSE){
  
  max_s = min(sqrt(p*log(n)), p)
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*p*4*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log.p=TRUE))
  if(debug){
    print(ss)
  }
  
  ts = ss
 
  
  
  if(bonferroni){
    ESACtreshinfo =  r_func_values(n,p)
    h = ESACtreshinfo[[2]]
    v = length(ESACtreshinfo[[1]])
    if(h <= v && v >1){
      corr = 3
    }else if(v>1){
      corr = 2
    }else{
      corr = 1
    }
    tol = tol/corr
    N = corr*N
  }
  toln = max(round(N*tol),1)
  
  
  res = .Call(cESAC_test_calibrate, as.integer(n), as.integer(p), as.integer(N), as.integer(toln),
              as, nu_as, as.integer(length(as)), as.integer(0), as.integer(ts), as.integer(fast), 
              as.integer(rescale_variance), as.integer(debug))
  if(bonferroni){
    ESACtreshinfo =  r_func_values(n,p)
    ttt = ESACtreshinfo[[1]]
    h = ESACtreshinfo[[2]]
    v = length(ttt)
    
    if(h <=length(ESACtreshinfo[[1]]) && v>1){
      
      res[[1]] = c(res[[1]][1], max(res[[1]][2:(h-1)] / ttt[2:(h-1)])*ttt[2:(h-1)], max(res[[1]][h:v] / ttt[h:v])*ttt[h:v] )
      res[[2]] = c(res[[2]][1], max(res[[2]][2:(h-1)] / ttt[2:(h-1)])*ttt[2:(h-1)], max(res[[2]][h:v] / ttt[h:v])*ttt[h:v] )
    }else if(v>1){
      res[[1]] = c(res[[1]][1], max(res[[1]][2:v] / ttt[2:v])*ttt[2:v])
      res[[2]] = c(res[[2]][1], max(res[[2]][2:v] / ttt[2:v])*ttt[2:v])
    }else{
      res[[1]] = max(res[[1]] / ttt)*ttt
      res[[2]] = max(res[[2]] / ttt)*ttt
    }
  }
  return(res)
}


#' @title Generates empirical penalty function \eqn{\gamma(t)} for the ESAC algorithm using Monte Carlo simulation
#' @description R wrapper for C function choosing the penalty function \eqn{\gamma(t)} by Monte Carlo simulation, as described in Appendix B in \insertCite{moen2023efficient;textual}{HDCD}.
#' @param n Number of observations
#' @param p Number time series
#' @param alpha Parameter for generating seeded intervals
#' @param K Parameter for generating seeded intervals
#' @param bonferroni If \code{TRUE}, a Bonferroni correction applied and the empirical penalty function \eqn{\gamma(t)} is chosen by simulating leading constants of \eqn{r(t)} through Monte Carlo simulation.
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @param tol False error probability tolerance
#' @param N Number of Monte Carlo samples used
#' @param tdf If NULL, samples are drawn from a Gaussian distribution. Otherwise, they are drawn from a \eqn{t}
#' distribution with \code{tdf} degrees of freedom.
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}
#' @param fast If \code{TRUE}, ESAC only tests for a change-point at the midpoint of each seeded interval
#' @return A list containing 
#'   \item{without_partial}{ a vector of values of \eqn{\gamma(t)} for \eqn{t \in \mathcal{T}} decreasing in \eqn{t}}
#'   \item{with_partial}{same as \code{without_partial}}
#'   \item{as}{vector of threshold values \eqn{a(t)} for \eqn{t \in \mathcal{T}} decreasing in \eqn{t}}
#'   \item{nu_as}{vector of conditional expectations \eqn{\nu_{a(t)}} of a thresholded Gaussian, for \eqn{t \in \mathcal{T}} decreasing in \eqn{t}}#' 
#' @examples
#' library(HDCD)
#' n = 50
#' p = 50
#' 
#' set.seed(100)
#' thresholds_emp = ESAC_calibrate(n,p, N=100, tol=1/100)
#' set.seed(100)
#' thresholds_emp_without_bonferroni = ESAC_calibrate(n,p, N=100, tol=1/100,bonferroni=FALSE)
#' thresholds_emp[[1]] # vector of \gamma(t) for t = p,...,1
#' thresholds_emp_without_bonferroni[[1]] # vector of \gamma(t) for t = p,...,1
#' 
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' # Adding a single sparse change-point:
#' X[1:5, 26:n] = X[1:5, 26:n] +2
#' 
#' res = ESAC(X, thresholds_test = thresholds_emp[[1]])
#' res$changepoints
#' @references
#' \insertAllCited{}
#' @export
ESAC_calibrate = function(n,p, alpha = 1.5, K = 5, N=1000, tol=0.001, bonferroni = TRUE, fast = FALSE,rescale_variance = TRUE,
                          tdf = NULL,debug=FALSE){
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
  
  max_s = min(sqrt(p*log(n)), p)
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*p*4*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log.p=TRUE))
  if(debug){
    print(ss)
  }
  
  ts = ss

  
  if(is.null(tdf)){
    tdf =0
  }
  if(debug){
    print("as:")
    print(as)
  }
  
  
  if(bonferroni){
    ESACtreshinfo =  r_func_values(n,p)
    h = ESACtreshinfo[[2]]
    v = length(ESACtreshinfo[[1]])
    if(h <= v && v >1){
      corr = 3
    }else if(v>1){
      corr = 2
    }else{
      corr = 1
    }
    tol = tol/corr
    N = corr*N
  }
  toln = max(round(N*tol),1)
  
  res= .Call(cESAC_calibrate, as.integer(n), as.integer(p), as.integer(N),
             as.integer(toln), as.integer(lens), as.integer(length(lens)),
             as.integer(K), as.numeric(as), as.numeric(nu_as), as.integer(length(as)),
             as.integer(0), as.integer(ts), as.integer(fast),
             as.integer(rescale_variance),as.integer(tdf),as.integer(debug))
  
  if(bonferroni){
    ESACtreshinfo =  r_func_values(n,p)
    ttt = ESACtreshinfo[[1]]
    h = length(ttt) +1
    v = length(ttt)
    
    if(h <=length(ESACtreshinfo[[1]]) && v>1){
      
      res[[1]] = c(res[[1]][1], max(res[[1]][2:(h-1)] / ttt[2:(h-1)])*ttt[2:(h-1)], max(res[[1]][h:v] / ttt[h:v])*ttt[h:v] )
      res[[2]] = c(res[[2]][1], max(res[[2]][2:(h-1)] / ttt[2:(h-1)])*ttt[2:(h-1)], max(res[[2]][h:v] / ttt[h:v])*ttt[h:v] )
    }else if(v>1){
      res[[1]] = c(res[[1]][1], max(res[[1]][2:v] / ttt[2:v])*ttt[2:v])
      res[[2]] = c(res[[2]][1], max(res[[2]][2:v] / ttt[2:v])*ttt[2:v])
    }else{
      res[[1]] = max(res[[1]] / ttt)*ttt
      res[[2]] = max(res[[2]] / ttt)*ttt
    }
  }
  return(res)
}


r_func_values = function(n,p){
  
  max_s = min(sqrt(p*log(n)), p)
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*p*4*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log.p=TRUE))

  
  
  thresholds = nu_as[]
  
  thresholds[2:length(ss)] = (ss[2:length(ss)]*log(exp(1)*4*p*log(n)/ss[2:length(ss)]^2)+ log(n^4))
  thresholds[1] =  (sqrt(p*log(n^4)) + log(n^4))
  
  twologn = min(p, floor(log(n)/2))
  tt = 1
  count = 0
  h = length(thresholds) +1
  while(tt <=twologn){
    h = h-1
    tt = 2*tt
  }
  

  
  return(list(thresholds, h))
}



rescale.variance <- function(x){
  p <- dim(x)[1]
  n <- dim(x)[2]
  for (j in 1:p){
    scale <- mad(diff(x[j,]))/sqrt(2)
    x[j,] <- x[j,] / scale
  }
  return(x)
}

rescale.variance.scales <- function(x){
  p <- dim(x)[1]
  n <- dim(x)[2]
  scales = rep(NA, p)
  for (j in 1:p){
    scales[j] <- mad(diff(x[j,]))/sqrt(2)
  }
  return(scales)
}