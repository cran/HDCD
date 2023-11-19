#' @useDynLib HDCD cPilliat_test
#' @useDynLib HDCD cPilliat_test_calibrate
#' @useDynLib HDCD sort_test
#' @useDynLib HDCD partialsum_test




#' @title Pilliat single change-point test
#' @description R wrapper function testing for a single change-point using the three test statistics in the multiple change point detection algorithm of \insertCite{pilliat_optimal_2022;textual}{HDCD}. See also Appendix E in \insertCite{moen2023efficient;textual}{HDCD}.
#' @param X Matrix of observations, where each row contains a time series
#' @param empirical If \code{TRUE}, detection thresholds are based on Monte Carlo simulation
#' @param threshold_d_const Leading constant for the analytical detection threshold for the dense statistic 
#' @param threshold_partial_const Leading constant for the analytical detection threshold for the partial sum statistic 
#' @param threshold_bj_const Leading constant for \eqn{p_0} when computing the detection threshold for the Berk-Jones statistic
#' @param threshold_dense Manually specified value of detection threshold for the dense statistic
#' @param thresholds_partial Vector of manually specified detection thresholds for the partial sum statistic, for sparsities/partial sums \eqn{t=1,2,4,\ldots,2^{\lfloor\log_2(p)\rfloor}}
#' @param thresholds_bj Vector of manually specified detection thresholds for the Berk-Jones statistic, order corresponding to \eqn{x=1,2,\ldots,x_0}
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @param tol If \code{empirical=TRUE}, \code{tol} is the false error probability tolerance
#' @param N If \code{empirical=TRUE}, \code{N} is the number of Monte Carlo samples used
#' @param rescale_variance If \code{TRUE}, each row of the data is re-scaled by a MAD estimate (see \code{\link{rescale_variance}})
#' @param fast If \code{TRUE}, only the mid-point of \eqn{(0,\ldots,n]} is tested for a change-point. Otherwise a test is performed at each candidate change-point poisition
#' @returns 1 if a change-point is detected, 0 otherwise
#' @examples 
#' library(HDCD)
#' n = 200
#' p = 200
#' 
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' Y = matrix(rnorm(n*p), ncol = n, nrow=p)
#' 
#' # Adding a single sparse change-point to X (and not Y):
#' X[1:5, 100:200] = X[1:5, 100:200] +1
#' 
#' # Vanilla Pilliat test:
#' resX = Pilliat_test(X)
#' resX
#' resY = Pilliat_test(Y)
#' resY
#' 
#' # Manually setting leading constants for the theoretical thresholds for the three 
#' # test statistics used
#' resX = Pilliat_test(X, 
#'                     threshold_d_const=4, 
#'                     threshold_bj_const=6, 
#'                     threshold_partial_const=4
#' )
#' resX 
#' resY = Pilliat_test(Y, 
#'                     threshold_d_const=4, 
#'                     threshold_bj_const=6, 
#'                     threshold_partial_const=4
#' )
#' resY
#' 
#' # Empirical choice of thresholds:
#' resX = Pilliat_test(X, empirical = TRUE, N = 100, tol = 1/100)
#' resX
#' resY = Pilliat_test(Y, empirical = TRUE, N = 100, tol = 1/100)
#' resY
#' 
#' # Manual empirical choice of thresholds (equivalent to the above)
#' thresholds_test_emp = Pilliat_test_calibrate(n,p, N=100, tol=1/100,bonferroni=TRUE)
#' resX = Pilliat_test(X, 
#'                     threshold_dense=thresholds_test_emp$threshold_dense, 
#'                     thresholds_bj = thresholds_test_emp$thresholds_bj, 
#'                     thresholds_partial = thresholds_test_emp$thresholds_partial
#' )
#' resX
#' resY = Pilliat_test(Y, 
#'                     threshold_dense=thresholds_test_emp$threshold_dense, 
#'                     thresholds_bj = thresholds_test_emp$thresholds_bj, 
#'                     thresholds_partial = thresholds_test_emp$thresholds_partial
#' )
#' resY
#' @references
#' \insertAllCited{}
#' @export
Pilliat_test = function(X, empirical =FALSE, N = 100, tol = 0.05,
                                    thresholds_partial=NULL, threshold_dense = NULL,
                                    thresholds_bj = NULL,
                                    threshold_d_const=4, threshold_bj_const=6, threshold_partial_const=4,
                                    rescale_variance = TRUE,fast = FALSE,
                                    debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  
  log2p = floor(log(p,base=2))+1
  
  xx = 1
  maxx = 1
  while(TRUE){
    logp0 = log(threshold_bj_const) - 2*log(pi) - 2*log(xx) - 2*log(n)
    qq = qbinom(p = logp0, size = p, prob = exp(log(2) + pnorm(xx,lower.tail = FALSE,log.p = TRUE)),
                lower.tail=FALSE,log.p = TRUE)

    if(qq==0){
      maxx = xx-1
      break
    }
    xx = xx+1
  }
  
  if(empirical){
    if(is.null(thresholds_partial) || is.null(threshold_dense) || is.null(threshold_dense)){
      res = Pilliat_test_calibrate(n,p, N=N, tol=tol,threshold_bj_const=threshold_bj_const,debug=debug)
      thresholds_partial = res$thresholds_partial
      thresholds_bj = res$thresholds_bj
      threshold_dense = res$threshold_dense
    }
  }else{
    
    # dense
    threshold_dense = threshold_d_const * (sqrt(p*log(2*n^2)) + log(2*n^2))
    
    #partial norm
    thresholds_partial = c()
    
    t = 1
    while(TRUE){
      thresholds_partial = c(thresholds_partial, t*log(2*exp(1)*p/t) + log(2*n^2))
      t = 2*t
      if(t>=p){
        break
      }
    }
    thresholds_partial = thresholds_partial*threshold_partial_const
    
    # berk jones
    
    thresholds_bj = c()
    xx = 1
    maxx = 1
    while(TRUE){
      logp0 = log(threshold_bj_const) - 2*log(pi) - 2*log(xx) - 2*log(n)
      qq = qbinom(p = logp0, size = p, prob = exp(log(2) + pnorm(xx,lower.tail = FALSE,log.p = TRUE)),
                  lower.tail=FALSE,log.p = TRUE)
      if(debug){
        print(xx)
        print(exp(log(2) + pnorm(xx,lower.tail = FALSE,log.p = TRUE)))
        print(qq)
      }
      if(qq==0){
        maxx = xx-1
        break
      }
      thresholds_bj = c(thresholds_bj,qq)
      xx = xx+1
    }
  }
  
  if(debug){
    print(threshold_dense)
    print(thresholds_partial)
    print(thresholds_bj)
  }
  

  res = .Call(cPilliat_test, as.numeric(X),as.integer(n),as.integer(p),
              as.numeric(thresholds_partial), as.numeric(threshold_dense), 
              as.integer(thresholds_bj),as.integer(maxx), as.integer(rescale_variance),
              as.integer(fast),as.integer(debug))

  return(res)
}


#' @title Generates detection thresholds for the Pilliat algorithm for testing for a single change-point using Monte Carlo simulation
#' @description R wrapper for function choosing detection thresholds for the Dense, Partial sum and Berk-Jones statistics in the multiple change-point detection algorithm of \insertCite{pilliat_optimal_2022;textual}{HDCD} for single change-point testing using Monte Carlo simulation. When \code{Bonferroni==TRUE}, the detection thresholds are chosen by simulating the leading constant in the theoretical detection thresholds given in \insertCite{pilliat_optimal_2022;textual}{HDCD}, similarly as described in Appendix B in \insertCite{moen2023efficient;textual}{HDCD} for ESAC. When \code{Bonferroni==TRUE}, the thresholds for the Berk-Jones statistic are theoretical and not chosen by Monte Carlo simulation.
#' @param n Number of observations
#' @param p Number time series
#' @param bonferroni If \code{TRUE}, a Bonferroni correction applied and the detection thresholds for each statistic is chosen by simulating the leading constant in the theoretical detection thresholds
#' @param threshold_bj_const Leading constant for \eqn{p_0} for the Berk-Jones statistic
#' @param debug If \code{TRUE}, diagnostic prints are provided during execution
#' @param tol False error probability tolerance
#' @param N Number of Monte Carlo samples used
#' @param rescale_variance If \code{TRUE}, each row of the data is rescaled by a MAD estimate
#' @param fast If \code{FALSE}, a change-point test is applied to each candidate change-point position in each interval. If FALSE, only the mid-point of each interval is considered
#' @return A list containing 
#'   \item{thresholds_partial}{vector of thresholds for the Partial Sum statistic (respectively for \eqn{t=1,2,4,\ldots,2^{\lfloor\log_2(p)\rfloor}} number of terms in the partial sum)}
#'   \item{threshold_dense}{threshold for the dense statistic}
#'   \item{thresholds_bj}{vector of thresholds for the Berk-Jones static (respectively for \eqn{x=1,2,\ldots,x_0})}
#' @examples 
#' library(HDCD)
#' n = 50
#' p = 50
#' 
#' set.seed(100)
#' thresholds_test_emp = Pilliat_test_calibrate(n,p, bonferroni=TRUE,N=100, tol=1/100)
#' set.seed(100)
#' thresholds_test_emp_without_bonferroni = Pilliat_test_calibrate(n,p, 
#'                                          bonferroni=FALSE,N=100, tol=1/100)
#' thresholds_test_emp # thresholds with bonferroni correction
#' thresholds_test_emp_without_bonferroni # thresholds without bonferroni correction
#' 
#' # Generating data
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' Y = matrix(rnorm(n*p), ncol = n, nrow=p)
#' 
#' # Adding a single sparse change-point to X (and not Y):
#' X[1:5, 25:50] = X[1:5, 25:50] +2
#' resX = Pilliat_test(X, 
#'                     threshold_dense=thresholds_test_emp$threshold_dense, 
#'                     thresholds_bj = thresholds_test_emp$thresholds_bj, 
#'                     thresholds_partial = thresholds_test_emp$thresholds_partial
#' )
#' resX
#' resY = Pilliat_test(Y, 
#'                     threshold_dense=thresholds_test_emp$threshold_dense, 
#'                     thresholds_bj = thresholds_test_emp$thresholds_bj, 
#'                     thresholds_partial = thresholds_test_emp$thresholds_partial
#' )
#' resY
#' @export
Pilliat_test_calibrate = function(n,p, N=100, tol=1/100,threshold_bj_const=6, bonferroni=TRUE,
                                              rescale_variance = TRUE,fast = FALSE,debug=FALSE){

  log2p = floor(log(p,base=2))+1
  if(bonferroni){
    tol = tol/3
    N = N*3
  }
  thresholds_bj = c()
  xx = 1
  maxx = 1
  delta = tol 
  
  while(TRUE){
    logp0 = log(6*delta) - 2*log(pi) - 2*log(xx) - 2*log(n)
    qq = qbinom(p = logp0, size = p, prob = exp(log(2) + pnorm(xx,lower.tail = FALSE,log.p = TRUE)),
                lower.tail=FALSE,log.p = TRUE)
    if(qq==0){
      maxx = xx-1
      break
    }
    thresholds_bj = c(thresholds_bj,qq)
    xx = xx+1
  }
  
  toln = max(round(N*tol),1)

  
  res= .Call(cPilliat_test_calibrate, as.integer(n), as.integer(p), as.integer(N), 
             as.integer(toln), as.integer(maxx), as.integer(log2p),as.integer(rescale_variance),
             as.integer(fast),as.integer(debug))
  
  if(bonferroni){
    pilliatthreshinfo = Pilliat_thresholds(n,p,tol)
    
    res$thresholds_bj = ceiling(max(res$thresholds_bj / thresholds_bj) * thresholds_bj)
    res$thresholds_partial = max(res$thresholds_partial / pilliatthreshinfo$thresholds_partial) * pilliatthreshinfo$thresholds_partial
  }
  return(res)
}





