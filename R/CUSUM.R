#' @useDynLib HDCD CUSUM_R
#' @useDynLib HDCD single_CUSUM_R



#' @title CUSUM transformation of a matrix
#' @description R wrapper for C function computing the CUSUM transformation of a matrix over an interval \eqn{(s,e]}. For compatibility with C indexing, the user should subtract \eqn{1} from both \eqn{s} and \eqn{e} when supplying the arguments to the function. If start and stop are not supplied, the CUSUM is computed over the full data, so \eqn{(s,e] = (0,n]}. In this case, \code{CUSUM} returns the same result as \code{cusum.transform} in the package \code{InspectChangepoint} \insertCite{inspectpackage}{HDCD}. 
#' @param X Matrix of observations, where each row contains a time series
#' @param start Starting point of interval over which the CUSUM should be computed, subtracted by one
#' @param stop Ending point of interval over which the CUSUM should be computed, subtracted by one
#' @return A matrix of CUSUM values. The \eqn{(i,j)}-th element corresponds to the CUSUM transformation of the \eqn{i}-th row of \eqn{X}, computed over the interval \eqn{(\code{start}+1,\code{end}+1]} and evaluated at position \eqn{\code{start}+1+j}, i.e. 
#' \eqn{\sqrt{\frac{e-v}{(e-s)(v-s)}}\sum_{t=s+1}^v X_{i,t} - \sqrt{\frac{v-s}{(e-s)(e-v)}}\sum_{t=v+1}^e X_{i,t}}, 
#' where \eqn{s = (\code{start}+1)}, \eqn{e = (\code{stop}+1)} and \eqn{v = \code{start}+1+j}.
#' @examples
#' n = 10
#' p = 10
#' set.seed(101)
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' # CUSUM over the full data (s,e] = (0,n]
#' X_cusum = CUSUM(X)
#' 
#' # CUSUM over (s,e] = (3,9]:
#' s = 3
#' e = 9
#' X_cusum = CUSUM(X, start = s-1, stop = e-1)
#' @references 
#' \insertAllCited{}
#' @export
CUSUM = function(X, start=NULL, stop=NULL){
  p = dim(X)[1]
  n = dim(X)[2]
  
  if(is.null(start)){
    start = -1
  }
  if(is.null(stop)){
    stop = n-1
  }
  
  res = .Call(CUSUM_R, X, as.integer(start), as.integer(stop), as.integer(p), as.integer(n))
  res = matrix(res[1:(p*(stop-start-1))], nrow = p, ncol = stop-start-1)
  return(res)
}


#' @title CUSUM transformation of matrix at a specific position
#' @description R wrapper for C function computing the CUSUM transformation of matrix over an interval \eqn{(s,e]} evaluated at a specific position. For compatibility with C indexing, the user should subtract \eqn{1} from \eqn{s}, \eqn{e} and \eqn{v} when supplying the arguments to the function. If start and stop are not supplied, the CUSUM is computed over the full data, so \eqn{(s,e] = (0,n]}. 
#' @param X Matrix of observations, where each row contains a time series
#' @param start Starting point of interval over which the CUSUM should be computed, subtracted by one
#' @param stop Ending point of interval over which the CUSUM should be computed, subtracted by one
#' @param pos Position at which the CUSUM should be evaluated, subtracted by one
#' @return A vector of CUSUM values, each corresponding to a row of the input matrix. The \eqn{i}-th element corresponds to the CUSUM transformation of the \eqn{i}-th row of \eqn{X}, computed over the interval \eqn{(\code{start}+1,\code{end}+1]} and evaluated at position \eqn{\code{pos}}, i.e. 
#' \eqn{\sqrt{\frac{e-v}{(e-s)(v-s)}}\sum_{t=s+1}^v X_{i,t} - \sqrt{\frac{v-s}{(e-s)(e-v)}}\sum_{t=v+1}^e X_{i,t}}, 
#' where \eqn{s = (\code{start}+1)}, \eqn{e = (\code{stop}+1)} and \eqn{v = \code{pos}+1}.
#' @examples
#' n = 10
#' p = 10
#' set.seed(101)
#' X = matrix(rnorm(n*p), ncol = n, nrow=p)
#' # CUSUM over the full data (s,e] = (0,n] evaluated at position v=4
#' position = 4
#' X_cusum_single = single_CUSUM(X,pos = position-1)
#' X_cusum_single
#' 
#' # verifying that this corresponds to the 4-th row of output of CUSUM():
#' X_cusum = CUSUM(X)
#' X_cusum[,4]
#' @export
single_CUSUM = function(X, start=NULL, stop=NULL, pos){
  p = dim(X)[1]
  n = dim(X)[2]
  if(is.null(start)){
    start = -1
  }
  if(is.null(stop)){
    stop = n-1
  }
  res = .Call(single_CUSUM_R, X, as.integer(start), as.integer(stop), as.integer(p), as.integer(pos),as.integer(n))
  
  
  return(res)
}