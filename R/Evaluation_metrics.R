#' @import mclust

#' @title Hausdorff distance between two sets
#' @description Computes the Hausdorff distance between two sets represented as vectors \code{v1} and \code{v2}. If \code{v1 == NULL} and \code{v2 != NULL}, then the largest distance between an element of \code{v1} and the set \eqn{\{1,n\}} is returned, and vice versa. If both vectors are \code{NULL}, \code{0} is returned. 
#' @param v1 Vector representing set 1
#' @param v2 Vector representing set 2
#' @param n Sample size (only relevant when either \code{v1} or \code{v2} is \code{NULL})
#' @returns The Hausdorff distance between \code{v1} and \code{v2}
#' @examples
#' library(HDCD)
#' n = 400
#' true_changepoints = c(50,100)
#' est_changepoints = c(51,110)
#' hausdorff(true_changepoints, est_changepoints,n)
#' hausdorff(true_changepoints, NULL,n)
#' hausdorff(NULL, est_changepoints,n)
#' hausdorff(NULL,NULL)
#' @export
hausdorff = function(v1, v2,n){
  if(is.null(v1)){
    if(!is.null(v2)){
      return(max(c(max(v2), max(n-v2))))
    }
    else{
      return(0)
    }
  }
  if(is.null(v2)){
    if(!is.null(v1)){
      return(max(c(max(v1), max(n-v1))))
    }
    else{
      return(0)
    }
  }

  max = 0
  for (v  in v1) {
    tmp = min(abs(v-v2))
    if(tmp>max){
      max=tmp
    }
  }
  for (v  in v2) {
    tmp = min(abs(v-v1))
    if(tmp>max){
      max=tmp
    }
  }
  return(max)
}


#' @title ARI
#' @description Computes the Adjusted Rand Index (ARI) of a vector of estimated change-points.
#' @param etas Vector of true change-points
#' @param eta_hats Vector of estimated change-points
#' @param n Sample size 
#' @returns The ARI
#' @examples
#' library(HDCD)
#' n = 400
#' true_changepoints = c(50,100)
#' est_changepoints = c(51,110)
#' ARI(true_changepoints, est_changepoints,n)
#' @export
ARI = function(etas, eta_hats, n){
  
  # create label for each data point
  seg_etas = rep(1, n)
  counter = 1
  if(!is.null(etas)){
    for (i in 1:length(etas)) {
      counter = counter+1
      seg_etas[(etas[i]+1):n] = counter
    }
    
  }
  
  seg_eta_hats = rep(1, n)
  counter = 1
  if(!is.null(eta_hats)){
    for (i in 1:length(eta_hats)) {
      counter = counter+1
      seg_eta_hats[(eta_hats[i]+1):n] = counter
    }
    
  }
  return(mclust::adjustedRandIndex(seg_etas, seg_eta_hats))
  
}

hausdorff_avg = function(v1, v2){
  dists = rep(NA, length(v1))
  for (i in 1:length(v1)) {
    v = v1[i]
    dists[i] = min(abs(v-v2))
    
  }
  
  return(mean(dists))
}