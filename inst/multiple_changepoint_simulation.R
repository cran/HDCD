# This code runs the simulations performed in Section 4.2 of Moen et al. (2024),	arXiv:2306.04702
# Note that the simulation study for the computational cost is in another file.
# Please note that the code for SUBSET must be downloaded from Github (https://github.com/SOTickle/SUBSET), 
#while the code for the method of Kaul et al (2021, Electronic Journal of Statistics) 
# must be obtained from the first author of that paper. 
# If SUBSET or the method of Kaul et al are not available,
# please set subset_included = FALSE and kaul_included = FALSE below.

## To recreate Table 3 in Section 4.2, uncomment lines 75 and 76

## Packages

library(InspectChangepoint)
library(hdbinseg)
library(doSNOW)
library(HDCD)
library(foreach)

# if running the code from Kaul, uncomment these if necessary:
#install.packages("MASS"); install.packages("rmutil")
#install.packages("parallel"); install.packages("doSNOW");
#install.packages("foreach");install.packages("InspectChangepoint")
#install.packages("pracma");

# source SUBSET_normal.R from https://github.com/SOTickle/SUBSET
subset_included = TRUE
subset_path = "... fill in/SUBSET ... "
kaul_included =TRUE
kaul_path = "... fill in ... "

# source code for the method of Kaul and Michailidis (2024, Statistica Sinica).
# Contact Abhishek <abhishek dot kaul at wsu.edu> for the source code
if(kaul_included){
  setwd(kaul_path)
  #unpack libraries
  source("unpack_libraries.R")
  source("K20/K20_bseg.R")
  source("alg1.R")
}

if(subset_included){
  source(sprintf("%s/SUBSET_normal.R", subset_path))
}




## Saving options
# maindir is the directory in which the results are stored, if save = TRUE
maindir = "... fill in ... "
dateandtime = gsub(" ", "--",as.character(Sys.time()))
dateandtime = gsub(":", ".", dateandtime)
savedir = file.path(maindir, dateandtime)
save = TRUE

if(save){
  # create subdirectory in maindir
  dir.create(savedir, showWarnings = FALSE)
  savedir = file.path(maindir, sprintf("%s/multi",dateandtime))
  dir.create(savedir, showWarnings = FALSE)

}


# Simulation options
N = 1000 #MC simulations
num_cores = 6 # cores to run in parallell
sparse_const = 4 # leading constant for sparse signal strength
dense_const = 4 # leading constant for dense signal strength
set.seed(1996)
rescale_variance = TRUE
Ncal = 1000
tol = 0.01#error tolerance for MC choice of thresholds


# select values of n and p:
ns = c(15,20)
ps = c(10,20)
#ns = c(100,200) # uncomment for having the same ns as in the article
#ps = c(100,1000,5000) #uncomment for having the same ps as in the article





## we now run calibrations, e.g. choosing thresholds/tuning parameters via MC simulations
## for the methods. these are stored in maindir/calibrates. if you have already 
## computed the mc thresholds, specify the path for the calibrates and set 
## loadcalibrates = TRUE below

loadcalibrates = FALSE
if(loadcalibrates){
  calibrates = readRDS(file=sprintf("%s/calibrates.Rda", maindir))
}else{
  cl <- makeCluster(num_cores,type="SOCK")
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = ((length(ns)*length(ps))), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  calibrates = foreach(z = 0:((length(ns)*length(ps))-1),.options.snow = opts) %dopar% {
  #for(z in 0:((length(ns)*length(ps))-1)){
    library(HDCD)
    library(InspectChangepoint)
    set.seed(300*z)
    rez = list()
    nind = floor(z/length(ps))+1
    pind = z%%length(ps)+1
    cc = ESAC_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K=4, alpha = 2, fast = TRUE,
                        rescale_variance = rescale_variance, debug=FALSE)
    cc4 = ESAC_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K=5, alpha = 1.5, fast = FALSE,
                         rescale_variance = rescale_variance, debug=FALSE)
    # cc2 = HDCD_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,debug=FALSE)
    lambda = sqrt(log(ps[pind]*log(ns[nind]))/2)
    cc2 = Inspect_calibrate(n = ns[nind], p = ps[pind], N=Ncal, tol=tol,lambda = lambda , alpha = 2, K = 4,eps=1e-10,
                            maxiter=10000,rescale_variance = rescale_variance, debug =FALSE)

    cc3 = Pilliat_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K = 4,
                            rescale_variance = rescale_variance, debug=FALSE)

    if(subset_included){
      source(sprintf("%s/main.R", subset_path))



      empirical_penalties = rep(NA, Ncal)
      if(rescale_variance){
        for (i in 1:Ncal) {
          mynulldata <- matrix(rnorm(ns[nind]*ps[pind],0,1),nrow=ps[pind],ncol=ns[nind],byrow=FALSE) # 5 variates with 1000 time points, no change
          empirical_penalties[i] = wbs_penaltyfinder(InspectChangepoint::rescale.variance(mynulldata), SUBSET.normal_penalty, 100)
        }
      }else{
        for (i in 1:Ncal) {
          mynulldata <- matrix(rnorm(ns[nind]*ps[pind],0,1),nrow=ps[pind],ncol=ns[nind],byrow=FALSE) # 5 variates with 1000 time points, no change
          empirical_penalties[i] = wbs_penaltyfinder(mynulldata, SUBSET.normal_penalty, 100)
        }
      }

      empirical_penalties = sort(empirical_penalties,decreasing = TRUE)
    }


    rez[[1]] = cc[[1]] #ESAC threshoolds
    rez[[2]] = cc[[2]] #ESAC thresholds
    rez[[3]] = cc2[[1]] #Inspect xi
    rez[[4]] = cc3[[1]] # pilliat
    rez[[5]] = cc3[[2]] # pilliat
    rez[[6]] = cc3[[3]] # pilliat
    rez[[7]] = ns[nind]
    rez[[8]] = ps[pind]
    if(subset_included){
      rez[[9]] = empirical_penalties[max(1,round(tol*Ncal))] # subset penalty

    }else{
      rez[[9]] = NULL
    }
    rez[[10]] = cc4[[1]] #ESAC threshoolds
    rez[[11]] = cc4[[2]] #ESAC thresholds
    rez
    #calibrates[[z+1]] = rez

  }
  close(pb)
  stopCluster(cl)
  saveRDS(calibrates, file=sprintf("%s/calibrates.Rda", savedir))
}


# here we create a function for generating synthetic data. There are seven different scenarios.
numconfig = 7
config = function(i,n,p){
  mus = matrix(0, nrow=p, ncol=n)
  etas = c()
  sparsities = c()
  sparsity = c()

  if(i==1){
    # no chgppts
  }
  else if(i==2){
    # 2 dense
    sparsity="Dense"
    etas = sort(sample(1:(n-1), 2))
    sparsities = sample(ceiling(sqrt(p*log(n))):p,length(etas), replace=TRUE)

    for (j in 1:length(etas)) {
      Delta = 1
      if(j==1){
        Delta = min(etas[1], etas[2]-etas[1])
      }else if(j<length(etas)){
        Delta = min(etas[j] - etas[j-1], etas[j+1] - etas[j])
      }else{
        Delta = min(etas[j]-etas[j-1], n -etas[j])
      }
      phi = dense_const/sqrt(Delta)*(p*log(n))^(1/4)
      k = sparsities[j]
      coords = sample(1:p, sparsities[j])
      diff = sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] + matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }

  }
  else if(i==3){
    # 2 sparse
    sparsity = "Sparse"
    etas = sort(sample(1:(n-1), 2))
    sparsities = sample(1:floor(sqrt(p*log(n))),length(etas), replace=TRUE )


    for (j in 1:length(etas)) {
      Delta = 1
      if(j==1){
        Delta = min(etas[1], etas[2]-etas[1])
      }else if(j<length(etas)){
        Delta = min(etas[j] - etas[j-1], etas[j+1] - etas[j])
      }else{
        Delta = min(etas[j]-etas[j-1], n -etas[j])
      }
      k = sparsities[j]
      phi = sparse_const/sqrt(Delta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
      coords = sample(1:p, sparsities[j])
      diff = sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }

  }
  else if(i==4){
    # 2 mixed
    sparsity="Mixed"
    etas = sort(sample(1:(n-1), 2))
    sparsities = sample(c(1,2), length(etas), replace=TRUE)
    for (j in 1:length(etas)) {
      Delta = 1
      if(j==1){
        Delta = min(etas[1], etas[2]-etas[1])
      }else if(j<length(etas)){
        Delta = min(etas[j] - etas[j-1], etas[j+1] - etas[j])
      }else{
        Delta = min(etas[j]-etas[j-1], n -etas[j])
      }
      if(sparsities[j]==1){
        # sparse
        sparsities[j] = sample(1:floor(sqrt(p*log(n))),1)
        k = sparsities[j]
        phi = sparse_const/sqrt(Delta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
      }else{
        #dense
        sparsities[j] = sample(ceiling(sqrt(p*log(n))):p,1)
        k = sparsities[j]
        phi = dense_const/sqrt(Delta)*(p*log(n))^(1/4)
      }
      coords = sample(1:p, sparsities[j])
      #diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }

  }
  else if(i==5){
    # 5 dense
    sparsity="Dense"
    etas = sort(sample(1:(n-1), 5))
    sparsities = sample(ceiling(sqrt(p*log(n))):p,length(etas), replace=TRUE)

    for (j in 1:length(etas)) {
      Delta = 1
      if(j==1){
        Delta = min(etas[1], etas[2]-etas[1])
      }else if(j<length(etas)){
        Delta = min(etas[j] - etas[j-1], etas[j+1] - etas[j])
      }else{
        Delta = min(etas[j]-etas[j-1], n -etas[j])
      }
      phi = dense_const/sqrt(Delta)*(p*log(n))^(1/4)
      k = sparsities[j]
      coords = sample(1:p, sparsities[j])
      diff = sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] + matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }
  }
  else if(i==6){
    # 5 sparse
    sparsity = "Sparse"
    etas = sort(sample(1:(n-1), 5))
    sparsities = sample(1:floor(sqrt(p*log(n))),length(etas), replace=TRUE )


    for (j in 1:length(etas)) {
      Delta = 1
      if(j==1){
        Delta = min(etas[1], etas[2]-etas[1])
      }else if(j<length(etas)){
        Delta = min(etas[j] - etas[j-1], etas[j+1] - etas[j])
      }else{
        Delta = min(etas[j]-etas[j-1], n -etas[j])
      }
      k = sparsities[j]
      phi = sparse_const/sqrt(Delta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
      coords = sample(1:p, sparsities[j])
      diff = sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }
  }
  else if(i==7){
    # 5 mixed
    sparsity="Mixed"
    etas = sort(sample(1:(n-1), 5))
    sparsities = sample(c(1,2), length(etas), replace=TRUE)
    for (j in 1:length(etas)) {
      Delta = 1
      if(j==1){
        Delta = min(etas[1], etas[2]-etas[1])
      }else if(j<length(etas)){
        Delta = min(etas[j] - etas[j-1], etas[j+1] - etas[j])
      }else{
        Delta = min(etas[j]-etas[j-1], n -etas[j])
      }
      if(sparsities[j]==1){
        # sparse
        sparsities[j] = sample(1:floor(sqrt(p*log(n))),1)
        k = sparsities[j]
        phi = sparse_const/sqrt(Delta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
      }else{
        #dense
        sparsities[j] = sample(ceiling(sqrt(p*log(n))):p,1)
        k = sparsities[j]
        phi = dense_const/sqrt(Delta)*(p*log(n))^(1/4)
      }
      coords = sample(1:p, sparsities[j])
      #diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }
  }
  return(list(etas, sparsities,mus,sparsity))
}

###### Simulation ########

# run simulation in parallell:
cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result = foreach(z = 1:N,.options.snow = opts) %dopar% {
  rez = list()
  set.seed(2*z)
  library(HDCD)
  library(hdbinseg)
  counter = 1
  if(kaul_included){
    setwd(kaul_path)
    #unpack libraries
    source("unpack_libraries.R")
    source("K20/K20_bseg.R")
    source("alg1.R")
  }
  if(subset_included){
    source(sprintf("%s/main.R", subset_path))
  }

  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]

      for (y in 1:numconfig) {
        noise = matrix(rnorm(n*p), nrow=p, ncol=n)
        conf = config(y, n, p)
        etas = conf[[1]]
        sparsities = conf[[2]]
        mus = conf[[3]]

        X = noise + mus


        rezi = list()
        rezi[["i"]] = i
        rezi[["j"]] = j
        rezi[["y"]] = y

       
        lambda = sqrt(log(p*log(n))/2)

          
        # running "fast" version of ESAC. Not included in the article:
        a = proc.time()
        res = ESAC(X[,], 1.5,1, empirical=TRUE,alpha = 2, K = 4, thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[2]],  fast =TRUE,
                    rescale_variance = rescale_variance, debug= FALSE)
        b=proc.time()
        if(res$changepointnumber==0){
          res$changepoints=NULL
        }
        rezi[["esac_fast_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["esac_fast_J"]]= res$changepointnumber
        rezi[["esac_fast_chgpts"]]= res$changepoints

        rezi[["esac_fast_hausd"]] = hausdorff(res$changepoints, etas,n)
        rezi[["esac_fast_J_error"]] = length(res$changepoints) - length(etas)
        if(res$changepointnumber==0){
          rezi[["esac_fast_ari"]] = ARI(etas, c(), n)
        }else{
          rezi[["esac_fast_ari"]] = ARI(etas, res$changepoints, n)
        }

        
        # Run standard version of ESAC:
        a = proc.time()
        res = ESAC (X[,], 1.5,1, empirical=TRUE,alpha = 1.5, K = 5, thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[10]], fast =FALSE,
                    rescale_variance = rescale_variance, debug= FALSE)
        b=proc.time()
        if(res$changepointnumber==0){
          res$changepoints=NULL
        }
        rezi[["esac_long_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["esac_long_J"]]= res$changepointnumber
        rezi[["esac_long_chgpts"]]= res$changepoints
        rezi[["esac_long_hausd"]] = hausdorff(res$changepoints, etas,n)
        rezi[["esac_long_J_error"]] = length(res$changepoints) - length(etas)
        if(res$changepointnumber==0){
          rezi[["esac_long_ari"]] = ARI(etas, c(), n)
        }else{
          rezi[["esac_long_ari"]] = ARI(etas, res$changepoints, n)
        }


        a = proc.time()
        res  = Pilliat(X[,], K = 4, alpha = 2, empirical = TRUE, threshold_dense = (calibrates[[j+(i-1)*length(ps)]])[[5]],
                       thresholds_partial = (calibrates[[j+(i-1)*length(ps)]])[[4]], thresholds_bj = (calibrates[[j+(i-1)*length(ps)]])[[6]],
                       rescale_variance = rescale_variance, debug =FALSE)
        b=proc.time()
        if(res$number_of_changepoints==0){
          res$changepoints=NULL
        }

        rezi[["pilliat_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["pilliat_J"]] = res$number_of_changepoints
        rezi[["pilliat_chgpts"]]= res$changepoints
        if(res$number_of_changepoints==0){
          rezi[["pilliat_hausd"]] = hausdorff(NULL, etas,n)
        }else{
          rezi[["pilliat_hausd"]] = hausdorff(res$changepoints, etas,n)
        }
        rezi[["pilliat_J_error"]] = length(res$changepoints) - length(etas)
        if(res$number_of_changepoints==0){
          rezi[["pilliat_ari"]] = ARI(etas, c(), n)
        }else{
          rezi[["pilliat_ari"]] = ARI(etas, res$changepoints, n)
        }
        #rezi[["pilliat_ari"]] = ARI(etas, res$changepoints, n)



        # run inspect:
        a = proc.time()
        res = Inspect(X[,], xi=(calibrates[[j+(i-1)*length(ps)]])[[3]], alpha = 2, K = 4,eps=1e-10,
                      lambda = lambda, maxiter=10000,
                      rescale_variance = rescale_variance, debug=FALSE)
        b=proc.time()
        if(res$changepointnumber==0){
          res$changepoints=NULL
        }
        rezi[["inspect_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["inspect_J"]]= res$changepointnumber
        rezi[["inspect_chgpts"]]= res$changepoints
        if(res$changepointnumber==0){
          rezi[["inspect_hausd"]] = hausdorff(NULL, etas,n)
        }else{
          rezi[["inspect_hausd"]] = hausdorff(res$changepoints, etas,n)
        }
        rezi[["inspect_J_error"]] = length(res$changepoints) - length(etas)
        if(res$changepointnumber==0){
          rezi[["inspect_ari"]] = ARI(etas, c(), n)
        }else{
          rezi[["inspect_ari"]] = ARI(etas, res$changepoints, n)
        }



        # subset
        if(subset_included){
          a = proc.time()
          if(rescale_variance){
            res = change_main(InspectChangepoint::rescale.variance(X), SUBSET.normal, 100,penalties= (calibrates[[j+(i-1)*length(ps)]])[[9]])
          }else{
            res = change_main(X, SUBSET.normal,100, penalties=(calibrates[[j+(i-1)*length(ps)]])[[9]])
          }
          b=proc.time()
          rezi[["subset_time"]] = (b-a)[1]+(b-a)[2]
          rezi[["subset_J"]]= length(res[[2]])
          rezi[["subset_chgpts"]]= sort(res[[2]])
          rezi[["subset_hausd"]] = hausdorff(sort(res[[2]]), etas,n)
          rezi[["subset_J_error"]] = length(res[[2]]) - length(etas)
          if(is.null(res[[2]])){
            rezi[["subset_ari"]] = ARI(etas, c(), n)
          }else{
            rezi[["subset_ari"]] = ARI(etas, sort(res[[2]]), n)
          }

        }

        # kaul (2024):
        if(kaul_included){

          if(rezi[["inspect_J"]] == 0){
            # inspect detected no changepoints
            rezi[["kaul24_time"]] = rezi[["inspect_time"]]
            rezi[["kaul24_J"]] = 0
            rezi[["kaul24_hausd"]] = rezi[["inspect_hausd"]]
            rezi[["kaul24_J_error"]] = rezi[["inspect_J_error"]]
            rezi[["kaul24_ari"]] = rezi[["inspect_ari"]]
            rezi[["kaul24_failed"]] = 0
          }else{
            # use inspect output as preliminary estimate:
            inspectpreliminary = function(x){
              aa = list()
              aa$t = rezi[["inspect_chgpts"]]
              return(aa)
            }
            a = proc.time()
            source("alg1.R")
            hh = tryCatch({
            if(rescale_variance){
              res = alg1(t(InspectChangepoint::rescale.variance(X[,])), prel.cp = inspectpreliminary)$tilde.t
            }else{
              res =alg1(t(X[,]), prel.cp = inspectpreliminary)$tilde.t
            }
            b=proc.time()

            rezi[["kaul24_time"]] = (b-a)[3] + rezi[["inspect_time"]]
            rezi[["kaul24_J"]]= length(res)
            rezi[["kaul24_chgpts"]]= sort(res)
            rezi[["kaul24_hausd"]] = hausdorff(sort(res), etas,n)
            rezi[["kaul24_J_error"]] = length(res) - length(etas)
            rezi[["kaul24_ari"]] = ARI(etas, sort(res), n)
            rezi[["kaul24_failed"]] = 0
            }, error = function(e) {
              return(-1)
            })
            if(hh==-1){
              # if alg1.R failed
              rezi[["kaul24_chgpts"]]= rezi[["inspect_chgpts"]]
              rezi[["kaul24_time"]] = rezi[["inspect_time"]]
              rezi[["kaul24_J"]] = rezi[["inspect_J"]]
              rezi[["kaul24_hausd"]] = rezi[["inspect_hausd"]]
              rezi[["kaul24_J_error"]] = rezi[["inspect_J_error"]]
              rezi[["kaul24_ari"]] = rezi[["inspect_ari"]]
              rezi[["kaul24_failed"]] = 1
            }
          }

        }

        # SBS
        a = proc.time()
        res = sbs.alg(X[,], cp.type = 1, thr = NULL, trim = NULL, height = NULL,
                      temporal = TRUE, scales = NULL, diag = FALSE, B = 100, q = 0.01,
                      do.parallel = 1)
        b=proc.time()
        if(length(res$ecp)==0){
          res$ecp = NULL
        }
        rezi[["sbs_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["sbs_J"]]= length(res$ecp)
        rezi[["sbs_chgpts"]]= sort(res$ecp)
        rezi[["sbs_hausd"]] = hausdorff(sort(res$ecp), etas,n)
        rezi[["sbs_J_error"]] = length(res$ecp) - length(etas)
        rezi[["sbs_ari"]] = ARI(etas, sort(res$ecp), n)
        if(p<=100){


          # DC
          a = proc.time()
          res = dcbs.alg(X[,], cp.type = 1, phi = -1, thr = NULL, trim = NULL,
                         height = NULL, temporal = TRUE, scales = NULL, diag = FALSE,
                         B = 100, q = 0.01, do.parallel = 1)
          b=proc.time()
          if(length(res$ecp)==0){
            res$ecp = NULL
          }
          rezi[["dc_time"]] = (b-a)[3]
          rezi[["dc_J"]]= length(res$ecp)
          rezi[["dc_chgpts"]]= sort(res$ecp)
          rezi[["dc_hausd"]] = hausdorff(sort(res$ecp), etas,n)
          rezi[["dc_J_error"]] = length(res$ecp) - length(etas)
          rezi[["dc_ari"]] = ARI(etas, sort(res$ecp), n)
        }else{

          rezi[["dc_time"]] = NA
          rezi[["dc_J"]]= NA
          rezi[["dc_chgpts"]]= NA
          rezi[["dc_hausd"]] = NA
          rezi[["dc_J_error"]] = NA
          rezi[["dc_ari"]] = NA
        }

        rezi[["hehe"]] = res


        rezi[["true_J"]] = length(etas)
        rezi[["true_etas"]] = etas
        rezi[["true_sparsities"]] = sparsities

        rez[[counter]] = rezi
        counter = counter+1



      }
    }
  }
  rez
}

close(pb)
stopCluster(cl)


# gather results into arrays
{
  esac_fast_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  esac_fast_Jerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  esac_fast_time = array(0, dim= c(length(ns), length(ps), numconfig))
  esac_fast_ari = array(0, dim= c(length(ns), length(ps), numconfig))

  esac_long_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  esac_long_Jerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  esac_long_time = array(0, dim= c(length(ns), length(ps), numconfig))
  esac_long_ari = array(0, dim= c(length(ns), length(ps), numconfig))

  inspect_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  inspect_Jerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  inspect_time = array(0, dim= c(length(ns), length(ps), numconfig))
  inspect_ari = array(0, dim= c(length(ns), length(ps), numconfig))

  pilliat_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  pilliat_Jerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  pilliat_time = array(0, dim= c(length(ns), length(ps), numconfig))
  pilliat_ari = array(0, dim= c(length(ns), length(ps), numconfig))

  if(subset_included){
    subset_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
    subset_Jerr = array(0, dim = c(length(ns), length(ps), numconfig) )
    subset_time = array(0, dim= c(length(ns), length(ps), numconfig))
    subset_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  }else{
    subset_hausd = array(NA, dim = c(length(ns), length(ps), numconfig) )
    subset_Jerr = array(NA, dim = c(length(ns), length(ps), numconfig) )
    subset_time = array(NA, dim= c(length(ns), length(ps), numconfig))
    subset_ari = array(NA, dim= c(length(ns), length(ps), numconfig))
  }

  sbs_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  sbs_Jerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  sbs_time = array(0, dim= c(length(ns), length(ps), numconfig))
  sbs_ari = array(0, dim= c(length(ns), length(ps), numconfig))

  dc_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  dc_Jerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  dc_time = array(0, dim= c(length(ns), length(ps), numconfig))
  dc_ari = array(0, dim= c(length(ns), length(ps), numconfig))

  if(kaul_included){
    kaul24_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
    kaul24_Jerr = array(0, dim = c(length(ns), length(ps), numconfig) )
    kaul24_time = array(0, dim= c(length(ns), length(ps), numconfig))
    kaul24_ari = array(0, dim= c(length(ns), length(ps), numconfig))
    kaul24_fail = array(0, dim= c(length(ns), length(ps), numconfig))
  }else{
    kaul24_hausd = array(NA, dim = c(length(ns), length(ps), numconfig) )
    kaul24_Jerr = array(NA, dim = c(length(ns), length(ps), numconfig) )
    kaul24_time = array(NA, dim= c(length(ns), length(ps), numconfig))
    kaul24_ari = array(NA, dim= c(length(ns), length(ps), numconfig))
    kaul24_fail = array(NA, dim= c(length(ns), length(ps), numconfig))
  }


  pilliat_hausd[,,1] = NA
  esac_long_hausd[,,1] = NA
  esac_fast_hausd[,,1] = NA
  inspect_hausd[,,1] = NA
  subset_hausd[,,1] = NA
  sbs_hausd[,,1] = NA
  dc_hausd[,,1] = NA

  for (z in 1:N) {
    list = result[[z]]
    len = length(list)

    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      i = sublist[["i"]]
      j = sublist[["j"]]

      esac_fast_Jerr[i,j,y] = esac_fast_Jerr[i,j,y] + abs(sublist[["esac_fast_J_error"]])/N
      esac_long_Jerr[i,j,y] = esac_long_Jerr[i,j,y] + abs(sublist[["esac_long_J_error"]])/N
      pilliat_Jerr[i,j,y] = pilliat_Jerr[i,j,y] + abs(sublist[["pilliat_J_error"]])/N
      inspect_Jerr[i,j,y] = inspect_Jerr[i,j,y] + abs(sublist[["inspect_J_error"]])/N
      sbs_Jerr[i,j,y] = sbs_Jerr[i,j,y] + abs(sublist[["sbs_J_error"]])/N
      dc_Jerr[i,j,y] = dc_Jerr[i,j,y] + abs(sublist[["dc_J_error"]])/N
      if(subset_included){
        subset_Jerr[i,j,y] = subset_Jerr[i,j,y] + abs(sublist[["subset_J_error"]])/N
      }
      if(kaul_included){
        kaul24_Jerr[i,j,y] = kaul24_Jerr[i,j,y] + abs(sublist[["kaul24_J_error"]])/N
      }

      esac_fast_time[i,j,y] = esac_fast_time[i,j,y] + sublist[["esac_fast_time"]]/N
      esac_long_time[i,j,y] = esac_long_time[i,j,y] + sublist[["esac_long_time"]]/N
      pilliat_time[i,j,y] = pilliat_time[i,j,y] + sublist[["pilliat_time"]]/N
      inspect_time[i,j,y] = inspect_time[i,j,y] + sublist[["inspect_time"]]/N
      sbs_time[i,j,y] = sbs_time[i,j,y] + sublist[["sbs_time"]]/N
      dc_time[i,j,y] = dc_time[i,j,y] + sublist[["dc_time"]]/N
      if(subset_included){
        subset_time[i,j,y] = subset_time[i,j,y] + abs(sublist[["subset_time"]])/N
      }
      if(kaul_included){
        kaul24_time[i,j,y] = kaul24_time[i,j,y] + abs(sublist[["kaul24_time"]])/N
      }

      esac_fast_ari[i,j,y] = esac_fast_ari[i,j,y] + sublist[["esac_fast_ari"]]/N
      esac_long_ari[i,j,y] = esac_long_ari[i,j,y] + sublist[["esac_long_ari"]]/N
      pilliat_ari[i,j,y] = pilliat_ari[i,j,y] + sublist[["pilliat_ari"]]/N
      inspect_ari[i,j,y] = inspect_ari[i,j,y] + sublist[["inspect_ari"]]/N
      sbs_ari[i,j,y] = sbs_ari[i,j,y] + sublist[["sbs_ari"]]/N
      dc_ari[i,j,y] = dc_ari[i,j,y] + sublist[["dc_ari"]]/N
      if(subset_included){
        subset_ari[i,j,y] = subset_ari[i,j,y] + abs(sublist[["subset_ari"]])/N
      }
      if(kaul_included){
        kaul24_ari[i,j,y] = kaul24_ari[i,j,y] + abs(sublist[["kaul24_ari"]])/N
      }


      if(y!= 1){
        esac_fast_hausd[i,j,y] = esac_fast_hausd[i,j,y] + sublist[["esac_fast_hausd"]]/N
        esac_long_hausd[i,j,y] = esac_long_hausd[i,j,y] + sublist[["esac_long_hausd"]]/N
        pilliat_hausd[i,j,y] = pilliat_hausd[i,j,y] + sublist[["pilliat_hausd"]]/N
        inspect_hausd[i,j,y] = inspect_hausd[i,j,y] + sublist[["inspect_hausd"]]/N
        sbs_hausd[i,j,y] = sbs_hausd[i,j,y] + sublist[["sbs_hausd"]]/N
        dc_hausd[i,j,y] = dc_hausd[i,j,y] + sublist[["dc_hausd"]]/N
        if(subset_included){
          subset_hausd[i,j,y] = subset_hausd[i,j,y] + abs(sublist[["subset_hausd"]])/N
        }
        if(kaul_included){
          kaul24_hausd[i,j,y] = kaul24_hausd[i,j,y] + abs(sublist[["kaul24_hausd"]])/N
        }
      }



    }

  }


}

if(save){
  saveRDS(result, file=sprintf("%s/result.RDA", savedir))
  saveRDS(esac_fast_hausd, file=sprintf("%s/esac_fast_haus.RDA", savedir))
  saveRDS(esac_fast_time, file=sprintf("%s/esac_fast_time.RDA", savedir))
  saveRDS(esac_fast_Jerr, file=sprintf("%s/esac_fast_Jerr.RDA", savedir))
  saveRDS(esac_fast_ari, file=sprintf("%s/esac_fast_ari.RDA", savedir))

  saveRDS(esac_long_hausd, file=sprintf("%s/esac_long_haus.RDA", savedir))
  saveRDS(esac_long_time, file=sprintf("%s/esac_long_time.RDA", savedir))
  saveRDS(esac_long_Jerr, file=sprintf("%s/esac_long_Jerr.RDA", savedir))
  saveRDS(esac_long_ari, file=sprintf("%s/esac_long_ari.RDA", savedir))

  saveRDS(inspect_hausd, file=sprintf("%s/inspect_haus.RDA", savedir))
  saveRDS(inspect_time, file=sprintf("%s/inspect_time.RDA", savedir))
  saveRDS(inspect_Jerr, file=sprintf("%s/inspect_Jerr.RDA", savedir))
  saveRDS(inspect_ari, file=sprintf("%s/inspect_ari.RDA", savedir))

  saveRDS(pilliat_hausd, file=sprintf("%s/pilliat_haus.RDA", savedir))
  saveRDS(pilliat_time, file=sprintf("%s/pilliat_time.RDA", savedir))
  saveRDS(pilliat_Jerr, file=sprintf("%s/pilliat_Jerr.RDA", savedir))
  saveRDS(pilliat_ari, file=sprintf("%s/pilliat_ari.RDA", savedir))

  if(subset_included){
    saveRDS(subset_hausd, file=sprintf("%s/subset_haus.RDA", savedir))
    saveRDS(subset_time, file=sprintf("%s/subset_time.RDA", savedir))
    saveRDS(subset_Jerr, file=sprintf("%s/subset_Jerr.RDA", savedir))
    saveRDS(subset_ari, file=sprintf("%s/subset_ari.RDA", savedir))
  }
  if(kaul_included){
    saveRDS(kaul24_hausd, file=sprintf("%s/kaul24_haus.RDA", savedir))
    saveRDS(kaul24_time, file=sprintf("%s/kaul24_time.RDA", savedir))
    saveRDS(kaul24_Jerr, file=sprintf("%s/kaul24_Jerr.RDA", savedir))
    saveRDS(kaul24_ari, file=sprintf("%s/kaul24_ari.RDA", savedir))
  }


  saveRDS(sbs_hausd, file=sprintf("%s/sbs_haus.RDA", savedir))
  saveRDS(sbs_time, file=sprintf("%s/sbs_time.RDA", savedir))
  saveRDS(sbs_Jerr, file=sprintf("%s/sbs_Jerr.RDA", savedir))
  saveRDS(sbs_ari, file=sprintf("%s/sbs_ari.RDA", savedir))

  saveRDS(dc_hausd, file=sprintf("%s/dc_haus.RDA", savedir))
  saveRDS(dc_time, file=sprintf("%s/dc_time.RDA", savedir))
  saveRDS(dc_Jerr, file=sprintf("%s/dc_Jerr.RDA", savedir))
  saveRDS(dc_ari, file=sprintf("%s/dc_ari.RDA", savedir))


  infofile<-file(sprintf("%s/parameters.txt", savedir))
  writeLines(c(sprintf("N = %d", N),
               sprintf("n = %s", paste(ns, collapse=",")),
               sprintf("p = %s", paste(ps, collapse=",")),
               sprintf("Sparse constant = %f", sparse_const),
               sprintf("Dense constant = %f", dense_const),
               sprintf("Rescale variance = %d", as.integer(rescale_variance))),
             infofile)
  close(infofile)
}



# creating tables. one for the hausdorff distance, and one for the estimation error of J
if(save){
  # output latex table
  printlines1 = c("%%REMEMBER to use package \\usepackage{rotating}!!",
                 " \\begin{table} \\centering",
                 "\\caption{Multiple changepoints, Hausdorff distance}",
                 "\\label{multihausdorff}",
                 "\\small",
                 "\\begin{adjustbox}{width=\\columnwidth}",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} cccc|ccccccc}",
                 "\\hline",
                 "\\multicolumn{4}{c|}{Parameters} & \\multicolumn{7}{c}{Hausdorff distance}  \\\\ \\hline ",
                 "$n$ & $p$ & Sparsity & J & \\text{ESAC} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} & \\text{Kaul et al} \\\\",
                 "\\hline \\")

  printlines2 = c("%%REMEMBER to use package \\usepackage{rotating}!!",
                  " \\begin{table} \\centering",
                  "\\caption{Multiple changepoints, estimation of $J$}",
                  "\\label{multijerror}",
                  "\\small",
                  "\\begin{adjustbox}{width=\\columnwidth}",
                  "\\begin{tabular}{@{\\extracolsep{1pt}} cccc|ccccccc}",
                  "\\hline",
                  "\\multicolumn{4}{c|}{Parameters} &\\multicolumn{7}{c}{$\\left | \\widehat{J}-J \\right |$} \\\\ \\hline ",
                  "$n$ & $p$ & Sparsity & J & \\text{ESAC} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} & \\text{Kaul et al} \\\\",
                  "\\hline \\")



  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:numconfig) {
        conf = config(y, n, p)
        etas = conf[[1]]

        sparsity = conf[[4]]
        if(is.null(sparsity)){
          sparsity="-"
        }

        
        string1 = sprintf("%d & %d & %s & %d ", n, p, sparsity, length(etas))
        string2 = sprintf("%d & %d & %s & %d ", n, p, sparsity, length(etas))

        if(y==1){
          for (t in 1:7) {
            string1 = sprintf("%s & - ", string1)
            #string2 = sprintf("%s & - ", string2)

          }
        }


        else{
          res = round(c(esac_long_hausd[i,j,y],  pilliat_hausd[i,j,y], inspect_hausd[i,j,y], sbs_hausd[i,j,y], subset_hausd[i,j,y],dc_hausd[i,j,y],kaul24_hausd[i,j,y]),digits=2)
          res[is.na(res)] = Inf
          minind = (res==min(na.omit(res)))

          for (t in 1:length(res)) {
            if(is.na(res[t])){
              string1 = sprintf("%s & - ", string1)
            }
            else if(minind[t]){
              string1 = sprintf("%s & \\textbf{%.2f} ", string1, res[t])
            }else{
              string1 = sprintf("%s & %.2f", string1, res[t])
            }
          }
        }



        res = round(c(esac_long_Jerr[i,j,y], pilliat_Jerr[i,j,y], inspect_Jerr[i,j,y], sbs_Jerr[i,j,y], subset_Jerr[i,j,y],dc_Jerr[i,j,y],kaul24_Jerr[i,j,y]),digits=2)
        res[is.na(res)] = Inf
        minind = (abs(res)==min(abs(res)))

        for (t in 1:length(res)) {
          if(minind[t]){
            string2 = sprintf("%s & \\textbf{%.2f} ", string2, res[t])
          }else{
            string2 = sprintf("%s & %.2f", string2, res[t])
          }
        }

        string1 = sprintf("%s \\\\",string1)
        string2 = sprintf("%s \\\\",string2)
        printlines1 = c(printlines1, string1)
        printlines2 = c(printlines2, string2)
      }
    }
  }

  printlines1 = c(printlines1, "\\hline \\multicolumn{4}{c|}{Average}")
  printlines2 = c(printlines2, "\\hline \\multicolumn{4}{c|}{Average}")

  res = round(c(mean(na.omit(as.vector(esac_long_hausd))), mean(na.omit(as.vector(pilliat_hausd))),mean(na.omit(as.vector(inspect_hausd))), mean(na.omit(as.vector(sbs_hausd))), mean(na.omit(as.vector(subset_hausd))), mean(na.omit(as.vector(dc_hausd))), mean(na.omit(as.vector(kaul24_hausd)))),digits=2)
  res[is.na(res)] = Inf
  minind = (res==min(res))
  string =""
  for (t in 1:length(res)) {
    if(minind[t]){
      string = sprintf("%s & \\textbf{%.2f} ", string, res[t])
    }else{
      string = sprintf("%s & %.2f", string, res[t])
    }
  }
  string = sprintf("%s \\\\", string)
  printlines1 = c(printlines1, string)

  res = round(c(mean(esac_long_Jerr), mean(pilliat_Jerr),mean(inspect_Jerr), mean(sbs_Jerr), mean(subset_Jerr), mean(dc_Jerr), mean(kaul24_Jerr)),digits=2)
  res[is.na(res)] = Inf
  minind = (res==min(res))
  string =""
  for (t in 1:length(res)) {
    if(minind[t]){
      string = sprintf("%s & \\textbf{%.2f} ", string, res[t])
    }else{
      string = sprintf("%s & %.2f", string, res[t])
    }
  }
  string = sprintf("%s \\\\", string)
  printlines2 = c(printlines2, string)

  printlines1 = c(printlines1, c("\\hline \\\\[-1.8ex]",
                               "\\end{tabular}",
                               "\\end{adjustbox}",
                               "\\end{table}"))
  printlines2 = c(printlines2, c("\\hline \\\\[-1.8ex]",
                                 "\\end{tabular}",
                                 "\\end{adjustbox}",
                                 "\\end{table}"))
  texfile<-file(sprintf("%s/table_hausd.tex", savedir))
  writeLines(printlines1, texfile)
  close(texfile)

  texfile<-file(sprintf("%s/table_Jerror.tex", savedir))
  writeLines(printlines2, texfile)
  close(texfile)

}


# smaller version of the table, where p = 5000 is omitted:
# run this only if length(ps)>=3
if(FALSE){
  # output latex table
  printlines1 = c("%%REMEMBER to use package \\usepackage{rotating}!!",
                  " \\begin{table} \\centering",
                  "\\caption{Multiple changepoints, Hausdorff distance}",
                  "\\label{multihausdorff}",
                  "\\small",
                  "\\begin{adjustbox}{width=\\columnwidth}",
                  "\\begin{tabular}{@{\\extracolsep{1pt}} cccc|ccccccc}",
                  "\\hline",
                  "\\multicolumn{4}{c|}{Parameters} & \\multicolumn{7}{c}{Hausdorff distance}  \\\\ \\hline ",
                  "$n$ & $p$ & Sparsity & J & \\text{ESAC} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} & \\text{Kaul et al} \\\\",
                  "\\hline \\")

  printlines2 = c("%%REMEMBER to use package \\usepackage{rotating}!!",
                  " \\begin{table} \\centering",
                  "\\caption{Multiple changepoints, estimation of $J$}",
                  "\\label{multijerror}",
                  "\\small",
                  "\\begin{adjustbox}{width=\\columnwidth}",
                  "\\begin{tabular}{@{\\extracolsep{1pt}} cccc|ccccccc}",
                  "\\hline",
                  "\\multicolumn{4}{c|}{Parameters} &\\multicolumn{7}{c}{$\\left | \\widehat{J}-J \\right |$} \\\\ \\hline ",
                  "$n$ & $p$ & Sparsity & J & \\text{ESAC} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} & \\text{Kaul et al} \\\\",
                  "\\hline \\")



  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:(length(ps)-1)){
      p = ps[j]
      for (y in 1:numconfig) {
        conf = config(y, n, p)
        etas = conf[[1]]
       
        sparsity = conf[[4]]
        if(is.null(sparsity)){
          sparsity="-"
        }

        
        string1 = sprintf("%d & %d & %s & %d ", n, p, sparsity, length(etas))
        string2 = sprintf("%d & %d & %s & %d ", n, p, sparsity, length(etas))

        if(y==1){
          for (t in 1:7) {
            string1 = sprintf("%s & - ", string1)
            #string2 = sprintf("%s & - ", string2)

          }
        }else{
          res = round(c(esac_long_hausd[i,j,y],  pilliat_hausd[i,j,y], inspect_hausd[i,j,y], sbs_hausd[i,j,y], subset_hausd[i,j,y],dc_hausd[i,j,y],kaul24_hausd[i,j,y]),digits=2)
          res[is.na(res)] = Inf
          minind = (res==min(na.omit(res)))

          for (t in 1:length(res)) {
            if(is.na(res[t])){
              string1 = sprintf("%s & - ", string1)
            }
            else if(minind[t]){
              string1 = sprintf("%s & \\textbf{%.2f} ", string1, res[t])
            }else{
              string1 = sprintf("%s & %.2f", string1, res[t])
            }
          }
        }



        res = round(c(esac_long_Jerr[i,j,y], pilliat_Jerr[i,j,y], inspect_Jerr[i,j,y], sbs_Jerr[i,j,y], subset_Jerr[i,j,y],dc_Jerr[i,j,y],kaul24_Jerr[i,j,y]),digits=2)
        res[is.na(res)] = Inf
        minind = (abs(res)==min(abs(res)))

        for (t in 1:length(res)) {
          if(minind[t]){
            string2 = sprintf("%s & \\textbf{%.2f} ", string2, res[t])
          }else{
            string2 = sprintf("%s & %.2f", string2, res[t])
          }
        }

      
        string1 = sprintf("%s \\\\",string1)
        string2 = sprintf("%s \\\\",string2)
        printlines1 = c(printlines1, string1)
        printlines2 = c(printlines2, string2)
      }
    }
  }

  printlines1 = c(printlines1, "\\hline \\multicolumn{4}{c|}{Average}")
  printlines2 = c(printlines2, "\\hline \\multicolumn{4}{c|}{Average}")

  res = round(c(mean(na.omit(as.vector(esac_long_hausd[,1:2,]))), mean(na.omit(as.vector(pilliat_hausd[,1:2,]))),mean(na.omit(as.vector(inspect_hausd[,1:2,]))), mean(na.omit(as.vector(sbs_hausd[,1:2,]))), mean(na.omit(as.vector(subset_hausd[,1:2,]))), mean(na.omit(as.vector(dc_hausd[,1:2,]))), mean(na.omit(as.vector(kaul24_hausd[,1:2,])))),digits=2)
  res[is.na(res)] = Inf
  minind = (res==min(res))
  string =""
  for (t in 1:length(res)) {
    if(minind[t]){
      string = sprintf("%s & \\textbf{%.2f} ", string, res[t])
    }else{
      string = sprintf("%s & %.2f", string, res[t])
    }
  }
  string = sprintf("%s \\\\", string)
  printlines1 = c(printlines1, string)

  res = round(c(mean(esac_long_Jerr[,1:2,]), mean(pilliat_Jerr[,1:2,]),mean(inspect_Jerr[,1:2,]), mean(sbs_Jerr[,1:2,]), mean(subset_Jerr[,1:2,]), mean(dc_Jerr[,1:2,]), mean(kaul24_Jerr[,1:2,])),digits=2)
  res[is.na(res)] = Inf
  minind = (res==min(res))
  string =""
  for (t in 1:length(res)) {
    if(minind[t]){
      string = sprintf("%s & \\textbf{%.2f} ", string, res[t])
    }else{
      string = sprintf("%s & %.2f", string, res[t])
    }
  }
  string = sprintf("%s \\\\", string)
  printlines2 = c(printlines2, string)

  printlines1 = c(printlines1, c("\\hline \\\\[-1.8ex]",
                                 "\\end{tabular}",
                                 "\\end{adjustbox}",
                                 "\\end{table}"))
  printlines2 = c(printlines2, c("\\hline \\\\[-1.8ex]",
                                 "\\end{tabular}",
                                 "\\end{adjustbox}",
                                 "\\end{table}"))
  texfile<-file(sprintf("%s/table_small_hausd.tex", savedir))
  writeLines(printlines1, texfile)
  close(texfile)

  texfile<-file(sprintf("%s/table_small_Jerror.tex", savedir))
  writeLines(printlines2, texfile)
  close(texfile)

}


# combined small table, as presented in the article:
# small table:
if(save){
  # output latex table
  printlines1 = c("%%REMEMBER to use package \\usepackage{rotating}!!",
                  " \\begin{table} \\centering",
                  "\\caption{Multiple changepoints, Hausdorff distance and estimation error of $J$}",
                  "\\label{multi_small_combined}",
                  "\\small",
                  "\\begin{adjustbox}{width=\\columnwidth}",
                  "\\begin{tabular}{@{\\extracolsep{1pt}} cccc|ccccccc}",
                  "\\hline",
                  "\\multicolumn{4}{c|}{Parameters} & \\multicolumn{7}{c}{Hausdorff distance ($|\\widehat{J}-J|$)}  \\\\ \\hline ",
                  "$n$ & $p$ & Sparsity & J & \\text{ESAC} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} & \\text{Kaul et al} \\\\",
                  "\\hline \\")



  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:(length(ps)-1)){
      p = ps[j]
      for (y in 1:numconfig) {
        conf = config(y, n, p)
        etas = conf[[1]]
       
        sparsity = conf[[4]]
        if(is.null(sparsity)){
          sparsity="-"
        }


        string1 = sprintf("%d & %d & %s & %d ", n, p, sparsity, length(etas))






        if(y==1){
          resmse = rep(NA,7)
          minindmse = rep(FALSE,7)
        }else{
          resmse = round(c(esac_long_hausd[i,j,y],  pilliat_hausd[i,j,y], inspect_hausd[i,j,y], sbs_hausd[i,j,y], subset_hausd[i,j,y],dc_hausd[i,j,y],kaul24_hausd[i,j,y]),digits=2)
          resmse[is.na(resmse)] = Inf
          minindmse = (resmse==min(na.omit(resmse)))
          resmse[is.infinite(resmse)] = NA
        }

        resJ = round(c(esac_long_Jerr[i,j,y], pilliat_Jerr[i,j,y], inspect_Jerr[i,j,y], sbs_Jerr[i,j,y], subset_Jerr[i,j,y],dc_Jerr[i,j,y],kaul24_Jerr[i,j,y]),digits=2)
        resJ[is.na(resJ)] = Inf
        minindJ = (abs(resJ)==min(abs(resJ)))
        resJ[is.infinite(resJ)] = NA

        for (t in 1:length(resmse)) {
          if(is.na(resmse[t])){
            string1 = sprintf("%s & - ", string1)
          }
          else if(minindmse[t]){
            string1 = sprintf("%s & \\textbf{%.2f} ", string1, resmse[t])
          }else{
            string1 = sprintf("%s & %.2f", string1, resmse[t])
          }

          if(is.na(resJ[t])){
            string1 = sprintf("%s (-) ", string1)
          }
          else if(minindJ[t]){
            string1 = sprintf("%s (\\textbf{%.2f}) ", string1, resJ[t])
          }else{
            string1 = sprintf("%s (%.2f)", string1, resJ[t])
          }

        }




       
        string1 = sprintf("%s \\\\",string1)
        printlines1 = c(printlines1, string1)
      }
    }
  }

  printlines1 = c(printlines1, "\\hline \\multicolumn{4}{c|}{Average}")


  resmse = round(c(mean(na.omit(as.vector(esac_long_hausd[,1:2,]))), mean(na.omit(as.vector(pilliat_hausd[,1:2,]))),mean(na.omit(as.vector(inspect_hausd[,1:2,]))), mean(na.omit(as.vector(sbs_hausd[,1:2,]))), mean(na.omit(as.vector(subset_hausd[,1:2,]))), mean(na.omit(as.vector(dc_hausd[,1:2,]))), mean(na.omit(as.vector(kaul24_hausd[,1:2,])))),digits=2)
  resmse[is.na(resmse)] = Inf
  minindmse = (resmse==min(resmse))


  dc_mod_Jerr = na.omit(as.vector(dc_Jerr[!is.infinite(dc_Jerr[,1:2,])]))
  resJ = round(c(mean(na.omit(esac_long_Jerr[,1:2,])), mean(na.omit(pilliat_Jerr[,1:2,])),mean(na.omit(inspect_Jerr[,1:2,])), mean(na.omit(sbs_Jerr[,1:2,])), mean(na.omit(subset_Jerr[,1:2,])), mean(dc_mod_Jerr), mean(na.omit(kaul24_Jerr[,1:2,]))),digits=2)
  resJ[is.na(resJ)] = Inf
  minindJ = (resJ==min(resJ))


  string =""
  for (t in 1:length(resmse)) {
    if(minindmse[t]){
      string = sprintf("%s & \\textbf{%.2f} ", string, resmse[t])
    }else{
      string = sprintf("%s & %.2f", string, resmse[t])
    }
    if(minindJ[t]){
      string = sprintf("%s (\\textbf{%.2f}) ", string, resJ[t])
    }else{
      string = sprintf("%s (%.2f)", string, resJ[t])
    }

  }
  string = sprintf("%s \\\\", string)
  printlines1 = c(printlines1, string)


  printlines1 = c(printlines1, c("\\hline \\\\[-1.8ex]",
                                 "\\end{tabular}",
                                 "\\end{adjustbox}",
                                 "\\end{table}"))

  texfile<-file(sprintf("%s/table_small_combined.tex", savedir))
  writeLines(printlines1, texfile)
  close(texfile)

}









