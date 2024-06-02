# This code runs the simulations performed in Section 4.2 of Moen et al. (2024),	arXiv:2306.04702
# NOTE: THIS CODE IS ONLY FOR TIMING THE METHODS; NOT FOR ASSESSING THEIR
#       STATISTICAL PERFORMANCE!
# To recreate the simulation from Section 4.2, uncomment lines 82,83, 608,609.
# Please note that the code for SUBSET must be downloaded from Github (https://github.com/SOTickle/SUBSET), 
# while the code for the method of Kaul et al (2021, Electronic Journal of Statistics) 
# must be obtained from the first author of that paper. 
# If SUBSET or the method of Kaul et al are not available,
# please set subset_included = FALSE and kaul_included = FALSE below.






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
subset_path = "... fill in/SUBSET ..."
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
  savedir = file.path(maindir, sprintf("%s/multi_timing",dateandtime))
  dir.create(savedir, showWarnings = FALSE)

}



# Simulation options
N = 24 #MC simulations, this is used in the paper
num_cores = 6 # cores to run in parallell
sparse_const = 4 # constant for sparse signal strength
dense_const = 4 #constant for dense signal strength
set.seed(1996)
rescale_variance = TRUE
tol = 0.01 #error tolerance for MC choice of thresholds



### First we vary n and let p be fixed
# select values of n and p:
ns = c(20,40)
ps = c(20)

#ns = seq(100,1000,by=100) #uncomment for the same values of n as in the paper
#ps = c(100) # uncomment for the same values of p as in the paper






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
      #diff = runif(k) *sample(c(-1,1), sparsities[j], replace=TRUE)
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
      #diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
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
      #diff = runif(k) *sample(c(-1,1), sparsities[j], replace=TRUE)
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
      #diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
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


      #for (y in 1:numconfig) {
      y = 4
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


      a = proc.time()
      res = ESAC (X[,], 1.5,1, empirical=FALSE,alpha = 1.5, K = 4, fast =FALSE,
                  rescale_variance = rescale_variance, debug= FALSE)
      b=proc.time()
      rezi[["esac_fast_time"]] = (b-a)[3]

      a = proc.time()
      res  = Pilliat(X, K = 2, alpha = 1.5,
                     rescale_variance = rescale_variance, test_all = TRUE,debug =FALSE)
      b=proc.time()

      rezi[["pilliat_time"]] = (b-a)[3]
      


      #inspect
      a = proc.time()
      res = Inspect(X[,], alpha = 1.5, K = 4,eps=1e-10, lambda = 4*sqrt(log(n*p)), xi = 4*sqrt(log(n*p)),
                    maxiter=10000,
                    rescale_variance = rescale_variance, debug=FALSE)
      b=proc.time()
      rezi[["inspect_time"]] = (b-a)[3]
      rezi[["inspect_J"]]= res$changepointnumber
      rezi[["inspect_chgpts"]]= res$changepoints
     




      # subset
      if(subset_included){
      a = proc.time()
      if(rescale_variance){
        res = change_main(InspectChangepoint::rescale.variance(X), SUBSET.normal, 100,penalties= c(4*log(n)))
      }else{
        res = change_main(X, SUBSET.normal,100, penalties= c(4*log(n)))
      }
      b=proc.time()
      rezi[["subset_time"]] = (b-a)[3]
     
      }


      # SBS
      a = proc.time()
      res = sbs.alg(X[,], cp.type = 1, thr = NULL, trim = NULL, height = NULL,
                    temporal = FALSE, scales = NULL, diag = FALSE, B = 100, q = 0.01,
                    do.parallel = 1)
      b=proc.time()
      if(length(res$ecp)==0){
        res$ecp = NULL
      }
      rezi[["sbs_time"]] = (b-a)[3]
     

      # kaul (2024):
      if(kaul_included){

        if(rezi[["inspect_J"]] == 0){
          # inspect detected no changepoints
          rezi[["kaul24_time"]] = rezi[["inspect_time"]]
        }else{
          # use inspect output as preliminary estimate:
          inspectpreliminary = function(x){
            aa = list()
            aa$t = rezi[["inspect_chgpts"]]
            return(aa)
          }

         
          rezult = tryCatch({
            a = 0
            b = 0
            if(rescale_variance){
              a = proc.time()
              res = alg1(t(InspectChangepoint::rescale.variance(X[,])), prel.cp = inspectpreliminary)$tilde.t
              b=proc.time()
            }else{
              a = proc.time()
              res =alg1(t(X[,]), prel.cp = inspectpreliminary)$tilde.t
              b=proc.time()
            }
            rezult = ((b-a)[3] + rezi[["inspect_time"]])

          }, error = function(e) {
            return(-1)
          })
          if(rezult==-1){
            # if alg1.R failed
            rezi[["kaul24_time"]] = rezi[["inspect_time"]]
          }else{
            rezi[["kaul24_time"]] = rezult
          }
        }

      }
      rezi[["inspect_time"]]
      rezi[["kaul24_time"]]

      # DC
      a = proc.time()
      res = dcbs.alg(X[,], cp.type = 1, phi = 0.5, thr = NULL, trim = NULL,
                     height = NULL, temporal = TRUE, scales = NULL, diag = FALSE,
                     B = 100, q = 0.01, do.parallel = 1)
      b=proc.time()
      if(length(res$ecp)==0){
        res$ecp = NULL
      }
      rezi[["dc_time"]] = (b-a)[3]
      
      rezi[["hehe"]] = res


      rezi[["true_K"]] = length(etas)
      rezi[["true_etas"]] = etas
      rezi[["true_sparsities"]] = sparsities

      rez[[counter]] = rezi
      counter = counter+1







    }
  }
  rez
}

close(pb)
stopCluster(cl)

# gather results into arrays
{

  esac_fast_time = array(0, dim= c(length(ns), length(ps), numconfig))


  inspect_time = array(0, dim= c(length(ns), length(ps), numconfig))


  pilliat_time = array(0, dim= c(length(ns), length(ps), numconfig))

  if(subset_included){

    subset_time = array(0, dim= c(length(ns), length(ps), numconfig))

  }else{

    subset_time = array(NA, dim= c(length(ns), length(ps), numconfig))
  }


  sbs_time = array(0, dim= c(length(ns), length(ps), numconfig))

  dc_time = array(0, dim= c(length(ns), length(ps), numconfig))

  if(kaul_included){

    kaul24_time = array(0, dim= c(length(ns), length(ps), numconfig))

  }else{

    kaul24_time = array(NA, dim= c(length(ns), length(ps), numconfig))

  }



  for (z in 1:N) {
    list = result[[z]]
    len = length(list)

    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      i = sublist[["i"]]
      j = sublist[["j"]]


      esac_fast_time[i,j,y] = esac_fast_time[i,j,y] + sublist[["esac_fast_time"]]/N
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



    }

  }


}

if(save){
  saveRDS(result, file=sprintf("%s/result.RDA", savedir))
  saveRDS(esac_fast_time, file=sprintf("%s/esac_fast_time.RDA", savedir))



  saveRDS(inspect_time, file=sprintf("%s/inspect_time.RDA", savedir))

  saveRDS(pilliat_time, file=sprintf("%s/pilliat_time.RDA", savedir))


  if(subset_included){
    saveRDS(subset_time, file=sprintf("%s/subset_time.RDA", savedir))

  }
  if(kaul_included){
    saveRDS(kaul24_time, file=sprintf("%s/kaul24_time.RDA", savedir))

  }

  saveRDS(sbs_time, file=sprintf("%s/sbs_time.RDA", savedir))


  saveRDS(dc_time, file=sprintf("%s/dc_time.RDA", savedir))


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


library(ggplot2)
library(reshape2)
times = cbind(ns, ESAC= 1000*esac_fast_time[,,4], Pilliat = 1000*pilliat_time[,,4], Inspect = 1000*inspect_time[,,4],
              SBS = 1000*sbs_time[,,4], SUBSET = 1000*subset_time[,,4], DC = 1000*dc_time[,,4], kaul = 1000*kaul24_time[,,4])

times = data.frame(log(times))
times$ns = ns
#times[,2:7] = log(times[,2:7])
colnames(times)[8] = "Kaul and Michailidis (2024)"
df <- melt(times ,  id.vars = 'ns', variable.name = 'Method')

p1 = ggplot(df, aes(ns,value))+ geom_line(mapping=aes(colour = Method)) + xlab("n") + ylab("Run time (ms), log scale") +
 ylim(c(0,10))
p1


























#### Now varying p: #####

# select values of n and p:
ps = c(20,40)
ns = c(20)

#ps = seq(100,1000,by=100) #uncomment for same value of p as in the paper
#ns = c(100) #uncomment for same values of n as in the paper





# run simulation in parallell:
cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result2 = foreach(z = 1:N,.options.snow = opts) %dopar% {
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


      y = 4
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


      a = proc.time()
      res = ESAC (X[,], 1.5,1, empirical=FALSE,alpha = 1.5, K = 4, fast =FALSE,
                  rescale_variance = rescale_variance, debug= FALSE)
      b=proc.time()
      rezi[["esac_fast_time"]] = (b-a)[3]

      a = proc.time()
      res  = Pilliat(X, K = 2, alpha = 1.5,
                     rescale_variance = rescale_variance, test_all = TRUE,debug =FALSE)
      b=proc.time()
     

      rezi[["pilliat_time"]] = (b-a)[3]
     

      #inspect
      a = proc.time()
      res = Inspect(X[,], alpha = 1.5, K = 4,eps=1e-10, lambda = 4*sqrt(log(n*p)), xi = 4*sqrt(log(n*p)),
                    maxiter=10000,
                    rescale_variance = rescale_variance, debug=FALSE)
      b=proc.time()
      rezi[["inspect_time"]] = (b-a)[3]
      rezi[["inspect_J"]]= res$changepointnumber
      rezi[["inspect_chgpts"]]= res$changepoints
      



      # subset
      if(subset_included){
        a = proc.time()
        if(rescale_variance){
          res = change_main(InspectChangepoint::rescale.variance(X), SUBSET.normal, 100,penalties= c(4*log(n)))
        }else{
          res = change_main(X, SUBSET.normal,100, penalties= c(4*log(n)))
        }
        b=proc.time()
        rezi[["subset_time"]] = (b-a)[3]
       
      }


      # SBS
      a = proc.time()
      res = sbs.alg(X[,], cp.type = 1, thr = NULL, trim = NULL, height = NULL,
                    temporal = FALSE, scales = NULL, diag = FALSE, B = 100, q = 0.01,
                    do.parallel = 1)
      b=proc.time()
      if(length(res$ecp)==0){
        res$ecp = NULL
      }
      rezi[["sbs_time"]] = (b-a)[3]
     
      # kaul (2024):
      if(kaul_included){

        if(rezi[["inspect_J"]] == 0){
          # inspect detected no changepoints
          rezi[["kaul24_time"]] = rezi[["inspect_time"]]
        }else{
          # use inspect output as preliminary estimate:
          inspectpreliminary = function(x){
            aa = list()
            aa$t = rezi[["inspect_chgpts"]]
            return(aa)
          }

         
          rezult = tryCatch({
            a = 0
            b = 0
            if(rescale_variance){
              a = proc.time()
              res = alg1(t(InspectChangepoint::rescale.variance(X[,])), prel.cp = inspectpreliminary)$tilde.t
              b=proc.time()
            }else{
              a = proc.time()
              res =alg1(t(X[,]), prel.cp = inspectpreliminary)$tilde.t
              b=proc.time()
            }
            rezult = ((b-a)[3] + rezi[["inspect_time"]])

          }, error = function(e) {
            return(-1)
          })
          if(rezult==-1){
            # if alg1.R failed
            rezi[["kaul24_time"]] = rezi[["inspect_time"]]
          }else{
            rezi[["kaul24_time"]] = rezult
          }
        }

      }
      rezi[["inspect_time"]]
      rezi[["kaul24_time"]]

      # DC
      a = proc.time()
      res = dcbs.alg(X[,], cp.type = 1, phi = 0.5, thr = NULL, trim = NULL,
                     height = NULL, temporal = TRUE, scales = NULL, diag = FALSE,
                     B = 100, q = 0.01, do.parallel = 1)
      b=proc.time()
      if(length(res$ecp)==0){
        res$ecp = NULL
      }
      rezi[["dc_time"]] = (b-a)[3]
      

      rezi[["hehe"]] = res


      rezi[["true_K"]] = length(etas)
      rezi[["true_etas"]] = etas
      rezi[["true_sparsities"]] = sparsities

      rez[[counter]] = rezi
      counter = counter+1





    }
  }
  rez
}

close(pb)
stopCluster(cl)


# gather results into arrays
{

  esac_fast_time2 = array(0, dim= c(length(ns), length(ps), numconfig))


  inspect_time2 = array(0, dim= c(length(ns), length(ps), numconfig))


  pilliat_time2 = array(0, dim= c(length(ns), length(ps), numconfig))

  if(subset_included){

    subset_time2 = array(0, dim= c(length(ns), length(ps), numconfig))

  }else{

    subset_time2 = array(NA, dim= c(length(ns), length(ps), numconfig))
  }


  sbs_time2 = array(0, dim= c(length(ns), length(ps), numconfig))

  dc_time2 = array(0, dim= c(length(ns), length(ps), numconfig))

  if(kaul_included){

    kaul24_time2 = array(0, dim= c(length(ns), length(ps), numconfig))

  }else{

    kaul24_time2 = array(NA, dim= c(length(ns), length(ps), numconfig))

  }



  for (z in 1:N) {
    list = result2[[z]]
    len = length(list)

    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      i = sublist[["i"]]
      j = sublist[["j"]]


      esac_fast_time2[i,j,y] = esac_fast_time2[i,j,y] + sublist[["esac_fast_time"]]/N
      pilliat_time2[i,j,y] = pilliat_time2[i,j,y] + sublist[["pilliat_time"]]/N
      inspect_time2[i,j,y] = inspect_time2[i,j,y] + sublist[["inspect_time"]]/N
      sbs_time2[i,j,y] = sbs_time2[i,j,y] + sublist[["sbs_time"]]/N
      dc_time2[i,j,y] = dc_time2[i,j,y] + sublist[["dc_time"]]/N
      if(subset_included){
        subset_time2[i,j,y] = subset_time2[i,j,y] + abs(sublist[["subset_time"]])/N
      }
      if(kaul_included){
        kaul24_time2[i,j,y] = kaul24_time2[i,j,y] + abs(sublist[["kaul24_time"]])/N
      }



    }

  }


}

if(save){
  saveRDS(result2, file=sprintf("%s/result2.RDA", savedir))
  saveRDS(esac_fast_time2, file=sprintf("%s/esac_fast_time2.RDA", savedir))




  saveRDS(inspect_time2, file=sprintf("%s/inspect_time2.RDA", savedir))

  saveRDS(pilliat_time2, file=sprintf("%s/pilliat_time2.RDA", savedir))


  if(subset_included){
    saveRDS(subset_time2, file=sprintf("%s/subset_time2.RDA", savedir))

  }
  if(kaul_included){
    saveRDS(kaul24_time2, file=sprintf("%s/kaul24_time2.RDA", savedir))

  }


  saveRDS(sbs_time2, file=sprintf("%s/sbs_time2.RDA", savedir))

  saveRDS(dc_time2, file=sprintf("%s/dc_time2.RDA", savedir))


  infofile<-file(sprintf("%s/parameters2.txt", savedir))
  writeLines(c(sprintf("N = %d", N),
               sprintf("n = %s", paste(ns, collapse=",")),
               sprintf("p = %s", paste(ps, collapse=",")),
               sprintf("Sparse constant = %f", sparse_const),
               sprintf("Dense constant = %f", dense_const),
               sprintf("Rescale variance = %d", as.integer(rescale_variance))),
             infofile)
  close(infofile)
}



library(reshape2)
times2 = cbind(ps, ESAC= 1000*esac_fast_time2[,,4], Pilliat = 1000*pilliat_time2[,,4], Inspect = 1000*inspect_time2[,,4],
              SBS = 1000*sbs_time2[,,4], SUBSET = 1000*subset_time2[,,4], DC = 1000*dc_time2[,,4], kaul = 1000*kaul24_time2[,,4])

times2 = data.frame(log(times2))
times2$ps = ps
colnames(times2)[8] = "Kaul and Michailidis (2024)"
#times[,2:7] = log(times[,2:7])
df <- melt(times2 ,  id.vars = 'ps', variable.name = 'Method')

p2 = ggplot(df, aes(ps,value))+ geom_line(mapping=aes(colour = Method)) + xlab("p") +
  theme(axis.title.y=element_blank())+
  ylim(c(0,10))
p2




# creating combined plot: (this is the plot in the paper)

library(patchwork)
p1+p2 + plot_layout(ncol = 2, guides = "collect")



ggsave(
  sprintf("%s/timings_multi.eps", savedir),
  plot = last_plot(),
  device = cairo_ps,
  scale = 1,
  width = 8,
  height = 3,
  dpi = 300,
  limitsize = TRUE,
  fallback_resolution=1000
)
















































