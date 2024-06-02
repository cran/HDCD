#### Simulation for multiple change-points
### Here we recreate the simulation study from Appendix D in Moen et al. (2024),	arXiv:2306.04702
##  To recreate the table in the paper, uncomment lines 37 and 38
library(doSNOW)
library(HDCD)
library(foreach)

maindir = "... fill in ... "
dateandtime = gsub(" ", "--",as.character(Sys.time()))
dateandtime = gsub(":", ".", dateandtime)
savedir = file.path(maindir, dateandtime)


save = TRUE

if(save){
  dir.create(savedir, showWarnings = FALSE)
  savedir = file.path(maindir, sprintf("%s/multi_ESAC_comparison",dateandtime))
  dir.create(savedir, showWarnings = FALSE)
  
}



N = 1000
num_cores = 6
sparse_const = 3.5
dense_const = 3.5
set.seed(1996)
rescale_variance = TRUE
Ncal = 1000
tol = 1/Ncal

ns = c(10,20)
ps = c(10,20)

#ns = c(200,500) #uncomment for using the same values of n as in the article
#ps = c(50,100,1000) #uncomment for using the same values of p as in the article



# calibrate ESAC thresholds:

cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = ((length(ns)*length(ps))), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
calibrates = foreach(z = 0:((length(ns)*length(ps))-1),.options.snow = opts) %dopar% {
  library(HDCD)
  set.seed(300*z)
  rez = list()
  nind = floor(z/length(ps))+1
  pind = z%%length(ps)+1
  cc = ESAC_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K=4, alpha = 2, fast = TRUE,
                      rescale_variance = rescale_variance, debug=FALSE)
 
  
  
  cc2 = ESAC_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K=4, alpha = 2, fast = FALSE,
                      rescale_variance = rescale_variance, debug=FALSE)
  
  cc3 = ESAC_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K=5, alpha = 1+1/2, fast = FALSE,
                       rescale_variance = rescale_variance, debug=FALSE)
  
  rez[[1]] = cc[[1]] #ESAC threshoolds
  rez[[2]] = cc[[2]] #ESAC thresholds
  rez[[3]] = cc2[[1]] 
  rez[[4]] = cc2[[2]] 
  rez[[5]] = cc3[[1]] 
  rez[[6]] = cc3[[2]] 
  
  rez[[7]] = ns[nind]
  rez[[8]] = ps[pind]
  rez
  #calibrates[[z+1]] = rez
  
}
close(pb)
stopCluster(cl) 


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


cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
#for(z in 1:N){
result = foreach(z = 1:N,.options.snow = opts) %dopar% {
  rez = list()
  set.seed(2*z)
  library(HDCD)
  library(hdbinseg)
  counter = 1

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
        
        
        
        a = proc.time()
        res = ESAC (X[,], 1.5,1, empirical=TRUE,alpha = 2, K = 4, thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[2]], fast =TRUE,
                    rescale_variance = rescale_variance, debug= FALSE)
        b=proc.time()
        rezi[["ESAC_fast_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["ESAC_fast_K"]]= res$changepointnumber
        rezi[["ESAC_fast_chgpts"]]= res$changepoints
        if(res$changepointnumber==0){
          rezi[["ESAC_fast_hausd"]] = n
        }else{
          rezi[["ESAC_fast_hausd"]] = hausdorff(res$changepoints, etas,n)
        }
        rezi[["ESAC_fast_K"]] = length(res$changepoints) - length(etas)

        
        
        a = proc.time()
        res = ESAC (X[,], 1.5,1, empirical=TRUE,alpha = 2, K = 4, thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[1]],  fast =TRUE,
                    rescale_variance = rescale_variance, debug= FALSE)
        b=proc.time()
        rezi[["ESAC_fast_np_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["ESAC_fast_np_K"]]= res$changepointnumber
        rezi[["ESAC_fast_np_chgpts"]]= res$changepoints
        if(res$changepointnumber==0){
          rezi[["ESAC_fast_np_hausd"]] = n
        }else{
          rezi[["ESAC_fast_np_hausd"]] = hausdorff(res$changepoints, etas,n)
        }
        rezi[["ESAC_fast_np_K"]] = length(res$changepoints) - length(etas)
        
        
        a = proc.time()
        res = ESAC (X[,], 1.5,1, empirical=TRUE,alpha = 2, K = 4, thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[3]],  fast =FALSE,
                    rescale_variance = rescale_variance, debug= FALSE)
        b=proc.time()
        rezi[["ESAC_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["ESAC_K"]]= res$changepointnumber
        rezi[["ESAC_chgpts"]]= res$changepoints
        if(res$changepointnumber==0){
          rezi[["ESAC_hausd"]] = n
        }else{
          rezi[["ESAC_hausd"]] = hausdorff(res$changepoints, etas,n)
        }
        rezi[["ESAC_K"]] = length(res$changepoints) - length(etas)

        a = proc.time()
        res = ESAC (X[,], 1.5,1, empirical=TRUE,alpha = 1.5, K = 5, thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[5]],  fast =FALSE,
                    rescale_variance = rescale_variance, debug= FALSE)
        b=proc.time()
        rezi[["ESAC_long_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["ESAC_long_K"]]= res$changepointnumber
        rezi[["ESAC_long_chgpts"]]= res$changepoints
        if(res$changepointnumber==0){
          rezi[["ESAC_long_hausd"]] = n
        }else{
          rezi[["ESAC_long_hausd"]] = hausdorff(res$changepoints, etas,n)
        }
        
        rezi[["ESAC_long_K"]] = length(res$changepoints) - length(etas)

        a = proc.time()
        res = ESAC (X[,], 1.5,1, empirical=TRUE,alpha = 2, K = 4, thresholds = (calibrates[[j+(i-1)*length(ps)]])[[3]], thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[3]],  fast =FALSE,
                    rescale_variance = rescale_variance, debug= FALSE)
        b=proc.time()
        rezi[["ESAC_emp_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["ESAC_emp_K"]]= res$changepointnumber
        rezi[["ESAC_emp_chgpts"]]= res$changepoints
        if(res$changepointnumber==0){
          rezi[["ESAC_emp_hausd"]] = n
        }else{
          rezi[["ESAC_emp_hausd"]] = hausdorff(res$changepoints, etas,n)
        }
        rezi[["ESAC_emp_K"]] = length(res$changepoints) - length(etas)

        a = proc.time()
        res = ESAC (X[,], 1.5,1, empirical=TRUE,alpha = 1.5, K = 5, thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[5]], thresholds = (calibrates[[j+(i-1)*length(ps)]])[[5]],  fast =FALSE,
                    rescale_variance = rescale_variance, debug= FALSE)
        b=proc.time()
        rezi[["ESAC_long_emp_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["ESAC_long_emp_K"]]= res$changepointnumber
        rezi[["ESAC_long_emp_chgpts"]]= res$changepoints
        if(res$changepointnumber==0){
          rezi[["ESAC_long_emp_hausd"]] = n
        }else{
          rezi[["ESAC_long_emp_hausd"]] = hausdorff(res$changepoints, etas,n)
        }
        rezi[["ESAC_long_emp_K"]] = length(res$changepoints) - length(etas)

        
        
        
        rezi[["true_K"]] = length(etas)
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

{
  ESAC_fast_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_fast_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_fast_time = array(0, dim= c(length(ns), length(ps), numconfig))
  ESAC_fast_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  ESAC_fast_np_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_fast_np_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_fast_np_time = array(0, dim= c(length(ns), length(ps), numconfig))
  ESAC_fast_np_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  ESAC_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_time = array(0, dim= c(length(ns), length(ps), numconfig))
  ESAC_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  ESAC_long_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_long_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_long_time = array(0, dim= c(length(ns), length(ps), numconfig))
  ESAC_long_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  ESAC_emp_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_emp_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_emp_time = array(0, dim= c(length(ns), length(ps), numconfig))
  ESAC_emp_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  ESAC_long_emp_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_long_emp_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  ESAC_long_emp_time = array(0, dim= c(length(ns), length(ps), numconfig))
  ESAC_long_emp_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  
  ESAC_fast_hausd[,,1] = NA
  ESAC_fast_np_hausd[,,1] = NA
  ESAC_hausd[,,1] = NA
  ESAC_emp_hausd[,,1] = NA
  ESAC_long_emp_hausd[,,1] = NA
  ESAC_long_hausd[,,1] = NA
  
  for (z in 1:N) {
    list = result[[z]]
    len = length(list)
    
    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      i = sublist[["i"]]
      j = sublist[["j"]]
      
      ESAC_fast_Kerr[i,j,y] = ESAC_fast_Kerr[i,j,y] + abs(sublist[["ESAC_fast_K"]])/N
      ESAC_fast_np_Kerr[i,j,y] = ESAC_fast_np_Kerr[i,j,y] + abs(sublist[["ESAC_fast_np_K"]])/N
      ESAC_Kerr[i,j,y] = ESAC_Kerr[i,j,y] + abs(sublist[["ESAC_K"]])/N
      ESAC_emp_Kerr[i,j,y] = ESAC_emp_Kerr[i,j,y] + abs(sublist[["ESAC_emp_K"]])/N
      ESAC_long_Kerr[i,j,y] = ESAC_long_Kerr[i,j,y] + abs(sublist[["ESAC_long_K"]])/N
      ESAC_long_emp_Kerr[i,j,y] = ESAC_long_emp_Kerr[i,j,y] + abs(sublist[["ESAC_long_emp_K"]])/N
      
      
      ESAC_fast_time[i,j,y] = ESAC_fast_time[i,j,y] + sublist[["ESAC_fast_time"]]/N
      ESAC_fast_np_time[i,j,y] = ESAC_fast_np_time[i,j,y] + sublist[["ESAC_fast_np_time"]]/N
      ESAC_time[i,j,y] = ESAC_time[i,j,y] + sublist[["ESAC_time"]]/N
      ESAC_emp_time[i,j,y] = ESAC_emp_time[i,j,y] + sublist[["ESAC_emp_time"]]/N
      ESAC_long_time[i,j,y] = ESAC_long_time[i,j,y] + sublist[["ESAC_long_time"]]/N
      ESAC_long_emp_time[i,j,y] = ESAC_long_emp_time[i,j,y] + sublist[["ESAC_long_emp_time"]]/N
      
      
      # ESAC_fast_ari[i,j,y] = ESAC_fast_ari[i,j,y] + sublist[["ESAC_fast_ari"]]/N
      # ESAC_fast_np_ari[i,j,y] = ESAC_fast_np_ari[i,j,y] + sublist[["ESAC_fast_np_ari"]]/N
      # ESAC_ari[i,j,y] = ESAC_ari[i,j,y] + sublist[["ESAC_ari"]]/N
      # ESAC_emp_ari[i,j,y] = ESAC_emp_ari[i,j,y] + sublist[["ESAC_emp_ari"]]/N
      # ESAC_long_ari[i,j,y] = ESAC_long_ari[i,j,y] + sublist[["ESAC_long_ari"]]/N
      # ESAC_long_emp_ari[i,j,y] = ESAC_long_emp_ari[i,j,y] + sublist[["ESAC_long_emp_ari"]]/N
      
      
     
      
      if(y!= 1){
        ESAC_fast_hausd[i,j,y] = ESAC_fast_hausd[i,j,y] + sublist[["ESAC_fast_hausd"]]/N
        ESAC_fast_np_hausd[i,j,y] = ESAC_fast_np_hausd[i,j,y] + sublist[["ESAC_fast_np_hausd"]]/N
        ESAC_hausd[i,j,y] = ESAC_hausd[i,j,y] + sublist[["ESAC_hausd"]]/N
        ESAC_emp_hausd[i,j,y] = ESAC_emp_hausd[i,j,y] + sublist[["ESAC_emp_hausd"]]/N
        ESAC_long_hausd[i,j,y] = ESAC_long_hausd[i,j,y] + sublist[["ESAC_long_hausd"]]/N
        ESAC_long_emp_hausd[i,j,y] = ESAC_long_emp_hausd[i,j,y] + sublist[["ESAC_long_emp_hausd"]]/N
      }
      
      
      
    }
    
  }
  
  
}

if(save){
  saveRDS(result, file=sprintf("%s/result.RDA", savedir))
  saveRDS(ESAC_fast_hausd, file=sprintf("%s/ESAC_fast_haus.RDA", savedir))
  saveRDS(ESAC_fast_time, file=sprintf("%s/ESAC_fast_time.RDA", savedir))
  saveRDS(ESAC_fast_Kerr, file=sprintf("%s/ESAC_fast_Kerr.RDA", savedir))
  saveRDS(ESAC_fast_ari, file=sprintf("%s/ESAC_fast_ari.RDA", savedir))
  
  saveRDS(ESAC_fast_np_hausd, file=sprintf("%s/ESAC_fast_np_haus.RDA", savedir))
  saveRDS(ESAC_fast_np_time, file=sprintf("%s/ESAC_fast_np_time.RDA", savedir))
  saveRDS(ESAC_fast_np_Kerr, file=sprintf("%s/ESAC_fast_np_Kerr.RDA", savedir))
  saveRDS(ESAC_fast_np_ari, file=sprintf("%s/ESAC_fast_np_ari.RDA", savedir))
  
  saveRDS(ESAC_hausd, file=sprintf("%s/ESAC_haus.RDA", savedir))
  saveRDS(ESAC_time, file=sprintf("%s/ESAC_time.RDA", savedir))
  saveRDS(ESAC_Kerr, file=sprintf("%s/ESAC_Kerr.RDA", savedir))
  saveRDS(ESAC_ari, file=sprintf("%s/ESAC_ari.RDA", savedir))
  
  saveRDS(ESAC_long_hausd, file=sprintf("%s/ESAC_long_haus.RDA", savedir))
  saveRDS(ESAC_long_time, file=sprintf("%s/ESAC_long_time.RDA", savedir))
  saveRDS(ESAC_long_Kerr, file=sprintf("%s/ESAC_long_Kerr.RDA", savedir))
  saveRDS(ESAC_long_ari, file=sprintf("%s/ESAC_long_ari.RDA", savedir))
  
  saveRDS(ESAC_emp_hausd, file=sprintf("%s/ESAC_emp_haus.RDA", savedir))
  saveRDS(ESAC_emp_time, file=sprintf("%s/ESAC_emp_time.RDA", savedir))
  saveRDS(ESAC_emp_Kerr, file=sprintf("%s/ESAC_emp_Kerr.RDA", savedir))
  saveRDS(ESAC_emp_ari, file=sprintf("%s/ESAC_emp_ari.RDA", savedir))
  
  saveRDS(ESAC_long_emp_hausd, file=sprintf("%s/ESAC_long_emp_haus.RDA", savedir))
  saveRDS(ESAC_long_emp_time, file=sprintf("%s/ESAC_long_emp_time.RDA", savedir))
  saveRDS(ESAC_long_emp_Kerr, file=sprintf("%s/ESAC_long_emp_Kerr.RDA", savedir))
  saveRDS(ESAC_long_emp_ari, file=sprintf("%s/ESAC_long_emp_ari.RDA", savedir))

  
  
  infofile<-file(sprintf("%s/parameters.txt", savedir))
  writeLines(c(sprintf("N = %d", N),
               sprintf("n = %d", n),
               sprintf("p = %d", p),
               sprintf("Sparse constant = %f", sparse_const), 
               sprintf("Dense constant = %f", dense_const), 
               sprintf("Rescale variance = %d", as.integer(rescale_variance))),
             infofile)
  close(infofile)
}

# creating table:
if(save){
  # output latex table
  printlines = c("%%REMEMBER to use package \\usepackage{rotating}!!",
                 " \\begin{table}[H] \\centering",
                 "\\caption{Multiple changepoints}",
                 "\\label{}",
                 "\\small",
                 "\\begin{adjustbox}{scale=0.55,center}",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} cccc|cccc|cccc|cccc|cccc}",
                 "\\hline", 
                 "\\multicolumn{4}{c|}{Parameters} & \\multicolumn{4}{c|}{Hausdorff distance} &\\multicolumn{4}{c|}{$\\left | \\widehat{K}-K \\right |$} &\\multicolumn{4}{c|}{ARI} &\\multicolumn{4}{c}{Time in miliseconds} \\\\ \\hline ",
                 "$n$ & $p$ & Sparsity & K & \\text{ESAC A} & \\text{ESAC B} & \\text{ ESAC C} & \\text{ESAC D}  & \\text{ESAC A} & \\text{ESAC B} & \\text{ESAC C} & \\text{ESAC D}  &\\text{ESAC A} & \\text{ESAC B} & \\text{ESAC C} & \\text{ESAC D} &\\text{ESAC A} & \\text{ESAC B} & \\text{ESAC C} & \\text{ESAC D} \\\\", 
                 "\\hline \\")
  
  
  for (i in 1:length(ns)) {
    
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:numconfig) {
        conf = config(y, n, p)
        etas = conf[[1]]
        #sparsities = conf[[2]]
        #mus = conf[[3]]
        sparsity = conf[[4]]
        if(is.null(sparsity)){
          sparsity="-"
        }
        
       
        string = sprintf("%d & %d & %s & %d ", n, p, sparsity, length(etas))
        
        
        if(y==1){
          for (t in 1:4) {
            string = sprintf("%s & - ", string)
            
          }
        }
        
        else{
          res = round(c(ESAC_fast_hausd[i,j,y], ESAC_fast_np_hausd[i,j,y], ESAC_hausd[i,j,y], ESAC_long_hausd[i,j,y]),digits=3)
          minind = (res==min(na.omit(res)))
          res = c(ESAC_fast_hausd[i,j,y], ESAC_fast_np_hausd[i,j,y], ESAC_hausd[i,j,y], ESAC_long_hausd[i,j,y])
          
          for (t in 1:length(res)) {
            if(is.na(res[t])){
              string = sprintf("%s & - ", string)
            }
            else if(minind[t]){
              string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
            }else{
              string = sprintf("%s & %.3f", string, res[t])
            }
          }
        }
        
        res = round(c(ESAC_fast_Kerr[i,j,y], ESAC_fast_np_Kerr[i,j,y], ESAC_Kerr[i,j,y], ESAC_long_Kerr[i,j,y]),digits=3)
        minind = (abs(res)==min(abs(res)))
        res = c(ESAC_fast_Kerr[i,j,y], ESAC_fast_np_Kerr[i,j,y], ESAC_Kerr[i,j,y], ESAC_long_Kerr[i,j,y])
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        res = round(c(ESAC_fast_ari[i,j,y], ESAC_fast_np_ari[i,j,y], ESAC_ari[i,j,y], ESAC_long_ari[i,j,y]),digits=3)
        minind = (abs(res)==max(abs(res)))
        res = c(ESAC_fast_ari[i,j,y], ESAC_fast_np_ari[i,j,y], ESAC_ari[i,j,y], ESAC_long_ari[i,j,y])
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        
        res = round(1000*c(ESAC_fast_time[i,j,y], ESAC_fast_np_time[i,j,y], ESAC_time[i,j,y], ESAC_long_time[i,j,y]),digits=3)
        minind = (res==min(res))
        res = 1000*c(ESAC_fast_time[i,j,y], ESAC_fast_np_time[i,j,y], ESAC_time[i,j,y], ESAC_long_time[i,j,y])
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        string = sprintf("%s \\\\", string)
        printlines = c(printlines, string)
      }
    }
  }
  
  printlines = c(printlines, c("\\hline \\\\[-1.8ex]",
                               "\\end{tabular}",
                               "\\end{adjustbox}",
                               "\\end{table}"))
  texfile<-file(sprintf("%s/table.tex", savedir))
  writeLines(printlines, texfile)
  close(texfile)
  
}

