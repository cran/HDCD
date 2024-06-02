# This code runs the simulations performed in Section 4.3 of Moen et al. (2024),	arXiv:2306.04702
# Please note that the code for SUBSET must be downloaded from Github, while the code for
# the method of Kaul et al (2021, Electronic Journal of Statistics) must be obtained from the first
# author of that paper. If SUBSET or the method of Kaul et al are not available,
# please set subset_included = FALSE and kaul_included = FALSE below.

# load libraries

#install.packages("InspectChangepoint", repos = "http:cran.us.r-project.org")
#install.packages("hdbinseg",repos = "http:cran.us.r-project.org")
#install.packages("/mn/sarpanitu/ansatte-u2/pamoen/project_inspect/simulations_nov_22/ESAC", repos=NULL, type="source")


library(doSNOW)
library(HDCD)
library(foreach)
library(MASS)
library(InspectChangepoint)
## same magnitude in all coords



# source SUBSET_normal.R from https://github.com/SOTickle/SUBSET
subset_included = TRUE
subset_path = "... fill in/SUBSET ... "
kaul_included =TRUE
kaul_path = "... fill in ... "



if(subset_included){
  source(sprintf("%s/SUBSET_normal.R", subset_path))
}


# source code for the method of Kaul et al (2021, EJS).
# Contact Abhishek <abhishek dot kaul at wsu.edu> for the source code
if(kaul_included){
  setwd(kaul_path)
  #unpack libraries
  source("unpack_libraries.R")
  source("K20/K20_bseg.R")
  source("alg1.R")
}



# saving options:
maindir = "... fill in ... "
dateandtime = gsub(" ", "--",as.character(Sys.time()))
dateandtime = gsub(":", ".", dateandtime)
savedir = file.path(maindir, dateandtime)

save = TRUE

if(save){
  dir.create(savedir, showWarnings = FALSE)
  savedir = file.path(maindir, sprintf("%s/single_misspec",dateandtime))
  dir.create(savedir, showWarnings = FALSE)

}



N = 50
num_cores = 12
sparse_const = 3
dense_const = 3
set.seed(1996)
rescale_variance = TRUE

even_spread = TRUE
ns = c(10)
ps = c(10)
#ns = c(200) # uncomment to get the same value of n as in the article
#ps = c(200) # uncomment to get the same value of p as in the article


n = ns
p = ps

etas = round(n/5)

simsbsN = 1000
simsbstoln = 1/simsbsN

#for sbs:
set.seed(1996)
pis = matrix(nrow = length(ns), ncol =length(ps))
for (i in 1:length(ns)) {
  for (j in 1:length(ps)) {
    pis[i,j] =  single_SBS_calibrate(ns[i],ps[j],simsbsN,simsbstoln,rescale_variance = rescale_variance,debug=FALSE)

  }
}

numconfig = 12




config = function(i,h,n,p){
  # h = 1 sparse
  # h = 2 dense
  
  noise = NULL
  mus = matrix(0, nrow=p, ncol=n)
  sparsities = c()
  sparsity = c()
  Delta = etas
  k = 1
  phi = 1

  if(h==1){
    sparsity = "Sparse"
    Delta = etas
    sparsities = sample(1:floor(sqrt(p*log(n))),1)
    k = sparsities[j]
    phi = sparse_const/sqrt(Delta)*sqrt(((k*log(exp(1)*p*log(n)/k^2)+ log(n))))

  }else{
    sparsity = "Dense"
    Delta = etas
    sparsities = sample((ceiling(sqrt(p*log(n)))):p,1)
    k = sparsities[j]
    phi = dense_const/sqrt(Delta)*(p*log(n))^(1/4)
  }

  coords = sample(1:p, sparsities)
  diff = sample(c(-1,1), sparsities, replace=TRUE)
  diff = diff/norm(diff, type="2")
  mus[coords, (etas+1):n] = mus[coords, (etas+1):n] +diff*phi

  if(i==1){
    # following model assumptions
    noise = matrix(rnorm(n*p), nrow=p, ncol=n)


  }

  if(i==2){
    #uniform
    noise = matrix(runif(n*p,min=-sqrt(3), max=sqrt(3)), nrow=p, ncol=n)


  }

  if(i==3){
    #student t
    noise = matrix(rt(n*p,df=3), nrow=p, ncol=n)/sqrt(3)


  }

  if(i==4){
    #heavy tailed
    noise = matrix(rt(n*p,df=10), nrow=p, ncol=n)*sqrt(8/10)

  }


  if(i==5){
    #cs, loc 1
    rho = 0.1
    noise = t(mvrnorm(n = n, rep(0,p), Sigma = outer(1:p,1:p,function(a,b){rho^(abs(a-b))}), tol = 1e-6))


  }

  if(i==6){
    #cs, loc 2
    rho = 0.4
    noise = t(mvrnorm(n = n, rep(0,p), Sigma = outer(1:p,1:p,function(a,b){rho^(abs(a-b))}), tol = 1e-6))


  }

  if(i==7){
    #cs 1
    rho = 0.1
    noise = t(mvrnorm(n, mu = rep(0,p), Sigma = matrix(rho,p,p)+diag(p)*(1-rho)))

  }

  if(i==8){
    #cs 2
    rho = 0.4
    noise = t(mvrnorm(n, mu = rep(0,p), Sigma = matrix(rho,p,p)+diag(p)*(1-rho)))

  }

  if(i==9){
    #temp 1
    rho = 0.1
    noise = matrix(0,p,n)
    noise[,1] = rnorm(p)
    for (j in 1:p){
      for (t in 2:n) noise[j,t] = noise[j,t-1]*sqrt(rho) + rnorm(1)*sqrt(1-rho)
    }


  }

  if(i==10){
    #temp 2
    rho = 0.4
    noise = matrix(0,p,n)
    noise[,1] = rnorm(p)
    for (j in 1:p){
      for (t in 2:n) noise[j,t] = noise[j,t-1]*sqrt(rho) + rnorm(1)*sqrt(1-rho)
    }


  }

  if(i==11){
    #async
    noise = matrix(rnorm(n*p), nrow=p, ncol=n)
    L = floor(Delta/2)

    for (zz in 1:length(coords)) {
      start = sample(seq(etas-L,etas+L ),1)
      coord = coords[zz]
      hh = mus[coord,n]
      mus[coord,] = 0
      mus[coord, (start+1):n] = hh
    }



  }

  if(i==12){
    #gradual
    noise = matrix(rnorm(n*p), nrow=p, ncol=n)
    L = floor(Delta/2)
    for (zz in 1:length(coords)) {
      start = sample(seq(etas-L,etas+L ),1)
      coord = coords[zz]
      hh = mus[coord,n]
      mus[coord,] = 0
      mus[coord, (etas+L+1):n] = hh
      start = etas - L+1
      stop = etas+L+1
      mus[coord, start:stop] = hh*(1:(stop-start+1))/(stop-start+1)
    }



  }

  X = mus + noise


  return(list(etas, sparsities,mus, X,sparsity,noise))
}



cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result_misspec = foreach(z = 1:N,.options.snow = opts) %dopar% {
  rez = list()
  set.seed(2*z)
  library(HDCD)
  library(hdbinseg)
  library(MASS)
  library(InspectChangepoint)

  if(subset_included){
    source(sprintf("%s/SUBSET_normal.R", subset_path))
  }
  if(kaul_included){
    setwd(kaul_path)
    #unpack libraries
    source("unpack_libraries.R")
    source("K20/K20_bseg.R")
    source("alg1.R")
  }

  counter = 1

  n = 1
  p = 1
  i = 1
  j=1
  n = ns[i]

  p = ps[j]

  for(h in 1:2){
    for (y in 1:numconfig) {
    
      conf = config(y, h,n, p)
      eta = conf[[1]]
      sparsities = conf[[2]]
      mus = conf[[3]]
      X = conf[[4]]
      sparsity = conf[[5]]
     

      rezi = list()
      
      rezi[["y"]] = y
      rezi[["h"]] = h


      # DC
      a = proc.time()
      res_dc = dcbs.alg(X[,],phi = -1,height=1, thr = 0,cp.type=1,temporal=FALSE )
      b=proc.time()

     
      rezi["dc_time"] = (b-a)[3]
      if(is.null(res_dc$ecp)){
        rezi["dc_res"] = round(n/2)
      }else{
        rezi["dc_res"]= res_dc$ecp
      }
      
      #then sbs
      a = proc.time()
      res_sbs = sbs.alg(X[,],height=1, thr = rep(sqrt(pis[i,j]),p),cp.type=1,temporal=FALSE )
      b=proc.time()

      rezi["sbs_time"] = (b-a)[3]
      if(!is.numeric(res_sbs$ecp)){
        if(!rescale_variance){
          ress = single_SBS(X[,],pis[i,j] )
        }else{
          ress = single_SBS(rescale.variance(X),pis[i,j] )
        }

        rezi["sbs_res"] = ress$pos

      }else{
        rezi["sbs_res"] = res_sbs$ecp
      }



      # now rescale variance
      X = matrix(rescale_variance(X)$X,byrow = FALSE, ncol = n, nrow = p)
      
      
      
      lambda = sqrt(log(p*log(n))/2)

      a = proc.time()
      res_inspect = single_Inspect(X[,], lambda)
      b=proc.time()
      rezi["inspect_time"] = (b-a)[3]
      rezi["inspect_res"]= res_inspect$pos

      # then ESAC
      a = proc.time()
      res_ESAC  = single_ESAC(X[,], 1.5,1, debug =FALSE)
      b=proc.time()
      rezi["ESAC_time"] = (b-a)[3]
      rezi["ESAC_res"] = res_ESAC$pos
      rezi["ESAC_s"] = res_ESAC$s

     



      # then subset
      if(subset_included){
        a = proc.time()
        res_subset  = SUBSET.normal(X[,])
        b=proc.time()
        rezi["subset_time"] = (b-a)[3]

        if(is.null(res_subset$cpt)){
          rezi["subset_res"] = round(n/2)
        }else{
          rezi["subset_res"] = res_subset$cpt
        }
      }

      if(kaul_included){
        res_kaul21 = K20(t(X))$cp
        if(length(res_kaul21)>1){
          res_kaul21 = 1
        }else if(res_kaul21 == "no change"){
          res_kaul21 = 1
        }

        rezi["kaul21_time"] =(b-a)[3]

        rezi["kaul21_res"] = res_kaul21
      }


      rezi["z"] = z
      rezi["h"] = h
      rezi["y"] = y
      rezi["k"] = sparsities
      rezi["eta"] = eta

      rez[[counter]] = rezi
      counter = counter+1



    }

  }
  rez
}
close(pb)
stopCluster(cl)

# gathering results into arrays:
{
  inspect_res = array(NA, dim = c(2,numconfig,N) )
  inspect_time = array(0, dim= c(2,numconfig))
  inspect_mse = array(0, dim= c(2,numconfig))
  ESAC_res = array(NA, dim = c(2,numconfig,N) )
  ESAC_s = array(NA, dim = c(2,numconfig,N) )
  ESAC_time = array(0, dim= c(2,numconfig))
  ESAC_mse = array(0, dim= c(2,numconfig))
  scan_res = array(NA, dim = c(2,numconfig,N) )
  scan_s = array(NA, dim = c(2,numconfig,N) )
  scan_time = array(0, dim= c(2,numconfig))
  scan_mse = array(NA, dim= c(2,numconfig))
  sbs_res = array(NA, dim = c(2,numconfig,N) )
  sbs_time = array(0, dim= c(2,numconfig))
  sbs_mse = array(0, dim= c(2,numconfig))
  subset_res = array(NA, dim = c(2,numconfig,N) )
  subset_time = array(0, dim= c(2,numconfig))
  subset_mse = array(0, dim= c(2,numconfig))
  dc_res = array(NA, dim = c(2,numconfig,N) )
  dc_time = array(0, dim= c(2,numconfig))
  dc_mse = array(0, dim= c(2,numconfig))
  kaul21_res = array(NA, dim = c(2,numconfig,N) )
  kaul21_time = array(0, dim= c(2,numconfig))
  kaul21_mse = array(0, dim= c(2,numconfig))

  for (z in 1:N) {
    list = result_misspec[[z]]
    len = length(list)

    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      h = sublist[["h"]]
      eta = sublist[["eta"]]

      inspect_res[h,y,z] = sublist["inspect_res"][[1]]
      inspect_time[h,y] = inspect_time[h,y] +  sublist["inspect_time"][[1]]/N
      inspect_mse[h,y] = inspect_mse[h,y] + (inspect_res[h,y,z] - eta)^2/N

      ESAC_res[h,y,z] = sublist["ESAC_res"][[1]]
      ESAC_time[h,y] = ESAC_time[h,y] + sublist["ESAC_time"][[1]]/N
      ESAC_mse[h,y] = ESAC_mse[h,y] + (ESAC_res[h,y,z] - eta)^2/N

      sbs_res[h,y,z] = sublist["sbs_res"][[1]]
      sbs_time[h,y] = sbs_time[h,y] +  sublist["sbs_time"][[1]]/N
      sbs_mse[h,y] = sbs_mse[h,y] + (sbs_res[h,y,z] - eta)^2/N

      if(subset_included){
        subset_res[h,y,z] = sublist["subset_res"][[1]]
        subset_time[h,y] = subset_time[h,y] +  sublist["subset_time"][[1]]/N
        subset_mse[h,y] = subset_mse[h,y] + (subset_res[h,y,z] - eta)^2/N
      }


      dc_res[h,y,z] = sublist["dc_res"][[1]]
      dc_time[h,y] = dc_time[h,y] + sublist["dc_time"][[1]]/N
      dc_mse[h,y] = dc_mse[h,y] + (dc_res[h,y,z] - eta)^2/N

      if(kaul_included){
        kaul21_res[h,y,z] = sublist["kaul21_res"][[1]]
        kaul21_time[h,y] = kaul21_time[h,y] + sublist["kaul21_time"][[1]]/N
        kaul21_mse[h,y] = kaul21_mse[h,y] + (kaul21_res[h,y,z] - eta)^2/N
      }


    }
  }
}





#saving:
if(save){
  saveRDS(inspect_mse, file=sprintf("%s/inspect_mse.RDA", savedir))
  saveRDS(inspect_time, file=sprintf("%s/inspect_time.RDA", savedir))
  saveRDS(inspect_res, file=sprintf("%s/inspect_res.RDA", savedir))
  saveRDS(ESAC_res, file=sprintf("%s/ESAC_res.RDA", savedir))
  saveRDS(ESAC_s, file=sprintf("%s/ESAC_s.RDA", savedir))
  saveRDS(ESAC_time, file=sprintf("%s/ESAC_time.RDA", savedir))
  saveRDS(ESAC_mse, file=sprintf("%s/ESAC_mse.RDA", savedir))
  saveRDS(scan_res, file=sprintf("%s/scan_res.RDA", savedir))
  saveRDS(scan_s, file=sprintf("%s/scan_s.RDA", savedir))
  saveRDS(scan_time, file=sprintf("%s/scan_time.RDA", savedir))
  saveRDS(scan_mse, file=sprintf("%s/scan_mse.RDA", savedir))
  saveRDS(sbs_mse, file=sprintf("%s/sbs_mse.RDA", savedir))
  saveRDS(sbs_res, file=sprintf("%s/sbs_res.RDA", savedir))
  saveRDS(sbs_time, file=sprintf("%s/sbs_time.RDA", savedir))
  saveRDS(subset_mse, file=sprintf("%s/subset_mse.RDA", savedir))
  saveRDS(subset_res, file=sprintf("%s/subset_res.RDA", savedir))
  saveRDS(subset_time, file=sprintf("%s/subset_time.RDA", savedir))
  saveRDS(dc_mse, file=sprintf("%s/dc_mse.RDA", savedir))
  saveRDS(dc_res, file=sprintf("%s/dc_res.RDA", savedir))
  saveRDS(dc_time, file=sprintf("%s/dc_time.RDA", savedir))
  saveRDS(kaul21_mse, file=sprintf("%s/kaul21_mse.RDA", savedir))
  saveRDS(kaul21_res, file=sprintf("%s/kaul21_res.RDA", savedir))
  saveRDS(kaul21_time, file=sprintf("%s/kaul21_time.RDA", savedir))

  infofile<-file(sprintf("%s/parameters.txt", savedir))
  writeLines(c(sprintf("N = %d", N),
               sprintf("ns = %s", paste(ns, sep=" ", collapse=" ")),
               sprintf("ps = %s", paste(ps, sep=" ", collapse=" ")),
               sprintf("Sparse constant = %f", sparse_const),
               sprintf("Dense constant = %f", dense_const),
               sprintf("Even spread = %d", even_spread)),
             infofile)
  close(infofile)
}


models = c("$M$",
           "$M_{\\text{Unif}}$",
           "$M_{t_{3}}$",
           "$M_{t_{10}}$",
           "$M_{\\text{cs, loc}}(\\rho = 0.1)$",
           "$M_{\\text{cs, loc}}(\\rho = 0.4)$",
           "$M_{\\text{cs}}(\\rho = 0.1)$",
           "$M_{\\text{cs}}(\\rho = 0.4)$",
           "$M_{\\text{AR}}(\\rho = 0.1)$",
           "$M_{\\text{AR}}(\\rho = 0.4)$",
           "$M_{\\text{async}}$",
           "$M_{\\text{grad}}$")

# creating table:
if(save){
  # output latex table
  printlines = c("\\begin{table}[H] \\centering",
                 "\\caption{Single change-point estimation under misspecified model}",
                 "\\label{tablesinglelocmisspec}",
                 "\\small",
                 "\\begin{adjustbox}{width=\\columnwidth}",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} cc|cccccc}",
                 "\\hline",
                 "\\multicolumn{2}{c|}{Parameters} & \\multicolumn{6}{c|}{MSE} \\\\ \\hline ",
                 "Model & Sparsity & \\text{ESAC} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET} & \\text{DC} & \\text{Kaul et. al 2021} \\\\",
                 "\\hline \\")



    for (y in 1:numconfig) {
      for(h in 1:2){
        if(h==1){
          dens = "Sparse"
        }else{
          dens = "Dense"
        }
        eta = round(0.2*n)
        
        string = sprintf("%s & %s", models[y], dens)

        res = round(c(ESAC_mse[h,y],inspect_mse[h,y], sbs_mse[h,y], subset_mse[h,y], dc_mse[h,y], kaul21_mse[h,y]),digits=1)
        minind = (res==min(res))

        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.1f} ", string, res[t])
          }else{
            string = sprintf("%s & %.1f", string, res[t])
          }
        }


        string = sprintf("%s \\\\", string)
        printlines = c(printlines, string)
      }
    }

 

  printlines = c(printlines, "\\hline \\multicolumn{2}{c|}{Average MSE}")

  res = round(c(mean(ESAC_mse),mean(inspect_mse), mean(sbs_mse), mean(subset_mse), mean(dc_mse), mean(kaul21_mse)),digits=1)
  res[is.na(res)] = Inf
  minind = (res==min(res))
  string =""
  for (t in 1:length(res)) {
    if(minind[t]){
      string = sprintf("%s & \\textbf{%.1f} ", string, res[t])
    }else{
      string = sprintf("%s & %.1f", string, res[t])
    }
  }
  string = sprintf("%s \\\\", string)
  printlines = c(printlines, string)


  printlines = c(printlines, c("\\hline \\\\[-1.8ex]",
                               "\\end{tabular}",
                               "\\end{adjustbox}",
                               "\\end{table}"))
  texfile<-file(sprintf("%s/table_even.tex", savedir))
  writeLines(printlines, texfile)
  close(texfile)

}
