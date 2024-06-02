# This code runs the simulations performed in Section 4.1 of Moen et al. (2024),	arXiv:2306.04702
# Please note that the code for SUBSET must be downloaded from Github (https://github.com/SOTickle/SUBSET), 
#while the code for the method of Kaul et al (2021, Electronic Journal of Statistics) 
# must be obtained from the first author of that paper. 
# If SUBSET or the method of Kaul et al are not available,
# please set subset_included = FALSE and kaul_included = FALSE below.


## specify if SUBSET and the method of Kaul et al (2021) are included:
subset_included = FALSE
kaul_included =FALSE

# if yes, specify the paths on your machine:
subset_path = "... fill in .../SUBSET"
kaul_path = "... fill in ..."

# Load libraries
library(doSNOW)
library(HDCD)
library(foreach)

# if running the code from Kaul, uncomment these if necessary:
#install.packages("MASS"); install.packages("rmutil")
#install.packages("parallel"); install.packages("doSNOW");
#install.packages("foreach");install.packages("InspectChangepoint")
#install.packages("pracma");


# source SUBSET_normal.R from https://github.com/SOTickle/SUBSET
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


## Saving options:

save = TRUE # if results should be saved to disk

# specify directory in which results should be saved:
maindir = "... fill in ... "
dateandtime = gsub(" ", "--",as.character(Sys.time()))
dateandtime = gsub(":", ".", dateandtime)
savedir = file.path(maindir, dateandtime)


# Creating subfolder with current time as name:
if(save){
  dir.create(savedir, showWarnings = FALSE)
  savedir = file.path(maindir, sprintf("%s/single",dateandtime))
  dir.create(savedir, showWarnings =FALSE )
  plotdir = file.path(maindir, sprintf("%s/single/plots",dateandtime))
  dir.create(plotdir, showWarnings =FALSE )

}



########### Parameter settings ###########

# Set your desired values of n and p here
ns = c(20)
ps = c(20)
#ns = c(200,500) # uncomment for the same values of n as in the article
#ps = c(100,1000,5000) # uncomment for the same values of p as in the article

kvals=4 #do not change this
totruns = length(ns)*length(ps)*kvals # do not change this


## Obtain threshold for SBS via MC simulation:
simsbsN = 1000 # number of MC simulations to obtain SBS threshold
simsbstoln = 1/simsbsN # quantile chosen for SBS; here 1/N
rescale_variance = TRUE
#for sbs:
set.seed(1996)
pis = matrix(nrow = length(ns), ncol =length(ps))
for (i in 1:length(ns)) {
  for (j in 1:length(ps)) {
    pis[i,j] =  single_SBS_calibrate(ns[i],ps[j],simsbsN,simsbstoln,rescale_variance = rescale_variance,debug=FALSE)

  }
}


# Set simulation parameters:
N = 1000 # number of MC simulations for the simulation study
num_cores = 6 # number of cores to use on the machine
sparse_const = 2.5 # leading constant for the signal strength in the sparse regime
dense_const = 2.5 # leading constant for the signal strength in the dense regime

# Choose if the change in the mean vector should be spread evenly among coordinates,
# or not
even_spread = TRUE



# this is for choosing the sparsity regime, do not change
kfunc = function(y, n, p){
  k=1
  if(y==1){
    k = 1
  }else if(y==2){
    k = ceiling(p^(1/3))
  }else if(y==3){
    k = ceiling(sqrt(p*log(n)))
  }else{
    k = p
  }
  return(k)
}


######### Simulation: #############

# Running loop in paralell
cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result = foreach(z = 1:N,.options.snow = opts) %dopar% {
#for (z in 1:N) {
  rez = list()
  counter = 1
  set.seed(z)
  library(HDCD)
  library(hdbinseg)
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


  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        rezi = list()
        # determine k:
        k = kfunc(y,n,p)

        mus = matrix(0, nrow=p, ncol=n)
        noise = matrix(rnorm(n*p), nrow=p, ncol=n)
        eta = round(0.2*n)


        diff = 0
        if(k<sqrt(p*log(n))){
          rootnorm = sparse_const/sqrt(eta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
          diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)

        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n))^(1/4)
          diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
        }
        if(!even_spread){
          diff = c(rnorm(k), rep(0,p-k))
          diff = diff/norm(diff, type="2")
          diff = diff*rootnorm
        }

        mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
        X = mus+noise


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
            ress = single_SBS(matrix(rescale_variance(X[,])$X,byrow = FALSE, ncol = n, nrow = p),pis[i,j] )
          }

          rezi["sbs_res"] = ress$pos

        }else{
          rezi["sbs_res"] = res_sbs$ecp
        }






        # now rescale variance manually, as the remaining methods to not
        # do this by default:

        X = matrix(rescale_variance(X)$X,byrow = FALSE, ncol = n, nrow = p)


        # now Inspect:
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

        # then Kaul et al (2021, Electronic Journal of Statistics)
        if(kaul_included){
          a = proc.time()
          res_kaul21 = K20(t(X))$cp
          b=proc.time()
          if(length(res_kaul21)>1){
            res_kaul21 = 1
          }else if(res_kaul21 == "no change"){
            res_kaul21 = 1
          }

          rezi["kaul21_time"] =(b-a)[3]

          rezi["kaul21_res"] = res_kaul21
        }

        # then SUBSET:
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

        # save run parameters:
        rezi["z"] = z
        rezi["i"] = i
        rezi["j"] = j
        rezi["y"] = y
        rezi["k"] = kfunc(y,n,p)
        rezi["eta"] = eta


        rez[[counter]] = rezi
        counter = counter+1

      }
    }
  }

  rez
}
close(pb)
stopCluster(cl)

# The results above are saved to a list. We unpack them into arrays here:
{
  inspect_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  inspect_time = array(0, dim= c(length(ns), length(ps), kvals))
  inspect_mse = array(0, dim= c(length(ns), length(ps), kvals))
  ESAC_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  ESAC_s = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  ESAC_time = array(0, dim= c(length(ns), length(ps), kvals))
  ESAC_mse = array(0, dim= c(length(ns), length(ps), kvals))
  scan_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  scan_s = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  scan_time = array(0, dim= c(length(ns), length(ps), kvals))
  scan_mse = array(NA, dim= c(length(ns), length(ps), kvals))
  sbs_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  sbs_time = array(0, dim= c(length(ns), length(ps), kvals))
  sbs_mse = array(0, dim= c(length(ns), length(ps), kvals))
  subset_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  subset_time = array(NA, dim= c(length(ns), length(ps), kvals))
  subset_mse = array(NA, dim= c(length(ns), length(ps), kvals))
  if(subset_included){
    subset_time = array(0, dim= c(length(ns), length(ps), kvals))
    subset_mse = array(0, dim= c(length(ns), length(ps), kvals))
  }
  dc_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  dc_time = array(0, dim= c(length(ns), length(ps), kvals))
  dc_mse = array(0, dim= c(length(ns), length(ps), kvals))

  kaul21_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  kaul21_time = array(NA, dim= c(length(ns), length(ps), kvals))
  kaul21_mse = array(NA, dim= c(length(ns), length(ps), kvals))
  if(kaul_included){
    kaul21_time = array(0, dim= c(length(ns), length(ps), kvals))
    kaul21_mse = array(0, dim= c(length(ns), length(ps), kvals))
  }

  for (z in 1:N) {
    list = result[[z]]
    len = length(list)

    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      i = sublist[["i"]]
      j = sublist[["j"]]
      eta = sublist[["eta"]]

      inspect_res[i,j,y,z] = sublist["inspect_res"][[1]]
      inspect_time[i,j,y] = inspect_time[i,j,y] +  sublist["inspect_time"][[1]]/N
      inspect_mse[i,j,y] = inspect_mse[i,j,y] + (inspect_res[i,j,y,z] - eta)^2/N

      ESAC_res[i,j,y,z] = sublist["ESAC_res"][[1]]
      ESAC_time[i,j,y] = ESAC_time[i,j,y] + sublist["ESAC_time"][[1]]/N
      ESAC_mse[i,j,y] = ESAC_mse[i,j,y] + (ESAC_res[i,j,y,z] - eta)^2/N

      sbs_res[i,j,y,z] = sublist["sbs_res"][[1]]
      sbs_time[i,j,y] = sbs_time[i,j,y] +  sublist["sbs_time"][[1]]/N
      sbs_mse[i,j,y] = sbs_mse[i,j,y] + (sbs_res[i,j,y,z] - eta)^2/N

      if(subset_included){
        subset_res[i,j,y,z] = sublist["subset_res"][[1]]
        subset_time[i,j,y] = subset_time[i,j,y] +  sublist["subset_time"][[1]]/N
        subset_mse[i,j,y] = subset_mse[i,j,y] + (subset_res[i,j,y,z] - eta)^2/N
      }

      dc_res[i,j,y,z] = sublist["dc_res"][[1]]
      dc_time[i,j,y] = dc_time[i,j,y] + sublist["dc_time"][[1]]/N
      dc_mse[i,j,y] = dc_mse[i,j,y] + (dc_res[i,j,y,z] - eta)^2/N

      if(kaul_included){
        kaul21_res[i,j,y,z] = sublist["kaul21_res"][[1]]
        kaul21_time[i,j,y] = kaul21_time[i,j,y] + sublist["kaul21_time"][[1]]/N
        kaul21_mse[i,j,y] = kaul21_mse[i,j,y] + (kaul21_res[i,j,y,z] - eta)^2/N
      }

    }
  }
}


# check the results:
inspect_mse
ESAC_mse
kaul21_mse
dc_mse
sbs_mse







#saving all simulation results:
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
  if(subset_included){
    saveRDS(subset_mse, file=sprintf("%s/subset_mse.RDA", savedir))
    saveRDS(subset_res, file=sprintf("%s/subset_res.RDA", savedir))
    saveRDS(subset_time, file=sprintf("%s/subset_time.RDA", savedir))
  }
  saveRDS(dc_mse, file=sprintf("%s/dc_mse.RDA", savedir))
  saveRDS(dc_res, file=sprintf("%s/dc_res.RDA", savedir))
  saveRDS(dc_time, file=sprintf("%s/dc_time.RDA", savedir))
  if(kaul_included){
    saveRDS(kaul21_mse, file=sprintf("%s/kaul21_mse.RDA", savedir))
    saveRDS(kaul21_res, file=sprintf("%s/kaul21_res.RDA", savedir))
    saveRDS(kaul21_time, file=sprintf("%s/kaul21_time.RDA", savedir))
  }

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


# Now we create two tables; the table containing the MSEs, and the table containing
# run times. These are saved to savedir, which is a subdirectory of maindir, which the
# user specifies

if(save){
  # output latex table, the first for MSE and the second for timings

  printlines1 = c("\\begin{table}[H] \\centering",
                 "\\caption{Single change-point estimation MSE}",
                 "\\label{tablesinglelocmse}",
                 "\\small",
                 "\\begin{adjustbox}{width=\\columnwidth}",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|cccccc}",
                 "\\hline",
                 "\\multicolumn{5}{c|}{Parameters} & \\multicolumn{6}{c}{Mean Squared Error}  \\\\ \\hline ",
                 "$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{ESAC} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET} & \\text{DC} & \\text{Kaul et al.} \\\\",
                 "\\hline \\")

  printlines2 = c("\\begin{table}[H] \\centering",
                  "\\caption{Single change-point estimation run time}",
                  "\\label{tablesinglelocruntime}",
                  "\\small",
                  "\\begin{adjustbox}{width=\\columnwidth}",
                  "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|cccccc}",
                  "\\hline",
                  "\\multicolumn{5}{c|}{Parameters}  &\\multicolumn{6}{c}{Time in milliseconds} \\\\ \\hline ",
                  "$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{ESAC} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET} & \\text{DC} & \\text{Kaul et al.} \\\\",
                  "\\hline \\")

  for (i in 1:length(ns)) {

    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        k = kfunc(y,n,p)
        eta = round(0.2*n)
        rootnorm = 0
        if(k<sqrt(p*log(n))){
          rootnorm = sparse_const/sqrt(eta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n))^(1/4)
        }
        string1 = sprintf("%d & %d & %d & %d & %.2f", n, p, k, eta, rootnorm)
        string2 = sprintf("%d & %d & %d & %d & %.2f", n, p, k, eta, rootnorm)
        # create row in table for MSEs:
        res = round(c(ESAC_mse[i,j,y],inspect_mse[i,j,y], sbs_mse[i,j,y], subset_mse[i,j,y], dc_mse[i,j,y], kaul21_mse[i,j,y]),digits=1)
        res[is.na(res)] = Inf
        minind = (res==min(res))

        for (t in 1:length(res)) {
          if(minind[t]){
            string1 = sprintf("%s & \\textbf{%.1f} ", string1, res[t])
          }else{
            string1 = sprintf("%s & %.1f", string1, res[t])
          }
        }

        string1 = sprintf("%s \\\\", string1)
        printlines1 = c(printlines1, string1)

        # now the run time
        res = round(c(ESAC_time[i,j,y],inspect_time[i,j,y], sbs_time[i,j,y], subset_time[i,j,y], dc_time[i,j,y], kaul21_time[i,j,y])*1000,digits=1)
        res[is.na(res)] = Inf
        minind = (res==min(res))

        for (t in 1:length(res)) {
          if(minind[t]){
            string2 = sprintf("%s & \\textbf{%.1f} ", string2, res[t])
          }else{
            string2 = sprintf("%s & %.1f", string2, res[t])
          }
        }
        string2 = sprintf("%s \\\\", string2)
        printlines2 = c(printlines2, string2)
      }
    }
  }
  printlines1 = c(printlines1, "\\hline \\multicolumn{5}{c|}{Average MSE}")
  printlines2 = c(printlines2, "\\hline \\multicolumn{5}{c|}{Average run time}")

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
  printlines1 = c(printlines1, string)

  res = round(1000*c(mean(ESAC_time),mean(inspect_time), mean(sbs_time), mean(subset_time), mean(dc_time), mean(kaul21_time)),digits=1)
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
  printlines2 = c(printlines2, string)

  printlines1 = c(printlines1, c("\\hline \\\\[-1.8ex]",
                               "\\end{tabular}",
                               "\\end{adjustbox}",
                               "\\end{table}"))
  printlines2 = c(printlines2, c("\\hline \\\\[-1.8ex]",
                                "\\end{tabular}",
                                "\\end{adjustbox}",
                                "\\end{table}"))
  texfile<-file(sprintf("%s/table_mse.tex", savedir))
  writeLines(printlines1, texfile)
  close(texfile)

  texfile<-file(sprintf("%s/table_timings.tex", savedir))
  writeLines(printlines2, texfile)
  close(texfile)

}



