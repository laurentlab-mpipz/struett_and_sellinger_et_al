Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    OMP_THREAD_LIMIT = "1"
    )

Sys.getenv(c(
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "OMP_THREAD_LIMIT"
))

options(error = quote({dump.frames(to.file=TRUE); q()}))

get_number_of_chromosomes_from_path_to_mhs <- function(
  path_to_mhs,
  snake_obj, 
  from_params=T) {
  if (from_params) {
      num_chr <- as.numeric(as.character(snake_obj@config$chromosome_number))
    } else {
      mhs <- read.table(path_to_mhs)
      num_chr <- length(unique(mhs[,1]))
    }
  return(num_chr)
}

sink(snakemake@log[["esmc_log"]])

unix::rlimit_as(7e10,1e11)

# save image for debugging
#save.image(); stop("dadada+++++++++++++++++++")
# setwd("../..")
# load(".RData")

# parameters
inference_model = snakemake@params[["inference_model"]]
path_to_package = snakemake@input[["esmc_package"]]
M=as.numeric(snakemake@wildcards[["sample_size"]]) # Haploid sequence number
NC=get_number_of_chromosomes_from_path_to_mhs(snakemake@input[["mhs"]], snakemake) # Number of chrosome/scaffold simulated
mu=as.numeric(snakemake@wildcards[["mutation_rate"]]) # Mutation rate
r=as.numeric(snakemake@wildcards[["recombination_rate"]]) # Recombination rate
n=snakemake@params[["number_of_hidden_states"]] # Number of hidden state
rho=r/mu
recent_sigma = snakemake@params[["recent_selfing_rate"]]
ancient_sigma = snakemake@params[["ancient_selfing_rate"]]

# install package
load_package <- tryCatch(
  {
    library(eSMC2)
  },
  error=function(cond) {
    message("need to install eSMC2 package")
    message("Here's the original error message:")
    message(cond)
    
    message("going to install eSMC2 via devtools")
    library(devtools)
    devtools::install_local(path_to_package)
    library(eSMC2)
    
    # Choose a return value in case of error
    return(NA)
  }
)
message("loaded eSMC2 package")


# save.image(); stop("dadada+++++++++++++++++++")

# run esmc on the multihetsep input
run_esmc <- function(path_to_mhs, M, mu, r, n, rho, NC, inference_model,
                     recent_sigma, ancient_sigma, outfilerds) {
  # read multihetsep file
  filepath = dirname(path_to_mhs)
  filename = basename(path_to_mhs)
  O_total=Get_real_data(filepath,M,filename,"\t")

  my_big_window <- 4
  my_pop_vect <- rep(2, floor(n/2)); stopifnot(sum(my_pop_vect) <= n); while (sum(my_pop_vect)!=n) my_pop_vect <- c(my_pop_vect, 1)

  total = list()
  s_ex = list()
  mu_ex = list()
  Tc = list()
  
  stopifnot(inference_model %in% c(
    "Free", "Constant", "GivenTransition", "OneTransition"))
  
  message("inference model ", inference_model)
  # calling the function
  if (inference_model == "Free") {
    test_=teSMC(n=n,rho=rho,model=list("Free"),O_total,maxit =20,
                BoxB=list(c(0.05,1)),BoxP=c(3,3),Boxr=list(c(1,1)),
                Boxs=list(c(0,0.99)),Constant_Pop=F,estimate="SF",NC=NC,
                Big_Window = my_big_window, pop_vect=my_pop_vect)
    total[[(1+length(total))]]=as.numeric(test_$Xi)
    s_ex[[(1+length(s_ex))]]=test_$sigma
    mu_ex[[(1+length(mu_ex))]]=mean(test_$mu)
    Tc[[(1+length(Tc))]]=test_$Tc
  } else if (inference_model == "Constant") {
    test_=teSMC(n=n,rho=rho,model=list("Constant"),O_total,maxit =20,
                BoxB=list(c(0.05,1)),BoxP=c(3,3),Boxr=list(c(1,1)),
                Boxs=list(c(0,0.99)),Constant_Pop=F,estimate="SF",NC=NC,
                Big_Window = my_big_window, pop_vect=my_pop_vect)
    total[[(1+length(total))]]=as.numeric(test_$Xi)
    s_ex[[(1+length(s_ex))]]=rep(test_$sigma,n)
    mu_ex[[(1+length(mu_ex))]]=mean(test_$mu)
    Tc[[(1+length(Tc))]]=test_$Tc
  } else if (inference_model == "GivenTransition") {
    current_selfing_value=recent_sigma # Give current selfing value 
    ancestral_selfing_value=ancient_sigma  # Give ancestral selfing value
    given= c(current_selfing_value,ancestral_selfing_value)
    
    test_=teSMC(n=n,rho=rho,model=list("Given Transition"),O_total,maxit =20,
                BoxB=list(c(0.05,1)),BoxP=c(3,3),Boxr=list(c(1,1)),
                Boxs=list(given),Constant_Pop=F,estimate="SF",NC=NC,
                Big_Window = my_big_window, pop_vect=my_pop_vect)
    total[[(1+length(total))]]=as.numeric(test_$Xi)
    s_ex[[(1+length(s_ex))]]=test_$sigma
    mu_ex[[(1+length(mu_ex))]]=mean(test_$mu)
    Tc[[(1+length(Tc))]]=test_$Tc    
  } else if (inference_model == "OneTransition") {
    test_=teSMC(n=n,rho=rho,model=list("One transition"),O_total,maxit =20,
                BoxB=list(c(0.05,1)),BoxP=c(3,3),Boxr=list(c(1,1)),
                Boxs=list(list(c(0.5,0.99),c(0,0.5))),Constant_Pop=F,
                estimate="SF",NC=NC, 
                Big_Window = my_big_window, pop_vect=my_pop_vect)
    total[[(1+length(total))]]=as.numeric(test_$Xi)
    s_ex[[(1+length(s_ex))]]=test_$sigma
    mu_ex[[(1+length(mu_ex))]]=mean(test_$mu)
    Tc[[(1+length(Tc))]]=test_$Tc    
  }

  result_list = list()
  result_list[[1]] = test_
  result_list[[2]] = total
  result_list[[3]] = s_ex
  result_list[[4]] = mu_ex
  result_list[[5]] = Tc

  # save full data
  saveRDS(result_list,file = outfilerds)
  cat("saved RDS"); cat("\n")

  return(result_list)
}

message("NC:", NC)

result <- run_esmc(
  path_to_mhs = snakemake@input[["mhs"]],
  M = M,
  mu = mu,
  r = r,
  n = n,
  rho = rho,
  NC = NC,
  inference_model = inference_model,
  recent_sigma = recent_sigma,
  ancient_sigma = ancient_sigma,
  outfilerds = snakemake@output[["full_data"]]
  )
cat("inference successful"); cat("\n")

sink()

