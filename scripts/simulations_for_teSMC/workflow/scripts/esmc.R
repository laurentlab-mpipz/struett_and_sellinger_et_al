
thread_limiter = "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1"
system2(thread_limiter)

sink(snakemake@log[["esmc_log"]])

# save image for debugging
# save.image(); stop()
# setwd("../..")
# load(".RData")

# Thibaut: Betterly we should only run this code on the mhs, because the code
# is not ready, yet, for the run on the tree sequences

# load source file of Thibaut's teSMC
source(snakemake@params[["esmc_package"]])

# run esmc on the multihetsep input
run_esmc <- function(path_mhs, nsam, mutation_rate, recombination_rate,
  sigma_init = 0.0, beta_prior = c(0.05, 1), sigma_prior = c(0, 0.999),
  estimate_rec = F, estimate_seed_bank = F, estimate_sigma = T, BW=T) {
  
  mu=mutation_rate # Mutation rate per position per generation 
  r=recombination_rate # recombination rate per position per generation 
  rho=r/mu # ratio recombination/mutation 
  M=nsam # Number of haplotypes
  simga_0=sigma_init # initially assumed selfing rate
  
  # Set boundaries
  BoxB=beta_prior #  min and max value of germination rate 
  Boxs=sigma_prior #  min and max value of selfing rate 
  
  ER=estimate_rec # TRUE to estimate recombination
  SB=estimate_seed_bank # TRUE to estimate seed bank
  SF=estimate_sigma # TRUE to estimate selfing

  # read data from file
  NC=1 # Number of analysed scaffold
  full_data=list()
  for (scaffold in 1:NC) {
    filename=path_mhs # Put character string of the multihetsep file name
    data=Get_real_data(NA,M,filename)
    full_data[[scaffold]]=data
  }
 
  if (NC == 1) {
    full_data = full_data[[1]]
  }
  
  result=teSMC(
    n=40,
    maxit=50,
    rho=rho,
    full_data,
    mu_r=mu,
    BoxB=BoxB,
    Boxs=Boxs,
    SB=SB,
    SF=SF,
    Rho=ER,
    Check=F,
    NC=NC,
    symbol_number=20,
    pop_vect=rep(2,20), # means per n states 1 sigma, sum of this vector must be equal to n, the number of states (= time bins)
    sigma=simga_0,
    BW=BW
  )
  
  return(result)
}

nsam = as.numeric(snakemake@wildcards[["sample_size"]])
mutation_rate = as.numeric(snakemake@wildcards[["mutation_rate"]])
recombination_rate = as.numeric(snakemake@wildcards[["recombination_rate"]])

result <- run_esmc(snakemake@input[["mhs"]], nsam, mutation_rate, recombination_rate)
saveRDS(result, file=snakemake@output[["data"]])

sink()
