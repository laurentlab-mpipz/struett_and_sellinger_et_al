
# make performance analysis with neuralnetwork (without pls)

library(reticulate); pd <- import("pandas")

source("scripts/rfunctions.R")

# read arguments from command line
args = commandArgs(trailingOnly=TRUE)

# explify number of requested pls
requested_pls <- as.numeric(as.character(args[6:length(args)]))

# explify number of tolerated simulations
num_tolerated_simulations <- as.numeric(as.character(args[4]))

# extract directory for distribution plots
directory_for_distribution_plots = paste0(strsplit(args[2],split = "/")[[1]][-length(strsplit(args[2],split = "/")[[1]])], collapse = "/")

# read in transformed data frames from txt files
df_abc = data.frame(read.table(args[2], header = T, nrows = 10000))
df_pod = data.frame(read.table(args[3], header = T, nrows = 10000))

# extract the sumstat composition
sumstats_abc <- obtain_columns_with_wanted_sumstats(df_abc, "LinearCombination")
sumstats_pod <- obtain_columns_with_wanted_sumstats(df_pod, "LinearCombination")

# extract the parameters
params_abc <- obtain_columns_with_wanted_sumstats(df_abc, "param")
params_pod <- obtain_columns_with_wanted_sumstats(df_pod, "param")

# define the numbers of pls to be tested
num_pls_components <- get_number_of_pls_components_for_analysis(requested_pls, ncol(sumstats_abc))

# extract pod number from output file path
nPOD <- strsplit(strsplit(args[1], "/")[[1]][3], "_")[[1]][2]

# read logit boundaries from yaml file
logit_boundaries <- get_logit_bounds_from_yaml(as.character(args[5]), colnames(params_abc))

# result list
list_of_pls_performance_analysis <- list()

# run loop for each number of pls components
for (x in 1:length(num_pls_components)) {
  
  nComp <- num_pls_components[x]
  message("number of components using for loclinear regression: ", nComp)
  
  # create directory for distribution plots
  directory_for_pls_distribution_plots <- paste(directory_for_distribution_plots, paste("pod", as.character(nPOD),sep = "_"), paste("pls", as.character(nComp), sep = "_"), sep="/")
  dir.create(file.path(directory_for_pls_distribution_plots), recursive = TRUE, showWarnings = FALSE)
  
  # reduce to number of wanted pls_components
  pls_sumstats_abc <- sumstats_abc[paste("LinearCombination", 0:(nComp-1), sep = "_")]
  pls_sumstats_pod <- sumstats_pod[paste("LinearCombination", 0:(nComp-1), sep = "_")]
  
  abc_result <- run_performance_analysis(params_abc, pls_sumstats_abc,
                                         params_pod, pls_sumstats_pod,
                                         "loclinear", num_tolerated_simulations,
                                         directory_for_pls_distribution_plots,
                                         logit_boundaries)
  abc_result$plsComp = nComp
  
  # save into list
  list_of_pls_performance_analysis[[x]] <- abc_result
  
  # clear environment
  rm(pls_sumstats_abc,
     pls_sumstats_pod,
     abc_result,
     nComp,
     directory_for_pls_distribution_plots
  )
}

# concatenate all evaluated analysis
pls_performance_analysis <- do.call("rbind.data.frame", list_of_pls_performance_analysis)

# save to pickle
pd$to_pickle(pls_performance_analysis, args[1])

