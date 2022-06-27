
# make performance analysis with neuralnetwork (without pls)

library(reticulate); pd <- import("pandas")

source("scripts/rfunctions.R")

# read arguments from command line
args = commandArgs(trailingOnly=TRUE)

# explify number of requested pls
requested_pls <- as.numeric(as.character(args[8:length(args)]))

# explify number of tolerated simulations
num_tolerated_simulations <- as.numeric(as.character(args[6]))

# extract directory for distribution plots
directory_for_distribution_plots = paste0(strsplit(args[4],split = "/")[[1]][-length(strsplit(args[3],split = "/")[[1]])], collapse = "/")

# read in transformed data frames from txt files
df_abc = data.frame(read.table(args[4], header = T))
df_pod = data.frame(read.table(args[5], header = T))

# extract the sumstat composition
sumstats_abc <- obtain_columns_with_wanted_sumstats(df_abc, "LinearCombination")
sumstats_pod <- obtain_columns_with_wanted_sumstats(df_pod, "LinearCombination")

# extract the parameters
params_abc_all <- obtain_columns_with_wanted_sumstats(df_abc, "param")
params_pod_all <- obtain_columns_with_wanted_sumstats(df_pod, "param")

# remove columns with fixed params
params_abc <- remove_columns_with_no_variance(params_abc_all)
params_pod <- params_pod_all[, colnames(params_pod_all) %in% colnames(params_abc)]

# define the numbers of pls to be tested; use version 2 which returns false if num_pls is not available; this is recommended if each num_pls is run in parallel (not in loop)
num_pls_components <- get_number_of_pls_components_for_analysis_2(requested_pls, ncol(sumstats_abc))

# extract pod number from output file path
nPOD <- strsplit(strsplit(args[1], "/")[[1]][3], "_")[[1]][2]

# # read logit boundaries from yaml file
# logit_boundaries <- get_logit_bounds_from_yaml(as.character(args[6]), colnames(params_abc))
# print(logit_boundaries)

# change to actual logit boundaries
logit_boundaries <- t(apply(params_abc, 2, range))
colnames(logit_boundaries) <- c("minimal", "maximal")

# result list
list_of_pls_performance_analysis <- list()
list_of_pls_point_estimates <- list()
list_of_pls_quants <- list()

# run loop for each number of pls components if wanted pls is available, else return empty dataframe
if (num_pls_components[1]) {
  for (x in 1:length(num_pls_components)) {
    
    nComp <- num_pls_components[x]
    message("number of components using for loclinear regression: ", nComp)
        
    # create directory for distribution plots
    directory_for_pls_distribution_plots <- paste(directory_for_distribution_plots, paste("pod", as.character(nPOD),sep = "_"), paste("pls", as.character(nComp), sep = "_"), sep="/")
    dir.create(file.path(directory_for_pls_distribution_plots), recursive = TRUE, showWarnings = FALSE)
    
    # reduce to number of wanted pls_components
    pls_sumstats_abc <- sumstats_abc[paste("LinearCombination", 0:(nComp-1), sep = "_")]
    pls_sumstats_pod <- sumstats_pod[paste("LinearCombination", 0:(nComp-1), sep = "_")]
    
    print("check stuff for debugging")
    print(x)
    print(num_pls_components)
    
    abc_result_list <- run_performance_analysis(params_abc, pls_sumstats_abc,
                                           params_pod, pls_sumstats_pod,
                                           "loclinear", num_tolerated_simulations,
                                           directory_for_pls_distribution_plots,
                                           logit_boundaries)
    abc_result <- abc_result_list[[1]]
    abc_result$plsComp = nComp
    
    abc_estims <- abc_result_list[[2]]
    abc_estims$plsComp = nComp
    
    abc_quants <- abc_result_list[[3]]
    abc_quants$plsComp = nComp
    
    # save into list
    list_of_pls_performance_analysis[[x]] <- abc_result
    list_of_pls_point_estimates[[x]] <- abc_estims
    list_of_pls_quants[[x]] <- abc_quants
    
    # clear environment
    rm(pls_sumstats_abc,
       pls_sumstats_pod,
       abc_result,
       nComp,
       directory_for_pls_distribution_plots
    )
  }
} else {
  # create a named data frame, such that pandas can save it
  # this will break the pipeline, if abc_result function is changed
  list_of_pls_performance_analysis[[1]] <-
    cbind.data.frame(point_estim=c("None"),
                     eval_mode=c("None"),
                     parameter=c("None"),
                     value=c("None"),
                     plsComp=c("None"))
  list_of_pls_point_estimates[[1]] <-
    cbind.data.frame(
      param_0=c("None"),
      param_1=c("None"),
      param_2=c("None"),
      param_3=c("None"),
      param_4=c("None"),
      param_5=c("None"),
      pointing=c("None"),
      plsComp=c("None")
    )
  
  list_of_pls_quants[[1]] <-
    cbind.data.frame(
      param_0=c("None"),
      param_1=c("None"),
      param_2=c("None"),
      param_3=c("None"),
      param_4=c("None"),
      param_5=c("None"),
      pointing=c("None"),
      plsComp=c("None")
    )
}

# concatenate all evaluated analysis
pls_performance_analysis <- do.call("rbind.data.frame", list_of_pls_performance_analysis)
pls_point_estims <- do.call(rbind.data.frame, list_of_pls_point_estimates)
pls_quants <- do.call(rbind.data.frame, list_of_pls_quants)

# save to pickle
pd$to_pickle(pls_performance_analysis, args[1])
pd$to_pickle(pls_point_estims, args[2])
pd$to_pickle(pls_quants, args[3])

