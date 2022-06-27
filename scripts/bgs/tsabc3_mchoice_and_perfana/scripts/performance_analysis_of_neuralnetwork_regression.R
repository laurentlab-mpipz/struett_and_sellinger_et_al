
# make performance analysis with neuralnetwork (without pls)

## for working inscript and debug
# setwd("/netscratch/dep_tsiantis/grp_laurent/struett/temp_tsabc/tsabc2")

library(reticulate); pd <- import("pandas")

source("scripts/rfunctions.R")

# read arguments from command line
args = commandArgs(trailingOnly=TRUE)
# args = args_inside

# extract the sumstat names composition to use
sscomp = get_list_of_sumstat_vectors(args[5], "NEXT")[[as.numeric(args[4])]]

# extract directory for distribution plots
directory_for_distribution_plots = paste0(strsplit(args[1],split = "/")[[1]][-length(strsplit(args[1],split = "/")[[1]])], collapse = "/")
directory_for_distribution_plots = paste0(directory_for_distribution_plots, "/", paste0(sscomp, collapse = "_") )
dir.create(file.path(directory_for_distribution_plots), recursive = TRUE, showWarnings = FALSE)

# read in data frames from pickle
df_abc = pd$read_pickle(args[2])
df_pod = pd$read_pickle(args[3])

# remove columns with no variance from abc simulations
df_abc <- remove_columns_with_no_variance(df_abc)
df_pod <- df_pod[which(colnames(df_pod) %in% colnames(df_abc))]

# extract the sumstat composition
sumstats_abc <- obtain_columns_with_wanted_sumstats(df_abc, sscomp)
sumstats_pod <- obtain_columns_with_wanted_sumstats(df_pod, sscomp)

# extract the parameters
params_abc_all <- obtain_columns_with_wanted_sumstats(df_abc, "param")
params_pod_all <- obtain_columns_with_wanted_sumstats(df_pod, "param")

# remove columns with fixed params
params_abc <- remove_columns_with_no_variance(params_pod_all)
params_pod <- params_pod_all[, colnames(params_pod_all) %in% colnames(params_abc)]

# read logit boundaries from yaml file
logit_boundaries <- get_logit_bounds_from_yaml(as.character(args[7]), colnames(params_abc))

# run the performance
abc_result <- run_performance_analysis(params_abc, sumstats_abc,
                                       params_pod, sumstats_pod,
                                       "neuralnet", as.numeric(args[6]),
                                       directory_for_distribution_plots,
                                       logit_boundaries)
abc_result$plsComp = "none"

pd$to_pickle(abc_result, args[1])