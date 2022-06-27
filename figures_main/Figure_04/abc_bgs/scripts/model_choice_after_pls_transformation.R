

library(reticulate)
library("MASS")

pd <- import("pandas")

# read arguments from command line
args = commandArgs(trailingOnly = TRUE)

# print(getwd())
# save.image(
#   "/netscratch/dep_tsiantis/grp_laurent/struett/temp_tsabc/mchoice/tsabc3/TEMP.RData"
# )
# stop("#######################################")
# load("TEMP.RData")


OUTFILE = args[1]
INFILE_ABC = args[2]
INFILE_ALT = args[3]
INFILE_POD = args[4]
NTOLERATED = args[5]
SSCOMP_NO = args[6]
SSSTRING = args[7]
PLS_COMP = args[8]

ntolerated = as.numeric(NTOLERATED)
sscomp_no = as.numeric(SSCOMP_NO)
nComp = as.numeric(PLS_COMP)

# read in data frames from pickle
df_abc = data.frame(read.table(INFILE_ABC, header = T))
df_alt = data.frame(read.table(INFILE_ALT, header = T))
df_pod = data.frame(read.table(INFILE_POD, header = T))

# # restrain the selfing rate of the alternative model
# message("here happens a manual restrain of the selfing rate in the alternative model to 0.5-1.0")
# df_alt <- df_alt[df_alt$param_3 >= 0.5 & df_alt$param_3<=1.0,]

source("scripts/rfunctions.R")
# extract the sumstat composition
all_sumstats_abc <-
  obtain_columns_with_wanted_sumstats(df_abc, "LinearCombination")
all_sumstats_alt <-
  obtain_columns_with_wanted_sumstats(df_alt, "LinearCombination")
all_sumstats_pod <-
  obtain_columns_with_wanted_sumstats(df_pod, "LinearCombination")

if (nComp <=  ncol(all_sumstats_abc)) {
  # extract number of pls
  sumstats_abc <-
    all_sumstats_abc[paste("LinearCombination", 0:(nComp - 1), sep = "_")]
  sumstats_alt <-
    all_sumstats_alt[paste("LinearCombination", 0:(nComp - 1), sep = "_")]
  sumstats_pod <-
    all_sumstats_pod[paste("LinearCombination", 0:(nComp - 1), sep = "_")]
  
  # model specification
  sumstats_abc$model = "tsigma"
  sumstats_alt$model = "csigma"
  sumstats_pod$model = "pod"
  
  # combine into 1 table
  df <- rbind.data.frame(sumstats_abc, sumstats_alt, sumstats_pod)
  model_vector = df$model
  df$model = NULL
  
  df_model_choice <- df[model_vector %in% c("tsigma", "csigma"), ]
  df_pods <- df[model_vector == "pod", ]
  model_vector <-
    model_vector[model_vector %in% c("tsigma", "csigma")]
  
  # model choice for all pods
  source("scripts/rfunctions/get_bayes_factor_for_each_pod.R")
  df_model_choice_result <- get_bayes_factor_for_each_pod(df_model_choice,
                                                          df_pods,
                                                          model_vector,
                                                          ntolerated)
  
  # combine results with pod param table
  df_model_choice_result <-
    cbind.data.frame(
      obtain_columns_with_wanted_sumstats(
        rbind.data.frame(df_pod, df_pod),
        c("ident", "ecex_time", "rand", "param")
      ),
      df_model_choice_result
    )
  print(str(df_model_choice_result))
  
  a <-
    sum(df_model_choice_result$bayes_factor[df_model_choice_result$regression_method ==
                                              "rejection"] > 1.6)
  b <-
    sum(df_model_choice_result$bayes_factor[df_model_choice_result$regression_method ==
                                              "mnlogistic"] > 1.6)
  b_2 <-
    sum(df_model_choice_result$bayes_factor[df_model_choice_result$regression_method ==
                                              "mnlogistic"] < 1 / 1.6)
  d <-
    sum(df_model_choice_result$bayes_factor[df_model_choice_result$regression_method ==
                                              "rejection"] > 1)
  cat("\n")
  message("sumstats: nComp/sscomp",
          paste(nComp, sscomp_no, sep = " ", collapse = " "))
  message("bayes factor > 1:\trejection\t", d / nrow(df_pod))
  message("bayes factor > 1.6:\trejection\t", a / nrow(df_pod))
  message("change to selfing prefered (k>1.6):\tmnlogistic\t",
          100 * b / nrow(df_pod),
          "%")
  message("constant selfing prefered (1/k<1.6):\tmnlogistic\t",
          100 * b_2 / nrow(df_pod),
          "%")
  message("no model preference:'tmnlogistic\t", 100 * (1 - (b + b_2) / nrow(df_pod)), "%")
  cat("\n")
  
  # save to output file
  pd$to_pickle(df_model_choice_result, OUTFILE)
  
  message(
    "finished model choice analysis for ",
    INFILE_POD,
    " sscomp: ",
    SSCOMP_NO,
    " n-pls: ",
    PLS_COMP
  )
} else {
  # save empty df to output file
  df_model_choice_result <-
    cbind.data.frame(obtain_columns_with_wanted_sumstats(
      rbind.data.frame(df_pod, df_pod),
      c("ident", "ecex_time", "rand", "param")
    ))
  df_model_choice_result$marginal_density_tsigma = NA
  df_model_choice_result$marginal_density_csigma = NA
  df_model_choice_result$bayes_factor = NA
  df_model_choice_result$regression_method = NA
  print(str(df_model_choice_result))
  
  pd$to_pickle(df_model_choice_result, OUTFILE)
  
  message("not as many pls components as requested..")
  message(
    "finished model choice analysis for ",
    INFILE_POD,
    " sscomp: ",
    SSCOMP_NO,
    " n-pls: ",
    PLS_COMP
  )
}
