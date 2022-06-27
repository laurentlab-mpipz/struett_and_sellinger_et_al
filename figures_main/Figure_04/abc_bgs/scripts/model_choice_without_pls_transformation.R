
# model choice after box-cox transformation and no further pls transformation

library(reticulate)
library("MASS")

pd <- import("pandas")

# read arguments from command line
args = commandArgs(trailingOnly=TRUE)

OUTFILE = args[1]
INFILE_ABC = args[2]
INFILE_ALT = args[3]
INFILE_POD = args[4]
NTOLERATED = args[5]
SSCOMP_NO = args[6]
SSSTRING = args[7]
PLS_COMP = args[8]

# for (i in 1:length(args)) {
#   message(i, ": ", args[i])
# }

# this is script is made for model choice without pls transformation
stopifnot(as.numeric(PLS_COMP) == 0)

ntolerated = as.numeric(NTOLERATED)
sscomp_no = as.numeric(SSCOMP_NO)
ntolerated = as.numeric(NTOLERATED)

# read in data frames from pickle
df_abc = pd$read_pickle(INFILE_ABC)
df_alt = pd$read_pickle(INFILE_ALT)
df_pod = pd$read_pickle(INFILE_POD)

# restrain the selfing rate of the alternative model
message("here happens a manual restrain of the selfing rate in the alternative model to 0.5-1.0")
df_alt <- df_alt[df_alt$param_3 >= 0.5 & df_alt$param_3<=1.0,]

# get names of sumstats classes from string
source("scripts/rfunctions/get_names_of_sumstat_classes_from_string.R")
which_sumstats = get_names_of_sumstat_classes_from_string(SSSTRING, sscomp_no)

source("scripts/rfunctions.R")
# extract the sumstat composition
sumstats_abc <- obtain_columns_with_wanted_sumstats(df_abc, which_sumstats)
sumstats_alt <- obtain_columns_with_wanted_sumstats(df_alt, which_sumstats)
sumstats_pod <- obtain_columns_with_wanted_sumstats(df_pod, which_sumstats)

# model specification
sumstats_abc$model = "tsigma"
sumstats_alt$model = "csigma"
sumstats_pod$model = "pod"

# reduce the pod table to colnames present in the abc table
cols_to_remain <- as.vector(Reduce(
  intersect, list(
    colnames(sumstats_abc),
    colnames(sumstats_alt),
    colnames(sumstats_pod)
    )
  ))

sumstats_pod <- sumstats_pod[cols_to_remain]

# combine into 1 table
df <- rbind.data.frame(sumstats_abc, sumstats_alt, sumstats_pod)
model_vector=df$model; df$model=NULL

# remove columns with no variance from abc simulations
df <- remove_columns_with_no_variance(df)

# transform stats
stats <- df

#force stats in [1,2]
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:length(stats)){
  myMax<-c(myMax, max(stats[,i]));
  myMin<-c(myMin, min(stats[,i]));
  stats[,i]<-1+(stats[,i]-myMin[i])/(myMax[i]-myMin[i]);
}

df_model_choice <- stats[model_vector %in% c("tsigma", "csigma"),]
df_pods <- stats[model_vector=="pod",]
model_vector <- model_vector[model_vector %in% c("tsigma", "csigma")]

# model choice for all pods
source("scripts/rfunctions/get_bayes_factor_for_each_pod.R")
df_model_choice_result <- get_bayes_factor_for_each_pod(
  df_model_choice,
  df_pods,
  model_vector,
  ntolerated
  )

# combine results with pod param table
df_model_choice_result <- 
  cbind.data.frame(obtain_columns_with_wanted_sumstats(
    rbind.data.frame(df_pod, df_pod),
    c("ident", "ecex_time", "rand", "param")),
    df_model_choice_result)
print(str(df_model_choice_result))

a <- sum(df_model_choice_result$bayes_factor[df_model_choice_result$regression_method=="rejection"]>1.6)
b <- sum(df_model_choice_result$bayes_factor[df_model_choice_result$regression_method=="mnlogistic"]>1.6)
b_2 <- sum(df_model_choice_result$bayes_factor[df_model_choice_result$regression_method=="mnlogistic"]<1/1.6)
d <- sum(df_model_choice_result$bayes_factor[df_model_choice_result$regression_method=="rejection"]>1)
cat("\n")
message("sumstats: ", paste0(which_sumstats, collapse = " "))
message("bayes factor > 1:\trejection\t", d/nrow(df_pod))
message("bayes factor > 1.6:\trejection\t", a/nrow(df_pod))
message("change to selfing prefered (k>1.6):\tmnlogistic\t", 100*b/nrow(df_pod), "%")
message("constant selfing prefered (1/k<1.6):\tmnlogistic\t", 100*b_2/nrow(df_pod), "%")
message("no model preference:'tmnlogistic\t", 100*(1-(b+b_2)/nrow(df_pod)), "%")
cat("\n")

# save to output file
pd$to_pickle(df_model_choice_result, OUTFILE)

message("finished model choice analysis for ", INFILE_POD, " sscomp: ",
        SSCOMP_NO, " n-pls: ", PLS_COMP)