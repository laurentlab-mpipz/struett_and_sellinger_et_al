
library(reticulate); pd <- import("pandas")

source("scripts/rfunctions.R")
message("loaded rfunctions")

# read arguments from command line
args = commandArgs(trailingOnly=TRUE)

# explicify output files
outfile_performance <- args[1]

# extract the sumstat names composition to use
sscomp = get_list_of_sumstat_vectors(args[5], "NEXT")[[as.numeric(args[4])]]

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
params_abc <- obtain_columns_with_wanted_sumstats(df_abc, "param")
params_pod <- obtain_columns_with_wanted_sumstats(df_pod, "param")
print(head(params_abc))

################################################################################
# find pls from script from abc toolbox
wanted_pls_nums <- unlist( sapply(as.integer(args[7:length(args)]), function(x) if(x <= ncol(sumstats_abc) ) return(as.character(x)) ) )
if (length(wanted_pls_nums) == 0 ) { wanted_pls_nums <- ncol(sumstats_abc) }
DO_PLS <- T
if (length(wanted_pls_nums) <= 1 ) {DO_PLS <- F}
message("I will make the pls transformation: ", DO_PLS)

abc_result_list  <- list()

# calculate the pls for each pls, separately
for (x in wanted_pls_nums) {
  
  #open File
  dir_prefix <- paste(strsplit(args[1], "/")[[1]][1:(length(strsplit(args[1], "/")[[1]])-1)], collapse = "/")
  dir_component_1 <- paste0(sscomp, collapse = "_")
  directory <- paste0(dir_prefix, "/pls/", dir_component_1, "_numComp_", x, "/", collapse = "");
  dir.create(file.path(directory), showWarnings = F, recursive = T);
  filename<-paste(sscomp, collapse = "_");
  numComp<-as.numeric(x);
  
  # reassign the stats
  stats<-sumstats_abc; params<-params_abc;
  print(dim(stats))
  
  # use only maximal 30,000 rows to calculate the pls
  if (nrow(stats)>30000) {
    stats<-stats[1:30000,]; params<-params[1:30000,]
  }
  
  #standardize the params
  for(i in 1:length(params)){params[,i]<-(params[,i]-mean(params[,i]))/sd(params[,i]);}
  
  #force stats in [1,2]
  myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
  for(i in 1:length(stats)){
    myMax<-c(myMax, max(stats[,i]));
    myMin<-c(myMin, min(stats[,i]));
    stats[,i]<-1+(stats[,i]-myMin[i])/(myMax[i]-myMin[i]);
  }
  
  #transform statistics via boxcox  
  library("MASS");	
  for(i in 1:length(stats)){		
    d<-cbind(stats[,i], params);
    mylm<-lm(as.formula(d), data=d)			
    myboxcox<-boxcox(mylm, lambda=seq(-50, 80, 1/10), plotit=T, interp=T, eps=1/50);	
    lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);			
    print(paste(names(stats)[i], myboxcox$x[myboxcox$y==max(myboxcox$y)]));
    myGM<-c(myGM, exp(mean(log(stats[,i]))));			
  }
  
  #standardize the BC-stats
  myBCMeans<-c(); myBCSDs<-c();
  for(i in 1:length(stats)){
    stats[,i]<-(stats[,i]^lambda[i] - 1)/(lambda[i]*myGM[i]^(lambda[i]-1));	
    myBCSDs<-c(myBCSDs, sd(stats[,i]));
    myBCMeans<-c(myBCMeans, mean(stats[,i]));		
    stats[,i]<-(stats[,i]-myBCMeans[i])/myBCSDs[i];
  }
  
  #perform pls
  library("pls");
  myPlsr<-plsr(as.matrix(params) ~ as.matrix(stats), scale=F, ncomp=numComp);
  
  if (DO_PLS) {
  #write pls to a file
  myPlsrDataFrame<-data.frame(comp1=myPlsr$loadings[,1]);
  for(i in 2:numComp) { myPlsrDataFrame<-cbind(myPlsrDataFrame, myPlsr$loadings[,i]); }
  write.table(cbind(names(stats), myMax, myMin, lambda, myGM, myBCMeans, myBCSDs, myPlsrDataFrame), file=paste(directory, "Routput_", filename, sep=""), col.names=F, row.names=F, sep="\t", quote=F);
  
  #make RMSE plot
  pdf(paste(directory, "RMSE_", filename, ".pdf", sep=""));
  plot(RMSEP(myPlsr));
  dev.off();
  }
  
  message("finished finding pls")
  if (!DO_PLS) message("using untransformed values")
  
  ## transform pls
  
  if (DO_PLS) {
  # file definitions
  def_txt <- paste(directory, "Routput_", filename, sep="");
  script <- "./scripts/transformer";
  sumstats_abc_non_transformed <- paste(directory, "sumstats_abc_non_transformed.txt", sep = "");
  sumstats_pod_non_transformed <- paste(directory, "sumstats_pod_non_transformed.txt", sep = "");
  sumstat_abc_transf <- paste(directory, "sumstats_abc_transformed.txt", sep = "");
  sumstat_pod_transf <- paste(directory, "sumstats_pod_transformed.txt", sep = "");
  
  write.table(sumstats_abc, file=sumstats_abc_non_transformed, col.names=T, row.names=F, sep="\t", quote=F);
  write.table(sumstats_pod, file=sumstats_pod_non_transformed, col.names=T, row.names=F, sep="\t", quote=F);
  
  # transformation of abc  
  system_call_value<-system(paste(script, def_txt, sumstats_abc_non_transformed, sumstat_abc_transf, "boxcox", sep=" "));
  message("abc transformation: system call value: ", system_call_value)
  
  # transformation of pod
  system_call_value<-system(paste(script, def_txt, sumstats_pod_non_transformed, sumstat_pod_transf, "boxcox", sep=" "));
  message("pod transformation: system call value: ", system_call_value)
  
  # load transformed values  
  sumstats_abc_transformed <- data.frame(read.table(sumstat_abc_transf, header = T, row.names = NULL));
  sumstats_pod_transformed <- data.frame(read.table(sumstat_pod_transf, header = T, row.names = NULL));
  } else {
    sumstats_abc_transformed <- sumstats_abc
    sumstats_pod_transformed <- sumstats_pod
  }
  
  ## abc performance
  directory_for_distribution_plots <- directory;
  abc_result <- run_performance_analysis(params_abc, sumstats_abc_transformed,
                                         params_pod, sumstats_pod_transformed,
                                         "loclinear", as.numeric(args[6]),
                                         directory_for_distribution_plots)
  abc_result$plsComp = as.numeric(x)
  message("got the abc_performance of pls ", x)
  
  # save everything into outputfile
  abc_result_list <- c(abc_result_list, list(abc_result))
  
  # delete the temporary files
  for (fn in c(sumstats_abc_non_transformed, sumstats_pod_non_transformed,
               sumstat_abc_transf, sumstat_pod_transf)) {
    if (file.exists(fn)) {file.remove(fn)}
  }
  message("finished loop pls no ", x)
}
message("calculations done, start formatting..")

# rbind the data frames containing the performance analysis for each result
abc_results_data_frame <- do.call("rbind.data.frame", abc_result_list)
# print(abc_results_data_frame)
message("rbind.data.frame done")

pd$to_pickle(abc_results_data_frame, outfile_performance)
message("finished all pls transformations")






