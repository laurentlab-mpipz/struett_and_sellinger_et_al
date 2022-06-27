library(abc)

get_bayes_factor_for_each_pod = function(df_model_choice,
                                         df_pods,
                                         model_vector,
                                         ntolerated) {
  stopifnot(nrow(df_model_choice) == length(model_vector))
  
  # if only single sumstat
  if (is.null(dim(df_pods)))
    df_pods <- matrix(df_pods)
  
  # downsample the model with more simulations
  n_abc = sum(model_vector == "tsigma")
  n_alt = sum(model_vector == "csigma")
  if (n_abc != n_alt) {
    df_abc = df_model_choice[model_vector == "tsigma", ]
    df_alt = df_model_choice[model_vector == "csigma", ]
    
    if (n_abc > n_alt) {
      df_abc <- df_abc[sample(n_abc, n_alt),]
    } else if (n_abc < n_alt) {
      df_alt <- df_alt[sample(n_alt, n_abc),]
    }
    
    df_abc$model = "tsigma"
    df_alt$model = "csigma"
    
    df_model_choice = rbind.data.frame(df_abc, df_alt)
    model_vector = df_model_choice$model
    df_model_choice$model = NULL
    
    df_model_choice = remove_columns_with_no_variance(df_model_choice)
    df_pods = df_pods[, colnames(df_pods) %in% colnames(df_model_choice)]
  }
  
  tol <- ntolerated / nrow(df_model_choice)
  
  # loop through pods and get bayes
  marginal_density_tsigma_mnlogistic <- numeric()
  marginal_density_csigma_mnlogistic <- numeric()
  marginal_density_tsigma_rejection <- numeric()
  marginal_density_csigma_rejection <- numeric()
  for (pod_index in 1:nrow(df_pods)) {
    a <- postpr(
      df_pods[pod_index,],
      model_vector,
      df_model_choice,
      tol = tol,
      method = "mnlogistic",
      corr = TRUE   # corr seems not to work
    )
    b <- summary(a)
    
    # if no regression done; the results structure is different
    if (is.null(b$mnlogistic)) {
      new_tsigma_mnlogistic <- NA
      new_csigma_mnlogistic <- NA
    } else {
      new_tsigma_mnlogistic <- b$mnlogistic$Prob["tsigma"]
      new_csigma_mnlogistic <- b$mnlogistic$Prob["csigma"]
    }
    
    if (is.null(b$rejection)) {
      new_tsigma_rejection <- b$Prob["tsigma"]
      new_csigma_rejection <- b$Prob["csigma"]
    } else {
      new_tsigma_rejection <- b$rejection$Prob["tsigma"]
      new_csigma_rejection <- b$rejection$Prob["csigma"]
    }
    
    marginal_density_tsigma_mnlogistic <-
      c(marginal_density_tsigma_mnlogistic,
        new_tsigma_mnlogistic)
    marginal_density_csigma_mnlogistic <-
      c(marginal_density_csigma_mnlogistic,
        new_csigma_mnlogistic)
    
    marginal_density_tsigma_rejection <-
      c(marginal_density_tsigma_rejection, new_tsigma_rejection)
    marginal_density_csigma_rejection <-
      c(marginal_density_csigma_rejection, new_csigma_rejection)
    
    message("row #",
            pod_index,
            " len ",
            length(marginal_density_tsigma_mnlogistic),
            "\n\n")
    rm(a)
  }
  
  bayes_rejection = marginal_density_tsigma_rejection / marginal_density_csigma_rejection
  bayes_mnlogistic = marginal_density_tsigma_mnlogistic / marginal_density_csigma_mnlogistic
  
  # saveRDS(
  #   object = list(
  #     a = marginal_density_tsigma_mnlogistic,
  #     b = marginal_density_csigma_mnlogistic,
  #     c = bayes_mnlogistic,
  #     d = marginal_density_tsigma_rejection,
  #     e = marginal_density_csigma_rejection,
  #     f = bayes_rejection
  #   ),
  #   file = "/netscratch/dep_tsiantis/grp_laurent/struett/temp_tsabc/mchoice/tsabc3/TEMP.RDS"
  # )
  # stop("#######################################################")
  # setwd(
  #   "~/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/struett/temp_tsabc/mchoice/tsabc3"
  # )
  # readRDS("TEMP.RDS")
  # rds_list = readRDS("TEMP.RDS")
  # marginal_density_tsigma_mnlogistic = rds_list[[1]]
  # marginal_density_csigma_mnlogistic = rds_list[[2]]
  # bayes_mnlogistic = rds_list[[3]]
  # marginal_density_tsigma_rejection = rds_list[[4]]
  # marginal_density_csigma_rejection = rds_list[[5]]
  # bayes_rejection = rds_list[[6]]
  
  df <-
    rbind.data.frame(
      cbind.data.frame(
        marginal_density_tsigma = marginal_density_tsigma_mnlogistic,
        marginal_density_csigma = marginal_density_csigma_mnlogistic,
        bayes_factor = bayes_mnlogistic,
        regression_method = "mnlogistic"
      ),
      cbind.data.frame(
        marginal_density_tsigma = marginal_density_tsigma_rejection,
        marginal_density_csigma = marginal_density_csigma_rejection,
        bayes_factor = bayes_rejection,
        regression_method = "rejection"
      )
    )
  
  return(df)
}
