# save.image(file = "My_Object.RData")
# stop("created image to start manual development")
# setwd(
#   "~/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/struett/temp_tsabc/a_thaliana/tsabc_at_estimates/"
# )
# load("My_Object.RData")


library(locfit)
library(reticulate)
library("MASS")
library(abc)


pd <- import("pandas")


# Obtain parameters for estimate
{
  OUTFILE_POSTERIORS  <- snakemake@output$posteriors
  INFILE_ABC  <- snakemake@input$abc
  INFILE_OBS <- snakemake@input$obs
  NTOLERATED <-
    as.numeric(snakemake@config$number_of_tolerated_simulations)
  SSCOMP_NO <- as.numeric(snakemake@wildcards$sscomp)
  PLS_COMP <- as.numeric(snakemake@wildcards$pls)
  PLOT_DIR <- snakemake@output$posterior_plots
}

# Read abc sim data and subset to pls
{
  # read in table
  df_abc <- data.frame(read.table(INFILE_ABC, header = T))
  
  
  # subset stats (remove parameters)
  all_sumstats_abc <-
    df_abc[, grep("LinearCombination", names(df_abc))]
  if (ncol(all_sumstats_abc) <  PLS_COMP) {
    message("  changed PLS set from ", PLS_COMP, " to ", ncol(all_sumstats_abc))
    PLS_COMP <- ncol(all_sumstats_abc)
  }
  sumstats_abc <-
    all_sumstats_abc[, paste("LinearCombination", 0:(PLS_COMP - 1), sep = "_")]
  
  
  # Get params and remove the columns that do not have variance
  params_abc_all <- df_abc[, grep("param", names(df_abc))]
  params_abc <- params_abc_all[, apply(params_abc_all, 2, var) != 0]
  
  
  message("  read and prepared simulations")
}


# Read observation data
{
  # read in table
  df_obs <- data.frame(read.table(INFILE_OBS, header = T))
  
  
  # subset stats (remove parameters)
  all_sumstats_obs <-
    df_obs[, grep("LinearCombination", names(df_obs))]
  sumstats_obs <-
    all_sumstats_obs[paste("LinearCombination", 0:(PLS_COMP - 1), sep = "_")]
  
  
  message("  read and prepared observations")
}


# Produce logit boundaries from actual priors
{
  logit_boundaries <- t(apply(params_abc, 2, range))
  colnames(logit_boundaries) <- c("minimal", "maximal")
}

# results lists
{
  prior_and_posteriors <- list()
}


{
  # parameter preparation
  my_tolerance <- NTOLERATED / nrow(sumstats_abc)
  quant_vec_for_ci <-
    c(
      0.005,
      0.025,
      0.05,
      0.1,
      0.25,
      0.375,
      0.45,
      0.5,
      0.55,
      0.625,
      0.75,
      0.9,
      0.95,
      0.975,
      0.995
    )
}


# make estimates with library abc
{
  # make an estimate for each set of sumstats in the table (i. e. the samples)
  for (i in 1:nrow(sumstats_obs)) {
    # extract summary statistic to estimate from
    my_target <- as.numeric(sumstats_obs[i, ])


    # create estimate
    flag <- TRUE
    tryCatch(
      expr = {
        abc_result <- abc(
          target = my_target,
          param = params_abc,
          sumstat = sumstats_abc,
          tol = my_tolerance,
          method = "loclinear",
          MaxNWts = 5000,
          transf = "logit",
          logit.bounds = logit_boundaries
        )
      },
      error = function(e) {
        message("* Caught an error on itertion ", i)
        flag <- FALSE
      }
    )
    if (!exists("abc_result"))
      next


    # clean the data from na's
    rejection_values <- na.omit(abc_result$unadj.values)
    adjusted_values <- na.omit(abc_result$adj.values)


    # save results into list
    prior_and_posteriors[[i]] <-
      list(rejection_values, adjusted_values)
    message("  saved priors and posteriors to list: ",
            i,
            " of ",
            nrow(sumstats_obs))


    # plot diagnostics
    tryCatch(
      expr = {
        dir.create(PLOT_DIR, showWarnings = FALSE)
        pdf(paste0(PLOT_DIR, "/diagnostic_plots_", i, ".pdf", collapse = ""))
        plot(abc_result, param = params_abc, ask = F)
        dev.off()
    },
      error = function(e) {
        message("* Caught a plotting error on itertion ", i)
      }
    )
  }
}


# make estimate on the mean of the statistics; i. e. that the sampling variance on the stats is reduced
{
  # calculate the mean of each sumstat to estimate from
  my_target <- apply(sumstats_obs, 2, mean)
  
  
  # create estimate
  rm("abc_result")
  tryCatch(
    expr = {
      abc_result <- abc(
        target = my_target,
        param = params_abc,
        sumstat = sumstats_abc,
        tol = my_tolerance,
        method = "loclinear",
        MaxNWts = 5000,
        transf = "logit",
        logit.bounds = logit_boundaries
      )
    },
    error = function(e) {
      message("* Caught an error on itertion ", i)
      flag <- FALSE
    }
  )
  
  
  # only obtain the results if there was no error
  if (exists("abc_result")) {
    # clean the data from na's
    rejection_values <- na.omit(abc_result$unadj.values)
    adjusted_values <- na.omit(abc_result$adj.values)
    
    
    # save results into list
    prior_and_posteriors[["meanEstim"]] <-
      list(rejection_values, adjusted_values)
    message("saved priors and posteriors to list: mean estimate")
    
    
    # plot diagnostics
    tryCatch(
      expr = {
        pdf(paste0(PLOT_DIR, "/diagnostic_plots_", "meanEstim", ".pdf", collapse = ""))
        plot(abc_result, param = params_abc, ask = F)
        dev.off()
    },
      error = function(e) {
        message("* Caught a plotting error on itertion ", "meanEstim")
      }
    )
  }
}


# save results list
{
  # add prior to the result list
  prior_and_posteriors[["prior"]] <- params_abc
  
  saveRDS(prior_and_posteriors, file = OUTFILE_POSTERIORS)
}
