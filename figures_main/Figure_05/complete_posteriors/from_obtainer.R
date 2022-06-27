run_performance_analysis = function(pp_abc, ss_abc, pp_pod, ss_pod, my_method,
                                    tolerance_level, directory_posteriors,
                                    logit_boundaries) {
  # calculate tolerance
  my_tolerance = tolerance_level/nrow(ss_abc)
  
  # dataframes to collect the point estimates
  df.estim.mean = setNames(data.frame(matrix(ncol = ncol(pp_pod), nrow = 0)), colnames(pp_pod))
  df.estim.mode = setNames(data.frame(matrix(ncol = ncol(pp_pod), nrow = 0)), colnames(pp_pod))
  df.estim.median = setNames(data.frame(matrix(ncol = ncol(pp_pod), nrow = 0)), colnames(pp_pod))
  df.estim.mean.inside = setNames(data.frame(matrix(ncol = ncol(pp_pod), nrow = 0)), colnames(pp_pod))
  df.estim.mode.inside = setNames(data.frame(matrix(ncol = ncol(pp_pod), nrow = 0)), colnames(pp_pod))
  df.estim.median.inside = setNames(data.frame(matrix(ncol = ncol(pp_pod), nrow = 0)), colnames(pp_pod))
  
  df.estim.quant = setNames(data.frame(matrix(ncol = (ncol(pp_pod)+1), nrow = 0)),
                                              colnames(c(colnames(pp_pod), "quant")))
  
  ################################################################################
  # loop through all single pods from a single pod set
  
  
  # dev quant for quants = c(0.99, 0.95, 0.9, 0.8, 0.5, 0.25, 0.1)
  quant_vec_for_ci <-
    c(0.005,
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
      0.995)
  
  cat("____________\n")
  for (i in 1:nrow(ss_pod)) {
    # inform about progress
    cat("\r", "start: ", i, " of ", nrow(ss_pod), " using ", my_method, " regression")

    # extract target summary statistics
    my_target = as.numeric(ss_pod[i,])
    target_parameters = as.numeric(pp_pod[i,])

    flag <- TRUE;
    tryCatch(
      expr = {
        abc_result = abc(
          target = my_target,
          param = pp_abc,
          sumstat = ss_abc,
          tol = my_tolerance,
          method = my_method,
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
    if (!exists("abc_result")) next
    # message("ended: ", i, " of ", nrow(ss_pod), " using ", my_method, " regression")
    
    # clean the data from na's
    adjusted_values = na.omit(abc_result$adj.values)
    rejection_values = na.omit(abc_result$unadj.values)

    # freeing ram
    rm(abc_result)
    
    row.estim.mean = numeric()
    row.estim.mode = numeric()
    row.estim.median = numeric()
    row.estim.mean.inside = numeric()
    row.estim.mode.inside = numeric()
    row.estim.median.inside = numeric()
    
    row.estim.quants = list()
    
    # loop through each parameter
    posterior_file = paste0(directory_posteriors, "/row_", i, ".pdf")
    pdf(posterior_file, width=ncol(pp_pod)*7)
    par(mfrow=c(1, ncol(pp_pod)))
    for (param_i in 1:ncol(pp_pod)) {
      param_name = colnames(pp_pod)[param_i]
      param_true = target_parameters[param_i]
      
      distr_prior = density(pp_abc[,param_i])
      distr_rejection = density(rejection_values[,param_i])

      tryCatch(
        expr = {
          distr_posterior = density(adjusted_values[,param_i])
        },
        error = function(e) {
          message("* Caught an error on itertion ", i, " for param ", param_i)
        }
      )
      if (!exists("distr_posterior")) {
        distr_posterior = distr_rejection
        adjusted_values = rejection_values
        message(" replaced adj by rej ")
      }

      min_parameter = min(c(distr_prior[["x"]], distr_rejection[['x']], distr_posterior[["x"]]), na.rm=T)
      max_density = max(c(distr_prior[["y"]], distr_rejection[['y']], distr_posterior[["y"]]), na.rm=T)
      
      # draw prior, posterior vs true value

      plot(distr_prior, lwd=5, col = "#181715", main=param_name,
           ylim=c(0, max_density))
      lines(distr_rejection, col = "#CCCCCA", lwd=2)
      lines(distr_posterior, col = "#C19C0F", lwd=5)
      abline(v=param_true, col = "#127049", lwd=5)
      legend("topleft", col=c("#181715", "#CCCCCA", "#C19C0F", "#127049"), pch = 15,
             legend=c("prior", "rejection_posterior", "adjusted_posterior",
                      "true_parameter"))

      # get point estimates and inside_ci
      this_adjusted_values = adjusted_values[,param_i]
      this_adjusted_values <- this_adjusted_values[!is.na(this_adjusted_values)]
      estim.mean = mean(this_adjusted_values)
      estim.mode = get_mode(this_adjusted_values)
      estim.median = median(this_adjusted_values)
      estim.mean.inside = get_confidence_for_bounding(this_adjusted_values, estim.mean)
      estim.mode.inside = get_confidence_for_bounding(this_adjusted_values, estim.mode)
      estim.median.inside = get_confidence_for_bounding(this_adjusted_values, estim.median)
      
      estim.quants = quantile(this_adjusted_values, probs = quant_vec_for_ci)
      rm(this_adjusted_values)

      # append point estimates and inside_ci to the lists
      row.estim.mean = c(row.estim.mean, estim.mean)
      row.estim.mode = c(row.estim.mode, estim.mode)
      row.estim.median = c(row.estim.median, estim.median)
      row.estim.mean.inside = c(row.estim.mean.inside, estim.mean.inside)
      row.estim.mode.inside = c(row.estim.mode.inside, estim.mode.inside)
      row.estim.median.inside = c(row.estim.median.inside, estim.median.inside)
      row.estim.quants[[length(row.estim.quants)+1]] = estim.quants
      
      rm(param_name, param_true, distr_prior, distr_rejection, distr_posterior,
         min_parameter, max_density, estim.mean, estim.mode, estim.median,
         estim.mean.inside, estim.mode.inside, estim.median.inside)
    }
    dev.off()
    
    # row.estim.quants contains the results for a single pod, but all params
    row.estim.quants <- do.call("cbind.data.frame", row.estim.quants)  # make it a data frame
    colnames(row.estim.quants) <- colnames(pp_pod)
    row.estim.quants$quant <- quant_vec_for_ci
    rownames(row.estim.quants) <- NULL
    
    row.estim.quants$pod <- i
        
    df.estim.mean = rbind.data.frame(df.estim.mean, row.estim.mean)
    df.estim.mode = rbind.data.frame(df.estim.mode, row.estim.mode)
    df.estim.median = rbind.data.frame(df.estim.median, row.estim.median)
    df.estim.mean.inside = rbind.data.frame(df.estim.mean.inside, row.estim.mean.inside)
    df.estim.mode.inside = rbind.data.frame(df.estim.mode.inside, row.estim.mode.inside)
    df.estim.median.inside = rbind.data.frame(df.estim.median.inside, row.estim.median.inside)

    df.estim.quant = rbind.data.frame(df.estim.quant, row.estim.quants)
    
    rm(row.estim.mean, row.estim.median, row.estim.mode,
       row.estim.mean.inside, row.estim.mode.inside, row.estim.median.inside,
       row.estim.quants)
  }
  message("made the pods")
  ################################################################################
  
  # get rb for each point estimate
  relative_bias_mean = data.frame(t(get_rb(df.estim.mean, params_pod)))
  relative_bias_mode = data.frame(t(get_rb(df.estim.mode, params_pod)))
  relative_bias_median = data.frame(t(get_rb(df.estim.median, params_pod)))
  message("made rb")
  
  # get rmse for each point estimate
  rmse_mean = data.frame(t(get_rmse(df.estim.mean, params_pod)))
  rmse_mode = data.frame(t(get_rmse(df.estim.mode, params_pod)))
  rmse_median = data.frame(t(get_rmse(df.estim.median, params_pod)))
  message("made rmse")
  
  # get factor x for each point estimate
  factorx_mean = data.frame(t(get_factorx(df.estim.mean, params_pod, 0.95)))
  factorx_mode = data.frame(t(get_factorx(df.estim.mode, params_pod, 0.95)))
  factorx_median = data.frame(t(get_factorx(df.estim.median, params_pod, 0.95)))
  message("made factorx")
  
  # how often within 95% confidence interval
  count_inside_ci95_mean = data.frame(t(get_count_inside_ci(df.estim.mean.inside, 0.95)));
  count_inside_ci95_mode = data.frame(t(get_count_inside_ci(df.estim.mode.inside, 0.95)));
  count_inside_ci95_median = data.frame(t(get_count_inside_ci(df.estim.median.inside, 0.95)));
  colnames(count_inside_ci95_mean) <- colnames(params_pod);
  colnames(count_inside_ci95_mode) <- colnames(params_pod);
  colnames(count_inside_ci95_median) <- colnames(params_pod);
  message("made cov95")
  
  # summarise into long format data frame
  relative_bias_mean$point_estim = factor("mean")
  relative_bias_median$point_estim = factor("median")
  relative_bias_mode$point_estim = factor("mode")
  rmse_mean$point_estim = factor("mean")
  rmse_median$point_estim = factor("median")
  rmse_mode$point_estim = factor("mode")
  factorx_mean$point_estim = factor("mean")
  factorx_median$point_estim = factor("median")
  factorx_mode$point_estim = factor("mode")
  count_inside_ci95_mean$point_estim = factor("mean")
  count_inside_ci95_mode$point_estim = factor("mode")
  count_inside_ci95_median$point_estim = factor("median")
  
  message("long point estim")
  
  relative_bias_mean$eval_mode = factor("rb")
  relative_bias_median$eval_mode = factor("rb")
  relative_bias_mode$eval_mode = factor("rb")
  rmse_mean$eval_mode = factor("rmse")
  rmse_median$eval_mode = factor("rmse")
  rmse_mode$eval_mode = factor("rmse")
  factorx_mean$eval_mode = factor("factorx")
  factorx_median$eval_mode = factor("factorx")
  factorx_mode$eval_mode = factor("factorx")
  count_inside_ci95_mean$eval_mode = factor("ci95")
  count_inside_ci95_mode$eval_mode = factor("ci95")
  count_inside_ci95_median$eval_mode = factor("ci95")
  
  message("long eval mode")
  
  performance = list()
  performance[[1]] <- rbind.data.frame(relative_bias_mean, relative_bias_median, relative_bias_mode,
                                  rmse_mean, rmse_median, rmse_mode, factorx_mean, factorx_median, factorx_mode,
                                  count_inside_ci95_mean, count_inside_ci95_median, count_inside_ci95_mode) %>%
    gather(., key="parameter", value="value", -point_estim, -eval_mode)

  # pack all point estimates into one data.frame
  colnames(df.estim.mean) <- colnames(params_pod)
  colnames(df.estim.mode) <- colnames(params_pod)
  colnames(df.estim.median) <- colnames(params_pod)
  df.estim.mean$pointing = "mean"
  df.estim.mode$pointing = "mode"
  df.estim.median$pointing = "median"
  
  performance[[2]] <- rbind.data.frame(df.estim.mean, df.estim.mode, df.estim.median)
  
  performance[[3]] <- df.estim.quant
  
  message("about to finish the abc performance evaluation")
  
  return(performance)
}