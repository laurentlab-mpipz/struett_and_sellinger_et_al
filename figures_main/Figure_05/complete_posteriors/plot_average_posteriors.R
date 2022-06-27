library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(wesanderson)


# prepare some parameters
{
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


# read posterior results list
{
  # find files
  infile_list <-
    list.files("./", ".RDS", recursive = T, full.names = T)
  
  
  # read in files
  result_list <- list()
  for (i in 1:length(infile_list)) {
    result_list[[infile_list[i]]] <- readRDS(infile_list[i])
  }
}


# prepare a data frame dummy; not actually needed, but good for overview
{
  # prepare a dataframe for final plotting
  quants_list <- list()
  quants_list[[1]] <- data.frame(
    pop = NA,
    sscomp = NA,
    pls = NA,
    distr = NA,
    param = NA,
    summarised = NA,
    "0.005" = NA,
    "0.025" = NA,
    "0.05" = NA,
    "0.1" = NA,
    "0.25" = NA,
    "0.375" = NA,
    "0.45" = NA,
    "0.5" = NA,
    "0.55" = NA,
    "0.625" = NA,
    "0.75" = NA,
    "0.9" = NA,
    "0.95" = NA,
    "0.975" = NA,
    "0.995" = NA
  )
  quant_index <- 2
}


# create the loop to go through each file
{
  for (result_name in names(result_list)) {
    # get the pop, pls and statset
    values <-
      strsplit(result_name, "/|_|\\.", perl = T)[[1]][c(5, 7, 9)]
    sscomp <- values[1] %>% as.numeric()
    pop <- values[2] %>% as.numeric()
    pls <- values[3] %>% as.numeric()
    
    
    # obtain qants
    all_estimates_of_this_combination <- result_list[[result_name]]
    
    
    # get rejection, adjusted, then remove from list, to prepare for later loop
    allpar_prior <- all_estimates_of_this_combination[["prior"]]
    allpar_posterior_adjusted_mean <-
      all_estimates_of_this_combination[["meanEstim"]][[2]]
    allpar_posterior_rejection_mean <-
      all_estimates_of_this_combination[["meanEstim"]][[1]]
    all_estimates_of_this_combination[["prior"]] <- NULL
    all_estimates_of_this_combination[["meanEstim"]] <- NULL
    
    
    # loop through the params of the mean posteriors
    for (par_ix in 1:ncol(allpar_prior)) {
      # obtain the parameter
      prior <- allpar_prior[, par_ix]
      posterior_adjusted_mean <-
        allpar_posterior_adjusted_mean[, par_ix]
      posterior_rejection_mean <-
        allpar_posterior_rejection_mean[, par_ix]
      
      
      # calculate quantiles
      quants_adj <-
        quantile(posterior_adjusted_mean, probs = quant_vec_for_ci)
      quants_rej <-
        quantile(posterior_rejection_mean, probs = quant_vec_for_ci)
      
      
      # add adjusted quantiles to the result list
      quants_list[[quant_index]] <- data.frame(
        pop = pop,
        sscomp = sscomp,
        pls = pls,
        distr = "adj",
        param = par_ix,
        summarised = "mean",
        "0.005" = quants_adj[1],
        "0.025" = quants_adj[2],
        "0.05" = quants_adj[3],
        "0.1" = quants_adj[4],
        "0.25" = quants_adj[5],
        "0.375" = quants_adj[6],
        "0.45" = quants_adj[7],
        "0.5" = quants_adj[8],
        "0.55" = quants_adj[9],
        "0.625" = quants_adj[10],
        "0.75" = quants_adj[11],
        "0.9" = quants_adj[12],
        "0.95" = quants_adj[13],
        "0.975" = quants_adj[14],
        "0.995" = quants_adj[15]
      )
      quant_index <- quant_index + 1
      
      
      # add unadj quantiles to the result list
      quants_list[[quant_index]] <- data.frame(
        pop = pop,
        sscomp = sscomp,
        pls = pls,
        distr = "rej",
        param = par_ix,
        summarised = "mean",
        "0.005" = quants_rej[1],
        "0.025" = quants_rej[2],
        "0.05" = quants_rej[3],
        "0.1" = quants_rej[4],
        "0.25" = quants_rej[5],
        "0.375" = quants_rej[6],
        "0.45" = quants_rej[7],
        "0.5" = quants_rej[8],
        "0.55" = quants_rej[9],
        "0.625" = quants_rej[10],
        "0.75" = quants_rej[11],
        "0.9" = quants_rej[12],
        "0.95" = quants_rej[13],
        "0.975" = quants_rej[14],
        "0.995" = quants_rej[15]
      )
      quant_index <- quant_index + 1
    }
    
    
    # loop through the remaining posteriors
    for (estim_index in 1:length(all_estimates_of_this_combination)) {
      cat(
        estim_index,
        " of ",
        length(all_estimates_of_this_combination),
        " of ",
        result_name,
        "\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r"
      )
      # obtain the parameter
      allpar_posterior_adjusted <-
        all_estimates_of_this_combination[[estim_index]][[2]]
      allpar_posterior_rejection <-
        all_estimates_of_this_combination[[estim_index]][[1]]
      
      
      # loop through the params
      for (par_ix in 1:ncol(allpar_prior)) {
        # obtain the parameter
        posterior_adjusted <- allpar_posterior_adjusted[, par_ix]
        posterior_rejection <- allpar_posterior_rejection[, par_ix]
        
        
        # calculate quantiles
        quants_adj <-
          quantile(posterior_adjusted, probs = quant_vec_for_ci)
        quants_rej <-
          quantile(posterior_rejection, probs = quant_vec_for_ci)
        
        
        # add adjusted quantiles to the result list
        quants_list[[quant_index]] <- data.frame(
          pop = pop,
          sscomp = sscomp,
          pls = pls,
          distr = "adj",
          param = par_ix,
          summarised = "point",
          "0.005" = quants_adj[1],
          "0.025" = quants_adj[2],
          "0.05" = quants_adj[3],
          "0.1" = quants_adj[4],
          "0.25" = quants_adj[5],
          "0.375" = quants_adj[6],
          "0.45" = quants_adj[7],
          "0.5" = quants_adj[8],
          "0.55" = quants_adj[9],
          "0.625" = quants_adj[10],
          "0.75" = quants_adj[11],
          "0.9" = quants_adj[12],
          "0.95" = quants_adj[13],
          "0.975" = quants_adj[14],
          "0.995" = quants_adj[15]
        )
        quant_index <- quant_index + 1
        
        
        # add unadj quantiles to the result list
        quants_list[[quant_index]] <- data.frame(
          pop = pop,
          sscomp = sscomp,
          pls = pls,
          distr = "rej",
          param = par_ix,
          summarised = "point",
          "0.005" = quants_rej[1],
          "0.025" = quants_rej[2],
          "0.05" = quants_rej[3],
          "0.1" = quants_rej[4],
          "0.25" = quants_rej[5],
          "0.375" = quants_rej[6],
          "0.45" = quants_rej[7],
          "0.5" = quants_rej[8],
          "0.55" = quants_rej[9],
          "0.625" = quants_rej[10],
          "0.75" = quants_rej[11],
          "0.9" = quants_rej[12],
          "0.95" = quants_rej[13],
          "0.975" = quants_rej[14],
          "0.995" = quants_rej[15]
        )
        quant_index <- quant_index + 1
      }
    }
  }
}



# fuse the list into the final result data frame
{
  df <- do.call("rbind.data.frame", quants_list) %>% na.omit()
  df$pop <- as.factor(df$pop)
  levels(df$pop) <- c("CEU", "IBnr", "Relicts")
  df$sscomp <- as.factor(df$sscomp)
  levels(df$sscomp) = c("SFS/LD", "TM_WIN", "SFS/LD/TM_WIN")
  df$param <- as.factor(df$param)
  levels(df$param) <-
    c("N_PRES",
      "N_ANC",
      "T_N",
      "SIGMA_PRES",
      "SIGMA_ANC",
      "T_SIGMA")
}


# explore by plotting
{
  # subset to wanted combination of points
  subdf <- df  %>%
    subset(summarised == "point")
  
  
  # plot
  subdf %>%
    # subset(sscomp == "SFS/LD/TM_WIN") %>%
    # subset(param %in% c("SIGMA_PRES")) %>%
    ggplot(aes(x = 1, y = X0.5, col = distr)) +
    geom_point(position = position_jitter(height = 0)) +
    facet_grid(param ~ sscomp + pop, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
}


# obtain the means
{
  subdf <- df  %>%
    subset(sscomp == "SFS/LD/TM_WIN") %>%
    subset(summarised == "point") %>%
    subset(distr == "adj") 
    # subset(param == "T_SIGMA")
  
  
  # mode
  for (this_pop in unique(subdf$pop)) {
    cat(this_pop, "\n")
    for (this_param in unique(subdf$param)) {
      this_df <- subdf %>%
        subset(pop == this_pop) %>%
        subset(param == this_param)
      cat("  ", this_param, " ", mean(this_df$X0.5), "\n")
    }
  }
  
  
  # lower
  for (this_pop in unique(subdf$pop)) {
    cat("lower: ", this_pop, "\n")
    for (this_param in unique(subdf$param)) {
      this_df <- subdf %>%
        subset(pop == this_pop) %>%
        subset(param == this_param)
      cat("  ", this_param, " ", mean(this_df$X0.25), "\n")
    }
  }
  
  # higher
  for (this_pop in unique(subdf$pop)) {
    cat("higher: ", this_pop, "\n")
    for (this_param in unique(subdf$param)) {
      this_df <- subdf %>%
        subset(pop == this_pop) %>%
        subset(param == this_param)
      cat("  ", this_param, " ", mean(this_df$X0.75), "\n")
    }
  }
}


# botain the values for the estimates on the mean stats
{
  subdf <- df  %>%
    subset(sscomp == "SFS/LD/TM_WIN") %>%
    subset(summarised == "mean") %>%
    subset(distr == "adj") 
  # subset(param == "T_SIGMA")
  
  
  # mode
  for (this_pop in unique(subdf$pop)) {
    cat(this_pop, "\n")
    for (this_param in unique(subdf$param)) {
      this_df <- subdf %>%
        subset(pop == this_pop) %>%
        subset(param == this_param)
      cat("  ", this_param, " ", mean(this_df$X0.5), "\n")
    }
  }
  
  
  # lower
  for (this_pop in unique(subdf$pop)) {
    cat("lower: ", this_pop, "\n")
    for (this_param in unique(subdf$param)) {
      this_df <- subdf %>%
        subset(pop == this_pop) %>%
        subset(param == this_param)
      cat("  ", this_param, " ", mean(this_df$X0.25), "\n")
    }
  }
  
  # higher
  for (this_pop in unique(subdf$pop)) {
    cat("higher: ", this_pop, "\n")
    for (this_param in unique(subdf$param)) {
      this_df <- subdf %>%
        subset(pop == this_pop) %>%
        subset(param == this_param)
      cat("  ", this_param, " ", mean(this_df$X0.75), "\n")
    }
  }
}

