# save.image(file = "My_Object.RData")
# stop(paste0("created image to start manual development:\n  ", getwd()))
# setwd(
#   "~/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/struett/temp_tsabc/a_thaliana/tsabc_at_estimates/"
# )
# load("My_Object.RData")


library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(wesanderson)


# functions to help the execution
{
  # return a numeric from a string that is actually a numeric preceded by any character
  substrX <- function(x) {
    return(as.numeric(substr(x, 2, nchar(x))))
  }
}


# Obtain parameters for estimate
{
  INFILE_POSTERIORS  <- snakemake@input$posteriors
  OUTFILE <- snakemake@output$posterior_plot
  SSCOMP_NO <- as.numeric(snakemake@wildcards$sscomp)
  PLS_COMP <- as.numeric(snakemake@wildcards$pls)
  OBS_NO <- as.numeric(snakemake@wildcards$obs)
}


# read the data
{
  df_posteriors <- readRDS(INFILE_POSTERIORS)
}


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
  sscomp <- SSCOMP_NO
  pls <- PLS_COMP
  pop <- OBS_NO
  
  
  # for plotting
  mcol <- wes_palette("Darjeeling1", 9, "continuous")
  mcol <- c(mcol[1:8], mcol[7:1])
}


# define quantiles to obtain
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


# obtain quantiles
{
  # obtain qants
  all_estimates_of_this_combination <- df_posteriors
  
  
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


# fuse the list into the final result data frame
{
  df <- do.call("rbind.data.frame", quants_list) %>% na.omit()
  df$pop <- as.factor(sapply(df$pop, function(x) {
    return(c("CEU", "IBnr", "Relicts")[as.numeric(x)])
  }))
  df$sscomp <- as.factor(sapply(df$sscomp, function(x) {
    return(c("SFS/LD", "TM_WIN", "SFS/LD/TM_WIN")[as.numeric(x)])
  }))
  df$param <- as.factor(sapply(df$param, function(x) {
    return(c(
      "N_PRES",
      "N_ANC",
      "T_N",
      "SIGMA_PRES",
      "SIGMA_ANC",
      "T_SIGMA"
    )[as.numeric(x)])
  }))
  
  # provide long format
  df <- df %>%
    pivot_longer(
      cols = -c(pop, sscomp, pls, distr, param, summarised),
      names_to = "quantile",
      values_to = "posterior",
      names_transform = list(quantile = substrX)
    )
  df$distr <- as.factor(df$distr)
  df$summarised <- as.factor(df$summarised)
  
  # reorder rows
  df <- df %>%
    arrange(quantile, posterior)
  
  
  # provide index for sampled
  df$sample_index_raw <-
    interaction(df$pop,
                df$sscomp,
                df$pls,
                df$distr,
                df$param,
                df$summarised,
                df$quantile)
  df$sample_index <- NA
  counter <- 0
  total <- length(unique(df$sample_index_raw))
  for (sir in unique(df$sample_index_raw)) {
    nsam <- df %>%
      subset(sample_index_raw == sir) %>%
      nrow()
    df$sample_index[df$sample_index_raw == sir] <- 1:nsam
    counter <- counter + 1
    cat(counter, " of ", total, "\r\r\r\r\r\r\r\r\r\r\r")
  }
  df$sample_index_raw <- NULL
  df$sample_index <- as.factor(df$sample_index)
}


# remove the mean estimate (estimate on the mean of all sumstats)
{
  df_mean <- df %>%
    subset(summarised == "mean")
  
  df <- df %>%
    subset(summarised == "point")
}


# plot average interquantile ranges
# 1) plot rejection posteriors
# 2) plot adjusted posteriors
{
  # 1) plot rejection posteriors
  subdf <- df %>%
    subset(distr == "rej")
  
  p_rej <- subdf %>%
    ggplot(aes(
      sample_index,
      posterior,
      col = as.factor(quantile),
      fill = sample_index
    )) +
    geom_point(show.legend = FALSE,
               shape = 16,
               size = 0.1) +
    # geom_line()+
    facet_grid(sscomp + pop ~ param, scales = "free") +
    scale_color_manual(values = mcol) +
    theme(
      legend.position = "none",
      # axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      aspect.ratio = 1
    ) +
    coord_flip()
}


{
  # 2) plot adjusted posteriors
  subdf <- df %>%
    subset(distr == "adj")
  
  p_adj <- subdf %>%
    ggplot(aes(
      sample_index,
      posterior,
      col = as.factor(quantile),
      fill = sample_index
    )) +
    geom_point(show.legend = FALSE,
               shape = 16,
               size = 0.1) +
    # geom_line()+
    facet_grid(sscomp + pop ~ param, scales = "free") +
    scale_color_manual(values = mcol) +
    theme(
      legend.position = "none",
      # axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      aspect.ratio = 1
    ) +
    coord_flip()
}


{
  mcol <-
    c(rep("gray", 7), wes_palette("Darjeeling1")[1], rep("gray", 7))
  
  df_mean$ismode <- df_mean$quantile == 0.5
  
  p_mean <- df_mean %>%
    subset(ismode == F) %>%
    ggplot(aes(
      x = distr,
      posterior,
      col = as.factor(quantile),
      fill = as.factor(quantile)
    )) +
    geom_point(
      show.legend = FALSE,
      shape = 25,
      size = 2,
      col = "gray",
      fill = "gray"
    ) +
    geom_point(
      data = df_mean %>% subset(ismode == T),
      show.legend = FALSE,
      shape = 25,
      size = 2,
      col = wes_palette("Darjeeling1")[1],
      fill = wes_palette("Darjeeling1")[1]
    ) +
    # geom_abline()+
    facet_grid(sscomp + pop ~ param, scales = "free") +
    # scale_color_manual(values = mcol) +
    # scale_fill_manual(values = mcol) +
    theme(
      legend.position = "none",
      # axis.text.x = element_blank(),
      # axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      aspect.ratio = 0.3,
      axis.line.y = element_blank(),
      strip.text.y = element_blank()
    ) +
    coord_flip()
  
  p_mean
}


p <-
  plot_grid(
    p_adj,
    p_rej,
    p_mean,
    ncol = 1,
    align = T,
    labels = c("loclinear", "rejection", "avg sumstat")
  )


{
  pdf(file = OUTFILE)
  show(p)
  dev.off()
  }
