# save.image(file = "My_Object.agg.plot.RData")
# stop("created image to start manual development")
# setwd(
#   "~/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/struett/temp_tsabc/a_thaliana/tsabc_at_estimates/"
# )
# load("My_Object.agg.plot.RData")


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
  INFILE_POSTERIORS_LIST  <- snakemake@input$agg_post
  OUTFILE <- snakemake@output$agg_sumplot
  OUTFILE_TABLE <- snakemake@output$table
}


# read the data
{
  df_posteriors_list <- list()
  for (i in 1:length(INFILE_POSTERIORS_LIST))
    df_posteriors_list[[i]] <- readRDS(INFILE_POSTERIORS_LIST[i])
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
  for (i in 1:length(df_posteriors_list)) {
    cat("analyzing experiment ",
        i,
        " of ",
        length(df_posteriors_list),
        "\n")
    
    
    # obtain pop, sscomp, pls
    infile_name <- INFILE_POSTERIORS_LIST[i]
    info_from_filename <-
      strsplit(infile_name, "_|/|\\.", perl = T)[[1]]
    sscomp <-
      as.numeric(info_from_filename[which(info_from_filename == "statset") + 1])
    pop <-
      as.numeric(info_from_filename[which(info_from_filename == "pop") + 1])
    pls <-
      as.numeric(info_from_filename[which(info_from_filename == "pls") + 1])
    
    
    # obtain qants
    all_estimates_of_this_combination <- df_posteriors_list[[i]]
    
    
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
}


# fuse the list into the final result data frame
{
  df <- do.call("rbind.data.frame", quants_list) %>% na.omit()
  df$pop <- factor(sapply(df$pop, function(x) {
    # + 1 because in the filename the index is zero-based
    return(c("CEU", "IBnr", "Relicts")[as.numeric(x) + 1])
  }),
  levels = c("CEU", "IBnr", "Relicts"))
  df$sscomp <- factor(sapply(df$sscomp, function(x) {
    return(c("SFS/LD", "TM_WIN", "SFS/LD/TM_WIN")[as.numeric(x)])
  }),
  levels = c("SFS/LD", "TM_WIN", "SFS/LD/TM_WIN"))
  df$param <- factor(
    sapply(df$param, function(x) {
      return(c(
        "N_PRES",
        "N_ANC",
        "T_N",
        "SIGMA_PRES",
        "SIGMA_ANC",
        "T_SIGMA"
      )[as.numeric(x)])
    }),
    levels = c(
      "N_PRES",
      "N_ANC",
      "T_N",
      "SIGMA_PRES",
      "SIGMA_ANC",
      "T_SIGMA"
    )
  )
  
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
      # strip.text.y = element_blank()
    ) +
    coord_flip()
  
  # p_mean
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


# create the table with the mode and two different CIs
{
  # depend the height of the plot on the number of populations
  n_pop <- length(unique(df$pop))
  n_pls <- length(unique(df$pls))
  n_sscomp <- length(unique(df$sscomp))
  
  pdf(
    file = OUTFILE,
    height = 7 * n_pop * n_pls * n_sscomp / 3,
    useDingbats = F
  )
  show(p)
  dev.off()
  }


{
  # get the values for the main table 1
  cat("", file = OUTFILE_TABLE)
  
  
  # estims from average posteriors
  sdf <- df %>%
    subset(distr == "adj") %>%
    subset(sscomp == "SFS/LD/TM_WIN") %>%
    subset(summarised == "point")
  
  
  # subset for each condition
  for (a in unique(df$param)) {
    cat("\n", file = OUTFILE_TABLE, append = T)
    for (b in sort(unique(df$pop))) {
      for (c in unique(df$pls)) {
        this_ssdf <- sdf %>%
          subset(param == a) %>%
          subset(pop == b) %>%
          subset(pls == c)
        
        
        # check for proper subsetting
        stopifnot(length(unique(this_ssdf$pop)) == 1)
        stopifnot(length(unique(this_ssdf$param)) == 1)
        stopifnot(length(unique(this_ssdf$pls)) == 1)
        
        
        # get values
        cat(
          paste0(c(
            b,
            paste0(c, " pls"),
            "of",
            levels(this_ssdf$sscomp)[unique(this_ssdf$sscomp)],
            a
          ), collapse = "\t"),
          
          
          # this is the result
          this_ssdf$posterior[this_ssdf$quantile == 0.5] %>% mean(),
          this_ssdf$posterior[this_ssdf$quantile == 0.25] %>% mean(),
          this_ssdf$posterior[this_ssdf$quantile == 0.75] %>% mean(),
          this_ssdf$posterior[this_ssdf$quantile == 0.025] %>% mean(),
          this_ssdf$posterior[this_ssdf$quantile == 0.975] %>% mean(),
          "\n",
          sep = "\t",
          file = OUTFILE_TABLE,
          append = T
        )
      }
    }
  }
}


{
  # get the values for the main table 1
  
  
  # estims from average posteriors
  sdf <- df_mean %>%
    subset(distr == "adj") %>%
    subset(sscomp == "SFS/LD/TM_WIN") %>%
    subset(summarised == "mean")
  
  
  # subset for each condition
  for (a in unique(df$param)) {
    cat("\n", file = OUTFILE_TABLE, append = T)
    for (b in sort(unique(df$pop))) {
      for (c in unique(df$pls)) {
        this_ssdf <- sdf %>%
          subset(param == a) %>%
          subset(pop == b) %>%
          subset(pls == c)
        
        
        # check for proper subsetting
        stopifnot(length(unique(this_ssdf$pop)) == 1)
        stopifnot(length(unique(this_ssdf$param)) == 1)
        stopifnot(length(unique(this_ssdf$pls)) == 1)
        
        
        # get values
        cat(
          paste0(c(
            b,
            paste0(c, " pls"),
            "of",
            levels(this_ssdf$sscomp)[unique(this_ssdf$sscomp)],
            a
          ), collapse = "\t"),
          
          
          # this is the result
          this_ssdf$posterior[this_ssdf$quantile == 0.5],
          this_ssdf$posterior[this_ssdf$quantile == 0.25],
          this_ssdf$posterior[this_ssdf$quantile == 0.75],
          this_ssdf$posterior[this_ssdf$quantile == 0.025],
          this_ssdf$posterior[this_ssdf$quantile == 0.975],
          "\n",
          sep = "\t",
          file = OUTFILE_TABLE,
          append = T
        )
      }
    }
  }
}
