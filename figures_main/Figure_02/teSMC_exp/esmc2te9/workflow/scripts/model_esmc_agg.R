
# save image for debugging
# save.image(); stop("\n_____________\n saved image\n=============\n\n")
# setwd("../..")
# load(".RData")
# load("agg_teSMC.RData")

message("The readout of the wildcards and the inference model from the file")
message("name unfortunately is very hard coded. Please try to avoid to change")
message("the rules provided input; or at least make sure that it works fine")

infiles = snakemake@input$full
outfile = snakemake@output$csv
mu <- as.numeric(as.character(snakemake@wildcards$mutation_rate))
wildcards <- snakemake@wildcards

# functions
get_csv_from_full_result <- function(full_result, mu, wildcards) {
  test_ = full_result[[1]] # not further used
  
  total = unlist(full_result[[2]])
  s_ex = full_result[[3]]
  mu_ex = unlist(full_result[[4]])
  Tc = unlist(full_result[[5]])
  
  Pop_ = mu_ex / mu
  test_ = total
  
  t = Tc * Pop_
  sigma = s_ex
  N = 10**(log10((test_) * 0.5 * Pop_))
  
  infmodel <- rev(strsplit(full_result$filename, "\\.")[[1]])[2]
  
  df <- cbind.data.frame(
    t, sigma, N, infmodel, as.character(full_result$filename)
  )
  names(df) <- c("t", "sigma", "pop_size", "inference_model", "file")
  
  # add column for each wildcard
  for (i in 1:length(wildcards)) {
    wname = names(wildcards)[i]
    if (wname != "") {
      # if is numeric, then make it numeric
      if (!is.na(as.numeric(wildcards[wname]))) {
        df[wname] = as.numeric(wildcards[wname])
      } else {
        df[wname] = wildcards[wname]
      }
      
    } 
  }
  
  return(df)
}

# proc
results_list <- lapply(infiles, function(f) {
  r = readRDS(f)
  r[["filename"]] = f
  return(r)
})

df <- do.call(rbind.data.frame, lapply(results_list, function(f) {
  get_csv_from_full_result(f, mu, wildcards)
}))

PLOT = FALSE
if (PLOT) {
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
  library(tidyverse)
  library(wesanderson)
  
  # log breaks and labels
  log_breaks = sapply(10**(-12:10), function(x) x*(1:9)) %>% as.vector() %>% unique() %>% sort()
  `%nin%` = Negate(`%in%`)
  log_labels = log_breaks; log_labels[log_labels %nin% 10**(-12:10)] = ""
  
  mtheme <- theme(
    aspect.ratio = 0.707/2,
    legend.position = "right"
  )
  
  mscale_x_continuous <- scale_x_continuous(
    breaks = log_breaks,
    labels = log_labels,
    trans = "log10",
    limits = c(1e4, NA),
    expand = expansion(mult = c(0, 0))
  )
  
  mlabs <- labs(x="generations ago")
  
  mscale_color_manual <- scale_color_manual(
    values = wes_palette("Darjeeling1")
  )
  
  d <- df %>%
    ggplot(aes(t+1, sigma, col=inference_model, fill=file))+
    geom_step()+
    mscale_x_continuous+
    mscale_color_manual+
    mlabs+labs(y="Selfing rate")+
    mtheme
  
  s <- df %>%
    ggplot(aes(t+1, pop_size, col=inference_model, fill=file))+
    geom_step()+
    mscale_x_continuous+
    scale_y_continuous(breaks = log_breaks, labels = log_labels, trans = "log10")+
    mscale_color_manual+
    mlabs+labs(y="Population size")+
    mtheme
  
  p <- plot_grid(d, s, ncol = 1, align = T, labels = "AUTO")+
    theme(legend.position = "bottom")
  
  show(p)
}

write.csv(x = df, file = outfile)
cat("successful aggregation\n")
