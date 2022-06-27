
library(tidyr)
library(tibble)

setwd("bgs_good_window/")
infiles = list.files(pattern = ".csv")
outfile = "point_estims.csv"



message("we define the threshold of 'sigma = 0.5' as the time of transition to selfing.")

df_list <- list(); counter <- 1; for (x in infiles) {
  df <- read.csv(x)
  df$fname <- x
  names(df)[names(df) == "X"] <- "rep"
  names(df)[2:(ncol(df)-1)] <- 1:40
  df <- df %>%
    pivot_longer(-c(rep, fname),
                 names_to = "time_window",
                 values_to = "xi") %>%
    separate(
      col = "fname",
      sep = "_+",
      into = c(
        "torm1",
        "torm2",
        "dtype",
        "optim_method",
        "optim_model",
        "torm3",
        "torm4",
        "torm5",
        "torm6",
        "dfe",
        "torm7",
        "tsigma",
        "torm8",
        "torm9",
        "masked",
        "torm10"
      )
    )
  
  # remove no content columns
  df[, grep("torm", names(df))] <- NULL
  
  # correct mask information
  df$masked <- df$masked == "masked"
  
  df$file <- x
  
  df_list[[counter]] <- df
  counter <- counter + 1
}; df <- do.call(rbind.data.frame, df_list)

# find the correct rescaled time based on the time window
rescaled_time_windows <- df %>%
  subset(dtype == "t")  # subset for time windows
df <- df %>%
  subset(dtype != "t")  # subset for population size and selfing

rescaled_time <- pbapply(df, 1, function(x) {
  rep_ <- x[["rep"]]
  dtype_ <- x[["dtype"]]
  optim_method_ <- x[["optim_method"]]
  optim_model_ <- x[["optim_model"]]
  dfe_ <- x[["dfe"]]
  tsigma_ <- x[["tsigma"]]
  masked_ <- x[["masked"]]
  time_window_ <- x[["time_window"]]
  xi_ <- x[["xi"]]
  
  rescaled_time <- (
    rescaled_time_windows %>%
      subset(rep == as.numeric(rep_)) %>%
      subset(optim_method == optim_method_) %>%
      subset(optim_model == optim_model_) %>%
      subset(dfe == dfe_) %>%
      subset(tsigma == tsigma_) %>%
      subset(masked == masked_) %>%
      subset(time_window == time_window_)
  )[["xi"]]
  stopifnot(length(rescaled_time) == 1)
  
  return(rescaled_time)
})
df$rescaled_time <- rescaled_time

# reformat data table
df$time_window <- as.numeric(df$time_window)
names(df)[names(df) == "xi"] <- "value"
df$value[df$dtype == "p"] <- 10**(df$value[df$dtype == "p"])  # Thibaut always talks about log10 values
df$optim_method <- factor(df$optim_method, levels = c("LH", "BW"))
df$tsigma <- as.numeric(df$tsigma)
df$t <- df$rescaled_time


# get single estim loop
uniq_files <- unique(df$file[df$dtype == "s"])

row_collector = list()
row_indexer = 1
for (u in uniq_files) {
  sdf_ <- df %>%
    subset(file == u)
  
  
  uniq_rep <- sdf_$rep %>% unique()
  
  for (urep in uniq_rep) {
    sdf <- sdf_ %>%
      subset(rep == urep)
    
    sdf <- sdf[order(sdf$t), ]
    
    if (length(which(sdf$value <= 0.5)) > 0) {
      this_tsigma <- min(sdf[which(sdf$value <= 0.5),]$t)
    } else {
      this_tsigma <- Inf
    }
    
    new_row <- tibble(
      rep = sdf$rep %>% unique(),
      dtype = sdf$dtype %>% unique(),
      optim_method = sdf$optim_method %>% unique(),
      optim_model = sdf$optim_model %>% unique(),
      dfe = sdf$dfe %>% unique(),
      tsigma = sdf$tsigma %>% unique(),
      masked = sdf$masked %>% unique(),
      file = sdf$file %>% unique(),
      tsigma_estim = this_tsigma
    )
    
    row_collector[[row_indexer]] = new_row
    row_indexer = row_indexer + 1
    
    rm(new_row, sdf, this_tsigma)
  }
  
}

tsigma_estim <- do.call("rbind.data.frame", row_collector)

write.csv(x=tsigma_estim, file = outfile)
