
library(tidyr)
library(tibble)

# wd <- getwd()
# cat("working dir:\n")
# cat(wd, "\n")
# save.image()
# stop(paste0(rep("#", 600), collapse = ""))
# setwd(dir = paste0("/Users/struett/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/",
#                    "struett/temp_mspts/generate_ts_data_07_manuscript_8hap/",
#                    collapse = "")
# )
# load(file = ".RData")

infile = snakemake@input$rds
outfile = snakemake@output$csv

message("we define the threshold of 'sigma = 0.5' as the time of transition to selfing.")

df <- readRDS(infile) %>%
  subset(inference_model %in% c("OneTransition", "Free"))
# df$t_sigma <- factor(df$t_sigma, levels = sort(unique(df$t_sigma)))

# get single estim loop
uniq_files <- unique(df$file)

row_collector = list()
row_indexer = 1
for (u in uniq_files) {
  sdf <- df %>%
    subset(file == u)
  
  sdf <- sdf[order(sdf$t), ]
  
  if (length(which(sdf$sigma <= 0.5)) > 0) {
    this_tsigma <- min(sdf[which(sdf$sigma <= 0.5),]$t)
  } else {
    this_tsigma <- Inf
  }
  
  new_row <- tibble(
    infmodel = sdf$inference_model %>% unique(),
    file = sdf$file %>% unique(),
    demography = sdf$demography %>% unique(),
    replication = sdf$replication %>% unique(),
    nsam = sdf$sample_size %>% unique(),
    mu = sdf$mutation_rate %>% unique(),
    r = sdf$recombination_rate %>% unique(),
    chromlen = sdf$chromosome_length %>% unique(),
    tsigma = sdf$t_sigma %>% unique(),
    tsigma_estim = this_tsigma
  )
  
  row_collector[[row_indexer]] = new_row
  row_indexer = row_indexer + 1
  
  rm(new_row, sdf, this_tsigma)
}

df <- do.call("rbind.data.frame", row_collector)

write.csv(x=df, file = outfile)
