
# save image for debugging
# save.image(); stop("\n_____________\n saved image\n=============\n\n")
# setwd("../..")
# load(".RData")
# load("agg_teSMC.RData")

message("The readout of the values and specifications from the given file")
message("name unfortunately is very hard coded. Please try to avoid to change")
message("the rules provided input; or at least make sure that it works fine")
message("the time of the selfing rate change is also pretty much hardcoded")
 
infiles = unlist(snakemake@input)
outfile = snakemake@output$rds
sigma_d = snakemake@config["selfing_rates_backward_in_time"][[1]]

get_transition_time_from_demogr_scenario <- function(file_name, sigma_d) {
  d = as.numeric(strsplit(strsplit(file_name, "/")[[1]][3], "_")[[1]][3])
  t = sigma_d[(d*4+1):(d*4+4)][4]
  return(c(d,t))
}

results_list <- lapply(infiles, function(f) {
  df = read.csv(f, row.names = 1)
  t_time_and_d_scenario = get_transition_time_from_demogr_scenario(f,sigma_d)
  d = t_time_and_d_scenario[1]
  t = t_time_and_d_scenario[2]
  df$demographic_scenario=d
  df$t_sigma=t
  return(df)
})

df <- do.call(rbind.data.frame, results_list)

saveRDS(df, outfile)

