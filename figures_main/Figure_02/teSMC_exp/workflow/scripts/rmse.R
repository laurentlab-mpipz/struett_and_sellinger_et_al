
# # save image for debugging
# save.image(); stop("dadada+++++++++++++++++++")
# setwd("../..")
# load(".RData")

library(tidyr)
library(cumSeg)

infile = snakemake@input$allinf
outfile = snakemake@output$rmse
infmodels = snakemake@config$teSMC$inference_models

sigma_d = snakemake@config["selfing_rates_backward_in_time"][[1]]
d = as.numeric(snakemake@wildcards$demography)

get_transition_time_from_demogr_scenario <- function(d, sigma_d) {
  t = sigma_d[(d*4+1):(d*4+4)][4]
  return(t)
}

t_sigma_true = get_transition_time_from_demogr_scenario(d, sigma_d)
demographic_inference <- read.csv(infile, row.names = 1)
wildcards_for_subsetting <- stringi::stri_remove_empty(names(snakemake@wildcards))

make_performance_analysis_for_tsigma <- function(
  tsigma_true, wildcard_names_for_subsetting, demographic_inference, wc_config) {

  # subset to single evolving sigma; this has implicitly already happened
  subdf <- demographic_inference
  for (wc in wildcard_names_for_subsetting) {
    if (is.na(as.numeric(wc_config[[wc]]))) { # this is for the non-numerics
      subdf <- subdf[as.character(subdf[[wc]]) == as.character(wc_config[[wc]]),]
    } else { # this is for numerics
      subdf <- subdf[as.numeric(subdf[[wc]]) == as.numeric(wc_config[[wc]]),]
    }
  }
  
  stopifnot(nrow(subdf) == nrow(demographic_inference))
  
  # define point estimate
  tr = c()
  tr_name = c()
  for (d in unique(subdf$inference_model)) {
    ddf <- subdf %>% subset(inference_model==d)
    t <- ddf$t
    sigma_t <- ddf$sigma

    if (d == "Constant") {
      tr = c(tr, max(t))
      tr_name = c(tr_name, d)
    } else if (d == "Free") {
      a = jumpoints(y = sigma_t, x = t, k = 1)
      fitted_sigma = a$fitted.values
      stopifnot(length(a$est.means)==2)
      tr = c(tr, t[min(which(fitted_sigma == a$est.means[2]))])
      tr_name = c(tr_name, d)
    } else if (d %in% c("GivenTransition", "OneTransition")) {
      tr = c(tr, t[min(which(sigma_t == unique(sigma_t)[2]))])
      tr_name = c(tr_name, d)
    } else {
      stop("something went wrong, there is an unknown ")
    }
  }
  
  # provide wildcards for dataframe
  df = data.frame(t_sigma_inferred = tr,
                  t_sigma_true = tsigma_true, 
                  infmodel = tr_name)
  
  for (wc in wildcard_names_for_subsetting) {
    if (is.na(as.numeric(wc_config[[wc]]))) { # this is for the non-numerics
      df[[wc]] = wc_config[[wc]]
    } else { # this is for numerics
      df[[wc]] = as.numeric(wc_config[[wc]])
    }
  }
  
  return(df)
  
}

results <- make_performance_analysis_for_tsigma(
  t_sigma_true, wildcards_for_subsetting, demographic_inference, snakemake@wildcards)

write.csv(x=results, file=outfile)
