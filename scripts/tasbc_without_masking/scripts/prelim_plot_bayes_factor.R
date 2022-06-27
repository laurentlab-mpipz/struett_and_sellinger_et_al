


library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(wesanderson)

getwd()

h = "~/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/struett/"
setwd(h)

w = "temp_tsabc/mchoice/tsabc3/data/model_choice/pod_0/sscomp_1/"
setwd(w)

pd <- reticulate::import("pandas")

df <- pd$read_pickle("model_choice_pls_0.table.gzip")

discretize_bayes_factor <- function(bf, do.factor = F) {
  if (bf <= 10 ** 0)
    dbf = "Negative"
  else if (bf <= 10 ** 0.5)
    dbf = "Barely worth mentioning"
  else if (bf <= 10 ** 1)
    dbf = "Substantial"
  else if (bf <= 10 ** 1.5)
    dbf = "Strong"
  else if (bf <= 10 ** 2)
    dbf = "Very strong"
  else if (bf > 10 ** 2)
    dbf =  "Decisive"
  else
    stop(paste0(c("unknown bayes: ", as.character(bf)), collapse = ""))
  
  dbf <- factor(
    dbf,
    levels = c(
      "Negative",
      "Barely worth mentioning",
      "Substantial",
      "Strong",
      "Very strong",
      "Decisive"
    )
  )
  
  if (do.factor) {
    return(dbf)
  } else {
    return(as.numeric(dbf)-1)
  }
}

df$discrete_bayes_factor <-
  sapply(df$bayes_factor, discretize_bayes_factor)

df %>%
  ggplot(aes(x = regression_method,
             y = discrete_bayes_factor,
             fill = regression_method)) +
  geom_point(aes(col = regression_method),
             position = position_jitter(width = 0.2, height = 0.1))+
  scale_y_continuous(breaks = 0:6)

df %>%
  ggplot(aes(x = regression_method,
             y = discrete_bayes_factor,
             fill = regression_method))+
  geom_violin(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.3, height = 0.05))




