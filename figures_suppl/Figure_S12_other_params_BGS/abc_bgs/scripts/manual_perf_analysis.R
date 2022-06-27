
library(reticulate)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(wesanderson)
message("loaded libraries")

setwd("~/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/struett/temp_tsabc/data_3_2/tsabc2/")
getwd()

pd <- reticulate::import("pandas")

df <- pd$read_pickle("data/table/performance/all_pod_performance.table.gzip",
                  compression = "infer")
df = na.omit(df)

df$value <- as.numeric(df$value)
df$plsComp <- as.numeric(df$plsComp)


setwd("~/Desktop/")

df$point_estim %>% unique()
df$point_estim = factor(df$point_estim, levels = c("mode", "median", "mean", "None"))
df$eval_mode %>% unique()
df$eval_mode = factor(df$eval_mode, levels = c("rb", "rmse", "ci95", "factorx", "None"))
df<-df[which(df$eval_mode!="None"),]
df$parameter=as.factor(df$parameter)
levels(df$parameter) <- c("N_r", "N_a", "t_N", "sigma_r", "sigma_a", "t_sigma")
df$plsComp = factor(as.character(df$plsComp),
                    levels = c(df$plsComp %>% unique() %>% sort()))
df$sscomp = as.factor(df$sscomp)
levels(df$sscomp) <- c(
  "sfs/ld",
  "tl_true",
  "tm_true",
  "tm_win",
  "sfs/ld/tm_win")

mcol = wes_palette("Darjeeling1", 5, "continuous")[c(1, 3:5)]

# add true times
mt = c(1000, 10000, 20000, 30000, 40000, 50000, 60000, 1000, 10000, 20000, 30000, 40000, 50000, 60000)
ms = c(0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9)
df$t_sigma_true = 0
df$s_a_true = 0

for (i in 0:(length(mt)-1)) {
  pod_no = as.character(i)
  df$t_sigma_true[df$t_sigma_true == pod_no] <- mt[(i+1)]
  df$s_a_true[df$s_a_true == pod_no] <- ms[(i+1)]
}

# p <- df %>%
#   subset(eval_mode %in% c("rb", "rmse")) %>%
#   subset(point_estim %in% c("mode", "median")) %>%
#   ggplot(aes(plsComp, value, fill=interaction(eval_mode, point_estim)))+
#   geom_col(position = position_dodge())+
#   scale_fill_manual(values = mcol)+
#   facet_grid(sscomp+eval_mode~parameter+pod, scales = "free_y")
# show(p)

df %>%
  subset(eval_mode %in% c("rb", "rmse")) %>% 
  subset(point_estim %in% c("mode", "median")) %>%
  subset(pod == "5") %>%
  subset(parameter %in% c("t_sigma", "t_N")) %>%
  # subset(!(sscomp %in% c("tm_true"))) %>%
  subset(sscomp %in% c("tl_true", "tm_win", "sfs/ld/tm_win")) %>%
  # nrow()
  
  ggplot(aes(plsComp, value, fill=interaction(eval_mode, point_estim)))+
  geom_col(position = position_dodge())+
  scale_fill_manual(values = mcol)+
  facet_grid(sscomp~eval_mode+parameter+pod, scales = "fixed")+
  theme()
  






