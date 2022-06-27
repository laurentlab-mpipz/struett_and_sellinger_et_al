setwd("~/Desktop/")

library(reticulate)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(wesanderson)

pd <- import("pandas")
'%notin%' <- Negate('%in%')

# read raw data and plot execution times
if (F) {
  # read raw data
  
  path = "/netscratch/dep_tsiantis/grp_laurent/struett/temp_tsabc/data_1/tsabc2/data/table"
  files = list.files(path=path, pattern = "pod_[0-9]+.table.gzip")
  full_paths = paste(path, files, sep="/")
  
  df_sim <- pd$read_pickle(paste0(path, "/", "params_and_sumstats.table.gzip"));
  df_sim$pod <- "sim"
  df_pod <- do.call("rbind.data.frame", lapply(full_paths, function(x){
    t = pd$read_pickle(x)
    t$path = x
    return(t)
  }))
  
  # get pod number
  a = df_pod %>%
    separate(col=path, into=c("torm", "pod_prelim"), sep="_pod_") %>%
    separate(col=pod_prelim, into=c("pod", "torm"), sep=".table.");
  a[,which(colnames(a) == "torm")] <- NULL;
  a$rand_1 <- NA;
  df_pod <- a; rm(a)
  
  # rbind df_sim and df_pod
  df <- rbind.data.frame(df_sim, df_pod);
  df$pod <- as.factor(df$pod);
  rm(df_pod, df_sim)
  
  # exec time
  p_pod = df %>%
    subset(pod != "sim") %>%
    ggplot(aes(x=pod, y=exec_time_0, fill=pod))+
    geom_violin(alpha = 0.6)+
    geom_jitter(shape = 20)+
    scale_fill_brewer(palette = "Set1")+
    theme(legend.position="bottom",
          aspect.ratio=0.707)
  p_sim = df %>%
    subset(pod == "sim") %>%
    ggplot(aes(x=pod, y=exec_time_0, fill=pod))+
    geom_jitter(alpha=0.3, shape = 20)+
    geom_violin(alpha = 0.6)+
    scale_fill_brewer(palette = "Set1")+
    scale_y_log10()+
    theme(legend.position="bottom")
  plot_grid(p_pod, p_sim, nrow=1, rel_widths=c(4, 1), align=T)+
    labs(caption=paste0("mean execution time: ", mean(df$exec_time_0)))
  
}

# read in performance analysis
path = "/netscratch/dep_tsiantis/grp_laurent/struett/temp_tsabc/data_1/tsabc2/data/performance"
pattern = ".gzip"
files <- list.files(path=path, pattern=pattern, recursive=T)
full_paths <- paste(path, files, sep="/")

df <- do.call("rbind.data.frame", lapply(full_paths, function(x){
  t = pd$read_pickle(x)
  t$path = x
  return(t)
}));
a <- df %>%
  separate(col=path, into=c("torm", "prelim"), sep="pod_") %>%
  separate(col=prelim, into=c("pod", "prelim"), sep="/perf_sscomp_") %>%
  separate(col=prelim, into=c("sscomp", "prelim"), sep="_method_") %>%
  separate(col=prelim, into=c("regression", "torm"), sep=".table.");
a[,which(names(a)=="torm")] = NULL
df = a; rm(a)

# modification of columns
df$value[df$eval_mode == "factorx"] = 1-1/df$value[df$eval_mode == "factorx"]
df$value[df$eval_mode == "ci95"] = 1-df$value[df$eval_mode == "ci95"]
df$sscomp = as.factor(as.character(df$sscomp))
levels(df$sscomp) = c("1\nsfs__ld",
                      "2\ntl_true",
                      "3\nsfs__ld__tl_true",
                      "4\ntm_true",
                      "5\ntm_win",
                      "6\nsfs__ld__tm_true",
                      "7\nsfs__ld__tm_win",
                      "8\ntm_true__tl_true",
                      "9\nsfs__ld__tm_true__tl_true")
df$parameter = as.factor(as.character(df$parameter))
levels(df$parameter) = c("pop_size_recent",
                         "pop_size_ancient",
                         "time_change_pop_size",
                         "sigma_recent",
                         "sigma_ancient",
                         "time_change_sigma")

p <- df %>%
  subset(point_estim == "mode") %>%
  # subset(eval_mode %in% c("rb", "rmse", "factorx")) %>%
  subset(parameter %notin% c("time_change_pop_size")) %>%
  subset(pod %notin% c(10:19)) %>%
  subset(parameter %in% c("time_change_sigma", "sigma_recent")) %>%
  subset(sscomp %in% c("1\nsfs__ld", "2\ntl_true", "5\ntm_win", "7\nsfs__ld__tm_win")) %>%
  ggplot(aes(parameter, value, fill=eval_mode))+
  geom_col(position = "dodge")+
  geom_hline(aes(yintercept=1), size=0.5, col = "gray")+
  facet_grid(sscomp ~ point_estim + pod)+
  # scale_y_continuous(limits = c(-0.5, 2.5))+
  # scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values = wes_palette(length(unique(df$parameter)), name = "Darjeeling1", type = "continuous"))+
  # scale_y_log10()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p

