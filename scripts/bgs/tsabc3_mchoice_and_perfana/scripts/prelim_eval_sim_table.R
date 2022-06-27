
library(reticulate)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
message("loaded libraries")

pd <- import("pandas")
message("loaded pandas")

df <- pd$read_pickle("data/table/params_and_sumstats.table.gzip")
message("loaded pickle")

# p_exec_violin <- ggplot(df, aes(x = 1, exec_time_0))+
#   geom_jitter(col = "darkblue", alpha = 0.3, shape = 20)+
#   geom_violin(fill = "hotpink2", alpha = 0.6)+
#   geom_point(aes(x = 1, y = mean(df$exec_time_0)), shape = 18, size = 5, col = "black", fill = "black")
# message("created glob")
# 
# print(p_exec_violin)
# message("done printing figure")

for (i in 1:ncol(df)) {
  cat(i);
  cat("\t");
  cat(names(df)[i]);
  cat("\t");
  cat(mean(df[,i]));
  cat("\n")
}