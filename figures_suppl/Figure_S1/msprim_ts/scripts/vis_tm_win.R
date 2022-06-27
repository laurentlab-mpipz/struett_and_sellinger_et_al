
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(reshape2)

pd <- reticulate::import("pandas")
np <- reticulate::import("numpy")

args <- commandArgs(trailingOnly=TRUE)
# setwd("..")
# outfile <- "data/vis/scenario_ts_sigma_0.95_ts_4.pdf"
# infiles <-
#   c("data/df_tl_true_scenario_ts_sigma_0.95_ts_4_.pickle.gzip",
#     "data/df_tm_true_scenario_ts_sigma_0.95_ts_4_.pickle.gzip",
#     "data/df_tm_win/df_tm_true_scenario_tNe_sigma_0.95_ts_0_.pickle.gzip"
#   )

outfile = args[1]
infiles = args[2:length(args)]

breaks = np$load("data/breaks/tmrca.bins.tmwin.0.tNe.50000.0.95.3.4e-09.20000000.CS.0.npy") %>% 
  as.numeric() %>% 
  format(trim = T, digits = 2, scientific = T) %>% 
  as.character()

message("outfile:\n\t", outfile)
message("infiles:\n\t", paste(infiles, collapse = "\n\t"))

# tl_true <- pd$read_pickle(infiles[1], "gzip")
# tm_true <- pd$read_pickle(infiles[1], "gzip")
tm_win <- pd$read_pickle(infiles[1], "gzip")

a <- apply(tm_win[,1:121], 2, mean)
m <- matrix(a, nrow = 11)

m <- m %>% melt()
m$Var1 <- factor(rep(breaks, 11), levels = breaks)
m$Var2 <- factor(rep(breaks, each=11), levels = breaks)
names(m)[3] <- "tp"
m$tp <- log10(m$tp+1)

p <- m %>%
  ggplot(aes(Var1, Var2, fill=tp))+
  geom_tile()+
  xlab("n-th segment")+
  ylab("(n+1)-th segment")+
  scale_fill_gradient(low="#ffffc7",
                      high="#7d0125",
                      limits=c(0, 0.3))+
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        # axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle=90, vjust = 0.5),
        legend.position = "right",
        legend.key.height = unit(0.5, "in"),
        legend.key.width = unit(0.05, "in"),
        panel.border = element_rect(colour = "black", fill=NA))

pdf(outfile)
show(p)
dev.off()









