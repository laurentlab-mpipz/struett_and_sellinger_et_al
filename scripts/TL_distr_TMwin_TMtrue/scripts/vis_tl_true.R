
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
# library(reshape2)
library(yaml)

pd <- reticulate::import("pandas")
np <- reticulate::import("numpy")

args <- commandArgs(trailingOnly=TRUE)
# outfile <- "data/vis/scenario_ts_sigma_0.95_ts_4.pdf"
# infiles <-
#   c("data/df_tl_true_scenario_ts_sigma_0.95_ts_4_.pickle.gzip",
#              "data/df_tm_true_scenario_ts_sigma_0.95_ts_4_.pickle.gzip"
#   )

outfile = args[1]
infiles = args[2:length(args)]

# breaks = np$load("data/breaks/tmrca.bins.wc.0.ts.50000.0.95.3.6e-08.20000000.CS.4.npy") %>% 
#   as.numeric() %>% 
#   format(trim = T, digits = 2, scientific = T) %>% 
#   as.character()

message("outfile:\n\t", outfile)
message("infiles:\n\t", paste(infiles, collapse = "\n\t"))

tl_true <- pd$read_pickle(infiles[1], "gzip")
names(tl_true)[1:2] <-c("TMRCA", "LENGTH")

# stop(paste("",
#            paste(rep("_", 80), collapse = ""),
#            paste(rep("+", 80), collapse = ""),
#            paste(rep("=", 80), collapse = ""),
#            sep = "\n"))

get_time_in_generations_from_configfile_and_infile_name <-
  function(my_configfile, my_filename) {
    my_Ne <- read_yaml(my_configfile)["population_size"][[1]] %>%
      gsub("_", "", .) %>%
      as.numeric()
    my_ts <- strsplit(my_filename, "_")[[1]][11] %>% as.numeric() * my_Ne
    return(my_ts)
  }

time_in_generations <- get_time_in_generations_from_configfile_and_infile_name(
  "config.yaml", infiles[1])

fill_rects = cbind.data.frame(
  xmin=c(0, 0),
  xmax=c(Inf, Inf),
  ymin=c(time_in_generations, 0),
  ymax=c(Inf, time_in_generations),
  fill=factor(c(1, 2))
)

message("create main plot")
p <- ggplot()+
  geom_rect(data=fill_rects, mapping=aes(xmin=xmin,
                                         xmax=xmax,
                                         ymin=ymin,
                                         ymax=ymax,
                                         fill=fill,
                                         color=NA),
            alpha=0.8)+
  geom_density_2d(data=tl_true, aes(LENGTH, TMRCA), bins = 12, size = 1.2, col="black")+
  scale_colour_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Dark2")+
  scale_x_continuous(
    limits = c(5e1, 3e5),
    trans = "log10",
    breaks = sort(unique(as.vector(sapply(1:10, function(x) {10**c(2,3,4,5)*x})))),
    labels = c("100", rep("", 8), "1,000", rep("", 8), "10,000", rep("", 8), "100,000", rep("", 8), "1,000,000"))+
  scale_y_continuous(
    limits = c(5e3, 1e6),
    trans = "log10",
    breaks = sort(unique(as.vector(sapply(1:10, function(x) {10**c(4,5)*x})))),
    labels = c("1e+04", rep("", 8), "1e+05", rep("", 8), "1e+06"))+
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.line = element_line(size = 1.1),
    axis.text.y = element_text(angle = 90, hjust = 0.5)
  )

message("create x marginal")
xplot <- ggplot(tl_true, aes(LENGTH))+
  geom_density(fill = "black")+
  scale_x_continuous(
    limits = c(4e1, 4e5),
    trans = "log10",
    breaks = sort(unique(as.vector(sapply(1:10, function(x) {10**c(1,2,3,4,5)*x})))),
    labels = c("10", rep("", 8), "100", rep("", 8), "1,000", rep("", 8), "10,000", rep("", 8), "100,000", rep("", 8), "1,000,000"))+
  theme(
    legend.position = "none",
    axis.line = element_line(size = 1.1),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks = element_blank()
  )

message("create y marginal")
yplot <- ggplot(tl_true, aes(TMRCA))+
  geom_density(fill = "black")+
  scale_x_continuous(
    limits = c(1e4, 1e6),
    trans = "log10",
    breaks = sort(unique(as.vector(sapply(1:10, function(x) {10**c(1,2,3,4,5)*x})))),
    labels = c("10", rep("", 8), "100", rep("", 8), "1,000", rep("", 8), "10,000", rep("", 8), "100,000", rep("", 8), "1,000,000"))+
  theme(
    legend.position = "none",
    axis.line = element_line(size = 1.1),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks = element_blank()
  )+
  coord_flip()

message("create composition")
ptot <- plot_grid(
  xplot, NULL, NULL,
  NULL, NULL, NULL,
  p, NULL, yplot,
  ncol = 3,
  align = "hv",
  rel_widths = c(2, -0.35, 1),
  rel_heights = c(1, -0.35, 2),
  scale = 1
)+
  theme(aspect.ratio = 1)

message("print to file..")
pdf(outfile, width = 5, height = 5)
plot(ptot)
dev.off()
message("printed to file")










