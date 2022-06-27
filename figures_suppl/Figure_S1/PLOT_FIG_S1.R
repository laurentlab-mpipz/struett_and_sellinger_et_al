library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(wesanderson)
# library(reshape2)
library(yaml)

logbreak <- sapply(10 ** (-10:10), function(x) {
  x * (1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10 ** (-10:10))
    return(as.character(log10(x)))
  else
    return("")
})


outfile <- "3model_compareison_tl_true.pdf"

# plot settings (TL_true)
mwidth <- 3.6*3
mheight <- mwidth/3
mcol <-
  c(wes_palette("Moonrise2")[2], wes_palette("Moonrise3")[3]) %>% rev()


reticulate::use_condaenv("/Users/abgushtdizi/Library/r-miniconda/envs/r-reticulate")
reticulate::use_python("/Users/abgushtdizi/Library/r-miniconda/envs/r-reticulate/bin/python")

pd <- reticulate::import("pandas")
np <- reticulate::import("numpy")

infiles = c("3model_subset_4Ne95sigma.gzip")


tl_true <- pd$read_pickle(infiles[1], "gzip")
names(tl_true)[1:2] <- c("TMRCA", "LENGTH")

tl_true$model[tl_true$model == "EX"] <- "Forward-Explicit"
tl_true$model[tl_true$model == "NA"] <- "Forward-Rescaled"
tl_true$model[tl_true$model == "CS"] <- "Coalescent-Rescaled"
tl_true$model <- factor(tl_true$model, levels = c(
  "Forward-Explicit",
  "Forward-Rescaled",
  "Coalescent-Rescaled"
))


get_time_in_generations_from_configfile_and_infile_name <-
  function(my_configfile, my_filename) {
    my_Ne <- read_yaml(my_configfile)["population_size"][[1]] %>%
      gsub("_", "", .) %>%
      as.numeric()
    
    my_ts <-
      strsplit(my_filename, "_")[[1]][12] %>% as.numeric() * my_Ne
    return(my_ts)
  }


time_in_generations <- 200000
  # get_time_in_generations_from_configfile_and_infile_name("msprim_ts/config.yaml", infiles[1])
fill_rects = cbind.data.frame(
  xmin = c(0, 0),
  xmax = c(Inf, Inf),
  ymin = c(time_in_generations, 0),
  ymax = c(Inf, time_in_generations),
  fill = factor(c(1, 2))
)

message("create main plot")
p <- ggplot() +
  geom_rect(
    data = fill_rects,
    mapping = aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill,
      color = NA
    ),
    alpha = 0.8
  ) +
  geom_density_2d(
    data = tl_true,
    aes(LENGTH, TMRCA),
    bins = 12,
    size = 1.2,
    col = "black"
  ) +
  facet_grid(. ~ model) +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Dark2") +
  # scale_fill_manual(values = mcol) +
  scale_x_continuous(
    limits = c(5e1, 3e5),
    trans = "log10",
    breaks = logbreak,
    labels = loglabel
  ) +
  scale_y_continuous(
    limits = c(5e3, 1e6),
    trans = "log10",
    breaks = sort(unique(as.vector(
      sapply(1:10, function(x) {
        10 ** c(4, 5) * x
      })
    ))),
    labels = c("1e+04", rep("", 8), "1e+05", rep("", 8), "1e+06")
  ) +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.line = element_line(size = 1.1),
    axis.text.y = element_text(angle = 90, hjust = 0.5)
  )



# message("print to file..")
# pdf(outfile,
#     width = mwidth,
#     height = mheight,
#     useDingbats = F)
# plot(p)
# dev.off()

message("print png")
png(
  paste0(outfile, ".png"),
  width = mwidth,
  height = mheight,
  units = "in",
  res = 600
)
plot(p)
dev.off()
message("printed to file")
