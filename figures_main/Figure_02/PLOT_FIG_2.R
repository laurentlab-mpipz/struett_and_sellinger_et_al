library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(wesanderson)
library(yaml)
library(stringr)

outfile <- "MAN_FIG_2.pdf"

this_working_dir <- paste0(
  "/Users/struett/Dropbox/phd/manuscripts/transitioning_to_selfing/",
  "figure_finals/Figure_02"
)

setwd(this_working_dir)

# plot specs
mcol <- c(
  wes_palette("Chevalier1")[2],
  wes_palette("Chevalier1")[1]
)
  wes_palette("Darjeeling1", 5, type = "continuous")[c(1:3, 5)]
logbreak <- sapply(10 ** (-10:10), function(x) {
  x * (1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10 ** (-10:10))
    return(as.character(log10(x)))
  else
    return("")
})
mlimits <- c(6e2, 3e5)

###
# constant
source("constant.R")

setwd(this_working_dir)

message(getwd())


###
# expansion
source("expansion.R")

setwd(this_working_dir)

message(getwd())


###
# crash
source("crash.R")

setwd(this_working_dir)

message(getwd())

plot_list <- list(
  p_tsigma_constant,
  p_tsigma_expansion_tNe10k,
  p_tsigma_expansion_tNe40k,
  p_tsigma_crash_tNe10k,
  p_tsigma_crash_tNe40k
)

linesize <- 0.5


plot_list <- lapply(plot_list, function(x) {
  return(
    x + theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.border = element_rect(color = "black", size = linesize),
      axis.line = element_blank(),
      axis.ticks = element_line(size = linesize * 0.8),
      axis.ticks.length = unit(0.015, "in"),
      plot.margin = margin(8, 8, 8, 8, unit = "pt"),
      legend.position = "none"
    )
  )
})

### multipanel
p <- plot_grid(
  plot_list[[1]],
  NULL,
  plot_list[[2]],
  plot_list[[3]],
  plot_list[[4]],
  plot_list[[5]],
  ncol = 2,
  # labels = c("A", "", "B", "C", "D", "E"),
  align = TRUE
)

mwidth <- 4.4
mheight <- 6.8

pdf(outfile,
    useDingbats = FALSE,
    height = mheight,
    width = mwidth)
show(p)
dev.off()

png(
  paste0(outfile, ".png"),
  width = mwidth,
  height = mheight,
  res = 600,
  units = "in"
)
show(p)
dev.off()
