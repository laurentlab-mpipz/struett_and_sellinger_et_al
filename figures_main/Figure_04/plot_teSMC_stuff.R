library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(wesanderson)
library(pbapply)

setwd("teSMC_bgs/")

# return mio if INF
maxInf <- function(x) {
  if (is.infinite(x)) {
    value <- 1e6
  } else {
    value <- x
  }
  return(value)
}


fnames_bw <-
  list.files(path = "bgs_good_window",
             pattern = "mat_*.csv",
             full.names = T)

csvFiles2df <- function(flist) {
  df <- do.call(rbind.data.frame, lapply(flist, function(x) {
    t <- read.csv(x)
    names(t) <- c("rep", 1:(ncol(t) - 1))
    t$fname <- strsplit(x, "/+", perl = TRUE)[[1]][2]
    return(t)
  }))
  
  return(df)
}

# read files with added filenames
df_bw <- csvFiles2df(fnames_bw)

#

###
# get estim vs true plots

###
# constant
thiswd <- getwd()
source("teSMC.bgs.R")
setwd(thiswd)

message(getwd())

# plot value setup
logbreak <- sapply(10 ** (-10:10), function(x) {
  x * (1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10 ** (-10:10))
    return(as.character(log10(x)))
  else
    return("")
})
msbreaks <- seq(0, 1, 0.1)
mslabels <- sapply(msbreaks, function(x) {
  if (x %in% c(0, 1))
    return(as.character(x))
  else
    return("")
})
mcol <- c(
  "gray",
  wes_palette("Chevalier1")[1]
)
xlim <- c(5e3, 2e5)
ylim <- xlim


df <- read.csv("bgs_good_window/point_estims.csv", row.names = 1)
df$dfe <- factor(df$dfe, levels = c(0, 1))

p_dfe0 <- df %>%
  subset(dfe == 0) %>%
  ggplot(aes(tsigma, tsigma_estim, col = masked, shape = dfe)) +
  geom_abline() +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.06,
                                    dodge.width = 0.0),
    shape = 16,
    size = 1.7,
    alpha = 0.8
  ) +
  # facet_grid(dfe ~ .) +
  scale_x_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel,
                     limits = xlim) +
  scale_y_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel,
                     limits = ylim) +
  scale_color_manual(values = mcol) +
  labs(
    x = expression("t"[sigma] ~ " [log"[10] ~ "]"),
    y = expression(hat("t"[sigma]) ~ " [log"[10] ~ "]"),
    col = "Masked"
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", size = 0.9),
    axis.line = element_blank(),
    axis.ticks.length = unit(-3, units = "pt")
  )

p_dfe1 <- df %>%
  subset(dfe == 1) %>%
  ggplot(aes(tsigma, tsigma_estim, col = masked, shape = dfe)) +
  geom_abline() +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.06,
                                    dodge.width = 0.0),
    shape = 16,
    size = 1.7,
    alpha = 0.8
  ) +
  # facet_grid(dfe ~ .) +
  scale_x_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel,
                     limits = xlim) +
  scale_y_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel,
                     limits = ylim) +
  scale_color_manual(values = mcol) +
  labs(
    x = expression("t"[sigma] ~ " [log"[10] ~ "]"),
    y = expression(hat("t"[sigma]) ~ " [log"[10] ~ "]"),
    col = "Masked"
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", size = 0.9),
    axis.line = element_blank(),
    axis.ticks.length = unit(-3, units = "pt")
  )


# working dir
setwd("../")
getwd()

p <- p_dfe0 

outfile = "teSMC.dfe0.bgs.pdf"
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


p <- p_dfe1

outfile = "teSMC.dfe1.bgs.pdf"
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

