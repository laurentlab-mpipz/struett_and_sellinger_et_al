
library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(wesanderson)
library(scales)
library(grid)

# for (i in 1:length(wes_palettes)) {
#   show(wes_palette(names(wes_palettes)[i]))
# }

logbreak <- sapply(10**(-10:10), function(x) {
  x*(1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10**(-10:10))
    return(as.character(log10(x)))
  else
    return("")
})

mcol = wes_palette("Darjeeling1")[1:2]
mcol = c(mcol[1], mcol[2])
show_col(mcol)
colfunc<-colorRampPalette(mcol)
show_col(colfunc(8))
mcol = colfunc(8)
mcol = c(mcol, rev(mcol)[2:length(mcol)])
show_col(mcol)

# setwd("~/Desktop/some_data_analysis/tsabc_tables_txt/")
infile = "quants.all_par.tsabc.bgs.dfe_0.csv"

pod_to_time = c(20000, 2000, 4000, 8000, 12000, 16000, 30000, 40000, 60000,
                80000, 100000, 200000, 1000, 3000, 5000, 6000, 7000, 9000,
                10000, 50000, 70000, 90000)

df <- fread(infile)

# Ne
df$param_0 %>% min()
df$param_0 %>% max()

# sigma_recent
df$param_1 %>% min()
df$param_1 %>% max()

# sigma_ancient
df$param_2 %>% min()
df$param_2 %>% max()

# tsigma
df$param_3 %>% min()
df$param_3 %>% max()

df$quant %>% unique() %>% length()

a = df %>%
  subset(quant == 0.5)

# Ne
hist(a$param_0, breaks="FD", col="gray")
hist(a$param_1, breaks="FD", col="gray")
hist(a$param_2, breaks="FD", col="gray")
hist(a$param_3, breaks="FD", col="gray")

mdf <- df %>%
  separate(col=file, sep="/", into=c("data", "point_estims", "pod_n", "pls", "end"))
mdf$data = NULL
mdf$point_estims = NULL

mdf <- mdf %>%
  separate(col=pod_n, sep="_", into=c("torm", "pod_n"))
mdf$torm = NULL
mdf$pod_n = as.numeric(mdf$pod_n)

# translate pod_n to tsigma
mdf$tsigma <- sapply(mdf$pod_n, function(x) {
  pod_index = x + 1
  return(pod_to_time[pod_index])
})

mdf %>%
  subset(sscomp == 5) %>%
  subset(plsComp == 5) %>%
  subset(quant %in% c(
    0.25, 0.75,
    0.5
  )) %>%
  # subset(pod==1) %>%
  ggplot(
    aes(
      x=tsigma,
      y=param_3,
      col=as.factor(quant),
      fill=as.factor(pod)
    ))+
  geom_line(show.legend = F)+
  # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
  geom_line(aes(tsigma, tsigma), col="black")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(
    trans = "log10",
    # limits = c(-2, 2)
  )+
  scale_color_manual(values = mcol[c(1, 8, 15)])+
  theme(
    aspect.ratio = 1
  )

r <- mdf %>%
  group_by(pod, quant, plsComp, sscomp) %>%
  summarise(param_3 = mean(param_3),
            tsigma = tsigma)

r %>%
  subset(sscomp == 1) %>%
  subset(plsComp==5) %>%
  ggplot(
    aes(
      x=tsigma,
      y=param_3,
      col=as.factor(quant)
    ))+
  geom_line(show.legend = F)+
  # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
  geom_line(aes(tsigma, tsigma), col="black")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(
    trans = "log10",
    # limits = c(-2, 2)
  )+
  scale_color_manual(values = mcol)+
  # scale_color_manual(values = mcol[c(1, 8, 15)])+
  theme(
    aspect.ratio = 1
  )

r %>%
  subset(tsigma == 20000) %>%
  subset(plsComp == 5) %>%
  subset(quant == 0.5) %>%
  subset(sscomp == 2)

# average per 100 pods
mdf <- df %>%
  separate(col=file, sep="/", into=c("data", "point_estims", "pod_n", "pls", "end"))
mdf$data = NULL
mdf$point_estims = NULL

mdf <- mdf %>%
  separate(col=pod_n, sep="_", into=c("torm", "pod_n"))
mdf$torm = NULL
mdf$pod_n = as.numeric(mdf$pod_n)

# translate pod_n to tsigma
mdf$tsigma <- sapply(mdf$pod_n, function(x) {
  pod_index = x + 1
  return(pod_to_time[pod_index])
})

mdf %>%  # to find the grouping factors
  subset(plsComp == 5) %>%
  subset(sscomp == 1) %>%
  subset(method == "loclinear") %>%
  subset(tsigma == 20000) %>%
  subset(pod_n == 0) %>%
  subset(quant == 0.005) %>%
  nrow()

data <- mdf %>%
  group_by(plsComp, sscomp, method, pod_n, quant, tsigma) %>%
  summarise(
    avquant_param_0 = mean(param_0),
    avquant_param_1 = mean(param_1),
    avquant_param_2 = mean(param_2),
    avquant_param_3 = mean(param_3),
    n = n(), 
  )

stopifnot(all(data$n == 100))

# rename sscomp, plsComp, params
data$sscomp <- as.factor(data$sscomp)
levels(data$sscomp) <- c("sfs/ld", "tl_true", "tm_true", "tm_win", "sfs/ld/tm_win")
data$plsComp <- as.factor(data$plsComp)
levels(data$plsComp) <- c("5-PLS", "10-PLS", "20-PLS", "30-PLS")
colnames(data)[colnames(data) == 'avquant_param_0'] <- 'pop_size'
colnames(data)[colnames(data) == 'avquant_param_1'] <- 'recent_selfing'
colnames(data)[colnames(data) == 'avquant_param_2'] <- 'ancestral_selfing'
colnames(data)[colnames(data) == 'avquant_param_3'] <- 'transitioning_time'


# all parameters to be plotted
{
  p_tsigma <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = transitioning_time,
               col = as.factor(quant))) +
    geom_line(show.legend = F) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    geom_line(aes(tsigma, tsigma), col = "black") +
    # geom_text() %>%
    facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      trans = "log10",
      limits = c(1e3, 5e5),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(aspect.ratio = 1)
  
  p_sigma_recent <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = recent_selfing,
               col = as.factor(quant))) +
    geom_line(show.legend = F) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    geom_line(aes(tsigma, 0.99), col="black")+
    # geom_text() %>%
    facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      # trans = "log10", breaks = logbreak, labels = loglabel,
      limits = c(0.5, 1),
      breaks = seq(0, 1, 0.1),
      labels = sapply(seq(0, 1, 0.1), function(x) {
        if (x %in% c(0, 1))
          return(as.character(x))
        else if (x == 0.5)
          return(".5")
        else
          return("")
      })
    ) +
    scale_color_manual(values = mcol) +
    theme(aspect.ratio = 1)
  
  p_sigma_ancient <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = ancestral_selfing,
               col = as.factor(quant))) +
    geom_line(show.legend = F) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    geom_line(aes(tsigma, 0.1), col="black")+
    # geom_text() %>%
    facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      # trans = "log10", breaks = logbreak, labels = loglabel,
      limits = c(0, .2),
      breaks = seq(0, 1, 0.01),
      labels = sapply(seq(0, 1, 0.01), function(x) {
        if (x %in% c(0, 1))
          return(as.character(x))
        else if (x == 0.1)
          return(".1")
        else if (x == 0.2)
          return(".2")
        else
          return("")
      })
    ) +
    scale_color_manual(values = mcol) +
    theme(aspect.ratio = 1)
  
  p_pop_size <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = pop_size,
               col = as.factor(quant))) +
    geom_line(show.legend = F) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    geom_line(aes(tsigma, 40000), col="black")+
    # geom_text() %>%
    facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      trans = "log10",
      limits = c(1e4, 2e5),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(aspect.ratio = 1)
  
  # legend plot
  data$quantile = as.factor(data$quant)
  p_legend <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = transitioning_time,
               col = quantile)) +
    geom_line(show.legend = T) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    geom_line(aes(tsigma, tsigma), col = "black") +
    # geom_text() %>%
    facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      trans = "log10",
      limits = c(1e3, 5e5),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(aspect.ratio = 1)
}

p_tot <- plot_grid(
  p_tsigma,
  p_pop_size,
  p_sigma_recent,
  p_sigma_ancient,
  align = T,
  labels = "auto",
  ncol = 2
)

mcol = wes_palette("Darjeeling1")[1:2]
mcol = c(mcol[1], mcol[2])
show_col(mcol)
colfunc<-colorRampPalette(mcol)
show_col(colfunc(8))
mcol = colfunc(8)
mcol = c(mcol, rev(mcol)[2:length(mcol)])
show_col(mcol)


msize = 2

logbreak <- sapply(10**(-10:10), function(x) {
  x*(1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10**(-10:10))
    return(as.character(log10(x)))
  else
    return("")
})

slogbreak <- logbreak
sloglabel <- loglabel

# the plots for the three stat compositions for t_sigma
if (F)
{
  logbreak <- slogbreak
  loglabel <- sloglabel
  
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "sfs/ld") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = transitioning_time,
               col = as.factor(quant))) +
    geom_line(aes(tsigma, tsigma), col = "black", size=msize) +
    geom_line(show.legend = F, size=msize) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    # geom_text() %>%
    # facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      trans = "log10",
      limits = c(1e3, 5e5),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(
      legend.position = "none",
      aspect.ratio = 1,
      panel.border = element_rect(colour = "black", size = msize),
      # text = element_text(size = 12),
      #   strip.background = element_blank(),
      # axis.title.y = element_blank(),
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8),
      # axis.ticks.y = element_blank(),
      # axis.ticks.x = element_line(size = 1),
      # axis.ticks.length.x = unit(-1, "lines"),
      # panel.spacing = unit(-0.9, "lines")
    )
  
  png("sfs_ld_rmu1.tsigma.param_estim.bgs.png", width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()


  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "tm_win") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = transitioning_time,
               col = as.factor(quant))) +
    geom_line(show.legend = F, size=msize) +
    geom_line(aes(tsigma, tsigma), col = "black", size=msize) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    # geom_text() %>%
    # facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      trans = "log10",
      limits = c(1e3, 5e5),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(
      legend.position = "none",
      aspect.ratio = 1,
      panel.border = element_rect(colour = "black", size = msize),
      # text = element_text(size = 12),
      #   strip.background = element_blank(),
      # axis.title.y = element_blank(),
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8),
      # axis.ticks.y = element_blank(),
      # axis.ticks.x = element_line(size = 1),
      # axis.ticks.length.x = unit(-1, "lines"),
      # panel.spacing = unit(-0.9, "lines")
    )
  png("tm_win_rmu1.tsigma.param_estim.bgs.png", width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()
  
  
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "sfs/ld/tm_win") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = transitioning_time,
               col = as.factor(quant))) +
    geom_line(show.legend = F, size=msize) +
    geom_line(aes(tsigma, tsigma), col = "black", size=msize) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    # geom_text() %>%
    # facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      trans = "log10",
      limits = c(1e3, 5e5),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(
      legend.position = "none",
      aspect.ratio = 1,
      panel.border = element_rect(colour = "black", size = msize),
      # text = element_text(size = 12),
      #   strip.background = element_blank(),
      # axis.title.y = element_blank(),
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8),
      # axis.ticks.y = element_blank(),
      # axis.ticks.x = element_line(size = 1),
      axis.ticks.length = unit(0.5, "lines"),
      # panel.spacing = unit(-0.9, "lines")
    )
  png("sfs_ld_tm_win_rmu1.tsigma.param_estim.bgs.png", width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()
  
}


# the plots for the three stat compositions for recent selfing
{
  which_param <- "recent_selfing"
  
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "sfs/ld") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = recent_selfing,
               col = as.factor(quant))) +
    geom_line(show.legend = F, size=msize) +
    geom_line(aes(tsigma, 0.99), col="black", size=msize)+
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      limits = c(0.5, 1.0),
      breaks = seq(0.5, 1.0, 0.1),
      labels = seq(0.5, 1.0, 0.1)
    ) +
    scale_color_manual(values = mcol) +
    theme(
      panel.border = element_rect(colour = "black", size = msize),
      legend.position = "none",
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8)
    )
  
  png(paste0("sfs_ld_rmu1.", which_param, ".param_estim.bgs.png"), width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()
  
  
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "tm_win") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = recent_selfing,
               col = as.factor(quant))) +
    geom_line(show.legend = F, size=msize) +
    geom_line(aes(tsigma, 0.99), col="black", size=msize)+
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      limits = c(0.5, 1.0),
      breaks = seq(0.5, 1.0, 0.1),
      labels = seq(0.5, 1.0, 0.1)
    ) +
    scale_color_manual(values = mcol) +
    theme(
      panel.border = element_rect(colour = "black", size = msize),
      legend.position = "none",
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8)
    )
  
  
  png(paste0("tm_win_rmu1.", which_param, ".param_estim.bgs.png"), width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()
  
  
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "sfs/ld/tm_win") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = recent_selfing,
               col = as.factor(quant))) +
    geom_line(show.legend = F, size=msize) +
    geom_line(aes(tsigma, 0.99), col="black", size=msize)+
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      limits = c(0.5, 1.0),
      breaks = seq(0.5, 1.0, 0.1),
      labels = seq(0.5, 1.0, 0.1)
    ) +
    scale_color_manual(values = mcol) +
    theme(
      panel.border = element_rect(colour = "black", size = msize),
      legend.position = "none",
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8)
    )
  
  
  png(paste0("sfs_ld_tm_win_rmu1.", which_param, ".param_estim.bgs.png"), width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()
  
}


# the plots for the three stat compositions for recent selfing
{
  which_param <- "ancestral_selfing"
  logbreak <- seq(0, 0.2, 0.05)
  loglabel <- logbreak
    
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "sfs/ld") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = ancestral_selfing,
               col = as.factor(quant))) +
    geom_line(show.legend = F, size=msize) +
    geom_line(aes(tsigma, 0.1), col="black", size=msize)+
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = slogbreak,
                       labels = sloglabel) +
    scale_y_continuous(
      limits = c(0, 0.2),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(
      panel.border = element_rect(colour = "black", size = msize),
      legend.position = "none",
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8)
    )
  
  png(paste0("sfs_ld_rmu1.", which_param, ".param_estim.bgs.png"), width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()
  
  
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "tm_win") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = ancestral_selfing,
               col = as.factor(quant))) +
    geom_line(show.legend = F, size=msize) +
    geom_line(aes(tsigma, 0.1), col="black", size=msize)+
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = slogbreak,
                       labels = sloglabel) +
    scale_y_continuous(
      limits = c(0, 0.2),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(
      panel.border = element_rect(colour = "black", size = msize),
      legend.position = "none",
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8)
    )
  
  
  png(paste0("tm_win_rmu1.", which_param, ".param_estim.bgs.png"), width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()
  
  
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "sfs/ld/tm_win") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = ancestral_selfing,
               col = as.factor(quant))) +
    geom_line(show.legend = F, size=msize) +
    geom_line(aes(tsigma, 0.1), col="black", size=msize)+
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = slogbreak,
                       labels = sloglabel) +
    scale_y_continuous(
      limits = c(0, 0.2),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(
      panel.border = element_rect(colour = "black", size = msize),
      legend.position = "none",
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8)
    )
  
  
  png(paste0("sfs_ld_tm_win_rmu1.", which_param, ".param_estim.bgs.png"), width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()
  
}

# the plots for the three stat compositions for population size
{
  which_param <- "pop_size"
  logbreak <- slogbreak
  loglabel <- sloglabel
  true_val <- 4e4
  
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "sfs/ld") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = pop_size,
               col = as.factor(quant))) +
    geom_line(show.legend = F, size=msize) +
    geom_line(aes(tsigma, true_val), col="black", size=msize)+
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = slogbreak,
                       labels = sloglabel) +
    scale_y_continuous(
      trans = "log10",
      limits = c(1e4, 2e5),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(
      panel.border = element_rect(colour = "black", size = msize),
      legend.position = "none",
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8)
    )
  
  png(paste0("sfs_ld_rmu1.", which_param, ".param_estim.bgs.png"), width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()
  
  
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "tm_win") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = pop_size,
               col = as.factor(quant))) +
    geom_line(show.legend = F, size=msize) +
    geom_line(aes(tsigma, true_val), col="black", size=msize)+
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = slogbreak,
                       labels = sloglabel) +
    scale_y_continuous(
      trans = "log10",
      limits = c(1e4, 2e5),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(
      panel.border = element_rect(colour = "black", size = msize),
      legend.position = "none",
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8)
    )
  
  
  png(paste0("tm_win_rmu1.", which_param, ".param_estim.bgs.png"), width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()
  
  
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "sfs/ld/tm_win") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = pop_size,
               col = as.factor(quant))) +
    geom_line(show.legend = F, size=msize) +
    geom_line(aes(tsigma, true_val), col="black", size=msize)+
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = slogbreak,
                       labels = sloglabel) +
    scale_y_continuous(
      trans = "log10",
      limits = c(1e4, 2e5),
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    theme(
      panel.border = element_rect(colour = "black", size = msize),
      legend.position = "none",
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(size = msize*0.8)
    )
  
  
  png(paste0("sfs_ld_tm_win_rmu1.", which_param, ".param_estim.bgs.png"), width = 7, height = 7, units = "in", res = 600)
  show(p1)
  dev.off()
  
}
