


if (FALSE) {
  source("./prelim_plot_mchoice.R")
  source("./prelim_plot_mchoice_5.R")
  source("./prelim_plot_mchoice_bgs.R")
}

mcol <- wes_palette("Cavalcanti1")[c(5, 2)]
mcol = c(mcol[1], "gray90", mcol[2])
show_col(mcol)
colfunc<-colorRampPalette(mcol)
show_col(colfunc(9))
mcol = colfunc(9)[c(4:9)]
mcol = colfunc(9)[c(1, 5:9)]
show_col(mcol)

data <- readRDS("quants.prepared.all_par.tsabc.ronmu1.RDS")

mcol = wes_palette("Darjeeling1")[1:2]
mcol = c(mcol[1], mcol[2])
show_col(mcol)
colfunc<-colorRampPalette(mcol)
show_col(colfunc(8))
mcol = colfunc(8)
mcol = c(mcol, rev(mcol)[2:length(mcol)])
show_col(mcol)

stopifnot(all(data$n == 100))

msize = 2

# pop_size
ylimits <- c(1e4, 2e5)
{
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "sfs/ld") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = pop_size,
               col = as.factor(quant))) +
    geom_line(aes(tsigma, 40000), col = "black", size=msize) +
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
      limits = ylimits,
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
  png("sfs_ld_rmu1.tsigma.pop_size.param_estim.png", width = 7, height = 7, units = "in", res = 100)
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
    geom_line(aes(tsigma, 40000), col = "black", size=msize) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    # geom_text() %>%
    # facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      trans = "log10",
      limits = ylimits,
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
  png("tm_win_rmu1.tsigma.pop_size.param_estim.png", width = 7, height = 7, units = "in", res = 100)
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
    geom_line(aes(tsigma, 40000), col = "black", size=msize) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    # geom_text() %>%
    # facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      trans = "log10",
      limits = ylimits,
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
  png("sfs_ld_tm_win_rmu1.tsigma.pop_size.param_estim.png", width = 7, height = 7, units = "in", res = 100)
  show(p1)
  dev.off()
  
}


# recent_selfing
ylimits <- c(0.5, 1); ybreaks = seq(0.5, 1, 0.1)
{
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "sfs/ld") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = recent_selfing,
               col = as.factor(quant))) +
    geom_line(aes(tsigma, 0.99), col = "black", size=msize) +
    geom_line(show.legend = F, size=msize) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    # geom_text() %>%
    # facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      # trans = "log10",
      limits = ylimits,
      breaks = ybreaks
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
  png("sfs_ld_rmu1.tsigma.recent_selfing.param_estim.png", width = 7, height = 7, units = "in", res = 100)
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
    geom_line(aes(tsigma, 0.99), col = "black", size=msize) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    # geom_text() %>%
    # facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      # trans = "log10",
      limits = ylimits,
      breaks = ybreaks
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
  png("tm_win_rmu1.tsigma.recent_selfing.param_estim.png", width = 7, height = 7, units = "in", res = 100)
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
    geom_line(aes(tsigma, 0.99), col = "black", size=msize) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    # geom_text() %>%
    # facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      # trans = "log10",
      limits = ylimits,
      breaks = ybreaks
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
  png("sfs_ld_tm_win_rmu1.tsigma.recent_selfing.param_estim.png", width = 7, height = 7, units = "in", res = 100)
  show(p1)
  dev.off()
  
}


# ancestral_selfing
ylimits <- c(0, 0.2); ybreaks = seq(0, 0.2, 0.02)
{
  p1 <- data %>%
    subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win")) %>%
    subset(sscomp == "sfs/ld") %>%
    subset(plsComp == "20-PLS") %>%
    ggplot(aes(x = tsigma,
               y = ancestral_selfing,
               col = as.factor(quant))) +
    geom_line(aes(tsigma, 0.1), col = "black", size=msize) +
    geom_line(show.legend = F, size=msize) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    # geom_text() %>%
    # facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      # trans = "log10",
      limits = ylimits,
      breaks = ybreaks
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
  png("sfs_ld_rmu1.tsigma.ancestral_selfing.param_estim.png", width = 7, height = 7, units = "in", res = 100)
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
    geom_line(aes(tsigma, 0.1), col = "black", size=msize) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    # geom_text() %>%
    # facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      # trans = "log10",
      limits = ylimits,
      breaks = ybreaks
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
  png("tm_win_rmu1.tsigma.ancestral_selfing.param_estim.png", width = 7, height = 7, units = "in", res = 100)
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
    geom_line(aes(tsigma, 0.1), col = "black", size=msize) +
    # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
    # geom_text() %>%
    # facet_grid(plsComp ~ sscomp) +
    scale_x_continuous(trans = "log10",
                       # limits = c(1e3, 5e5),
                       breaks = logbreak,
                       labels = loglabel) +
    scale_y_continuous(
      # trans = "log10",
      limits = ylimits,
      breaks = ybreaks
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
  png("sfs_ld_tm_win_rmu1.tsigma.ancestral_selfing.param_estim.png", width = 7, height = 7, units = "in", res = 100)
  show(p1)
  dev.off()
  
}

