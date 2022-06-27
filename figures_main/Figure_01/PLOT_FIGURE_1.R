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

# plot settings (TL_true)
time_point <- 4
mwidth <- 4.4
mheight <- mwidth
mcol <-
  c(wes_palette("Moonrise2")[2], wes_palette("Moonrise3")[3]) %>% rev()

reticulate::use_condaenv("/Users/abgushtdizi/Library/r-miniconda/envs/r-reticulate")
reticulate::use_python("/Users/abgushtdizi/Library/r-miniconda/envs/r-reticulate/bin/python")

pd <- reticulate::import("pandas")
np <- reticulate::import("numpy")

args <- commandArgs(trailingOnly = TRUE)
outfile = args[1]
infiles = args[2:length(args)]

if (T) {
  outfile <- paste0("scenario_tNe_sigma_0.95_ts_", time_point, ".pdf")
  infiles <-
    c(
      paste0(
        "msprim_ts/data/df_tl_true/df_tl_true_scenario_tNe_sigma_0.95_ts_",
        time_point,
        "_.pickle.gzip"
      )
    )
  
}
message("outfile:\n\t", outfile)
message("infiles:\n\t", paste(infiles, collapse = "\n\t"))
{
  tl_true <- pd$read_pickle(infiles[1], "gzip")
  names(tl_true)[1:2] <- c("TMRCA", "LENGTH")
  
  
  get_time_in_generations_from_configfile_and_infile_name <-
    function(my_configfile, my_filename) {
      my_Ne <- read_yaml(my_configfile)["population_size"][[1]] %>%
        gsub("_", "", .) %>%
        as.numeric()
      
      my_ts <-
        strsplit(my_filename, "_")[[1]][12] %>% as.numeric() * my_Ne
      return(my_ts)
    }
  
  
  time_in_generations <-
    get_time_in_generations_from_configfile_and_infile_name("msprim_ts/config.yaml", infiles[1])
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
  
  
  message("create x marginal")
  xplot <- ggplot(tl_true, aes(LENGTH)) +
    geom_density(fill = "black") +
    scale_x_continuous(
      limits = c(4e1, 4e5),
      trans = "log10",
      breaks = sort(unique(as.vector(
        sapply(1:10, function(x) {
          10 ** c(1, 2, 3, 4, 5) * x
        })
      ))),
      labels = c(
        "10",
        rep("", 8),
        "100",
        rep("", 8),
        "1,000",
        rep("", 8),
        "10,000",
        rep("", 8),
        "100,000",
        rep("", 8),
        "1,000,000"
      )
    ) +
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
  yplot <- ggplot(tl_true, aes(TMRCA)) +
    geom_density(fill = "black") +
    scale_x_continuous(
      limits = c(1e4, 1e6),
      trans = "log10",
      breaks = sort(unique(as.vector(
        sapply(1:10, function(x) {
          10 ** c(1, 2, 3, 4, 5) * x
        })
      ))),
      labels = c(
        "10",
        rep("", 8),
        "100",
        rep("", 8),
        "1,000",
        rep("", 8),
        "10,000",
        rep("", 8),
        "100,000",
        rep("", 8),
        "1,000,000"
      )
    ) +
    theme(
      legend.position = "none",
      axis.line = element_line(size = 1.1),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks = element_blank()
    ) +
    coord_flip()
  
  message("create composition")
  ptot <- plot_grid(
    xplot,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    p,
    NULL,
    yplot,
    ncol = 3,
    align = "hv",
    rel_widths = c(2,-0.35, 1),
    rel_heights = c(1,-0.35, 2),
    scale = 1
  ) +
    theme(aspect.ratio = 1)
  
  message("print to file..")
  pdf(outfile, width = mwidth, height = mheight, useDingbats = F)
  plot(ptot)
  dev.off()
  
  png(paste0(outfile, ".png"), width = 4.4, height = 4.4, units = "in", res = 300)
  plot(ptot)
  dev.off()
  message("printed to file")
}

if (T) {
  outfile <- paste0("scenario_ts_sigma_0.95_ts_", time_point, ".pdf")
  infiles <-
    c(
      paste0(
        "msprim_ts/data/df_tl_true/df_tl_true_scenario_ts_sigma_0.95_ts_",
        time_point,
        "_.pickle.gzip"
      )
    )
  
}
message("outfile:\n\t", outfile)
message("infiles:\n\t", paste(infiles, collapse = "\n\t"))
{
  tl_true <- pd$read_pickle(infiles[1], "gzip")
  names(tl_true)[1:2] <- c("TMRCA", "LENGTH")
  
  
  get_time_in_generations_from_configfile_and_infile_name <-
    function(my_configfile, my_filename) {
      my_Ne <- read_yaml(my_configfile)["population_size"][[1]] %>%
        gsub("_", "", .) %>%
        as.numeric()
      
      my_ts <-
        strsplit(my_filename, "_")[[1]][12] %>% as.numeric() * my_Ne
      return(my_ts)
    }
  
  
  time_in_generations <-
    get_time_in_generations_from_configfile_and_infile_name("msprim_ts/config.yaml", infiles[1])
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
  
  
  message("create x marginal")
  xplot <- ggplot(tl_true, aes(LENGTH)) +
    geom_density(fill = "black") +
    scale_x_continuous(
      limits = c(4e1, 4e5),
      trans = "log10",
      breaks = sort(unique(as.vector(
        sapply(1:10, function(x) {
          10 ** c(1, 2, 3, 4, 5) * x
        })
      ))),
      labels = c(
        "10",
        rep("", 8),
        "100",
        rep("", 8),
        "1,000",
        rep("", 8),
        "10,000",
        rep("", 8),
        "100,000",
        rep("", 8),
        "1,000,000"
      )
    ) +
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
  yplot <- ggplot(tl_true, aes(TMRCA)) +
    geom_density(fill = "black") +
    scale_x_continuous(
      limits = c(1e4, 1e6),
      trans = "log10",
      breaks = sort(unique(as.vector(
        sapply(1:10, function(x) {
          10 ** c(1, 2, 3, 4, 5) * x
        })
      ))),
      labels = c(
        "10",
        rep("", 8),
        "100",
        rep("", 8),
        "1,000",
        rep("", 8),
        "10,000",
        rep("", 8),
        "100,000",
        rep("", 8),
        "1,000,000"
      )
    ) +
    theme(
      legend.position = "none",
      axis.line = element_line(size = 1.1),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks = element_blank()
    ) +
    coord_flip()
  
  message("create composition")
  ptot2 <- plot_grid(
    xplot,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    p,
    NULL,
    yplot,
    ncol = 3,
    align = "hv",
    rel_widths = c(2,-0.35, 1),
    rel_heights = c(1,-0.35, 2),
    scale = 1
  ) +
    theme(aspect.ratio = 1)
  
  message("print to file..")
  pdf(outfile,
      width = mwidth,
      height = mheight,
      useDingbats = FALSE)
  plot(ptot2)
  dev.off()
  
  png(paste0(outfile, ".png"), width = 4.4, height = 4.4, units = "in", res = 300)
  plot(ptot2)
  dev.off()
  message("printed to file")
}


# create plot along the sequence
tl_true$POS <- cumsum(tl_true$LENGT)
fill_rects$xmin <- rep(-Inf, nrow(fill_rects))
mdata <- tl_true %>%
  subset(POS <= 5e6)

{
  pseq <- ggplot() +
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
    geom_step(data = mdata, aes(POS, TMRCA), col = "black") +
    scale_colour_brewer() +
    # geom_hline(yintercept = max(fill_rects$ymin),
    #            col = mcol[2],
    #            size = 2) +
    scale_y_continuous(
      limits = c(1e4, 1e6),
      trans = "log10",
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_fill_brewer(palette = "Dark2") +
    # scale_fill_manual(values = mcol) +
    theme(
      aspect.ratio = 1 / 3,
      legend.position = "none",
      axis.line = element_line(size = 1.1),
      axis.text.y = element_text(angle = 90, hjust = 0.5)
    )
  
  pdf("pseq.pdf", width = mwidth, height = mheight, useDingbats = F)
  plot(pseq)
  dev.off()
  
  png(paste0("pseq.pdf", ".png"), width = 8.8, height = 4.4, units = "in", res = 300)
  plot(pseq)
  dev.off()
  }




# tm win; fig D, E
library(reshape2)

pd <- reticulate::import("pandas")
np <- reticulate::import("numpy")

if (T) {
  outfile <- paste0("tm_true_scenario_tNe_sigma_0.95_ts_", time_point, ".pdf")
  infiles <-
    c(
      paste0(
        "msprim_ts/data/df_tm_true/df_tm_true_scenario_tNe_sigma_0.95_ts_",
        time_point,
        "_.pickle.gzip"
      )
    )
  
}

breaks = np$load("msprim_ts/data/breaks/tmrca.bins.wc.0.ts.50000.0.95.3.4e-09.20000000.CS.4.npy") %>% 
  as.numeric() %>% 
  format(trim = T, digits = 2, scientific = T) %>% 
  as.character()

{
  message("outfile:\n\t", outfile)
  message("infiles:\n\t", paste(infiles, collapse = "\n\t"))
  
  # tl_true <- pd$read_pickle(infiles[1], "gzip")
  tm_true <- pd$read_pickle(infiles[1], "gzip")
  
  a <- apply(tm_true[,1:441], 2, mean)
  m <- matrix(a, nrow = 21)
  
  m <- m %>% melt()
  m$Var1 <- factor(rep(breaks, 21), levels = breaks)
  m$Var2 <- factor(rep(breaks, each=21), levels = breaks)
  names(m)[3] <- "tp"
  
  stopifnot((m %>%
               subset(Var1 == Inf) %>%
               subset(tp != 0) %>%
               nrow()) == 0)
  
  stopifnot((m %>%
               subset(Var2 == Inf) %>%
               subset(tp != 0) %>%
               nrow()) == 0)
  
  m <- m %>% subset(Var1 != Inf) %>% subset(Var2 != Inf)
  
  ptp_popsize <- m %>%
    ggplot(aes(Var1, Var2, fill=tp))+
    geom_tile()+
    xlab("n-th segment")+
    ylab("(n+1)-th segment")+
    scale_fill_gradient(
      low="#ffffc7", high="#7d0125",
      # low=wes_palette("IsleofDogs1")[5], high=wes_palette("IsleofDogs1")[2],
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
  
  pdf(outfile, width = mwidth, height = mheight, useDingbats = F)
  plot(ptp_popsize)
  dev.off()
  
  png(paste0(outfile, ".png"), width = 4.4, height = 4.4, units = "in", res = 300)
  plot(ptp_popsize)
  dev.off()
  }


if (T) {
  outfile <- paste0("tm_true_scenario_ts_sigma_0.95_ts_", time_point, ".pdf")
  infiles <-
    c(
      paste0(
        "msprim_ts/data/df_tm_true/df_tm_true_scenario_ts_sigma_0.95_ts_",
        time_point,
        "_.pickle.gzip"
      )
    )
  
}
{
  message("outfile:\n\t", outfile)
  message("infiles:\n\t", paste(infiles, collapse = "\n\t"))
  
  # tl_true <- pd$read_pickle(infiles[1], "gzip")
  tm_true <- pd$read_pickle(infiles[1], "gzip")
  
  a <- apply(tm_true[,1:441], 2, mean)
  m <- matrix(a, nrow = 21)
  
  m <- m %>% melt()
  m$Var1 <- factor(rep(breaks, 21), levels = breaks)
  m$Var2 <- factor(rep(breaks, each=21), levels = breaks)
  names(m)[3] <- "tp"
  
  stopifnot((m %>%
               subset(Var1 == Inf) %>%
               subset(tp != 0) %>%
               nrow()) == 0)
  
  stopifnot((m %>%
               subset(Var2 == Inf) %>%
               subset(tp != 0) %>%
               nrow()) == 0)
  
  m <- m %>% subset(Var1 != Inf) %>% subset(Var2 != Inf)
  
  ptp_sigma <- m %>%
    ggplot(aes(Var1, Var2, fill=tp))+
    geom_tile()+
    xlab("n-th segment")+
    ylab("(n+1)-th segment")+
    scale_fill_gradient(
      low="#ffffc7", high="#7d0125",
      # low=wes_palette("IsleofDogs1")[5], high=wes_palette("IsleofDogs1")[2],
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
  
  pdf(outfile, width = mwidth, height = mheight, useDingbats = F)
  plot(ptp_sigma)
  dev.off()
  
  png(paste0(outfile, ".png"), width = 4.4, height = 4.4, units = "in", res = 300)
  plot(ptp_sigma)
  dev.off()
}


# complete
row1 <- plot_grid(ptot, ptot2, nrow = 1, align = T)
row2 <- plot_grid(pseq)
row3 <- plot_grid(ptm_popsize, ptp_sigma, row=2, align = T)

my_final_figure <- plot_grid(row1, row2, row3, ncol=1, align = T,
                             rel_heights = c(2, 1, 2))

pdf("my_final_figure.pdf", width = 8.8, height = 9, useDingbats = F)
show(my_final_figure)
dev.off()









