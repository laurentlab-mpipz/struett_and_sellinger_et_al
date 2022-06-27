


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

mu_r_one_data$rho_on_theta = 1
mu_r_five_data$rho_on_theta = 5
mu_r_one_bgs_data$rho_on_theta = "1w/bgs"

df <- tibble(
  rbind(
    mu_r_five_data,
    mu_r_one_data,
    mu_r_one_bgs_data
  )
)

df$rho_on_theta <- factor(df$rho_on_theta, levels=c("1", "5", "1w/bgs"))
    
df %>%
  ggplot(aes(
        x = param_5, y = prop * 100, fill = discrete_bayes_factor
      )) +
      geom_area() +
      geom_hline(
        yintercept = 5,
        col = "white",
        size = 0.1
      ) +
      geom_hline(
        yintercept = 20,
        col = "white",
        size = 0.1
      ) +
      geom_hline(
        yintercept = 50,
        col = "white",
        size = 0.1
      ) +
      geom_hline(
        yintercept = 80,
        col = "white",
        size = 0.1
      ) +
      geom_hline(
        yintercept = 95,
        col = "white",
        size = 0.1
      ) +
      geom_vline(
        xintercept = 50000,
        col = "white",
        size = 0.1
      ) +
      # facet_grid(sscomp~regression_method+pls)+
      facet_grid(rho_on_theta + sscomp ~ pls) +
      scale_x_continuous(
        trans = "log",
        breaks = logbreak,
        labels = loglabel
      ) +
      scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
      scale_fill_manual(values = mcol) +
      theme(
        # axis.text.x = element_text(angle = 60, hjust = 1),
        aspect.ratio = 1,
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        # axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(-0.1, "lines"),
        panel.spacing = unit(-0.9, "lines")
      ) +
      labs(x = "time after transition [log10]",
           fill = "model support",
           main = "dada")
    
p1 <- df %>%
  subset(rho_on_theta == 1) %>%
  subset(sscomp == "sfs/ld") %>%
  subset(pls == "full") %>%
  ggplot(aes(x = param_5, y = prop * 100, fill = discrete_bayes_factor)) +
  geom_area() +
  geom_hline(yintercept = 5,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 20,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 50,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 80,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 95,
             col = "white",
             size = 0.1) +
  geom_vline(xintercept = 50000,
             col = "white",
             size = 0.1) +
  scale_x_continuous(trans = "log",
                     breaks = logbreak[logbreak >= 1e3],
                     labels = loglabel[logbreak >= 1e3]) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    legend.position = "none",
    aspect.ratio = 1,
    # text = element_text(size = 12),
  #   strip.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(size = 1),
    axis.ticks.length.x = unit(-1, "lines"),
    # panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support",
       main = "dada")
p1

png("sfs_ld_rmu1.png", width = 7, height = 7, units = "in", res = 100)
show(p1)
dev.off()


p1 <- df %>%
  subset(rho_on_theta == 1) %>%
  subset(sscomp == "tm_win") %>%
  subset(pls == "20-pls") %>%
  ggplot(aes(x = param_5, y = prop * 100, fill = discrete_bayes_factor)) +
  geom_area() +
  geom_hline(yintercept = 5,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 20,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 50,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 80,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 95,
             col = "white",
             size = 0.1) +
  geom_vline(xintercept = 50000,
             col = "white",
             size = 0.1) +
  scale_x_continuous(trans = "log",
                     breaks = logbreak[logbreak >= 1e3],
                     labels = loglabel[logbreak >= 1e3]) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    legend.position = "none",
    aspect.ratio = 1,
    # text = element_text(size = 12),
    #   strip.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(size = 1),
    axis.ticks.length.x = unit(-1, "lines"),
    # panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support",
       main = "dada")
p1

png("tm_win_rmu1.png", width = 7, height = 7, units = "in", res = 100)
show(p1)
dev.off()


p1 <- df %>%
  subset(rho_on_theta == 1) %>%
  subset(sscomp %in% c("sfs/ld/tm_win")) %>%
  subset(pls == "20-pls") %>%
  ggplot(aes(x = param_5, y = prop * 100, fill = discrete_bayes_factor)) +
  geom_area() +
  geom_hline(yintercept = 5,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 20,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 50,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 80,
             col = "white",
             size = 0.1) +
  geom_hline(yintercept = 95,
             col = "white",
             size = 0.1) +
  geom_vline(xintercept = 50000,
             col = "white",
             size = 0.1) +
  scale_x_continuous(trans = "log",
                     breaks = logbreak[logbreak >= 1e3],
                     labels = loglabel[logbreak >= 1e3]) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    legend.position = "none",
    aspect.ratio = 1,
    # text = element_text(size = 12),
    #   strip.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(size = 1),
    axis.ticks.length.x = unit(-1, "lines"),
    # panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support",
       main = "dada")
p1

png("sfs_ld_tm_win_rmu1.png", width = 7, height = 7, units = "in", res = 100)
show(p1)
dev.off()


## find the upper time estimate
a <- df %>%
  subset(rho_on_theta == 1) %>%
  subset(sscomp %in% c("sfs/ld/tm_win")) %>%
  subset(pls == "20-pls") %>%
  subset(!discrete_bayes_factor %in% c("Negative", "Barely worth mentioning")) %>%
  select(!matches("N")) %>%
  select(!matches("pls")) %>%
  select(!matches("sscomp")) %>%
  select(!matches("pod")) %>%
  select(!matches("sscomp")) %>%
  group_by(param_5) %>%
  summarise(
    support = sum(prop)
  )

# get max point that is above threshold
thres <- 0.8
time.max.low <- a$param_5[1]
prop.max.low <- a$support[1]
for (i in 2:nrow(a)) {
  time.here <- a$param_5[i]
  prop.here <- a$support[i]
  
  if (prop.here < thres) {
    break
  } else {
    time.max.low <- time.here
    prop.max.low <- prop.here
  }
}
time.min.high <- a$param_5[nrow(a)]
prop.min.high <- a$support[nrow(a)]
for (i in (nrow(a)-1):1) {
  time.here <- a$param_5[i]
  prop.here <- a$support[i]
  
  if (prop.here < thres) {
    time.min.high <- time.here
    prop.min.high <- prop.here
  } else {
    break
  }
}

x1 <- time.max.low; x2 <- time.min.high
y1 <- prop.max.low; y2 <- prop.min.high


plot(a$param_5, a$support, log="x", type = "b")
abline(h = 0.8, col="red")
abline(h = y1, v=x1)
abline(h = y2, v=x2)


last_time <- 0
reg.fun <- function(x) return(-1.5e-5 * x + 1.57)
for (x in seq(50000, 60000, .001)) {
  if (reg.fun(x)<thres) break
  else last_time <- x
}
cat("last time point that is at least 0.8 good: ", last_time)
last_time_Ne <- last_time / 20200



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
png("sfs_ld_rmu1.tsigma.param_estim.png", width = 7, height = 7, units = "in", res = 100)
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
png("tm_win_rmu1.tsigma.param_estim.png", width = 7, height = 7, units = "in", res = 100)
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
png("sfs_ld_tm_win_rmu1.tsigma.param_estim.png", width = 7, height = 7, units = "in", res = 100)
show(p1)
dev.off()

