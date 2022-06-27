# old plots


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
                     breaks = logbreak,
                     labels = loglabel) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    # axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none",
    aspect.ratio = 1,
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    # axis.text = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "lines"),
    panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support",
       main = "dada")

p2 <- df %>%
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
                     breaks = logbreak,
                     labels = loglabel) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    # axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none",
    aspect.ratio = 1,
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    # axis.text = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "lines"),
    panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support",
       main = "dada")

p3 <- df %>%
  subset(rho_on_theta == 1) %>%
  subset(sscomp == "sfs/ld/tm_win") %>%
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
                     breaks = logbreak,
                     labels = loglabel) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    # axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none",
    aspect.ratio = 1,
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    # axis.text = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "lines"),
    panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support",
       main = "dada")

p4 <- df %>%
  subset(rho_on_theta == 5) %>%
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
                     breaks = logbreak,
                     labels = loglabel) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    # axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none",
    aspect.ratio = 1,
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    # axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "lines"),
    panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support")

p5 <- df %>%
  subset(rho_on_theta == 5) %>%
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
                     breaks = logbreak,
                     labels = loglabel) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    # axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none",
    aspect.ratio = 1,
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    # axis.text = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "lines"),
    panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support",
       main = "dada")

p6 <- df %>%
  subset(rho_on_theta == 5) %>%
  subset(sscomp == "sfs/ld/tm_win") %>%
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
                     breaks = logbreak,
                     labels = loglabel) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    # axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none",
    aspect.ratio = 1,
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    # axis.text = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "lines"),
    panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support",
       main = "dada")

p7 <- df %>%
  subset(rho_on_theta == "1w/bgs") %>%
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
                     breaks = logbreak,
                     labels = loglabel) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    # axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none",
    aspect.ratio = 1,
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    # axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "lines"),
    panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support")

p8 <- df %>%
  subset(rho_on_theta == "1w/bgs") %>%
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
                     breaks = logbreak,
                     labels = loglabel) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    # axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none",
    aspect.ratio = 1,
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    # axis.text = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "lines"),
    panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support",
       main = "dada")

p9 <- df %>%
  subset(rho_on_theta == "1w/bgs") %>%
  subset(sscomp == "sfs/ld/tm_win") %>%
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
                     breaks = logbreak,
                     labels = loglabel) +
  scale_y_continuous(breaks = c(5, 95), labels = c("0%", "100%")) +
  scale_fill_manual(values = mcol) +
  theme(
    # axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none",
    aspect.ratio = 1,
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    # axis.text = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "lines"),
    panel.spacing = unit(-0.9, "lines")
  ) +
  labs(x = "time after transition [log10]",
       fill = "model support",
       main = "dada")



p <- plot_grid(p1, p2, p3,
               NULL, NULL, NULL,
               p4, p5, p6,
               NULL, NULL, NULL,
               p7, p8, p9,
               ncol = 3,
               align = T,
               labels = c(
                 "a", "b", "c",
                 "", "", "",
                 "d", "e", "f",
                 "", "", "",
                 "g", "h", "i"
               ),
               rel_heights = c(1, -0.1, 1, -0.1, 1)
)

title <- ggdraw() + 
  draw_label("Model choice 'transition to selfing' vs 'constant recombination'",
             fontface='bold')

p <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))

pdf("tsabc_model_choice_incl_bgs.pdf",
    width = 12,
    height = 13)
show(p)
dev.off()





