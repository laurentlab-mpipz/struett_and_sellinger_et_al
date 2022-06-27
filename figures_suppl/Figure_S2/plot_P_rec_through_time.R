library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(wesanderson)
# library(reshape2)
library(yaml)


message("this is a complete manual plot")


m1 <- "Constant outcrossing"
m2 <- "Outcrossing to selfing"
m3 <- "Selfing to outcrossing"

t1 <- 0.35
t2 <- 1 - t1
shift <- 0.005

df <- rbind.data.frame(
  cbind.data.frame(
    time = c(0, 1),
    p = c(0, 1) + shift,
    model = m1
  ),
  cbind.data.frame(
    time = c(0, t1, 1),
    p = c(0, t1, 0.5),
    model = m2
  ),
  cbind.data.frame(
    time = c(0, t2, 1),
    p = c(0, t2 - 0.5, 0.5),
    model = m3
  )
)

df$model <- factor(df$model, levels = c(m1, m3, m2))
mcol <- wes_palette("Darjeeling1")[c(1:3)]

p <- df %>% ggplot(aes(time, p, col = model)) +
  geom_vline(xintercept = t1, col = mcol[3], linetype="dashed", alpha = 0.4) +
  geom_vline(xintercept = t2, col = mcol[2], linetype="dashed", alpha = 0.4) +
  geom_line(size=1) +
  scale_y_continuous(limits = c(0, 1 + shift),
                     labels = 0,
                     breaks = 0) +
  scale_x_continuous(limits = c(0, 1),
                     labels = 0,
                     breaks = 0) +
  scale_color_manual(values = mcol) +
  theme(
    aspect.ratio = 0.707,
    legend.position = "bottom",
    axis.line = element_line(size = 1.1, arrow = grid::arrow(length = unit(0.3, "cm"))),
    # axis.ticks = element_blank(),
    axis.text = element_blank()
  )

mwidth <- 4.4
mheight <- mwidth

pdf("prec.pdf", width = mwidth, height = mheight, useDingbats = F)
plot(p)
dev.off()

png(paste0("prec.pdf", ".png"), width = 8.8, height = 4.4, units = "in", res = 300)
plot(p)
dev.off()

