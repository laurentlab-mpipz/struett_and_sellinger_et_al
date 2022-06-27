
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(wesanderson)

infile = snakemake@input$rds
outfile = snakemake@output$pdf

mcol = wes_palette("Darjeeling1")[c(1,2,3,5)]

df <- readRDS(infile)
df$t_sigma <- factor(df$t_sigma, levels = sort(unique(df$t_sigma)))

s <- ggplot(df, aes(t+1, sigma, col=inference_model, fill=file))+
  geom_step()+
  geom_vline(aes(xintercept=as.numeric(as.character(t_sigma))))+
  scale_x_continuous(trans = "log", limits = c(6000,NA))+
  scale_color_manual(values = mcol)+
  facet_grid(t_sigma~.)+
  theme(
    aspect.ratio = 0.707/2,
    legend.position = "bottom"
  )

d <- ggplot(df, aes(t+1, pop_size, col=inference_model, fill=file))+
  geom_hline(yintercept = 40000)+
  geom_hline(yintercept = 85000)+
  geom_step()+
  geom_vline(xintercept = 1e5)+
  scale_x_continuous(trans = "log", limits = c(6000,NA))+
  scale_y_continuous(trans = "log", limits = c(1.5e4, NA))+
  scale_color_manual(values = mcol)+
  facet_grid(t_sigma~.)+
  theme(
    aspect.ratio = 0.707/2,
    legend.position = "bottom"
  )

p <- plot_grid(s, d, nrow = 1, align = T, labels = "AUTO")

variable_hight = (0.5 + length(unique(df$t_sigma))) * 21/13.5

pdf(outfile, width = 12, height = variable_hight)
show(p)
dev.off()


