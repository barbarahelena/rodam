#### RODAM aldo-renin-arr

## Libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

## Output folder
resultsfolder <- "results/reninaldo"
dir.create(resultsfolder, showWarnings = FALSE)

## Load dataset
df_new <- rio:: import("data/clinicaldata.RDS")

(plaldo <- ggplot(data = df_new, aes(x = Site, y = Aldo, fill = Site)) +
        geom_violin()+
        geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
        theme_Publication() + 
        scale_fill_manual(values = pal_futurama()(4)[3:4], guide = "none") + 
        labs(title = "Aldosterone levels", y = "Aldosterone concentration (pg/ml)", x = "") +
        stat_compare_means(method = "wilcox.test", label.y = 450))
ggsave(plaldo, filename = "results/reninaldo/aldo.pdf", device = "pdf", width = 4, height = 5)
ggsave(plaldo, filename = "results/reninaldo/aldo.svg", device = "svg", width = 4, height = 5)

(plrenin <- ggplot(data = df_new, aes(x = Site, y = Renin, fill = Site)) +
        geom_violin()+
        geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
        theme_Publication() + 
        scale_fill_manual(values = pal_futurama()(4)[3:4], guide = "none") + 
        labs(title = "Renin levels", y = "Renin concentration (pg/ml)", x = "") +
        stat_compare_means(method = "wilcox.test", label.y = 100))
ggsave(plrenin, filename = "results/reninaldo/renin.pdf", device = "pdf", width = 4, height = 5)
ggsave(plrenin, filename = "results/reninaldo/renin.svg", device = "svg", width = 4, height = 5)

(plarr <- ggplot(data = df_new, aes(x = Site, y = ARR, fill = Site)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_manual(values = pal_futurama()(4)[3:4], guide = "none") + 
    labs(title = "Aldosterone-renin ratio", y = "Aldosterone-renin ratio", x = "") +
    stat_compare_means(method = "wilcox.test", label.y = 90))
ggsave(plarr, filename = "results/reninaldo/arr.pdf", device = "pdf", width = 4, height = 5)
ggsave(plarr, filename = "results/reninaldo/arr.svg", device = "svg", width = 4, height = 5)

## Ggarrange
ggarrange(plaldo, plrenin, plarr, labels = c("A", "B", "C"),
          nrow = 1)
ggsave("results/reninaldo/reninaldo.pdf", width = 12, height = 5)
ggsave("results/reninaldo/reninaldo.svg", width = 12, height = 5)
