## Diversity metrics

## Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ggsci)
library(phyloseq)

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

## Load data
phydata <- readRDS("data/phyloseq_sampledata.RDS")
df_new <- rio:: import("data/clinicaldata.RDS")
tab <- as.data.frame(t(as(phydata@otu_table, 'matrix')))
tab_matrix <- t(as(phydata@otu_table, 'matrix'))
# counts <- sample_sums(phydata@otu_table)
# counts # samples should all sum up to 14932

## Output folder
resultsfolder <- "results/alphadiversity"
dir.create(resultsfolder, showWarnings = FALSE)

# Diversity metrics
## Shannon plot urban-rural
shannon <- vegan::diversity(tab, index = 'shannon')
df_shan <- data.frame(ID = names(shannon), shannon = shannon)
df_shan <- left_join(df_shan, df_new, by = "ID")
(plshan <- ggplot(data = df_shan, aes(x = Site, y = shannon, fill = Site)) +
    geom_violin() +
    scale_fill_manual(values = pal_cosmic()(4)[2:4], guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    stat_compare_means(label.y = 5.5) +
    labs(title = "Shannon index", y = "Shannon index", x="") + 
    theme_Publication())
ggsave(plshan, filename = "results/alphadiversity/shannon.svg", width = 6, height = 5)
ggsave(plshan, filename = "results/alphadiversity/shannon.pdf", width = 6, height = 5)

(plshan <- ggplot(data = df_shan, aes(x = SodiumInt, y = shannon, color = SodiumInt)) +
        geom_jitter() +
        geom_smooth(method = "lm")+
        # scale_fill_manual(values = pal_cosmic()(4)[2:4], guide = "none") +
        # geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
        # stat_compare_means(label.y = 5.5) +
        labs(title = "Shannon index", y = "Shannon index", x="") + 
        theme_Publication())

## Species richness Urban-rural
specrich <- specnumber(tab)
dfspec <- data.frame(ID = names(specrich), richness = specrich)
dfspec <- left_join(dfspec, df_new, by = "ID")

(plrich <- ggplot(data = dfspec, aes(x = Site, y = richness, fill = Site)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_manual(values = pal_cosmic()(4)[2:4], guide = "none") + 
    labs(title = "Richness", y = "Number of ASVs", x = "") +
    stat_compare_means(label.y = 800))
ggsave(plrich, filename = "results/alphadiversity/richness.pdf", width = 4, height = 5)
ggsave(plrich, filename = "results/alphadiversity/richness.svg", width = 4, height = 5)

(plrich <- ggplot(data = dfspec, aes(x = SodiumInt, y = richness, color = SodiumInt)) +
        geom_jitter()+
        geom_smooth(method = "lm") +
        stat_cor()+
        # geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
        theme_Publication() + 
        # scale_fill_manual(values = pal_cosmic()(4)[2:4], guide = "none") + 
        labs(title = "Richness", y = "Number of ASVs", x = "") +
        stat_compare_means(label.y = 610))

## Faith's PD urban-rural
faith <- picante::pd(samp = tab_matrix, tree = phydata@phy_tree)
dffai <- as.data.frame(faith)
dffai$ID <- rownames(faith)
dffai <- left_join(dffai, df_new, by = "ID")

(plfaith <- ggplot(data = dffai, aes(x = Site, y = PD, fill = Site)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_manual(values = pal_cosmic()(4)[2:4], guide = "none") + 
    labs(title = "Faith's PD", y = "Faith's phylogenetic diversity", x = "") +
    stat_compare_means(label.y = 80))
ggsave(plfaith, filename = "results/alphadiversity/faiths.pdf", device = "pdf", width = 4, height = 5)
ggsave(plfaith, filename = "results/alphadiversity/faiths.svg", device = "svg", width = 4, height = 5)

(plrich <- ggplot(data = dffai, aes(x = SodiumInt, y = PD, color = SodiumInt)) +
        geom_jitter()+
        geom_smooth(method = "lm") +
        stat_cor()+
        theme_Publication() + 
        labs(title = "Faiths", y = "Number of ASVs", x = "") +
        stat_compare_means(label.y = 80))

## Ggarrange urban-rural
pl_total <- ggarrange(plshan, plrich, plfaith, 
                      # labels = c("A", "B", "C"), 
                      nrow =1)
ggsave(pl_total, filename = "results/alphadiversity/alphadivplots.pdf", width = 11, height = 5)
ggsave(pl_total, filename = "results/alphadiversity/alphadivplots.svg", width = 11, height = 5)

