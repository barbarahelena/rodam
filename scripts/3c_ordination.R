## Diversity metrics

## Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(doParallel)
registerDoParallel(6)

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

#### Load data ####
phydata <- readRDS("data/phyloseq_sampledata.RDS")
df_new <- rio::import("data/clinicaldata.RDS")
tab <- as.data.frame(t(as(phydata@otu_table, 'matrix')))
tab_matrix <- t(as(phydata@otu_table, 'matrix'))

#### Output folder ####
resultsfolder <- "results/ordination"
dir.create(resultsfolder, showWarnings = FALSE)

#### Bray-Curtis distance: male-female ####
print('Bray-Curtis distance total dataset')
bray <- vegan::vegdist(tab, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
head(expl_variance_bray)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$ID <- rownames(dbray)
dbray <- left_join(dbray, df_new, by = 'ID') # add metadata / covariates

print('PERMANOVA..')
set.seed(1234)
dfanova <- df_new %>%
    slice(match(sample_names(phydata), ID)) # distance matrix and metadata must have the same sample order
all(dfanova$ID == sample_names(phydata)) # TRUE
dim(dfanova)
res1 <- adonis2(bray ~ Site, data = dfanova) # PERMANOVA
print(res1)

print('plotting..')
(braycurt <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    geom_point(aes(color = Site), size = 1, alpha = 0.7) +
    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
    scale_color_manual(values = pal_cosmic()(4)[c(2:4)]) +
    scale_fill_manual(values = pal_cosmic()(4)[c(2:4)], guide = "none") +
    theme_Publication() +
    labs(color = "", fill = "", title = "PCoA Bray-Curtis distance") +
    stat_ellipse(geom = "polygon", aes(color = Site, fill = Site), type = "norm", 
                 alpha = 0.1, linewidth = 1.0) + 
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res1$`Pr(>F)`, ", r2 = ", format(round(res1$R2[1],3), nsmall = 3))))
ggsave(braycurt, filename = "results/ordination/PCoA_BrayCurtis.pdf", device = "pdf", width = 5, height = 5)
ggsave(braycurt, filename = "results/ordination/PCoA_BrayCurtis.svg", device = "svg", width = 5, height = 5)

#### Weighted UniFrac ####
print('Weighted UniFrac')
wunifrac <- UniFrac(phydata, normalized = T, weighted = T, parallel = T)
pcoord <- ape::pcoa(wunifrac, correction = "cailliez")
expl_var_wu <- pcoord$values$Rel_corr_eig * 100
head(expl_var_wu)
dfpc <- pcoord$vectors[, c('Axis.1', 'Axis.2')] # get PCoA coordinates
dfpc <- as.data.frame(dfpc)
dfpc$ID <- rownames(dfpc)
dfpc <- left_join(dfpc, df_new, by = 'ID') # add metadata / covariates

print('PERMANOVA..')
res2 <- adonis2(wunifrac ~ Site, data = dfanova) # PERMANOVA
print(res2)

print('plotting..')
(unifracpl <- dfpc %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    geom_point(aes(color = Site), size = 1, alpha = 0.7) +
    theme_Publication() +
    stat_ellipse(geom = "polygon", aes(color = Site, fill = Site), type = "norm", alpha = 0.1,
                 linewidth = 1.0) + 
    scale_color_manual(values = pal_cosmic()(4)[c(2:4)]) +
    scale_fill_manual(values = pal_cosmic()(4)[c(2:4)], guide = "none") +
    labs(title = 'PCoA Weighted UniFrac', color = "", fill = "",
         x = paste0('PCo1 (', round(expl_var_wu[1], digits = 1),'%)'),
         y = paste0('PCo2 (', round(expl_var_wu[2], digits = 1),'%)')) +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res1$`Pr(>F)`, ", r2 = ", format(round(res2$R2[1],3), nsmall = 3))))
ggsave(unifracpl, filename = "results/ordination/PCoA_WeightedUnifrac.pdf", device = "pdf", width = 5, height = 5)
ggsave(unifracpl, filename = "results/ordination/PCoA_WeightedUnifrac.svg", device = "svg", width = 5, height = 5)

(unifrac_ab <- dfpc %>% 
        ggplot(aes(Axis.1, Axis.2)) +
        geom_point(aes(color = FecalSample_AB), size = 1, alpha = 0.7) +
        xlab(paste0('PCo1 (', round(expl_var_wu[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_var_wu[2], digits = 1),'%)')) +
        theme_Publication() +
        scale_color_manual(values = pal_cosmic()(3)) +
        ggtitle('PCoA Weighted UniFrac') +
        guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
        labs(color = "Antibiotics use") +
        stat_ellipse(aes(color = FecalSample_AB), type = "norm"))
ggsave(unifrac_ab, filename = "results/ordination/PCoA_WeightedUnifrac_AB.pdf", device = "pdf", width = 5, height = 5)
ggsave(unifrac_ab, filename = "results/ordination/PCoA_WeightedUnifrac_AB.svg", device = "svg", width = 5, height = 5)
