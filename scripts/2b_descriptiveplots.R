#### RODAM diet descriptives

## Libraries
library(dplyr)
library(ggsci)
library(ggplot2)
library(forcats)
library(ggpubr)
library(mixOmics)

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
resultsfolder <- "results/diet"
dir.create(resultsfolder, showWarnings = FALSE)

## Load dataset
df_new <- readRDS("data/clinicaldata.RDS")
vanish <- readRDS("data/ids_vanishing.RDS")
blossom <- readRDS("data/ids_blossom.RDS")
df_new <- left_join(df_new, left_join(vanish, blossom, by = "ID"), by = "ID")

## PCA diet
df_diet <- df_new %>% dplyr::select(ID, Site, vanishing, blossom, TotalCalories, Fibre, Proteins, 
                             Carbohydrates, Fat, SodiumInt) %>% 
    filter(!is.na(TotalCalories)) %>% 
    filter(TotalCalories < 7500)
df_diet2 <- df_diet %>% dplyr::select(-ID, -Site, -vanishing, -blossom, -TotalCalories) %>% 
    mutate(across(everything(.), scale))
matdiet <- as.matrix(df_diet2)
#matdiet <- scale(matdiet)
tunediet <- tune.pca(matdiet, ncomp = 5, scale = TRUE)
plot(tunediet)
pc <- mixOmics::pca(matdiet, ncomp = 4)
pcs <- as.data.frame(pc$variates$X)
pcs <- pcs %>% mutate(ID = df_diet$ID, Site = df_diet$Site, vanishing = df_diet$vanishing,
                      blossom = df_diet$blossom
                      )
expvar_diet <- pc$prop_expl_var$X[1:4]
loadings <- as.data.frame(pc$loadings$X)
loadings$Variables <- rownames(loadings)
# loadings <- loadings %>% filter(Variables %in% c("TotalCalories", "SodiumInt")) %>% 
#     mutate(Variables = fct_recode(Variables, "Sodium" = "SodiumInt"))
(pcadiet <- pcs %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = Site), size = 1, alpha = 1.0) +
        xlab(paste0('PC1 (', round(expvar_diet[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_diet[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = Site, fill = Site), linewidth = 1.0,
                     alpha = 0.1, type = "norm")+
        scale_color_manual(values = pal_cosmic()(4)[c(2,3,4)]) +
        scale_fill_manual(values = pal_cosmic()(4)[c(2,3,4)], guide = "none") +
        labs(color = "", title = "PCA diet")+
        geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), 
                     arrow = arrow(length = unit(1/2, "picas")),
                     color = "black", linewidth = 0.9) +
        annotate("text", x = (loadings$PC1*9), y = (loadings$PC2*9),
                 label = loadings$Variables)
)
ggsave(pcadiet, filename = "results/diet/PCA_diet_loading_3.pdf", device = "pdf", 
       width = 5, height = 6)
ggsave(pcadiet, filename = "results/diet/PCA_diet_loading_3.svg", device = "svg", 
       width = 5, height = 6)

(pcadiet <- pcs %>% 
        ggplot(aes(PC3, PC4)) +
        geom_point(aes(color = Site), size = 1, alpha = 1.0) +
        xlab(paste0('PC3 (', round(expvar_diet[3]*100, digits = 1),'%)')) +
        ylab(paste0('PC4 (', round(expvar_diet[4]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = Site, fill = Site), linewidth = 1.0,
                     alpha = 0.1, type = "norm")+
        scale_color_manual(values = pal_cosmic()(4)[c(2,3,4)]) +
        scale_fill_manual(values = pal_cosmic()(4)[c(2,3,4)], guide = "none") +
        labs(color = "", title = "PCA diet")+
        geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC3*2), yend = (PC4*2)), 
                     arrow = arrow(length = unit(1/2, "picas")),
                     color = "black", linewidth = 0.9) +
        annotate("text", x = (loadings$PC3*3), y = (loadings$PC4*3),
                 label = loadings$Variables)
)


df <- left_join(df_new, pcs, by = c("ID", "Site")) %>% 
    dplyr::select(everything(.), DietPC1=PC1, DietPC2=PC2)
saveRDS(df, "data/clinicaldata_pcdiet.RDS")

comp <- list(c("Rural Ghana", "Urban Ghana"), c("Urban Ghana", "Amsterdam"), 
             c("Rural Ghana", "Amsterdam"))

pl1 <- ggplot(df_diet, aes(x=Site, y=TotalCalories))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_cosmic()(4)[c(2,3,4)], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'Kcal', title = "Total calories")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl1

pl2 <- ggplot(df_diet, aes(x=Site, y=Fibre))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_cosmic()(4)[c(2,3,4)], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Fibre")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl2

pl3 <- ggplot(df_diet, aes(x=Site, y=Proteins))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_cosmic()(4)[c(2,3,4)], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Proteins")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl3

pl4 <- ggplot(df_diet, aes(x=Site, y=Carbohydrates))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_cosmic()(4)[c(2,3,4)], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Carbohydrates")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl4

pl5 <- ggplot(df_diet, aes(x=Site, y=Fat))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_cosmic()(4)[c(2,3,4)], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Fat")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl5

pl6 <- ggplot(df_diet, aes(x=Site, y=SodiumInt))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_cosmic()(4)[c(2,3,4)], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'mg', title = "Sodium")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl6

ggarrange(ggarrange(pl1, pl2, pl3, pl4, pl5, pl6, nrow = 2, ncol = 4,
          labels = c(LETTERS[1:6])), 
          ggarrange(pcadiet, NULL, labels = c("G", ""), widths = c(1.5, 0.5)),
            nrow = 2, heights = c(2.0, 1.5))
ggsave("results/diet/dietaryintake.pdf", width = 12, height = 17)
ggsave("results/diet/dietaryintake.svg", width = 12, height = 17)



(pcadietvanish <- pcs %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = vanishing), size = 1, alpha = 0.75) +
        xlab(paste0('PC1 (', round(expvar_diet[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_diet[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = vanishing, fill = vanishing), linewidth = 1.0,
                     alpha = 0.1, type = "norm")+
        scale_color_manual(values = c("grey70",pal_cosmic()(7)[6])) +
        scale_fill_manual(values = c("grey70",pal_cosmic()(7)[6]), guide = "none") +
        labs(color = "", title = "PCA diet - vanish")+
        geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*2), yend = (PC2*2)), 
                     arrow = arrow(length = unit(1/2, "picas")),
                     color = "black", linewidth = 0.9) +
        annotate("text", x = (loadings$PC1*3), y = (loadings$PC2*3),
                 label = loadings$Variables)
)

comp <- list(c("controls", "vanishing"))
pl1 <- ggplot(df_diet, aes(x=vanishing, y=TotalCalories))+
    geom_violin(aes(fill=vanishing))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[6]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'Kcal', title = "Total calories")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl1

pl2 <- ggplot(df_diet, aes(x=vanishing, y=Fibre))+
    geom_violin(aes(fill=vanishing))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[6]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Fibre")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl2

pl3 <- ggplot(df_diet, aes(x=vanishing, y=Proteins))+
    geom_violin(aes(fill=vanishing))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[6]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Proteins")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl3

pl4 <- ggplot(df_diet, aes(x=vanishing, y=Carbohydrates))+
    geom_violin(aes(fill=vanishing))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[6]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Carbohydrates")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl4

pl5 <- ggplot(df_diet, aes(x=vanishing, y=Fat))+
    geom_violin(aes(fill=vanishing))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[6]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Fat")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl5

pl6 <- ggplot(df_diet, aes(x=vanishing, y=SodiumInt))+
    geom_violin(aes(fill=vanishing))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[6]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'mg', title = "Sodium")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl6

ggarrange(ggarrange(pl1, pl2, pl3, pl4, pl5, pl6, nrow = 2, ncol = 4,
                    labels = c(LETTERS[1:6])), 
          ggarrange(pcadietvanish, NULL, labels = c("G", ""), widths = c(1.5, 0.5)),
          nrow = 2, heights = c(2.0, 1.5))
ggsave("results/diet/vanish_diet.pdf", width = 12, height = 17)
ggsave("results/diet/vanish_diet.svg", width = 12, height = 17)


(pcadietblossom <- pcs %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = blossom), size = 1, alpha = 1.0) +
        xlab(paste0('PC1 (', round(expvar_diet[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_diet[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = blossom, fill = blossom), linewidth = 1.0,
                     alpha = 0.1, type = "norm")+
        scale_color_manual(values = c("grey70",pal_cosmic()(7)[7])) +
        scale_fill_manual(values = c("grey70",pal_cosmic()(7)[7]), guide = "none") +
        labs(color = "", title = "PCA diet - blossom")+
        geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*2), yend = (PC2*2)), 
                     arrow = arrow(length = unit(1/2, "picas")),
                     color = "black", linewidth = 0.9) +
        annotate("text", x = (loadings$PC1*3), y = (loadings$PC2*3),
                 label = loadings$Variables)
)

comp <- list(c("controls", "blossom"))
pl1 <- ggplot(df_diet, aes(x=blossom, y=TotalCalories))+
    geom_violin(aes(fill=blossom))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[7]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'Kcal', title = "Total calories")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl1

pl2 <- ggplot(df_diet, aes(x=blossom, y=Fibre))+
    geom_violin(aes(fill=blossom))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[7]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Fibre")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl2

pl3 <- ggplot(df_diet, aes(x=blossom, y=Proteins))+
    geom_violin(aes(fill=blossom))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[7]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Proteins")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl3

pl4 <- ggplot(df_diet, aes(x=blossom, y=Carbohydrates))+
    geom_violin(aes(fill=blossom))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[7]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Carbohydrates")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl4

pl5 <- ggplot(df_diet, aes(x=blossom, y=Fat))+
    geom_violin(aes(fill=blossom))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[7]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Fat")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl5

pl6 <- ggplot(df_diet, aes(x=blossom, y=SodiumInt))+
    geom_violin(aes(fill=blossom))+
    scale_fill_manual(values = c("grey70",pal_cosmic()(7)[7]), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'mg', title = "Sodium")+
    ggpubr::stat_compare_means(method = "t.test", comparisons = comp,
                               tip.length = 0, label="p.signif", hide.ns = TRUE)
pl6

ggarrange(ggarrange(pl1, pl2, pl3, pl4, pl5, pl6, nrow = 2, ncol = 4,
                    labels = c(LETTERS[1:6])), 
          ggarrange(pcadietblossom, NULL, labels = c("G", ""), widths = c(1.5, 0.5)),
          nrow = 2, heights = c(2.0, 1.5))
ggsave("results/diet/blossom_diet.pdf", width = 12, height = 17)
ggsave("results/diet/blossom_diet.svg", width = 12, height = 17)


ggarrange(pcadietvanish, pcadietblossom, labels = c("A", "B"), nrow = 1)
ggsave("results/diet/vanishblossom_diet.pdf", width = 10, height = 5)
ggsave("results/diet/vanishblossom_diet.svg", width = 10, height = 5)
