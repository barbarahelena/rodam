#### RODAM diet descriptives

## Libraries
library(dplyr)
library(ggsci)
library(ggplot2)
library(forcats)
library(ggpubr)

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
df_new <- rio:: import("data/clinicaldata.RDS")

## PCA diet
df_diet <- df_new %>% dplyr::select(ID, Site, TotalCalories, Fibre, Proteins, 
                             Carbohydrates, Fat, SodiumInt) %>% 
    filter(!is.na(TotalCalories)) %>% 
    filter(TotalCalories < 7500)
df_diet2 <- df_diet %>% dplyr::select(-ID, -Site, -TotalCalories) %>% mutate(across(everything(.), scale))
matdiet <- as.matrix(df_diet2)
#matdiet <- scale(matdiet)
tunediet <- tune.pca(matdiet, ncomp = 5, scale = TRUE)
plot(tunediet)
pc <- mixOmics::pca(matdiet, ncomp = 4)
pcs <- as.data.frame(pc$variates$X)
pcs <- pcs %>% mutate(ID = df_diet$ID, Site = df_diet$Site)
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


df <- left_join(df_new, pc2, by = c("ID", "Site")) %>% dplyr::select(everything(.), DietPC1=PC1, DietPC2=PC2)
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


# pl7 <- ggplot(df_diet, aes(x=Site, y=(AlcoholIntake+0.01)))+
#     geom_violin(aes(fill=Site))+
#     scale_fill_manual(values = pal_cosmic()(4)[c(2,3,4)], guide = FALSE)+
#     geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
#     theme_Publication()+
#     scale_y_log10()+
#     theme(legend.position = 'none')+
#     labs(x='', y = 'log10(grams)', title = "Alcohol")+
#     ggpubr::stat_compare_means()
# pl7

ggarrange(pl1, pl2, pl3, pl4, pl5, pl6, nrow = 2, ncol = 3,
          labels = c(LETTERS[1:6]))
ggsave("results/diet/dietaryintake.pdf", width = 15, height = 12)
ggsave("results/diet/dietaryintake.svg", width = 15, height = 12)

