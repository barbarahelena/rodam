## Process results machine learning
options(scipen=999)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(stringr)
library(forcats)

colfam <- list(
    "Bacteroidaceae" = "#2F509E",
    "Bifidobacteriaceae" = "#DC9445",
    "Clostridiaceae" = "#BE4A90",
    "Enterobacteriaceae" = "#8C57A2",
    "Lachnospiraceae" = "#CD7467",
    "Lactobacillaceae" = "#2F509E",
    "Oscillospiraceae" = "#77567D",
    "Prevotellaceae" = "#4CA56E",
    "Ruminococcaceae" = "#CF4E9C",
    "Unknown family" = "#58545E"
)


source("scripts/functions.R")

path_true <- 'rural_urban/output_XGB_class_rural_urban_2024_02_05__21-24-45'
path_permuted <- 'rural_urban/output_XGB_class_rural_urban_2024_02_05__21-37-26_PERMUTED'
data_path <- 'rural_urban/input_data'
labels <- c("Urban", "Rural")

# plot_feature_importance_class(path_true, 20)
pl1 <- plot_feature_importance_color_microbiome(path_true, 20)
# plot_features_tests_class(data_path, path_true, top_n=20, labels)
# plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)

grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
plarr1 <- ggarrange(svg_grob, pl1, nrow = 1, widths = c(1.0, 1.3))
plarr1b <- annotate_figure(plarr1, top = text_grob("Rural Ghana - Urban Ghana", color = "black", face = "bold", size = 14))

path_true <- 'urban_ams/output_XGB_class_urban_ams_2024_02_21__14-01-19'
path_permuted <- 'urban_ams/output_XGB_class_urban_ams_2024_02_21__14-10-44_PERMUTED'
data_path <- 'urban_ams/input_data'
labels <- c("Amsterdam", "Urban Ghana")

# plot_feature_importance_class(path_true, 20)
pl2 <- plot_feature_importance_color_microbiome(path_true, 20)
# plot_features_tests_class(data_path, path_true, top_n=20, labels)
# plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)

grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
plarr2 <- ggarrange(svg_grob, pl2, nrow = 1, widths = c(1.0, 1.3))
plarr2b <- annotate_figure(plarr2, top = text_grob("Urban Ghana - Amsterdam", color = "black", face = "bold", size = 14))


## Open data
mb <- readRDS("data/phyloseq_sampledata.RDS")
feat1 <- rio::import('rural_urban/output_XGB_class_rural_urban_2024_02_05__21-24-45/feature_importance.txt') %>% slice(1:10)
feat2 <- rio::import('urban_ams/output_XGB_class_urban_ams_2024_02_21__14-01-19/feature_importance.txt') %>% slice(1:10)
feats <- rbind(feat1, feat2) %>% filter(!duplicated(FeatName))
dim(feats)
mat <- t(as(mb@otu_table, "matrix"))
mat <- mat[,feats$FeatName]
dim(mat)
head(mat)[1:5,1:5]
tax <- readRDS("data/taxtable.RDS")
tax <- tax %>% filter(ASV %in% colnames(mat)) %>% arrange(Family) %>% 
    mutate(Family = ifelse(is.na(Family), "Unknown family", Family),
           Family = as.factor(Family),
           Family = fct_inorder(Family),
           ASV = as.factor(ASV),
           ASV = fct_inorder(ASV),
           Tax = make.unique(Tax),
           Tax = factor(Tax, levels = Tax),
           Tax = fct_recode(Tax,"Bifidobacterium brev/kashiw/long"="Bifidobacterium breve/kashiwanohense/longum"))
mat <- mat[,match(tax$ASV,colnames(mat))]
dfmb <- as.data.frame(mat)
dfmb$ID <- rownames(mat)

clin <- readRDS("data/clinicaldata.RDS")
df <- left_join(clin, dfmb, by = "ID")
head(df)[1:5,1:5]

dflong <- df %>% select(ID, Site, 68:ncol(.)) %>% pivot_longer(., 3:20, names_to = "ASV") %>% 
    mutate(Tax = factor(ASV, levels = ASV, labels = make.unique(as.character(tax$Tax[match(ASV, tax$ASV)]))),
           Tax = fct_reorder(Tax, match(ASV, tax$ASV)),
           value = (value / 15000 ) *100 
    )

tax2 <- tax %>% mutate(Family = case_when(duplicated(Family) ~ "",
                                          .default = Family),
                       model = case_when(ASV %in% feat1$FeatName & ASV %in% feat2$FeatName ~ "both",
                                        ASV %in% feat1$FeatName ~ "rural-urban",
                                         ASV %in% feat2$FeatName ~ "urban-Ams",
                                         ),
                       model = as.factor(model)
)
(bar <- ggplot(tax, aes(x = Tax, y = 0, fill = Family)) +
        geom_tile() +
        geom_text(data = tax2, aes(label = Family, y = 0.5), vjust = 0, hjust = 0, angle = 20,
                  size = rel(3), check_overlap = TRUE, position = "identity")+
        coord_cartesian(clip = "off", ylim = c(0,3))+
        theme_void() +
        scale_fill_manual(values = colfam, guide = "none"))

(bar2 <- ggplot(tax2, aes(x = Tax, y = 0, fill = model)) +
        labs(x = "", y = "")+
        geom_tile(alpha = 0.7) +
        coord_cartesian(clip = "off", ylim = c(-0.5,0.5))+
        theme_Publication()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = "right",
              plot.margin = unit(c(0,0,0,0),"mm"),
              legend.title = element_blank(),
              legend.key.size = unit(0.4, "cm"),
              legend.spacing = unit(0, "cm")) +
        scale_fill_manual(values = pal_cosmic()(5)[c(5,2,3)]))

(plot <- ggplot(data = dflong, aes(x = Tax, y = value + 0.00001, 
                                fill = Site, group = interaction(ASV,Site))) +
                    geom_boxplot(width = 0.5, outlier.size = 0.8, position = position_dodge(0.6)) +
                    theme_Publication() +
                    scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
                    labs(y = "Abundance in %", x = "", fill = "", title = "") +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    scale_y_log10()+
                    theme(axis.text.x = element_blank(),
                          legend.direction = "vertical",
                          legend.position = "right"))
(comppl <- plot %>% insert_bottom(bar2, height = .05) %>% insert_top(bar, height = .18))
(compl_anno <- annotate_figure(as.ggplot(comppl), 
                               top = text_grob("Abundance of best predictors", 
                                               color = "black", face = "bold", size = 14)))

## Assemble figure
ggarrange(plarr1b, plarr2b, compl_anno, nrow = 3, labels = c("A", "B", "C"))
ggsave("results/ml_abundance/ml_abundance.pdf", width = 17, height = 22)
ggsave("results/ml_abundance/ml_abundance.svg", width = 17, height = 22)
