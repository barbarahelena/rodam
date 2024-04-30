## Process results machine learning

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(stringr)
library(forcats)
library(aplot)
library(ComplexHeatmap)
library(ggplotify)

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
                legend.key.size= unit(0.4, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

colfam <- list(
    "Bifidobacteriaceae" = "#DC9445",
    "Christensenellaceae" = "#358DB9",
    "Clostridiaceae" = "#BE4A90",
    "Prevotellaceae" = "#4CA56E",
    "Lactobacillaceae" = "#2F509E",
    "Lachnospiraceae" = "#CD7467",
    "Oscillospiraceae" = "#77567D",
    "Rikenellaceae" = "#2E2A2B",
    "Streptococcaceae" = "#4F7CB2",
    "Succinivibrionaceae" = "#84A29C",
    "Unknown family" = "#58545E"
)

# cols <- colorRampPalette(c(pal_cosmic(palette = "hallmarks_light")(10)))

draw_corrplot <- function(df){
    rescor <-   rcorr(as.matrix(df), type="spearman")
    corplot <- corrplot(rescor$r,  type = "upper", tl.col = "black", tl.cex = 0.6, cl.cex = 0.6, number.cex = 0.5,
                        order = 'hclust', hclust.method="ward.D", tl.srt = 45, insig="blank", sig.level = 0.05,
                        p.mat=rescor$P, method="color", mar=c(0,0,1,0), addCoef.col = "black", diag = F, 
                        addgrid.col = "grey")
}

## Process ML results
source("scripts/functions.R")
### Rural
path_true <- 'rural_urban_dich/output_XGB_class_rural_urban_dich_2024_02_27__16-57-54'
pl1 <- plot_feature_importance_color_microbiome(path_true, 10)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
plarr1 <- ggarrange(svg_grob, pl1, nrow = 1, widths = c(1.0, 1.3))
plarr1b <- annotate_figure(plarr1, top = text_grob("Rural - Urban Ghana", color = "black", face = "bold", size = 14))
### Urban
path_true <- 'urban_ams_dich/output_XGB_class_urban_ams_2024_02_27__22-38-49'
pl2 <- plot_feature_importance_color_microbiome(path_true, 10)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
plarr2 <- ggarrange(svg_grob, pl2, nrow = 1, widths = c(1.0, 1.3))
plarr2b <- annotate_figure(plarr2, top = text_grob("Urban Ghana - Amsterdam", color = "black", face = "bold", size = 14))

## Make prevalence plot
mb <- readRDS("data/phyloseq_sampledata.RDS")
feat1 <- rio::import('rural_urban_dich/output_XGB_class_rural_urban_dich_2024_02_27__16-57-54/feature_importance.txt') %>% slice(1:10)
feat2 <- rio::import('urban_ams_dich/output_XGB_class_urban_ams_2024_02_27__22-38-49/feature_importance.txt') %>% slice(1:10)
feats <- rbind(feat1, feat2) %>% filter(!duplicated(FeatName))
dim(feats)
mat <- t(as(mb@otu_table, "matrix"))
mat <- mat[,feats$FeatName]
mat <- ifelse(mat > 0, 1, 0)
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
vanish <- df %>% select(Site, 68:ncol(.)) %>% 
    pivot_longer(., 2:ncol(.), names_to = "ASV", values_to = "abundance") %>% 
    group_by(Site, ASV) %>% summarise(mean = mean(abundance)) %>% 
    pivot_wider(., id_cols = "ASV", names_from = "Site", values_from = "mean") %>% 
    filter(`Amsterdam` < `Rural Ghana` & `Urban Ghana` < `Rural Ghana`)
blossum <- df %>% select(Site, 68:ncol(.)) %>% 
    pivot_longer(., 2:ncol(.), names_to = "ASV", values_to = "abundance") %>% 
    group_by(Site, ASV) %>% summarise(mean = mean(abundance)) %>% 
    pivot_wider(., id_cols = "ASV", names_from = "Site", values_from = "mean") %>% 
    filter(`Rural Ghana` < `Amsterdam` & `Urban Ghana` < `Amsterdam`)
df <- df %>% mutate(across(c(68:ncol(.)), as.factor)) %>% 
                mutate(across(c(68:ncol(.)), ~fct_recode(.x, "Present"="1", "Absent"="0")))
    
head(df)[1:5,1:5]


dflong <- df %>% select(ID, Site, 68:ncol(.)) %>% pivot_longer(., 3:21, names_to = "ASV")
dfsum <- dflong %>% group_by(Site, ASV) %>% summarise(percentage = mean(value == "Present")*100) %>% 
    mutate(Tax = factor(ASV, levels = ASV, labels = make.unique(as.character(tax$Tax[match(ASV, tax$ASV)]))),
           Tax = fct_reorder(Tax, match(ASV, tax$ASV)),
           )

tax2 <- tax %>% mutate(Family = case_when(duplicated(Family) ~ "",
                                          .default = Family),
                       vanishblossum = case_when(
                           ASV %in% vanish$ASV ~ "vanish",
                           ASV %in% blossum$ASV ~ "blossum"
                       ),
                       vanishblossum = as.factor(vanishblossum)
                       )
(bar <- ggplot(tax, aes(x = Tax, y = 0, fill = Family)) +
    geom_tile() +
    geom_text(data = tax2, aes(label = Family, y = 0.5), vjust = 0, hjust = 0, angle = 20,
              size = rel(3), check_overlap = TRUE, position = "identity")+
    coord_cartesian(clip = "off", ylim = c(0,3))+
    theme_void() +
    scale_fill_manual(values = colfam, guide = "none"))

(bar2 <- ggplot(tax2, aes(x = Tax, y = 0, fill = vanishblossum)) +
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
        scale_fill_manual(values = pal_cosmic()(7)[c(7,6)]))

(plot <- ggplot(data = dfsum, aes(x = Tax, group = Site, fill = Site, y = percentage)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") + 
    theme_Publication() +
    scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
    labs(y = "Presence %", x = "", fill = "") +
    theme(axis.text.x = element_blank(),
            legend.direction = "vertical",
          legend.position = "right"))

# options("aplot_guides" = "keep")
(comppl <- plot %>% insert_bottom(bar2, height = .05) %>% insert_top(bar, height = .18))
(compl_anno <- annotate_figure(as.ggplot(comppl), 
                               top = text_grob("Prevalence of best predictors", 
                                               color = "black", face = "bold", size = 14)))
# ggsave("results/vanishblossum/barplot.pdf", width = 12, height = 6)

## Heatmap
mat <- t(as(mb@otu_table, "matrix"))
mat <- ifelse(mat > 0, 1, 0)
dim(mat)
head(mat)[1:5,1:5]
tk <- apply(mat, 2, function(x) sum(x > 0) > (0.3*length(x)))
dfmb <- mat[,tk]
dfmb <- as.data.frame(dfmb)
dfmb$ID <- rownames(dfmb)
df <- left_join(clin, dfmb, by = "ID")
rownames(df) <- df$ID
mbdf <- df %>% select(ID, Site, 68:ncol(.)) %>% 
    arrange(Site)
ordercol <- rev(order(colSums(mbdf[,4:ncol(mbdf)])))
matdf <- as.matrix(mbdf[4:ncol(mbdf)])
rownames(matdf) <- mbdf$ID
numb <- summary(mbdf$Site)
# pdf("results/heatmap.pdf", width = 15, height = 12)
hm <- Heatmap(matdf, name = "presence", 
        col = c("white", "#BE4A90"),
        row_split = factor(c(rep(names(numb[1]), numb[[1]]),
                             rep(names(numb[2]), numb[[2]]),
                             rep(names(numb[3]), numb[[3]])
        ), 
        levels = c(names(numb)[1:3])),
        cluster_row_slices = TRUE,
        row_order = 1:nrow(matdf),
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_order = ordercol
)
# dev.off()

(pl_total <- ggarrange(as.ggplot(hm), plarr1b, plarr2b, compl_anno, 
                       nrow = 4, labels = c("A", "B", "C", "D"),
                       heights = c(0.5, 0.7, 0.7, 0.9)))
ggsave("results/vanishblossom/presenceabsence.pdf", height = 25, width = 14)
# ggsave("results/vanishblossom/presenceabsence.svg", height = 25, width = 14)

# ## Data Model Rural - urban
# mat <- t(as(mb@otu_table, "matrix"))
# mat <- mat[,feat1$FeatName]
# mat <- ifelse(mat > 0, 1, 0)
# colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
# dfmb <- as.data.frame(mat)
# colnames(dfmb) <- case_when(colnames(dfmb) == "Bifidobacterium breve/kashiwanohense/longum" ~
#                                 "Bifidobacterium brev/kashiw/long",
#                             .default = colnames(dfmb))
# dfmb$ID <- rownames(mat)
# 
# df <- left_join(clin, dfmb, by = "ID") %>% 
#     mutate(across(c(68:ncol(.)), as.factor)) %>% 
#     mutate(across(c(68:ncol(.)), ~fct_recode(.x, "Present"="1", "Absent"="0")))
# 
# dflong <- df %>% select(ID, Site, 67:ncol(.)) %>% pivot_longer(., 3:ncol(.), names_to = "ASV")
# dfsum <- dflong %>% group_by(Site, ASV) %>% summarise(percentage = mean(value == "Present")*100)
# 
# (pl1 <- ggplot(data = dfsum, aes(x = ASV, group = Site, fill = Site, y = percentage)) +
#     geom_bar(stat = "identity", position = "dodge") + 
#     theme_Publication() +
#     scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
#     labs(y = "Presence %", x = "", fill = "", title = "Predictors rural - urban Ghana") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)))
# ggsave(pl1, filename = "results/vanishingmicrobes_grouped_1.pdf", width = 10, height = 5, )
# 
# ## Data model urban-ams
# mat <- t(as(mb@otu_table, "matrix"))
# mat <- mat[,feat2$FeatName]
# mat <- ifelse(mat > 0, 1, 0)
# colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
# dfmb <- as.data.frame(mat)
# dfmb$ID <- rownames(mat)
# 
# df <- left_join(clin, dfmb, by = "ID") %>% 
#     mutate(across(c(56:65), as.factor)) %>% 
#     mutate(across(c(56:65), ~fct_recode(.x, "Present"="1", "Absent"="0")))
# 
# dflong <- df %>% select(ID, Site, 68:ncol(.)) %>% pivot_longer(., 3:ncol(.), names_to = "ASV")
# dfsum <- dflong %>% group_by(Site, ASV) %>% summarise(percentage = mean(value == "Present")*100)
# 
# pl2 <- ggplot(data = dfsum, aes(x = ASV, group = Site, fill = Site, y = percentage)) +
#     geom_bar(stat = "identity", position = "dodge") + 
#     theme_Publication() +
#     scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
#     labs(y = "Presence %", x = "", fill = "", title = "Predictors urban Ghana - Amsterdam") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggsave(pl2, filename = "results/vanishingmicrobes_grouped_2.pdf", width = 10, height = 5)
# 
# ggarrange(pl1, pl2, nrow = 2, 
#           common.legend = TRUE, 
#           legend = "bottom",
#           labels = c("A", "B"))
# ggsave("results/vanishingmicrobes_total.pdf", width = 10, height = 10)
# 
# 
# ## Total
# mb <- readRDS("data/phyloseq_sampledata.RDS")
# feat1 <- rio::import('rural_urban_dich/output_XGB_class_rural_urban_dich_2024_02_27__16-57-54/feature_importance.txt') %>% slice(1:20)
# feat2 <- rio::import('urban_ams_dich/output_XGB_class_urban_ams_2024_02_27__22-38-49/feature_importance.txt') %>% slice(1:20)
# feats <- rbind(feat1, feat2) %>% filter(!duplicated(FeatName))
# dim(feats)
# mat <- t(as(mb@otu_table, "matrix"))
# mat <- mat[,feats$FeatName]
# dim(mat)
# head(mat)[1:5,1:5]
# tax <- readRDS("data/taxtable.RDS")
# colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
# dfmb <- as.data.frame(mat)
# dfmb$ID <- rownames(mat)
# clin <- readRDS("data/clinicaldata.RDS")
# df <- left_join(clin, dfmb, by = "ID")
# dflong <- df %>% select(ID, Site, 68:ncol(.)) %>% pivot_longer(., 3:40, names_to = "ASV") %>% 
#     mutate(value = log10(value + 1))
# dfwide <- pivot_wider(dflong, id_cols = 1, names_from = "ASV", values_from = "value")
# dfwide <- dfwide[,-1]
# names(dfwide)
# 
# pdf("results/corrplot.pdf", width = 15, height = 15)
#     draw_corrplot(dfwide)
# dev.off()
# 
# 
# ## spirochetes
# clin <- readRDS("data/clinicaldata.RDS")
# mb <- readRDS("data/phyloseq_sampledata.RDS")
# mb <- subset_taxa(mb, Family=="Spirochaetaceae")
# mat <- t(as(mb@otu_table, "matrix"))
# tax <- readRDS("data/taxtable.RDS")
# colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
# dfmb <- as.data.frame(mat)
# dfmb$ID <- rownames(mat)
# df <- left_join(clin, dfmb, by = "ID")
# dflong <- df %>% select(ID, Site, 67:ncol(.)) %>% pivot_longer(., 3:ncol(.), names_to = "ASV") %>% 
#     mutate(value = (value / 15000 ) *100 )
# dflong %>% group_by(ASV) %>% summarise(mean = mean(value))
# dfwide <- pivot_wider(dflong, id_cols = 1, names_from = "ASV", values_from = "value")
# (pl3 <- ggplot(data = dflong, aes(x = ASV, y = value + 0.00001, 
#                                   fill = Site, group = interaction(ASV,Site))) +
#         geom_boxplot(width = 0.3) +
#         theme_Publication() +
#         scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
#         labs(y = "Abundance in %", x = "", fill = "", title = "Spirochaetaceae") +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#         scale_y_log10())
# ggsave(pl3, filename = "results/spirochaetaceae_abundance.pdf", width = 7, height = 5)
# 
# mb <- readRDS("data/phyloseq_sampledata.RDS")
# mb <- subset_taxa(mb, Family=="Spirochaetaceae")
# mat <- t(as(mb@otu_table, "matrix"))
# mat <- ifelse(mat > 0, 1, 0)
# tax <- readRDS("data/taxtable.RDS")
# colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
# dfmb <- as.data.frame(mat)
# dfmb$ID <- rownames(mat)
# 
# df <- left_join(clin, dfmb, by = "ID") %>% 
#     mutate(across(c(67:ncol(.)), as.factor)) %>% 
#     mutate(across(c(67:ncol(.)), ~fct_recode(.x, "Present"="1", "Absent"="0")))
# 
# head(df)[1:5,1:5]
# colnames(df)
# 
# dflong <- df %>% select(ID, Site, 67:ncol(.)) %>% pivot_longer(., 3:ncol(.), names_to = "ASV")
# dftotal <- dflong
# dfsumtot <- dflong %>% group_by(ASV) %>% summarise(percentage = mean(value == "Present")*100) %>% print(.)
# mean(dfsumtot$percentage)
# dfsum <- dflong %>% group_by(Site, ASV) %>% summarise(percentage = mean(value == "Present")*100)
# dfwide <- dfsum %>% pivot_wider(., id_cols = "Site", names_from = "ASV", values_from = "percentage")
# 
# (pl3 <- ggplot(data = dfsum, aes(x = ASV, group = Site, fill = Site, y = percentage)) +
#         geom_bar(stat = "identity", position = "dodge") + 
#         theme_Publication() +
#         scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
#         labs(y = "Presence %", x = "", fill = "", title = "Spirochaetaceae") +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1)))
# ggsave(pl3, filename = "results/spirochaetaceae_prevalence.pdf", width = 10, height = 5)

