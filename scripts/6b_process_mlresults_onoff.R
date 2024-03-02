## Process results machine learning

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(stringr)
library(forcats)

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

draw_corrplot <- function(df){
    rescor <-   rcorr(as.matrix(df), type="spearman")
    corplot <- corrplot(rescor$r,  type = "upper", tl.col = "black", tl.cex = 0.6, cl.cex = 0.6, number.cex = 0.5,
                        order = 'hclust', hclust.method="ward.D", tl.srt = 45, insig="blank", sig.level = 0.05,
                        p.mat=rescor$P, method="color", mar=c(0,0,1,0), addCoef.col = "black", diag = F, 
                        addgrid.col = "grey")
    # ord <- dimnames(corplot)[[1]]
    # df_ord <- df[,match(ord, names(df))]
    # rescor2 <- rcorr(as.matrix(df_ord), type="spearman")
    # mycol <- ifelse((rescor2$P < 0.05), 1,0)
    # mycol <- rescor2$r * mycol
    # mycol <- ifelse(mycol==0, NA, mycol)
    # corrplot::corrplot(mycol,  type = "upper", tl.col = "black", tl.cex = 0.6, cl.cex = 0.6, number.cex = 0.6,
    #                    order = 'original', tl.srt = 45, insig="blank", sig.level = 0.05,
    #                    p.mat=rescor2$P, method="color", mar=c(0,0,1,0), addCoef.col = "black", diag = F, 
    #                    addgrid.col = "grey")
    
}

source("scripts/functions.R")

path_true <- 'rural_urban_dich/output_XGB_class_rural_urban_dich_2024_02_27__16-57-54'
path_permuted <- 'rodam/rural_urban_dich/output_XGB_class_rural_urban_dich_2024_02_27__17-22-42_PERMUTED'
data_path <- 'rural_urban_dich/input_data'
labels <- c("Urban", "Rural")

plot_feature_importance_class(path_true, 20)
plot_feature_importance_color_microbiome(path_true, 20)

path_true <- 'urban_ams_dich/output_XGB_class_urban_ams_2024_02_27__22-38-49'
path_permuted <- 'urban_ams_dich/output_XGB_class_urban_ams_2024_02_27__22-55-49_PERMUTED'
data_path <- 'urban_ams_dich/input_data'
labels <- c("Amsterdam", "Urban Ghana")

plot_feature_importance_class(path_true, 20)
plot_feature_importance_color_microbiome(path_true, 20)

## Data
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
colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
dfmb <- as.data.frame(mat)
dfmb$ID <- rownames(mat)

clin <- readRDS("data/clinicaldata.RDS")
df <- left_join(clin, dfmb, by = "ID") %>% 
    mutate(across(c(56:74), as.factor)) %>% 
    mutate(across(c(56:74), ~fct_recode(.x, "Present"="1", "Absent"="0")))
    
head(df)[1:5,1:5]
colnames(df)

dflong <- df %>% select(ID, Site, 56:74) %>% pivot_longer(., 3:21, names_to = "ASV")
dftotal <- dflong
dfsum <- dflong %>% group_by(Site, ASV) %>% summarise(percentage = mean(value == "Present")*100)
dfwide <- dfsum %>% pivot_wider(., id_cols = "Site", names_from = "ASV", values_from = "percentage")

plots <- list()
for(a in 2:20){
    dfwide$var <- dfwide[[a]]
    varname <- colnames(dfwide)[a]
    print(varname)
    pl <- ggplot(data = dfwide, aes(x = Site, y = var, fill = Site)) +
        geom_bar(stat = "identity") +
        theme_Publication() +
        scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
        labs(title = varname, y = "Presence %", x = "", fill = "") +
        theme(plot.title = element_text(size = rel(0.7)),
              axis.text.x = element_blank())
    plots[[a-1]] <- pl
    dfwide$var <- NULL
}
ggarrange(plotlist = plots, nrow = 4, ncol = 5, 
          common.legend = TRUE, legend = "bottom")
ggsave("results/vanishingmicrobes.pdf", width = 15, height = 20)

ggplot(data = dfsum, aes(x = ASV, group = Site, fill = Site, y = percentage)) +
    geom_bar(stat = "identity", position = "dodge") + 
    theme_Publication() +
    scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
    labs(y = "Presence %", x = "", fill = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/vanishingmicrobes_grouped.pdf", width = 20, height = 10)


## Data Model Rural - urban
mat <- t(as(mb@otu_table, "matrix"))
mat <- mat[,feat1$FeatName]
mat <- ifelse(mat > 0, 1, 0)
colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
dfmb <- as.data.frame(mat)
colnames(dfmb) <- case_when(colnames(dfmb) == "Bifidobacterium breve/kashiwanohense/longum" ~
                                "Bifidobacterium breve",
                            .default = colnames(dfmb))
dfmb$ID <- rownames(mat)

df <- left_join(clin, dfmb, by = "ID") %>% 
    mutate(across(c(56:65), as.factor)) %>% 
    mutate(across(c(56:65), ~fct_recode(.x, "Present"="1", "Absent"="0")))

dflong <- df %>% select(ID, Site, 56:65) %>% pivot_longer(., 3:12, names_to = "ASV")
dfsum <- dflong %>% group_by(Site, ASV) %>% summarise(percentage = mean(value == "Present")*100)

(pl1 <- ggplot(data = dfsum, aes(x = ASV, group = Site, fill = Site, y = percentage)) +
    geom_bar(stat = "identity", position = "dodge") + 
    theme_Publication() +
    scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
    labs(y = "Presence %", x = "", fill = "", title = "Predictors rural - urban Ghana") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(pl1, filename = "results/vanishingmicrobes_grouped_1.pdf", width = 10, height = 5, )

## Data model urban-ams
mat <- t(as(mb@otu_table, "matrix"))
mat <- mat[,feat2$FeatName]
mat <- ifelse(mat > 0, 1, 0)
colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
dfmb <- as.data.frame(mat)
dfmb$ID <- rownames(mat)

df <- left_join(clin, dfmb, by = "ID") %>% 
    mutate(across(c(56:65), as.factor)) %>% 
    mutate(across(c(56:65), ~fct_recode(.x, "Present"="1", "Absent"="0")))

dflong <- df %>% select(ID, Site, 56:65) %>% pivot_longer(., 3:12, names_to = "ASV")
dfsum <- dflong %>% group_by(Site, ASV) %>% summarise(percentage = mean(value == "Present")*100)

pl2 <- ggplot(data = dfsum, aes(x = ASV, group = Site, fill = Site, y = percentage)) +
    geom_bar(stat = "identity", position = "dodge") + 
    theme_Publication() +
    scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
    labs(y = "Presence %", x = "", fill = "", title = "Predictors urban Ghana - Amsterdam") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(pl2, filename = "results/vanishingmicrobes_grouped_2.pdf", width = 10, height = 5)

ggarrange(pl1, pl2, nrow = 2, 
          common.legend = TRUE, 
          legend = "bottom",
          labels = c("A", "B"))
ggsave("results/vanishingmicrobes_total.pdf", width = 10, height = 10)


## Relative abundances of vanishing microbes
## Data
mb <- readRDS("data/phyloseq_sampledata.RDS")
feat1 <- rio::import('rural_urban_dich/output_XGB_class_rural_urban_dich_2024_02_27__16-57-54/feature_importance.txt') %>% slice(1:10)
feat2 <- rio::import('urban_ams_dich/output_XGB_class_urban_ams_2024_02_27__22-38-49/feature_importance.txt') %>% slice(1:10)
feats <- rbind(feat1, feat2) %>% filter(!duplicated(FeatName))
dim(feats)
mat <- t(as(mb@otu_table, "matrix"))
mat <- mat[,feats$FeatName]
dim(mat)
head(mat)[1:5,1:5]
tax <- readRDS("data/taxtable.RDS")
colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
dfmb <- as.data.frame(mat)
dfmb$ID <- rownames(mat)

clin <- readRDS("data/clinicaldata.RDS")
df <- left_join(clin, dfmb, by = "ID")

## Data Model Rural - urban
mat <- t(as(mb@otu_table, "matrix"))
mat <- mat[,feat1$FeatName]
colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
dfmb <- as.data.frame(mat)
colnames(dfmb) <- case_when(colnames(dfmb) == "Bifidobacterium breve/kashiwanohense/longum" ~
                                "Bifidobacterium breve",
                            .default = colnames(dfmb))
dfmb$ID <- rownames(mat)
df <- left_join(clin, dfmb, by = "ID")
dflong <- df %>% select(ID, Site, 56:65) %>% pivot_longer(., 3:12, names_to = "ASV") %>% 
    mutate(value = (value / 15000 ) *100 )
df1 <- dflong

options(scipen=100)
(pl1 <- ggplot(data = dflong, aes(x = ASV, y = value + 0.00001, 
                                  fill = Site, group = interaction(ASV,Site))) +
        geom_boxplot(width = 0.3) +
        theme_Publication() +
        scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
        labs(y = "Abundance in counts", x = "", fill = "", title = "Predictors rural - urban Ghana") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_log10())
# ggsave(pl1, filename = "results/vanishingmicrobes_grouped_1.pdf", width = 10, height = 5, )

## Data model urban-ams
mat <- t(as(mb@otu_table, "matrix"))
mat <- mat[,feat2$FeatName]
colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
dfmb <- as.data.frame(mat)
dfmb$ID <- rownames(mat)

df <- left_join(clin, dfmb, by = "ID")
dflong <- df %>% select(ID, Site, 56:65) %>% pivot_longer(., 3:12, names_to = "ASV") %>% 
    mutate(value = (value / 15000 ) *100 )
df2 <- dflong

options(scipen=100)
(pl2 <- ggplot(data = dflong, aes(x = ASV, y = value + 0.00001, 
                                  fill = Site, group = interaction(ASV,Site))) +
        geom_boxplot(width = 0.3) +
        theme_Publication() +
        scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
        labs(y = "Abundance in %", x = "", fill = "", title = "Predictors urban Ghana - Amsterdam") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_log10())

# ggsave(pl2, filename = "results/vanishingmicrobes_grouped_2.pdf", width = 10, height = 5)

ggarrange(pl1, pl2, nrow = 2, 
          common.legend = TRUE, 
          legend = "bottom",
          labels = c("A", "B"))
ggsave("results/vanishingmicrobes_abundance.pdf", width = 10, height = 10)


## Total
mb <- readRDS("data/phyloseq_sampledata.RDS")
feat1 <- rio::import('rural_urban_dich/output_XGB_class_rural_urban_dich_2024_02_27__16-57-54/feature_importance.txt') %>% slice(1:10)
feat2 <- rio::import('urban_ams_dich/output_XGB_class_urban_ams_2024_02_27__22-38-49/feature_importance.txt') %>% slice(1:10)
feats <- rbind(feat1, feat2) %>% filter(!duplicated(FeatName))
dim(feats)
mat <- t(as(mb@otu_table, "matrix"))
mat <- mat[,feats$FeatName]
dim(mat)
head(mat)[1:5,1:5]
tax <- readRDS("data/taxtable.RDS")
colnames(mat) <- make.unique(tax$Tax[match(colnames(mat), tax$ASV)])
dfmb <- as.data.frame(mat)
dfmb$ID <- rownames(mat)
clin <- readRDS("data/clinicaldata.RDS")
df <- left_join(clin, dfmb, by = "ID")
dflong <- df %>% select(ID, Site, 56:74) %>% pivot_longer(., 3:21, names_to = "ASV") %>% 
    mutate(value = log10(value + 1))
dfwide <- pivot_wider(dflong, id_cols = 1, names_from = "ASV", values_from = "value")
dfwide <- dfwide[,-1]
names(dfwide)

pdf("results/corrplot.pdf", width = 15, height = 15)
    draw_corrplot(dfwide)
dev.off()
