## Process results machine learning

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(stringr)
library(forcats)
library(ComplexHeatmap)

## Data
mb <- readRDS("data/phyloseq_sampledata.RDS")
mat <- t(as(mb@otu_table, "matrix"))
mat <- ifelse(mat > 0, 1, 0)
dim(mat)
head(mat)[1:5,1:5]
tk <- apply(mat, 2, function(x) sum(x > 0) > (0.3*length(x)))
dfmb <- mat[,tk]
dfmb <- as.data.frame(dfmb)
dfmb$ID <- rownames(dfmb)
clin <- readRDS("data/clinicaldata.RDS")
df <- left_join(clin, dfmb, by = "ID")
rownames(df) <- df$ID
mbdf <- df %>% select(ID, Site, 67:ncol(.)) %>% 
    arrange(Site)
ordercol <- rev(order(colSums(mbdf[,4:ncol(mbdf)])))
matdf <- as.matrix(mbdf[4:ncol(mbdf)])
rownames(matdf) <- mbdf$ID
numb <- summary(mbdf$Site)
pdf("results/heatmap.pdf", width = 15, height = 12)
    Heatmap(matdf, name = "presence", 
            col = c("white", "dodgerblue4"),
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
dev.off()

