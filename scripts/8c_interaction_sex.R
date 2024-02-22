## Interactions sex

library(rio)
library(phyloseq)
library(dplyr)
library(tidyr)
library(forcats)
library(patchwork)
library(compositions)

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
resultsfolder <- "results/interactions_sex"
dir.create(resultsfolder, showWarnings = FALSE)

interactions <- function(df, dfname){
    ## Variable selection
    dfsub <- df %>% select(Site, Sex, tail(names(.), 20))
    
    ## Models
    res_ia <- c()
    for (i in c((ncol(dfsub)-19):ncol(dfsub))){
        dfsub$asv <- NULL    
        dfsub$asv <- dfsub[[i]]                   
        m0 <- lm(asv ~ Site + Site*Sex, data = dfsub)
        
        tax <- colnames(dfsub)[i]
        m0 <- tidy(m0, conf.int=T)[4,]
        
        resRow <- cbind(tax, m0$term, exp(m0$estimate), exp(m0$conf.low), exp(m0$conf.high), m0$p.value)
        colnames(resRow) <- c("Metabolite", "interaction",
                              "m0-est", "m0-l95", "m0-u95", "m0-p")
        res_ia <- rbind(res_ia, resRow)
        dfsub$met <- NULL 
    }
    
    resia <- as.data.frame(res_ia)
    afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
    afronden5 <- function(x) return(as.numeric(format(round(x, 5),5)))
    resia2 <- resia %>% 
        mutate_at(c(3:6), as.character) %>% 
        mutate_at(c(3:6), as.numeric) %>% 
        mutate_at(c(3:5), afronden2) %>% 
        mutate_at(c(6), afronden5) %>% 
        mutate(
            `m0-q` = p.adjust(`m0-p`, 'fdr')
        )
    openxlsx::write.xlsx(resia2, file.path("results/interactions_sex", str_c(dfname,"_interactions.xlsx")))
}

## Opening RODAM file
rodam_mb <- readRDS("data/phyloseq_sampledata.RDS")
rodamclin <- readRDS("data/clinicaldata.RDS")
rodamotu <- as.data.frame(t(as(rodam_mb@otu_table, "matrix")))
tax <- readRDS("data/tax_table.RDS")
best_pred <- rio::import('rural_urban/output_XGB_class_rural_urban_2024_01_14__00-47-28/feature_importance.txt')

## Preparation for models
best_pred <- best_pred %>% arrange(-RelFeatImp) %>% slice(1:20)
dfmb <- rodamotu %>% select(best_pred$FeatName)
dfmb2 <- as.data.frame(clr(dfmb+1))
dfmb2$ID <- rownames(dfmb2)
df <- left_join(rodamclin, dfmb2, by='ID')
dim(df)
colnames(df)[(ncol(df)-19):ncol(df)] <- make.unique(tax$Tax[match(colnames(df)[(ncol(df)-19):ncol(df)], tax$ASV)])
dfmb2 <- df[(ncol(df)-19):ncol(df)]

interactions(df, "sex")

