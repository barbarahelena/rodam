## Linear models
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
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
resultsfolder <- "results/glm_clr"
dir.create(resultsfolder, showWarnings = FALSE)

## Functions
theme_Empty <- function() {
    theme_minimal() +
        theme(axis.title = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
              axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), 
              axis.ticks.y = element_blank(),
              plot.margin=unit(c(4,0,0,0), "mm")
        )
}

distribution <- function(df, col){
    var <- df[[col]]
    pl <- ggplot(df) + 
        geom_density(aes(var), fill = pal_cosmic()(6)[6], color = pal_cosmic()(6)[6]) +      
        theme_Empty()
}

afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
afronden3 <- function(x) return(as.numeric(format(round(x, 3),3)))

## Opening RODAM file
rodam_mb <- readRDS("data/phyloseq_sampledata.RDS")
rodamclin <- readRDS("data/clinicaldata_pcdiet.RDS")
dietids <- c("G0421", "G2545")
rodamclin <- rodamclin %>% mutate(
    TotalCalories = case_when(ID %in% dietids ~ NA, .default = TotalCalories),
    Proteins = case_when(ID %in% dietids ~ NA, .default = Proteins),
    Fat = case_when(ID %in% dietids ~ NA, .default = Fat),
    Fibre = case_when(ID %in% dietids ~ NA, .default = Fibre),
    Carbohydrates = case_when(ID %in% dietids ~ NA, .default = Carbohydrates),
    SodiumInt = case_when(ID %in% dietids ~ NA, .default = SodiumInt),
    AlcoholIntake = log(AlcoholIntake+0.01)
)
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

df_rur <- df %>% filter(Site == "Rural Ghana")
df_urb <- df %>% filter(Site == "Urban Ghana")

## Regression models
restotal <- c()
for (i in c(55:59)) {
    df$asv <- df[[i]]
    taxasv <- colnames(df)[i]
    for(j in c(32,34,35,37)){
        df$diet <- scale(df[[j]])
        m0 <- lm(asv ~ diet, data = df)
        m1 <- lm(asv ~ diet + Age + Sex + BMI + DM + BristolScale, data = df)
        dietname <- colnames(df)[j]
        m0 <- tidy(m0, conf.int=T)[2,]
        m1 <- tidy(m1, conf.int=T)[2,]
        resRow <- cbind(taxasv, dietname, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                        m1$estimate, m1$conf.low, m1$conf.high, m1$p.value)
        colnames(resRow) <- c("ASV", "Diet",
                              "m0-est", "m0-l95", "m0-u95", "m0-p", 
                              "m1-est", "m1-l95", "m1-u95", "m1-p")
        restotal <- rbind(restotal, resRow)
        df$diet <- NULL
    }
    df$asv <- NULL
}

restotal <- as.data.frame(restotal)
restotal2 <- restotal %>% 
    mutate_at(c(3:10), as.character) %>% 
    mutate_at(c(3:10), as.numeric) %>% 
    mutate_at(c(3:5, 7:9), afronden2) %>% 
    mutate(
        `m0-q` = p.adjust(`m0-p`, 'fdr'),
        `m1-q` = p.adjust(`m1-p`, 'fdr')
    ) %>% 
    mutate_at(c(6,10), afronden3)
openxlsx::write.xlsx(restotal2, "results/glm_clr/strata_urb.xlsx")

restotal3 <- restotal2 %>% 
    pivot_longer(c(3:12), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0", "1"), 
                          labels = c("Unadjusted", "Adjusted")),
           ASV = factor(ASV, levels = colnames(df)[(ncol(df)-19):ncol(df)]),
           ASV = fct_rev(ASV), 
           Diet = factor(Diet, levels = colnames(df)[c(32,34,35,37)]),
           sigq = case_when(q < 0.05 ~ paste0("q<0.05"), q >= 0.05 ~ paste0("not sig")),
           sigq = as.factor(sigq))

(pltotaal <- ggplot(restotal3, aes(x=ASV, y=est, color=model, shape = sigq)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.7)) +
        scale_shape_manual(values = c(21,19))+
        # scale_y_continuous(breaks = c(-2:8))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.7)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "ASVs with differential abundance between geographical locations",
             x = "", y = "Difference (log10-transformed counts) for SD increase in macronutrient group", shape = "", color = "") +
        scale_color_cosmic() +
        coord_flip()+
        facet_wrap(~Diet))
ggsave("results/glm_clr/lm_urb.pdf", width = 9, height = 8)
ggsave("results/glm_clr/lm_urb.svg", width = 9, height = 8)
