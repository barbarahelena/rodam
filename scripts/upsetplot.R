## Process results machine learning

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(stringr)
library(forcats)
library(ggvenn)
library(UpSetR)

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

names(df)[56:74]
bugs <- c("Weissella spp.", "Succinivibrio spp.", "Prevotella_9 spp.", "Prevotella_9 spp..1",
          "Streptococcus spp.", "Lachnospiraceae spp.", "Lachnospiraceae spp..1", 
          "Lachnospiraceae UCG-004 spp.", "Lachnospiraceae spp..2", "Lactobacillales spp.",
          "Cellulosilyticum spp.", "Clostridium sensu stricto 1 spp.", "UCG-005 spp.")
dfvenn <- df %>% select(bugs) %>% mutate(across(everything(.), 
                                                ~ case_when(.x == "Absent" ~ 1, .x == "Present" ~ 0)))
                                          
upset(as.data.frame(dfvenn), nsets=19, sets.bar.color = "#56B4E9", order.by = "freq",
      empty.intersections = TRUE)

dfbugs <- df %>% 
    select(ID, Site, bugs) %>% 
    mutate(
        across(bugs, ~ case_when(.x == "Absent" ~ 1,
                                          .x == "Present" ~ 0)),
        across(bugs, as.logical),
        sum = rowSums(across(bugs)),
        vanishing = case_when(
            sum >= 10 ~ TRUE,
            .default = FALSE
        )
    )
summary(dfbugs$vanishing)
ggplot(dfbugs, aes(x = sum, fill = vanishing)) +
    geom_histogram(binwidth = 1) +
    theme_Publication() +
    ggsci::scale_fill_jama(guide = "none") +
    facet_wrap(~Site)
    

dfclin <- left_join(clin, dfbugs, by = c("ID", "Site")) %>% 
    mutate(CRP=log(CRP),
           Obesity = case_when(BMI >=30 ~ "Yes", BMI < 30 ~ "No"),
           vanishing = case_when(
               vanishing == TRUE ~ "absent",
               vanishing == FALSE ~ "present"
           ),
           vanishing = fct_relevel(as.factor(vanishing), "absent", after=1L))
names(dfclin)
dfclin1 <- dfclin %>% select(Age, BMI, SBP, DBP, HbA1c, GFR, Site, vanishing)
units <- data.frame(
    var = names(dfclin1)[1:6],
    units = c("Years", "kg/m2", "mmHg", "mmHg","mmol/mol", "ml/min")
)

plot_lin <- list()
for(a in c(1:6)) {
    dfclin1$var <- dfclin1[[a]]
    varname <- names(dfclin1)[a]
    unit <- units$unit[match(varname, units$var)]
    pl <- ggplot(data = dfclin1, aes(x = vanishing, y = var)) +
        geom_violin(aes(alpha = vanishing), fill = pal_cosmic()(6)[6]) +
        geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
        theme_Publication() +
        stat_compare_means(comparisons = list(c("present", "absent")), label = "p.signif", 
                           hide.ns = TRUE, tip.length = 0, method = "wilcox.test") +
        scale_alpha_manual(values = c(0.7, 1.0), guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs(y = unit, title = varname, x = "")
    print(pl)
    plot_lin[[a]] <- pl
    dfclin1$var <- NULL
}

plots_cont <- ggarrange(plotlist = plot_lin, nrow = 2, ncol = 3)

dfclin2 <- dfclin %>% select(Obesity, Hypertension=HT, Diabetes=DM, Hyperchol, 
                             Site, vanishing) %>% 
    mutate(across(1:4, ~case_when(.x == "No" ~ 0,
                                  .x == "Yes" ~ 1)),
           across(1:4, ~ as.numeric(as.character(.x))))

plot_fact <- list()
res <- c()
sum <- c()
for(a in c(1:4)) {
    dfclin2$var <- dfclin2[[a]]
    varname <- names(dfclin2)[a]
    dfsum <- dfclin2 %>% group_by(vanishing) %>% reframe(
        prevalence = sum(var, na.rm = TRUE)
        )
    counts <- dfclin2 %>% group_by(vanishing) %>% count(.)
    dfsum <- left_join(dfsum, counts, by = "vanishing") %>% mutate(prev = prevalence / n * 100,
                                                                   varr = varname)
    sum <- rbind(sum, dfsum)
    print(dfsum)
    pval <- chisq.test(dfsum$prevalence)
    pvaltab <- data.frame(varr = varname, .y. = "prev", 
                          group1 = levels(dfclin2$vanishing)[1], group2 = levels(dfclin2$vanishing)[2],
                          p.value = pval$p.value) %>% 
        mutate(p.signif = case_when(
            p.value < 0.0001 ~ "****",
            p.value < 0.001 ~ "***",
            p.value < 0.01 ~ "**",
            p.value < 0.05 ~ "*",
            .default=""
        ),
        y.position = max(dfsum$prev + 4))
    print(pvaltab)
    res <- rbind(res, pvaltab)
    pl <- ggplot(data = dfsum, aes(x = vanishing, y = prev)) +
        geom_bar(aes(alpha = vanishing), stat = "identity", fill = pal_cosmic()(6)[6]) +
        theme_Publication() +
        stat_pvalue_manual(pvaltab, label = "p.signif", tip.length = 0) +
        scale_alpha_manual(values = c(0.7, 1.0), guide = "none") +
        labs(y = "prevalence in % subjects", title = varname, x = "") +
        scale_y_continuous(limits = c(0,65))
    print(pl)
    plot_fact[[a]] <- pl
    dfclin2$var <- NULL
}
plots_factors <- ggarrange(plotlist = plot_fact, nrow = 1)

# res1 <- res %>% mutate(
#     group1 = str_c(group1, ".", varr),
#     group2 = str_c(group2, ".", varr)
# )

# ggplot(data = sum, aes(x = interaction(vanishing, varr), y = prev, group = varr)) +
#     geom_bar(aes(alpha = vanishing), stat = "identity",
#              fill = pal_cosmic()(6)[6], position = "dodge") +
#     theme_Publication() +
#     scale_alpha_manual(values = c(0.7,1.0), guide = "none") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(x = "", y = "prevalence in % of subjects") +
#     stat_pvalue_manual(data = res1,
#                        label = "p.signif", tip.length = 0)

ggplot(data = sum, aes(x = varr, y = prev, group = vanishing)) +
    geom_bar(aes(alpha = vanishing), stat = "identity",
             fill = pal_cosmic()(6)[6], position = "dodge") +
    theme_Publication() +
    scale_alpha_manual(values = c(0.7,1.0), guide = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "prevalence in % of subjects") +
    scale_y_continuous(limits = c(0, 65))
    
ggarrange(plots_cont, plots_factors, nrow = 2, heights = c(2,1))
ggsave("results/prevalence_health.pdf", width = 8, height = 12)

vanishingID <- dfclin %>% filter(vanishing == "absent") %>% select(ID, vanishing)
saveRDS(vanishingID, "data/ids_vanishing.RDS")
