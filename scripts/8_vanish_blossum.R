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
df <- left_join(clin, dfmb, by = "ID") 
group <- df %>% dplyr::select(Site, 69:ncol(.)) %>% 
    pivot_longer(., 2:ncol(.), names_to = "ASV", values_to = "abundance") %>% 
    group_by(Site, ASV) %>% summarise(mean = mean(abundance)) %>% 
    pivot_wider(., id_cols = "ASV", names_from = "Site", values_from = "mean") %>% 
    filter(`Amsterdam` < `Rural Ghana` & `Urban Ghana` < `Rural Ghana`)
group2 <- df %>% dplyr::select(Site, 69:ncol(.)) %>% 
    pivot_longer(., 2:ncol(.), names_to = "ASV", values_to = "abundance") %>% 
    group_by(Site, ASV) %>% summarise(mean = mean(abundance)) %>% 
    pivot_wider(., id_cols = "ASV", names_from = "Site", values_from = "mean") %>% 
    filter(`Rural Ghana` < `Amsterdam` & `Urban Ghana` < `Amsterdam`)
df <- df %>%
    mutate(across(c(69:ncol(.)), as.factor)) %>% 
    mutate(across(c(69:ncol(.)), ~fct_recode(.x, "Present"="1", "Absent"="0")))
head(df)[1:5,1:5]
colnames(df)

bugs <- group$ASV
len <- length(bugs)
dfvenn <- df %>% dplyr::select(all_of(bugs)) %>% mutate(across(everything(.), 
                                                ~ case_when(.x == "Absent" ~ 0, 
                                                            .x == "Present" ~ 1)))
                                          
# upset(as.data.frame(dfvenn), nsets=13, sets.bar.color = "#56B4E9", order.by = "freq")

dfbugs <- df %>% 
    dplyr::select(ID, Site, bugs) %>% 
    mutate(
        across(bugs, ~ case_when(.x == "Absent" ~ 1,
                                          .x == "Present" ~ 0)),
        across(bugs, as.logical),
        sum = rowSums(across(bugs)),
        vanishing = case_when(
            sum < 10 ~ FALSE,
            .default = TRUE
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
               vanishing == TRUE ~ "vanishing",
               vanishing == FALSE ~ "controls"
           ),
           vanishing = fct_relevel(as.factor(vanishing), "vanishing", after=1L))
names(dfclin)
dfclin1 <- dfclin %>% dplyr::select(Age, BMI, SBP, DBP, HbA1c, GFR, Site, vanishing)
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
        stat_compare_means(comparisons = list(c("vanishing", "controls")), 
                           label = "p.signif", 
                           hide.ns = TRUE, tip.length = 0, method = "wilcox.test") +
        scale_alpha_manual(values = c(0.7, 1.0), guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        labs(y = unit, title = varname, x = "")
    print(pl)
    plot_lin[[a]] <- pl
    dfclin1$var <- NULL
}

(plots_cont1 <- ggarrange(plotlist = plot_lin, nrow = 1, ncol = 6))

dfclin2 <- dfclin %>% dplyr::select(Obesity, Hypertension=HT, Diabetes=DM, Hyperchol, 
                             Site, vanishing) %>% 
    mutate(across(1:4, ~case_when(.x == "No" ~ 0,
                                  .x == "Yes" ~ 1)),
           across(1:4, ~ as.numeric(as.character(.x))))

# plot_fact <- list()
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
#     pl <- ggplot(data = dfsum, aes(x = vanishing, y = prev)) +
#         geom_bar(aes(alpha = vanishing), stat = "identity", fill = pal_cosmic()(6)[6]) +
#         theme_Publication() +
#         stat_pvalue_manual(pvaltab, label = "p.signif", tip.length = 0) +
#         scale_alpha_manual(values = c(0.7, 1.0), guide = "none") +
#         labs(y = "prevalence in % subjects", title = varname, x = "") +
#         scale_y_continuous(limits = c(0,65))
#     print(pl)
#     plot_fact[[a]] <- pl
    dfclin2$var <- NULL
}
# plots_factors <- ggarrange(plotlist = plot_fact, nrow = 1)

(plots_factors1 <- ggplot(data = sum, aes(x = varr, y = prev)) +
    geom_bar(aes(alpha = vanishing, group = vanishing), stat = "identity",
             fill = pal_cosmic()(6)[6], position = "dodge") +
    theme_Publication() +
    scale_alpha_manual(values = c(0.7,1.0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_signif(y_position = c(20, 35, 63, 40),
                xmin = c(0.7,1.7, 2.7,3.7), 
                xmax = c(1.3,2.3,3.3,4.3),
                annotation=c(rev(res$p.signif)), tip_length=0) +
    labs(x = "", y = "prevalence in % of subjects", alpha = "") +
    scale_y_continuous(limits = c(0, 65)))
    
(pl_vanish <- ggarrange(plots_cont1, plots_factors1, nrow = 2, heights = c(2,1)))
ggsave("results/vanishblossom/vanishing_health.pdf", width = 8, height = 12)

vanishingID <- dfclin %>% dplyr::select(ID, vanishing)
saveRDS(vanishingID, "data/ids_vanishing.RDS")


#### blossom ####
df <- df %>%
    mutate(across(c(69:ncol(.)), as.factor)) %>% 
    mutate(across(c(69:ncol(.)), ~fct_recode(.x, "Present"="1", "Absent"="0")))
head(df)[1:5,1:5]
colnames(df)

bugs2 <- group2$ASV
len <- length(bugs2)
dfvenn <- df %>% dplyr::select(bugs2) %>% mutate(across(everything(.), 
                                                ~ case_when(.x == "Absent" ~ 0, .x == "Present" ~ 1)))

# upset(as.data.frame(dfvenn), nsets=13, sets.bar.color = "#56B4E9", order.by = "freq")

dfbugs <- df %>% 
    dplyr::select(ID, Site, bugs2) %>% 
    mutate(
        across(bugs2, ~ case_when(.x == "Absent" ~ 0,
                                 .x == "Present" ~ 1)),
        across(bugs2, as.logical),
        sum = rowSums(across(bugs2)),
        blossom = case_when(
            sum < 4 ~ FALSE,
            .default = TRUE
        )
    )
summary(dfbugs$blossom)
ggplot(dfbugs, aes(x = sum, fill = blossom)) +
    geom_histogram(binwidth = 1) +
    theme_Publication() +
    ggsci::scale_fill_jama(guide = "none") +
    facet_wrap(~Site)


dfclin <- left_join(clin, dfbugs, by = c("ID", "Site")) %>% 
    mutate(CRP=log(CRP),
           Obesity = case_when(BMI >=30 ~ "Yes", BMI < 30 ~ "No"),
           blossom = case_when(
               blossom == FALSE ~ "controls",
               blossom == TRUE ~ "blossom"
           ),
           blossom = fct_relevel(as.factor(blossom), "blossom", after=1L))
names(dfclin)
dfclin1 <- dfclin %>% dplyr::select(Age, BMI, SBP, DBP, HbA1c, GFR, Site, blossom)
units <- data.frame(
    var = names(dfclin1)[1:6],
    units = c("Years", "kg/m2", "mmHg", "mmHg","mmol/mol", "ml/min")
)

plot_lin <- list()
for(a in c(1:6)) {
    dfclin1$var <- dfclin1[[a]]
    varname <- names(dfclin1)[a]
    unit <- units$unit[match(varname, units$var)]
    pl <- ggplot(data = dfclin1, aes(x = blossom, y = var)) +
        geom_violin(aes(alpha = blossom), fill = pal_cosmic()(7)[7]) +
        geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
        theme_Publication() +
        stat_compare_means(comparisons = list(c("blossom", "controls")), 
                           label = "p.signif", 
                           hide.ns = TRUE, tip.length = 0, method = "wilcox.test") +
        scale_alpha_manual(values = c(0.7, 1.0), guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        labs(y = unit, title = varname, x = "")
    print(pl)
    plot_lin[[a]] <- pl
    dfclin1$var <- NULL
}

plots_cont2 <- ggarrange(plotlist = plot_lin, nrow = 1)

dfclin2 <- dfclin %>% dplyr::select(Obesity, Hypertension=HT, Diabetes=DM, Hyperchol, 
                             Site, blossom) %>% 
    mutate(across(1:4, ~case_when(.x == "No" ~ 0,
                                  .x == "Yes" ~ 1)),
           across(1:4, ~ as.numeric(as.character(.x))))

# plot_fact <- list()
res <- c()
sum <- c()
for(a in c(1:4)) {
    dfclin2$var <- dfclin2[[a]]
    varname <- names(dfclin2)[a]
    dfsum <- dfclin2 %>% group_by(blossom) %>% reframe(
        prevalence = sum(var, na.rm = TRUE)
    )
    counts <- dfclin2 %>% group_by(blossom) %>% count(.)
    dfsum <- left_join(dfsum, counts, by = "blossom") %>% mutate(prev = prevalence / n * 100,
                                                                   varr = varname)
    sum <- rbind(sum, dfsum)
    print(dfsum)
    pval <- chisq.test(dfsum$prevalence)
    pvaltab <- data.frame(varr = varname, .y. = "prev",
                          group1 = levels(dfclin2$blossom)[1], group2 = levels(dfclin2$blossom)[2],
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
#     pl <- ggplot(data = dfsum, aes(x = blossom, y = prev)) +
#         geom_bar(aes(alpha = blossom), stat = "identity", fill = pal_cosmic()(7)[7]) +
#         theme_Publication() +
#         stat_pvalue_manual(pvaltab, label = "p.signif", tip.length = 0) +
#         scale_alpha_manual(values = c(0.7, 1.0), guide = "none") +
#         labs(y = "prevalence in % subjects", title = varname, x = "") +
#         scale_y_continuous(limits = c(0,65))
#     print(pl)
#     plot_fact[[a]] <- pl
    dfclin2$var <- NULL
}
# plots_factors <- ggarrange(plotlist = plot_fact, nrow = 1)

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

(plots_factors2 <- ggplot(data = sum, aes(x = varr, y = prev)) +
                    geom_bar(aes(alpha = blossom, group = blossom), stat = "identity",
                             fill = pal_cosmic()(7)[7], position = "dodge") +
                    theme_Publication() +
                    scale_alpha_manual(values = c(0.7,1.0)) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    geom_signif(y_position = c(35, 60, 40),
                                xmin = c(1.7, 2.7,3.7), 
                                xmax = c(2.3,3.3,4.3),
                                annotation=c(rev(res$p.signif[1:3])), tip_length=0) +
                    labs(x = "", y = "prevalence in % of subjects", alpha = "") +
                    scale_y_continuous(limits = c(0, 65)))

# (pl_blossom <- ggarrange(plots_cont, ggarrange(plots_factors, NULL), nrow = 2, heights = c(1,1)))
# ggsave("results/blossom_health.pdf", width = 10, height = 10)

blossomID <- dfclin %>% dplyr::select(ID, blossom)
saveRDS(blossomID, "data/ids_blossom.RDS")

ggarrange(plots_cont1, plots_cont2, ggarrange(plots_factors1, plots_factors2),
          nrow = 3, labels = c("A", "", "B"), heights = c(1,1,1.3))

ggsave("results/vanishblossom/vanish_blossom.pdf", width = 13, height = 16)
ggsave("results/vanishblossom/vanish_blossom.svg", width = 13, height = 16)

