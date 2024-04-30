## Explained variance
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(gghalves)
library(ggsci)
library(ggbeeswarm2)
library(tidyverse)
library(ggpubr)
theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                #family = 'Helvetica'
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.line.y = element_line(colour="black"),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line.x = element_line(colour="black"),
                axis.ticks.x = element_line(),
                axis.ticks.y = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(5,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(face = "italic", size=rel(0.6))
        ))
} 

# Making table
ev_list <- list()
groups <- c("kcal", "fibre", "fat", "carbohydrates", "proteins", "sodium",
            "age", "bmi", "physicalactivity", "bristol",
            "ldl", "gfr", "crp")
for(g in groups){
    li <- list.files(path = file.path("covariatemodels", g))
    a <- str_detect(li, regex(g, ignore_case = T))
    b <- !(str_detect(li, 'PERMUTED'))
    ev_list[[g]] <- rio::import(file.path("covariatemodels", g, li[which(a&b)], 'model_results_per_iteration.txt'))
}

df <- data.frame()
for (i in c(1:13)) {
    group <- str_split(names(ev_list[i]), pattern = '_', 2, simplify = T)[,1]
    mean_ev <- mean(ev_list[[i]]$`Explained Variance`)
    sd_ev <- sd(ev_list[[i]]$`Explained Variance`)
    row <- cbind(group, mean_ev, sd_ev)
    df<- rbind(df, row)
}

for(g in groups){
    li <- list.files(path = file.path("covariatemodels", g))
    a <- str_detect(li, regex(g, ignore_case = T))
    b <- !(str_detect(li, 'PERMUTED'))
    ev_list[[g]] <- rio::import(file.path("covariatemodels", g, li[which(a&b)], 'model_results_per_iteration.txt'))
}

df3 <- df %>% 
    mutate(
        mean_ev = as.numeric(mean_ev),
        sd_ev = as.numeric(sd_ev),
        group = as.factor(group)
    )
gem <- df3 %>% mutate(
    mean_ev = format(round(mean_ev*100, 2), nsmall = 2)
)
gem <- gem %>% arrange(desc(as.numeric(mean_ev)))

for(a in c(1:13)){
    ev_list[[a]]$outcome <- names(ev_list)[a]
}

compl <- bind_rows(ev_list)

df4 <- compl %>% select(`Explained Variance`, outcome) %>% 
    mutate(
        expvar = as.numeric(`Explained Variance`),
        outcome = as.factor(outcome)
    )

dfmean <- df4 %>% group_by(outcome) %>% summarise(mean = mean(expvar)*100) %>% 
    arrange(desc(mean))
df4 <- df4 %>% mutate(outcome2 = fct_reorder(outcome, .x = expvar,
                                             .fun = median, .desc = TRUE),
                      outcome2 = fct_recode(outcome2, "carbs"="carbohydrates"),
                    outcome2 = fct_recode(outcome2, "activity"="physicalactivity")
)

(pl3 <- df4 %>% filter(outcome %in% c("sodium", "proteins", "fat", "carbohydrates", "kcal", "fibre")) %>% 
        ggplot(., aes(x = outcome2, fill = outcome2, 
                        y = expvar*100, groups = outcome2))+
        annotate("rect", xmin=-Inf, xmax=Inf,
                 ymin=-10, ymax=0, alpha=0.5, fill="grey")+
        geom_hline(yintercept = 0, linetype = "dashed")+
        geom_half_violin(side = "r", scale = "width", alpha = 0.75, position = position_dodge(0.75)) +
        geom_half_boxplot(side = "r", fill = "white",
                     outlier.shape = NA, width = 0.2, position = position_dodge(0.75)) +
        geom_beeswarm(aes(color = outcome2), side = -1L, method = "swarm", alpha = 0.5, spacing = 0.8,
                      corral = "none", corral.width = 0.4) +
        labs(title = 'Macronutrients',
             x = '', y = 'Explained variance (%)', linetype = '', size = '') +
        scale_y_continuous(limits = c(-10, 30), n.breaks = 8)+
        scale_color_manual(guide = "none", values = c(pal_cosmic()(4)[2:4], rep("grey50",3))) +
        scale_fill_manual(guide = "none", values = c(pal_cosmic()(4)[2:4], rep("grey50",3))) +
        theme_Publication())

gem2 <- gem %>% filter(group %in% c("sodium", "proteins", "fat", 
                                    "carbohydrates", "kcal", "fibre")) %>% 
    mutate(mean_ev = as.numeric(mean_ev),
            mean_ev = case_when(mean_ev < 0 ~ NA, 
                                .default = mean_ev))
df5 <- df4 %>% filter(outcome2 %in% c("sodium", "proteins", "fat", 
                                      "carbs", "kcal", "fibre")) %>% droplevels(.)

for(i in 1:3){
    pl3 <- pl3 + annotation_custom(
        grob = textGrob(label = str_c(as.character(gem2$mean_ev[i]),'%'),
                        hjust = 0, gp = gpar(cex = 0.8)),
        xmin = levels(df5$outcome2)[i],      # Vertical position of the textGrob
        xmax = levels(df5$outcome2)[i],
        ymin = max(df5$expvar[which(levels(df5$outcome2)[i] == df5$outcome2)])*100 + 2,
        ymax = max(df5$expvar[which(levels(df5$outcome2)[i] == df5$outcome2)])*100 + 2
        )
}
pl3
# ggsave("results/240310_expvar_diet_all.pdf", width = 7, height = 5)
# ggsave("results/240310_expvar_diet_all.svg", width = 7, height = 5)


(pl4 <- df4 %>% filter(!outcome %in% c("sodium", "proteins", "fat", 
                                       "carbohydrates", "kcal", "fibre")) %>% 
        ggplot(., aes(x = outcome2, fill = outcome2, 
                      y = expvar*100, groups = outcome2))+
        annotate("rect", xmin=-Inf, xmax=Inf,
                 ymin=-10, ymax=0, alpha=0.5, fill="grey")+
        geom_hline(yintercept = 0, linetype = "dashed")+
        geom_half_violin(side = "r", scale = "width", alpha = 0.75, position = position_dodge(0.75)) +
        geom_half_boxplot(side = "r", fill = "white",
                          outlier.shape = NA, width = 0.2, position = position_dodge(0.75)) +
        geom_beeswarm(aes(color = outcome2), side = -1L, method = "swarm", alpha = 0.5, spacing = 0.8,
                      corral = "none", corral.width = 0.4) +
        labs(title = 'Host factors - continuous',
             x = '', y = 'Explained variance (%)', linetype = '', size = '') +
        scale_y_continuous(limits = c(-10, 30), n.breaks = 8)+
        scale_color_manual(guide = "none", values = c(pal_cosmic()(7)[5:7], rep("grey50",4))) +
        scale_fill_manual(guide = "none", values = c(pal_cosmic()(7)[5:7], rep("grey50",4))) +
        theme_Publication())

gem2 <- gem %>% filter(!group %in% c("sodium", "proteins", "fat", 
                                    "carbohydrates", "kcal", "fibre")) %>% 
    mutate(mean_ev = as.numeric(mean_ev),
           mean_ev = case_when(mean_ev < 0 ~ NA, 
                               .default = mean_ev))
df5 <- df4 %>% filter(!outcome2 %in% c("sodium", "proteins", "fat", 
                                      "carbs", "kcal", "fibre")) %>% droplevels(.)

for(i in 1:3){
    pl4 <- pl4 + annotation_custom(
        grob = textGrob(label = str_c(as.character(gem2$mean_ev[i]),'%'),
                        hjust = 0, gp = gpar(cex = 0.8)),
        xmin = levels(df5$outcome2)[i],      # Vertical position of the textGrob
        xmax = levels(df5$outcome2)[i],
        ymin = max(df5$expvar[which(levels(df5$outcome2)[i] == df5$outcome2)])*100 + 2,
        ymax = max(df5$expvar[which(levels(df5$outcome2)[i] == df5$outcome2)])*100 + 2
    )
}
pl4

ggarrange(pl3, pl4, labels = c("A", "B"), nrow = 1)
# ggsave("results/230315_expvar_allfactors.pdf", width = 12, height = 6)
# ggsave("results/230315_expvar_allfactors.svg", width = 12, height = 6)


## Classification models
# Making table
ev_list <- list()
groups <- c("women", "diabetes", "hyperchol", "hypertension", "bpmed", "probiotics",
            "occupationmanual", "active")
for(g in groups){
    li <- list.files(path = file.path("covariatemodels", g))
    a <- str_detect(li, regex(g, ignore_case = T))
    b <- !(str_detect(li, 'PERMUTED'))
    ev_list[[g]] <- rio::import(file.path("covariatemodels", g, li[which(a&b)], 'model_results_per_iteration.txt'))
}

df <- data.frame()
for (i in c(1:8)) {
    group <- str_split(names(ev_list[i]), pattern = '_', 2, simplify = T)[,1]
    mean_auc <- mean(ev_list[[i]]$ROC_AUC_scores)
    sd_auc <- sd(ev_list[[i]]$ROC_AUC_scores)
    row <- cbind(group, mean_auc, sd_auc)
    df<- rbind(df, row)
}

for(g in groups){
    li <- list.files(path = file.path("covariatemodels", g))
    a <- str_detect(li, regex(g, ignore_case = T))
    b <- !(str_detect(li, 'PERMUTED'))
    ev_list[[g]] <- rio::import(file.path("covariatemodels", g, li[which(a&b)], 'model_results_per_iteration.txt'))
}

df3 <- df %>% 
    mutate(
        mean_auc = as.numeric(mean_auc),
        sd_auc = as.numeric(sd_auc),
        group = as.factor(group)
    )
gem <- df3 %>% mutate(
    mean_auc = format(round(mean_auc, 2), nsmall = 2)
)
gem <- gem %>% arrange(desc(as.numeric(mean_auc)))

for(a in c(1:8)){
    ev_list[[a]]$outcome <- names(ev_list)[a]
}

compl <- bind_rows(ev_list)

df4 <- compl %>% select(ROC_AUC_scores, outcome) %>% 
    mutate(
        auc = as.numeric(ROC_AUC_scores),
        outcome = as.factor(outcome)
    )

dfmean <- df4 %>% group_by(outcome) %>% summarise(mean = mean(auc)*100) %>% 
    arrange(desc(mean))
df4 <- df4 %>% mutate(outcome2 = fct_reorder(outcome, .x = auc,
                                             .fun = median, .desc = TRUE),
                      outcome2 = fct_recode(outcome2, "sex" = "women")
)

(pl5 <- df4 %>% filter(!outcome2 %in% c("active", "diabetes", "hyperchol")) %>% 
        ggplot(., aes(x = outcome2, fill = outcome2, 
                      y = auc, groups = outcome2))+
        annotate("rect", xmin=-Inf, xmax=Inf,
                 ymin=0.40, ymax=0.50, alpha=0.5, fill="grey")+
        geom_hline(yintercept = 0.50, linetype = "dashed")+
        geom_half_violin(side = "r", scale = "width", alpha = 0.75, position = position_dodge(0.75)) +
        geom_half_boxplot(side = "r", fill = "white",
                          outlier.shape = NA, width = 0.2, position = position_dodge(0.75)) +
        geom_beeswarm(aes(color = outcome2), side = -1L, method = "swarm", 
                      alpha = 0.5, spacing = 0.8,
                      corral = "none", corral.width = 0.4) +
        labs(title = 'Host factors - binary',
             x = '', y = 'Area under the curve (AUC)', linetype = '', size = '') +
        scale_y_continuous(limits = c(0.40, 0.65))+
        scale_color_manual(guide = "none", values = c(pal_cosmic()(10)[c(9:10)],
                                                      "#CD7467", "#84A29C", 
                                                      "grey50")) +
        scale_fill_manual(guide = "none", values = c(pal_cosmic()(10)[c(9:10)], 
                                                     "#CD7467", "#84A29C", 
                                                     "grey50")) +
        theme_Publication())

gem2 <- gem %>% filter(!group %in% c("diabetes", "bpmed", "active", "hyperchol")) %>% 
    mutate(mean_auc = as.numeric(mean_auc),
           mean_auc = case_when(mean_auc < 0.51 ~ NA, 
                               .default = mean_auc))
df5 <- df4 %>% filter(!outcome2 %in% c("diabetes", "active", "hyperchol")) %>% droplevels(.)

for(i in 1:4){
    pl5 <- pl5 + annotation_custom(
        grob = textGrob(label = str_c(as.character(gem2$mean_auc[i])),
                        hjust = 0, gp = gpar(cex = 0.8)),
        xmin = levels(df5$outcome2)[i],      # Vertical position of the textGrob
        xmax = levels(df5$outcome2)[i],
        ymin = max(df5$auc[which(levels(df5$outcome2)[i] == df5$outcome2)]) + 0.02,
        ymax = max(df5$auc[which(levels(df5$outcome2)[i] == df5$outcome2)]) + 0.02
    )
}
pl5

ggarrange(pl3, pl4, pl5, labels = c("A", "B", "C"), nrow = 2, ncol = 2)
ggsave("results/covariatemodels/expvar_allfactors.pdf", width = 12, height = 12)
ggsave("results/covariatemodels/expvar_allfactors.svg", width = 12, height = 12)
