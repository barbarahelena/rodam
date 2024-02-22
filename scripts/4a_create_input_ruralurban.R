# create XGB input files for aldosterone and renin models

library(dplyr)
library(phyloseq)
rm(list=ls())

# make data for machine learning XGB classification models

# writes input data files for XGB models as tab-delimited 
# subject ids and feature ids are written as separate tab-delimited files
# write X data / predictors
write_data <- function(x, data_path){
    x <- as.matrix(x)
    if(any(is.na(x))){
        cat('There are missing values in the input data!\n')
    }
    write.table(x, file.path(data_path, 'X_data.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(colnames(x), file.path(data_path,'feat_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(rownames(x), file.path(data_path,'subject_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
}

# write y / predicted outcome
write_y <- function(x, name_y, data_path){
    if(missing(name_y)){
        cat('\n\nYou need to provide a name for the y data file!\n')
    }
    if(!name_y %in% c('y_binary.txt', 'y_reg.txt')){
        cat('\nThe file name is not compatible with XGBeast!\n' )
    }
    if(any(is.na(x))){
        cat('\nThere are missing values in the outcome data!\n')
    }
    write.table(x, file = file.path(data_path, name_y), row.names = F, col.names = F, sep = '\t', quote = F)
}


## Open ADC dataframe
df <- readRDS('data/clinicaldata.RDS')
head(df)
any(is.na(df$Site)) # FALSE
df$Site <- case_when(
                df$Site=="Rural Ghana" ~ 0,
                df$Site=="Urban Ghana" ~ 1)
summary(as.factor(df$Site))
df <- df %>% dplyr::select(ID, Site)

## Open RDS file with OTU table
mb <- readRDS('data/phyloseq_sampledata.RDS')
otu <- t(as(mb@otu_table, "matrix"))
tk <- apply(otu, 2, function(x) sum(x > 10) > (0.30*length(x)))
mbdf <- otu[,tk]
# dim(mbdf)
# gghistogram(mbdf[,6], bins = 30)
mbdf <- as.data.frame(mbdf)
dim(mbdf)

# rownames(df) <- df$ID
# mbdf$ID <- rownames(mbdf)
# dfmb <- left_join(mbdf, df, by = "ID")
# rowSums(otu)
# dfmb %>% group_by(Site) %>% summarise(number = median(ASV_220))
# ggplot(data = dfmb, aes(x = as.factor(Site), y = ASV_220, fill = as.factor(Site))) +
#     geom_violin() +
#     geom_boxplot(fill = "white", width = 0.1) +
#     scale_y_log10()+
#     scale_fill_manual(values = pal_futurama()(4)[3:4]) +
#     stat_compare_means()+
#     theme_Publication()

# Put clinical data and microbiome data in same sequence of IDs
df <- df[match(rownames(mbdf), df$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(df$ID == rownames(mbdf)) # TRUE
df$ID
rownames(mbdf)

# make input data
path <- 'rural_urban'
dir.create(path)
dir.create("rural_urban/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(df$Site)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

