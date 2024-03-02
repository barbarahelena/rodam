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
df <- readRDS('data/clinicaldata.RDS') %>% filter(Site %in% c("Urban Ghana", "Amsterdam"))
head(df)
any(is.na(df$Site)) # FALSE
df$Site <- case_when(
                df$Site=="Urban Ghana" ~ 0,
                df$Site=="Amsterdam" ~ 1)
summary(as.factor(df$Site))
df <- df %>% dplyr::select(ID, Site)

## Open RDS file with OTU table
mb <- readRDS('data/phyloseq_sampledata.RDS')
mb <- prune_samples(sample_names(mb) %in% df$ID, mb)
otu <- t(as(mb@otu_table, "matrix"))
otu_dich <- ifelse(otu > 0, 1, 0)
tk <- apply(otu_dich, 2, function(x) sum(x > 0) > (0.05*length(x)))
mbdf <- otu_dich[,tk]
dim(mbdf)
mbdf <- as.data.frame(mbdf)
dim(mbdf)

# Put clinical data and microbiome data in same sequence of IDs
df <- df[match(rownames(mbdf), df$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(df$ID == rownames(mbdf)) # TRUE
df$ID
rownames(mbdf)

# make input data
path <- 'urban_ams_dich'
dir.create(path)
dir.create("urban_ams_dich/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(df$Site)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))



## Open RDS file with OTU table
mb <- readRDS('data/phyloseq_sampledata.RDS')
mb <- prune_samples(sample_names(mb) %in% df$ID, mb)
mb <- tax_glom(mb, "Genus")
tax <- as(mb@tax_table, "matrix")
saveRDS(tax, "data/taxtable_genus.RDS")
otu <- t(as(mb@otu_table, "matrix"))
otu_dich <- ifelse(otu > 0, 1, 0)
# tk <- apply(otu_dich, 2, function(x) sum(x > 0) > (0.05*length(x)))
mbdf <- otu_dich
dim(mbdf)
mbdf <- as.data.frame(mbdf)
dim(mbdf)

# Put clinical data and microbiome data in same sequence of IDs
df <- df[match(rownames(mbdf), df$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(df$ID == rownames(mbdf)) # TRUE
df$ID
rownames(mbdf)

# make input data
path <- 'urban_ams_dich_genus'
dir.create(path)
dir.create("urban_ams_dich_genus/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(df$Site)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))



## Open RDS file with OTU table
mb <- readRDS('data/phyloseq_sampledata.RDS')
mb <- prune_samples(sample_names(mb) %in% df$ID, mb)
mb <- tax_glom(mb, "Species")
tax <- as(mb@tax_table, "matrix")
saveRDS(tax, "data/taxtable_species.RDS")
otu <- t(as(mb@otu_table, "matrix"))
otu_dich <- ifelse(otu > 0, 1, 0)
# tk <- apply(otu_dich, 2, function(x) sum(x > 0) > (0.05*length(x)))
mbdf <- otu_dich
dim(mbdf)
mbdf <- as.data.frame(mbdf)
dim(mbdf)

# Put clinical data and microbiome data in same sequence of IDs
df <- df[match(rownames(mbdf), df$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(df$ID == rownames(mbdf)) # TRUE
df$ID
rownames(mbdf)

# make input data
path <- 'urban_ams_dich_species'
dir.create(path)
dir.create("urban_ams_dich_species/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(df$Site)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))
