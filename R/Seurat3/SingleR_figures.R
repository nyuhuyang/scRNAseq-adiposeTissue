# conda activate r4.0
library(Seurat)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
pred = readRDS("output/Adipose_25_20201027_pred.rds")
object = readRDS(file = "data/Adipose_25_20201027.rds")

singlerDF = data.frame("label.fine" = pred$pruned.labels,
                       row.names = rownames(object@meta.data))
singler_mt <- read.csv("../SingleR/data/BlueprintEncodeImmuneCell.csv",
                       stringsAsFactors = F)
singler_mt = singler_mt[singler_mt$label.fine %in% singlerDF$label.fine,]
singlerDF$label.main = plyr::mapvalues(singlerDF$label.fine,
                                          from = singler_mt$label.fine,
                                          to =singler_mt$label.main)
table(singlerDF$label.main)
##############################
# process color scheme
##############################
object <- AddMetaData(object = object,metadata = singlerDF["label.main"])
object <- AddMetaColor(object = object, label= "label.main", colors = Singler.colors)

object <- AddMetaData(object = object,metadata = singlerDF["label.fine"])
object <- AddMetaColor(object = object, label= "label.fine", colors = Singler.colors)


lapply(c("label.fine","label.main"), function(group.by)
    TSNEPlot.1(object = object, label = T, label.repel = T,group.by = group.by,
        cols =  ExtractMetaColor(object,group.by),
        no.legend = T,
        pt.size = 1,label.size = ifelse(group.by=="label.fine",3,5),
        do.print = T,do.return = F,width=7, height=7,
        title = paste(group.by,"by Blueprint + Encode + Monaco Immune")))
