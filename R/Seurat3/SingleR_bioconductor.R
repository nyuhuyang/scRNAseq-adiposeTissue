#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0
library(SingleR)
library(SingleCellExperiment)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

# ====== load single cell =============
object = readRDS(file = "data/Adipose_25_20201027.rds")

sce <- SingleCellExperiment(list(logcounts=object[["RNA"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# ====== load reference =============
blue_encode <- BlueprintEncodeData()
remove = grepl("CD4|CD8|Tregs|B-cells",blue_encode$label.fine)
blue_encode = blue_encode[,!remove]

immue_exp <- DatabaseImmuneCellExpressionData()

common <- Reduce(intersect, list(rownames(sce),
                                 rownames(blue_encode),
                                 rownames(immue_exp)))
length(common)
combine_ref = do.call("cbind", list(blue_encode[common,],
                                    immue_exp[common,]))
table(combine_ref$label.fine)
system.time(trained <- trainSingleR(ref = combine_ref,
                                    labels=combine_ref$label.fine))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/Adipose_25_20201027_pred.rds")
