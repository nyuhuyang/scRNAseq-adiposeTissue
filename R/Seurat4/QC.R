########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","ggplot2","GEOquery","Matrix"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
        }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# GSE129363 Single Cell RNASeq profiling of stromal vascular fraction from Subcutaneous and visceral adipose tissue
getGEOSuppFiles("GSE129363",baseDir = "data") 
getGEOSuppFiles("GSE135775",baseDir = "data") 
getGEOSuppFiles("GSE135776",baseDir = "data") 
getGEOSuppFiles("GSE136229",baseDir = "data") 



#======== GSE129363 =======
GSE129363_counts <- Read10X.1(data.dir = "data/GSE129363/",
                       barcodes.fileName = 'GSE129363_Discovery_Cohort_barcodes.tsv',
                       features.fileName = 'GSE129363_Discovery_Cohort_genes.tsv.gz',
                       matrix.fileName = 'GSE129363_Discovery_Cohort_matrix.mtx')
gsub(".*-","",colnames(GSE129363_counts)) %>% as.character %>% as.integer %>% sort %>% table

GSE129363_data = read.table(gzfile("data/GSE129363/GSE129363_SVF_Normalised_Data.txt.gz"),header = T)
colnames(GSE129363_data) %<>% gsub("\\.","-",.)
data <- as.matrix(GSE129363_data) %>% as("dgCMatrix")
GSE129363_meta.data = read.table(gzfile("data/GSE129363/GSE129363_Discovery_Cohort_CellAnnotation.txt.gz"),
                            header = T)
GSE129363_meta.data$SampleName %<>% gsub("Ate-","AT-",.)
CellID = GSE129363_meta.data$CellID
rownames(GSE129363_meta.data) = CellID
GSE129363_counts = GSE129363_counts[,CellID]
GSE129363_data = GSE129363_data[,CellID]

table(CellID == colnames(GSE129363_data))
table(CellID == colnames(GSE129363_counts))

GSE129363_meta.data = GSE129363_meta.data[colnames(GSE129363),]
GSE129363 <- CreateSeuratObject(counts = GSE129363_counts,
                                meta.data = GSE129363_meta.data,
                                min.cells = 0,
                                min.features = 0)
GSE129363[["RNA"]]@data = data
#======== GSE136229 =======
GSE136229_data = read.table(gzfile("data/GSE136229/GSE136229_PMD18_Normalised_Exp.txt.gz"),header = T)
colnames(GSE136229_data) %<>% gsub("\\.","-",.)


# read and select mitochondial genes
mito = "^MT-"
(mito.features <- grep(pattern = mito, x = rownames(GSE129363), value = TRUE))
GSE129363[["percent.mt"]] <- PercentageFeatureSet(object = GSE129363, pattern = mito)

g1 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = GSE129363, features = features, ncol = 3, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=12),legend.position="none")
})
save(g1,file= paste0(path,"g1","_25_",gsub("-","",Sys.Date()),".Rda"))

