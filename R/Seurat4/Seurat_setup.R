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
GSE129363$orig.ident = GSE129363$SampleName
Idents(GSE129363) ="orig.ident"
#======== GSE136229 =======
GSE136229_data = read.table(gzfile("data/GSE136229/GSE136229_PMD18_Normalised_Exp.txt.gz"),header = T)
colnames(GSE136229_data) %<>% gsub("\\.","-",.)

#======== GSE135775 =======
GSE135775_data = read.table(gzfile("data/GSE135775/GSE135775_matrixRPKM_MSC.txt.gz"),header = T)

#======== GSE135776 =======
GSE135776_data1 = read.table(gzfile("data/GSE135776/GSE135776_matrixRPKM_Adipocyte_SVF.txt.gz"),header = T)
GSE135776_data2 = read.table(gzfile("data/GSE135776/GSE135776_matrixRPKM_AdiposeTissue.txt.gz"),header = T)


# read and select mitochondial genes
mito = "^MT-"
(mito.features <- grep(pattern = mito, x = rownames(GSE129363), value = TRUE))
GSE129363[["percent.mt"]] <- PercentageFeatureSet(object = GSE129363, pattern = mito)


g1 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = GSE129363, features = features, ncol = 1, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=8,angle = 90),legend.position="none")
})
save(g1,file= paste0(path,"g1","_25_",gsub("-","",Sys.Date()),".Rda"))


jpeg(paste0(path,"S1_QC.jpeg"), units="in", width=10, height=7,res=600)
ggarrange(g1[[1]], g1[[2]],g1[[3]],ncol = 3, nrow = 1)
dev.off()

# =========  Normalization ===================
object <- GSE129363
DefaultAssay(object) = "RNA"
object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = head(VariableFeatures(object), 10),
                     repel = TRUE)

jpeg(paste0(path,"S1_VariableFeaturePlot.jpeg"), units="in", width=10, height=7,res=600)
plot2
dev.off()
#======1.7 run UMAP =========================
object %<>% ScaleData
npcs <- 100
object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))

jpeg(paste0(path,"S1_ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = npcs)
dev.off()

# ========== Determine the 'dimensionality' of the dataset ==============
object %<>% JackStraw(num.replicate = 20,dims = npcs)
object %<>% ScoreJackStraw(dims = 1:npcs)
a <- seq(1,100, by = 10)
b <- a+9
save.path = paste0(path,"JackStrawPlots/")
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

for(i in seq_along(a)){
        jpeg(paste0(save.path,"JackStrawPlot_",i,"_", a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
        print(JackStrawPlot(object, dims = a[i]:min(b[i],100)))
        Progress(i,length(a))
        dev.off()
}

npcs = 23 # 50
object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))

object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.6)

tmp = object@reductions$tsne@cell.embeddings[,"tSNE_2"]
object@reductions$tsne@cell.embeddings[,"tSNE_2"] = object@reductions$tsne@cell.embeddings[,"tSNE_1"]
object@reductions$tsne@cell.embeddings[,"tSNE_1"] = tmp
object@reductions$tsne@cell.embeddings = -object@reductions$tsne@cell.embeddings

TSNEPlot.1(object, group.by="orig.ident",pt.size = 0.5,label = T,
           label.repel = T,alpha = 0.9, cols = Singler.colors,
           no.legend = T,label.size = 4, repel = T, title = "No Integration with 25 samples with 23 PCA",
           do.print = T, do.return = F)
TSNEPlot.1(object, pt.size = 1,label = T,
           label.repel = F,alpha = 1, 
           no.legend = T,label.size = 8, repel = T, 
           do.print = T, do.return = F)
saveRDS(object, file = "data/Adipose_25_20201027.rds")

# serial resolution and generate seurat
object = readRDS(file = "data/Adipose_25_20201027.rds")
DefaultAssay(object) = "RNA"
resolutions = c(seq(0.001,0.009, by = 0.001),seq(0.01,0.09, by = 0.01),seq(0.1,5, by = 0.1))
save.path = paste0(path,"serial_resolutions/")
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

for(i in 1:length(resolutions)){
        object %<>% FindClusters(resolution = resolutions[i])
        Idents(object) = paste0("RNA_snn_res.",resolutions[i])
        TSNEPlot.1(object, group.by=paste0("RNA_snn_res.",resolutions[i]),pt.size = 0.3,label = T,
                   label.repel = T,alpha = 0.9,
                   do.return = F,
                   no.legend = T,label.size = 4, repel = T, 
                   title = paste("res =",resolutions[i]),
                   do.print = T, save.path = save.path)
        Progress(i,length(resolutions))
}


save.path = paste0(path,"FeaturePlots/")
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)


features <- FilterGenes(object,c("CD3D","CD4","CD8A",
                                 "MS4A7","CD14","FCGR1A",
                                 "GNLY","CD19","CD38"))

features <- FilterGenes(object,c("CD34","CFD","KRT18",
                                 "PECAM1","VWF","CLDN5",
                                 "PTPRC","SRGN","LAPTM5"))
FeaturePlot.1(object,features = features, pt.size = 0.005, cols = c("gray", "red"),
              alpha = 1,reduction = "tsne",
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 3,
              units = "in",width=9, height=9, no.legend = T, save.path = save.path)
