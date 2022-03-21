library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)
Oral_ST_P13Tumor=Load10X_Spatial("/fs/home/hanya/Project/Carcinogenesis/Spatial_processed/Data/P13Tumor_A1",filename="filtered_feature_bc_matrix.h5")
range(Oral_ST_P13Tumor$nFeature_Spatial)
range(Oral_ST_P13Tumor$nCount_Spatial)
Oral_ST_P13Tumor$orig.ident="P13_Tumor"
Oral_ST_P13Tumor <- PercentageFeatureSet(Oral_ST_P13Tumor, "^MT-", col.name = "percent_mito")
Oral_ST_P13Tumor = Oral_ST_P13Tumor[, Oral_ST_P13Tumor$nFeature_Spatial > 200 & Oral_ST_P13Tumor$percent_mito < 25 ]
Oral_ST_P13Tumor <- Oral_ST_P13Tumor[!grepl("^mt-", rownames(Oral_ST_P13Tumor)), ]
#saveRDS(Oral_ST_P13Tumor,file="Oral_ST_P13Tumor_seurat.rds")

Oral_ST_P13Normal=Load10X_Spatial("/fs/home/hanya/Project/Carcinogenesis/Spatial_processed/Data/P13NormalOLK_B1",filename="filtered_feature_bc_matrix.h5")
Oral_ST_P13Normal$orig.ident="P13_NormalOLK"
Oral_ST_P13Normal <- PercentageFeatureSet(Oral_ST_P13Normal, "^MT-", col.name = "percent_mito")
Oral_ST_P13Normal = Oral_ST_P13Normal[, Oral_ST_P13Normal$nFeature_Spatial > 200 & Oral_ST_P13Normal$percent_mito < 25 ]
Oral_ST_P13Normal <- Oral_ST_P13Normal[!grepl("^mt-", rownames(Oral_ST_P13Normal)), ]
#saveRDS(Oral_ST_P13Normal,file="Oral_ST_P13Normal_seurat.rds")


#merge ST datas from one patient
setwd("/fs/home/hanya/Project/Carcinogenesis/Spatial_processed/")

Oral_ST_P13 <- merge(Oral_ST_P13Normal, Oral_ST_P13Tumor)
Oral_ST_P13 <- PercentageFeatureSet(Oral_ST_P13, "^MT-", col.name = "percent_mito")
Oral_ST_P13@project.name="Oral_ST_P13"

pdf("Oral_ST_P13_QC.pdf",width = 12,height = 8)
SpatialFeaturePlot(Oral_ST_P13, features = c("nCount_Spatial", "nFeature_Spatial","percent_mito"))
dev.off()

Oral_ST_P13 = Oral_ST_P13[, Oral_ST_P13$nFeature_Spatial > 200 & Oral_ST_P13$percent_mito < 25 ]
Oral_ST_P13 <- Oral_ST_P13[!grepl("^mt-", rownames(Oral_ST_P13)), ]

pdf("Oral_ST_P13_after_QC.pdf",width = 12,height = 12)
SpatialFeaturePlot(Oral_ST_P13, features = c("nCount_Spatial", "nFeature_Spatial","percent_mito"))
dev.off()

#Integrated two ST data by SCTransform
Oral_ST_P13 <- Oral_ST_P13[!grepl("^MT-", rownames(Oral_ST_P13)), ]
Oral_ST_P13 <- SCTransform(Oral_ST_P13, assay = "Spatial", verbose = TRUE, method = "poisson")


Oral_ST_P13 <- RunPCA(Oral_ST_P13, assay = "SCT", verbose = FALSE)
Oral_ST_P13 <- FindNeighbors(Oral_ST_P13, reduction = "pca", dims = 1:30)
Oral_ST_P13 <- FindClusters(Oral_ST_P13, verbose = FALSE)
Oral_ST_P13 <- RunUMAP(Oral_ST_P13, reduction = "pca", dims = 1:30)

pdf("Oral_ST_P13_Seurat_cluster.pdf",width = 12,height = 5)
DimPlot(Oral_ST_P13, reduction = "umap", group.by = c("ident", "orig.ident"))
dev.off()


#Quite often there are strong batch effects between different ST sections, so it may be a good idea to integrate the data across sections.
Oral_ST_P13_SCT = list(P13_Normal = Oral_ST_P13Normal, P13_Tumor = Oral_ST_P13Tumor)
# run SCT on both datasets
Oral_ST_P13_SCT = lapply(Oral_ST_P13_SCT, SCTransform, assay = "Spatial", method = "poisson")

options(future.globals.maxSize = 10000 * 1024^2) 
st.features = SelectIntegrationFeatures(Oral_ST_P13_SCT, nfeatures = 3000, verbose = FALSE)
Oral_ST_P13_SCT <- PrepSCTIntegration(object.list = Oral_ST_P13_SCT, anchor.features = st.features,verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = Oral_ST_P13_SCT, normalization.method = "SCT",verbose = FALSE, anchor.features = st.features)
Oral_ST_P13_Int <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

Oral_ST_P13_Int <- RunPCA(Oral_ST_P13_Int, verbose = FALSE)
Oral_ST_P13_Int <- FindNeighbors(Oral_ST_P13_Int, dims = 1:30)
Oral_ST_P13_Int <- FindClusters(Oral_ST_P13_Int, verbose = FALSE)
Oral_ST_P13_Int <- RunUMAP(Oral_ST_P13_Int, dims = 1:30)
pdf("Oral_ST_P13_Seurat_cluster.pdf",width = 12,height = 5)
DimPlot(Oral_ST_P13_Int, reduction = "umap", group.by = c("ident", "orig.ident"))
dev.off()

saveRDS(Oral_ST_P13_Int,file="Oral_ST_P13_Integrated_Seurat.rds")