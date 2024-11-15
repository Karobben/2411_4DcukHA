library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(ggplot2)


# read the Counts
S1.data <- Read10X(data.dir = "all_GEX/duck_G1_GEX/outs/filtered_feature_bc_matrix")
S1 <- CreateSeuratObject(counts = S1.data, project = "duck_G1")
S2.data <- Read10X(data.dir = "all_GEX/duck_Group3_GEX/outs/filtered_feature_bc_matrix")
S2 <- CreateSeuratObject(counts = S2.data, project = "duck_Group3")
S3.data <- Read10X(data.dir = "all_GEX/miniHA3_G1_GEX/outs/filtered_feature_bc_matrix")
S3 <- CreateSeuratObject(counts = S3.data, project = "miniHA3_G1")
S4.data <- Read10X(data.dir = "all_GEX/miniHA3_G3_GEX/outs/filtered_feature_bc_matrix")
S4 <- CreateSeuratObject(counts = S4.data, project = "miniHA3_G3")

# merge the all
HA_5 <- merge(S1, y = c(S2, S3, S4), add.cell.ids = c('duck_G1', 'duck_Group3', 'miniHA3_G1', 'miniHA3_G3'), project = "duck_HA")

# check the mitochondria determinant
HA_5[["percent.mt"]] <- PercentageFeatureSet(HA_5, pattern = "^MT-")
VlnPlot(HA_5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(HA_5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HA_5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# HA_5 <- NormalizeData(HA_5, normalization.method = "LogNormalize", scale.factor = 10000)
HA_5 <- NormalizeData(HA_5)
HA_5 <- FindVariableFeatures(HA_5, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(HA_5), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(HA_5)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggsave('plot/variable_features.png')


all.genes <- rownames(HA_5)
HA_5 <- ScaleData(HA_5, features = all.genes)

HA_5 <- RunPCA(HA_5, features = VariableFeatures(object = HA_5))
# Examine and visualize PCA results a few different ways
print(HA_5[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(HA_5, dims = 1:2, reduction = "pca")
DimPlot(HA_5, reduction = "pca")
ggsave('plot/PCA_DimPlot.png')

HA_5 <- FindNeighbors(HA_5, dims = 1:10)
HA_5 <- FindClusters(HA_5, resolution = 0.5)
# Look at cluster IDs of the first 5 cells

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
HA_5 <- RunUMAP(HA_5, dims = 1:10)
DimPlot(HA_5, reduction = "umap")
ggsave('plot/UMAP_DimPlot.png')

DimPlot(HA_5, reduction = "umap", split.by = 'orig.ident', ncol = 3)
ggsave('plot/UMAP_DimPlot_sample.png', width = 7, height = 5)

HA_5 <- RunTSNE(HA_5, dims = 1:10, check_duplicates = FALSE)
DimPlot(HA_5, reduction = "tsne")
ggsave('plot/Tsne_DimPlot.png')

DimPlot(HA_5, reduction = "tsne", split.by = 'orig.ident', ncol = 3)
ggsave('plot/Tsne_DimPlot_sample.png', width = 7, height = 5)

HA_5.markers <- FindAllMarkers(JoinLayers(HA_5), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

HA_5.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(HA_5, features = top20$gene) + NoLegend()

ggsave('plot/FeaturesHeatMap.png', width = 20, height = 20)

VlnPlot(HA_5, features = c("CD27", "CD38", "SDC1")) 
ggsave('plot/Vln_CDgenes3.png', width = 12, height = 3)

c9_list = data.frame(row.names(HA_5@meta.data)[HA_5@meta.data$seurat_clusters ==9])

write.table(c9_list, 'result/c9_id.txt', quote =F, row.names = F)
write.csv(HA_5@meta.data, "result/GEX_anno.csv")

saveRDS(HA_5, file = "result/AllGEX.rds")
# HA_5 <- readRDS("result/AllGEX.rds")
######################################################
## Check cells
######################################################

TB <-read.table('data/AB_screening.tsv')
TB$ID <- paste(TB$V1, data.frame(str_split_fixed(TB$V2, '_', 2))[[1]], sep = '_')

HA_5@meta.data$CellName = row.names(HA_5@meta.data)


HA_5@meta.data$Screening = NA
HA_5@meta.data$Screening <- TB$V3[match(HA_5@meta.data$CellName, TB$ID)]
DimPlot(HA_5, reduction = "tsne",   cells.highlight = 'Screening')


Plot = data.frame(HA_5@reductions$umap@cell.embeddings)
Plot <- cbind(Plot, HA_5@meta.data)

ggplot(Plot, aes(umap_1, umap_2, color =  factor(seurat_clusters))) +
  geom_point(alpha = .2) +
  geom_point(data = Plot[!is.na(Plot$Screening),], aes(shape = Screening), color = 'black') +
  theme_bw()
ggsave('plot/Umap_DimPlot_vdj_positive.png')



Plot = data.frame(HA_5@reductions$tsne@cell.embeddings)
Plot <- cbind(Plot, HA_5@meta.data)

ggplot(Plot, aes(tSNE_1, tSNE_2, color =  factor(seurat_clusters))) +
  geom_point(alpha = .2) +
  geom_point(data = Plot[!is.na(Plot$Screening),], aes(shape = Screening), color = 'black') +
  theme_bw()
ggsave('plot/Tsne_DimPlot_vdj_positive.png')


table(Plot[!is.na(Plot$Screening), c("Screening", 'seurat_clusters')])


# Diff from 2 groups
Screened <- subset(HA_5, subset = CellName %in% TB$ID )
Screened <- JoinLayers(Screened)

Screened@meta.data


Screened.features <- FindMarkers(Screened,
                    group.by = 'Screening',
                    ident.1 = "Positive",
                    min.pct = 0.25,
                    logfc.threshold = 0.25,
                    ident.2 = "Negative"
                    )

Screened.features %>%
    group_by() %>%
    top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(Screened, features = top20$gene) + NoLegend()
ggsave('plot/FeaturesHeatMap_Screened.png')

VlnPlot(Screened, features = "TMEM255B", group.by = "Screening") 
ggsave('plot/Vln_TMEM255B.png')

ggsave('plot/Vln_TMEM255B.png', width = 3, height = 3)

VlnPlot(Screened, features = c("CD27", "CD38", "SDC1"), group.by = "Screening") 
ggsave('plot/Vln_CDgenes3_screen.png', width = 12, height = 3)

head(Screened.features[order(abs(Screened.features$avg_log2FC), decreasing = T),], 30)


###################
## Based on mini HA
###################

library(reshape2)

# read HA_5 
HA_5 <- readRDS("result/AllGEX.rds")
# read productive VDJ 
VDJ_TB <- read.table("result/Duck_numbring_class.tsv", sep = '\t', header = T)


HA_5[["CellName"]] <- colnames(HA_5)
miniGroup <- c("miniHA3_G1", "miniHA3_G3")
HA_5[["SType"]] <- "G"
HA_5@meta.data$SType[HA_5@meta.data$orig.ident %in% miniGroup] <- 'mini'

HA_5.VDJPositive <- subset(HA_5, subset = CellName %in% VDJ_TB$cell_id )

TB_clusCnt.Posi <- data.frame(t(as.data.frame.matrix(table(HA_5.VDJPositive@meta.data[c('SType', 'seurat_clusters')]))))
TB_clusCnt.all <- data.frame((table(HA_5@meta.data[c('seurat_clusters')])))

TB_clus_count <- cbind(TB_clusCnt.Posi, all = TB_clusCnt.all$Freq)
TB_clus_count$Class <- row.names(TB_clus_count)
TB_clus_count$all <- TB_clus_count$all - TB_clus_count$G - TB_clus_count$mini
TB_clus_count.m <- melt(TB_clus_count, id.var='Class')


head(TB_clus_count.m)

ggplot(TB_clus_count.m, aes(Class, value, fill = variable)) + geom_bar(stat = 'identity', position = 'fill') +
  theme_bw()
ggsave('plot/ClusterInVDJProductive.png', w = 7, h = 4.5)



































## find a group of cell
TB <- read.csv("result/All_AB_filter.tsv", sep = '\t')
tmp <- data.frame(ID1 = c("Mich", "mini1", "mini2"), ID2 = rev(unique(TB$sample)))

C_lst <- paste(tmp$ID1[match(TB$sample, tmp$ID2)], as.data.frame(str_split_fixed(TB$sequence_id, '_', 2))[[1]], sep = '_')
C_lst <- unique(C_lst)

HA_5@meta.data$VDJ_productive = "others"
HA_5@meta.data$VDJ_productive[rownames(HA_5@meta.data) %in% unique(C_lst)] = 'positive'


## All_results
tmp1 <- read.csv('result/VDJ_anno_Mich15_HA_VDJ.tsv', sep = '\t')
tmp2 <- read.csv('result/VDJ_anno_miniHA_1_VDJ.tsv', sep = '\t')
tmp3 <- read.csv('result/VDJ_anno_miniHA_2_VDJ.tsv', sep = '\t')

tmp1$Sampel = "Mich"
tmp2$Sampel = "mini1"
tmp3$Sampel = "mini2"

tmp_all <- rbind(tmp1, tmp2, tmp3)
tmp_all$Cell = paste(tmp_all$Sampel, as.data.frame(str_split_fixed(tmp_all$sequence_id, '_', 2))[[1]], sep = '_')

HA_5@meta.data$in_VDJ = F
HA_5@meta.data$in_VDJ[rownames(HA_5@meta.data) %in% tmp_all$Cell] = T

saveRDS(HA_5, file = "result/AllGEX.rds")

DimPlot(HA_5, reduction = "tsne", label = T)
DimPlot(HA_5, reduction = "tsne", group.by = 'in_VDJ')
DimPlot(HA_5, reduction = "tsne", group.by = 'VDJ_productive')

P1 <- DimPlot(HA_5, reduction = "umap", label = T)
P2 <- DimPlot(HA_5, reduction = "umap", group.by = 'in_VDJ')
P3 <- DimPlot(subset(HA_5, cells = tmp_all$Cell), reduction = "umap", group.by = 'VDJ_productive')
P4 <- ggplot(HA_5@meta.data, aes(seurat_clusters, fill = VDJ_productive)) + geom_bar(positio= 'fill') + 
  theme_bw()

P1/P4|P2/P3
ggsave('plot/GEX_VDJ_Scatter_bar_umap.png', w = 10.5, h = 7.87)


P1 <- DimPlot(HA_5, reduction = "tsne", label = T)
P2 <- DimPlot(HA_5, reduction = "tsne", group.by = 'in_VDJ')
P3 <- DimPlot(subset(HA_5, cells = tmp_all$Cell), reduction = "tsne", group.by = 'VDJ_productive')
P4 <- ggplot(HA_5@meta.data, aes(seurat_clusters, fill = VDJ_productive)) + geom_bar(positio= 'fill') + 
  theme_bw()

P1/P4|P2/P3
ggsave('plot/GEX_VDJ_Scatter_bar_tsne.png', w = 10.5, h = 7.87)




# Subset of the HA_5 onlye 

HA_5_subset <- subset(HA_5, cells = tmp_all$Cell)

DimPlot(HA_5_subset, reduction = "umap", label = T)
DimPlot(HA_5_subset, reduction = "umap", group.by = 'VDJ_productive')

DimPlot(HA_5_subset, reduction = "tsne", label = T)
DimPlot(HA_5_subset, reduction = "tsne", group.by = 'VDJ_productive')


all.genes <- rownames(HA_5_subset)
HA_5_subset <- ScaleData(HA_5_subset, features = all.genes)

HA_5_subset <- RunPCA(HA_5_subset, features = VariableFeatures(object = HA_5_subset))
# Examine and visualize PCA results a few different ways
DimPlot(HA_5_subset, reduction = "pca")
DimPlot(HA_5_subset, reduction = "pca", group.by = 'VDJ_productive' )

HA_5_subset <- FindNeighbors(HA_5_subset, dims = 1:10)
HA_5_subset <- FindClusters(HA_5_subset, resolution = 0.5)
# Look at cluster IDs of the first 5 cells

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
HA_5_subset <- JoinLayers(HA_5_subset)

HA_5_subset <- RunUMAP(HA_5_subset, dims = 1:10)
HA_5_subset <- RunTSNE(HA_5_subset, dims = 1:10, check_duplicates = FALSE)

P1 <- DimPlot(HA_5_subset, reduction = "tsne", label = T)
P2 <- DimPlot(HA_5_subset, reduction = "tsne", group.by = 'in_VDJ')
P3 <- DimPlot(HA_5_subset, reduction = "tsne", group.by = 'VDJ_productive')
P4 <- ggplot(HA_5_subset@meta.data, aes(seurat_clusters, fill = VDJ_productive)) + geom_bar(positio= 'fill') + 
  theme_bw()

P1/P4|P2/P3
ggsave('plot/GEX_VDJ_SubClass_Scatter_bar_tsne.png', w = 10.5, h = 7.87)

P1 <- DimPlot(HA_5_subset, reduction = "umap", label = T)
P2 <- DimPlot(HA_5_subset, reduction = "umap", group.by = 'in_VDJ')
P3 <- DimPlot(HA_5_subset, reduction = "umap", group.by = 'VDJ_productive')
P4 <- ggplot(HA_5_subset@meta.data, aes(seurat_clusters, fill = VDJ_productive)) + geom_bar(positio= 'fill') + 
  theme_bw()

P1/P4|P2/P3
ggsave('plot/GEX_VDJ_SubClass_Scatter_bar_umap.png', w = 10.5, h = 7.87)



# features and annotation
M_VDJ_P <- FindMarkers(HA_5_subset, ident.1 = 'positive', ident.2 = NULL, group.by = 'VDJ_productive' )
M_VDJ_C5_umap <- FindMarkers(HA_5_subset, ident.1 = 5, ident.2 = NULL)

L1 = row.names(M_VDJ_P[abs(M_VDJ_P$avg_log2FC)>=2,])
L2 = row.names(M_VDJ_C5_umap[abs(M_VDJ_C5_umap$avg_log2FC)>=2,])
L12 = table(c(L1, L2)) 
L12_Gene = names(L12)[L12==2]
write.table(L12_Gene, 'result/Class_VDJ_intersect_Genes.txt', quote = F, row.names = F)


library(clusterProfiler)
library(org.Ss.eg.db)

gene.list <- bitr(row.names(M_VDJ_P), fromType = "SYMBOL",
        toType = 'ENTREZID',
        OrgDb = org.Ss.eg.db)

gene.list$log2FC <- M_VDJ_P$avg_log2FC[match(gene.list$SYMBOL, rownames(M_VDJ_P))]

geneList <- gene.list$log2FC
names(geneList) <- gene.list$ENTREZID
geneList <- sort(geneList, decreasing = T)


search_kegg_organism('pig', by='scientific_name')


kk <- enrichKEGG(gene         = names(geneList)[abs(geneList) >=2],
                 organism     = 'ssc',
                 pvalueCutoff = 0.05)

kk_tb <- kk@result[kk@result$pvalue<= 0.05,]
kk_tb <- kk_tb[kk_tb$qvalue<= 0.05,]

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'ssc',
               pvalueCutoff = 0.05,
               verbose      = FALSE)


library(enrichplot)

head(kk2@result,2)

i = 1
gseaplot2(kk2, geneSetID = i, title = kk2$Description[i], color = 'salmon', pvalue_table = TRUE)

ggsave('GSEM_KK.png', w = 5.46, h = 3.23)


GOgse_CC <- gseGO(geneList     = geneList,
              OrgDb        = org.Ss.eg.db,
              ont          = 'CC',
              minGSSize    = 10,
              maxGSSize    = 500,
              eps = 1e-10,
              pvalueCutoff = 0.1,
              verbose      = FALSE)

i = 7
gseaplot2(GOgse_CC, geneSetID = i, title = GOgse_CC$Description[i], color = 'salmon', pvalue_table = TRUE)
