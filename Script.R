#install.packages("Seurat")
#install.packages("igraph")
#install.packages("tidyverse")
#install.packages('BiocManager')
#BiocManager::install('glmGamPoi')
#BiocManager::install('presto')
#BiocManager::install('HGNChelper')

library(Seurat)
library(igraph)
library(tidyverse)
library('glmGamPoi')
library('presto')


#Upload the data 
indx <- list.files(pattern = "*_matrix.mtx")
for (i in 1:length(indx)) { 
  assign(paste("dat",i,sep=""),ReadMtx(
    mtx = paste(gsub("_matrix.mtx","",indx[i]),"_matrix.mtx",sep=""), 
    features = paste(gsub("_matrix.mtx","",indx[i]),"_features.tsv",sep=""),
    paste(gsub("_matrix.mtx","",indx[i]),"_barcodes.tsv",sep=""),feature.column = 1))
} 

files <- mget(paste0("dat", 1:9))
for (i in 1:length(files)) {
  colnames(files[[i]]) <- paste(gsub("_matrix.mtx","",indx[i]),colnames(files[[i]]),sep="")
}
rm(list = paste0("dat", 1:9))

#Creating Seurat object
for (i in 1:length(files)) { 
  assign(paste("dat",i,sep=""),
         CreateSeuratObject(files[[i]]))
} 

seurat_list <- mget(paste0("dat", 1:9))
core <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])
core

rm(list = paste0("dat", 1:9))

#QC

#Adding mitochondrial percent score
core[["percent.mt"]] <- PercentageFeatureSet(object = core, pattern = "^MT-")
#sanity check
head(core@meta.data, 10)
dim(core)

#QC plots
VlnPlot(core, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        pt.size = 0) +
  theme(
    legend.position = "none",
  )

dim(core)
core <- subset(core, nFeature_RNA > 500 & nFeature_RNA < 4000 & 
                 nCount_RNA> 1000 & nCount_RNA< 10000 & percent.mt < 5)
dim(core)

#Normalization the data
core <- SCTransform(
  core,
  verbose=TRUE,
  return.only.var.genes=T
)
dim(core@assays$SCT)
core <- RunPCA(core)

DimHeatmap(core, dims = 1:12, cells = 500, balanced = TRUE)
ElbowPlot(core,ndims=100)

npcs <- 40

core <- RunUMAP(core, 
                dims = 1:npcs, 
                verbose = TRUE, 
                seed.use = 42)

DimPlot(core, 
        reduction = "umap", 
        pt.size = 0.3)

#Data integration
core <- IntegrateLayers(object = core, method = "CCAIntegration", normalization.method = "SCT", new.reduction = "integrated.cca",
                        verbose = T)
core[["RNA"]] <- JoinLayers(core[["RNA"]])
core@assays$RNA$counts[1:4,1:4]
Layers(core[["RNA"]])
DefaultLayer(core[["RNA"]]) <- "counts"
core[["RNA"]] <- JoinLayers(core[["RNA"]])

#Clustering
#Neighbors
core <- FindNeighbors(core, reduction = "integrated.cca", dims = 1:30, k.param = 20)
#Cluster
core <- FindClusters(core, resolution = 0.5)
#UMAP
core <- RunUMAP(core, dims = 1:30, reduction = "integrated.cca")
core <- PrepSCTFindMarkers(core)

#Gene of interest
data.markers <- FindAllMarkers(
  core,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

Idents(core) <- "seurat_clusters"

goi <- data.markers %>% 
  as_tibble() %>% 
  group_by(cluster) %>% 
  slice_max(n = 3, order_by = avg_log2FC) %>% 
  pull(gene) %>% 
  unique()

DotPlot(core, features = goi, cluster.idents = TRUE) + coord_flip()

DimPlot(core, reduction = "umap", group.by = "seurat_clusters")

#Check the clusters
cluster_var <- 'seurat_clusters'
clusters <- unique(core@meta.data[[cluster_var]])
clusters <- clusters[order(clusters)]
df <- data.frame()
for(i in 1:length(clusters)){
  cur_df <- as.data.frame(core@meta.data %>% subset(seurat_clusters == clusters[i]) %>% .$orig.ident %>% table())
  
  cur_df$Freq <- cur_df$Freq * 1/(sum(cur_df$Freq))
  
  cur_df$cluster <- clusters[i]
  df <- rbind(df, cur_df)
  print(i)
}
df_counts <- df %>%
  group_by(cluster) %>%
  summarise(
    Freq = unique(Freq),
    n_cells = n()
  )
df1 <- df %>%
  group_by(cluster) %>%
  mutate(n_cells = n()) %>%
  ungroup()
p <- ggplot(df1, aes(y=Freq, x=cluster, fill=.)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  ylab('normalized proportion') +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
  )
print(p)

DimPlot(core, reduction = "umap", 
        group.by = c("orig.ident", "seurat_clusters")) & theme(legend.position = "none")


saveRDS(core,"scRNA_GBM_clustered.RDS")
