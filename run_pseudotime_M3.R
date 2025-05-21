setwd("/home/mdiaz/scRNA-seq_LungProject")

library(batchelor)
library(terra)
library(ggrastr)
library(BiocGenerics)
library(DelayedArray)
library(DelayedMatrixStats)
library(limma)
library(lme4)
library(S4Vectors)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(HDF5Array)
library(monocle3)
library(Seurat)
library(SeuratDisk)
library(Seurat)
library(zellkonverter)
library(scater)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Cargar .h5ad limpiamente
sce <- readH5AD("filtrado_celltypist_tipos_interes_cleaned.h5ad")

seurat_obj <- as.Seurat(sce, counts = "X", data = "X")

head(seurat_obj@meta.data)


# Normalización de datos
seurat_obj <- NormalizeData(object = seurat_obj, verbose = FALSE)

# Selección de características variables
seurat_obj <- FindVariableFeatures(object = seurat_obj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')

# Escalado de datos
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

# Análisis de componentes principales (PCA)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

# Encontrar vecinos para la construcción del gráfico KNN
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)

# Reducción de dimensionalidad con UMAP
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30, verbose = FALSE)

# Visualización inicial
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

DefaultAssay(seurat_obj) <- "originalexp"

# Asegurar que la ranura 'counts' exista (necesaria para convertir a cell_data_set)
if (is.null(seurat_obj[["originalexp"]]@counts)) {
  seurat_obj[["originalexp"]]@counts <- seurat_obj[["originalexp"]]@data
}


# Convertir el objeto Seurat a un cell_data_set de Monocle3
cds <- as.cell_data_set(seurat_obj)

colData(cds)$final_label <- seurat_obj@meta.data$final_label


# Clustering en Monocle3 (sin partición específica)
cds <- cluster_cells(cds)

# Visualización básica de las células sin grafo de trayectoria
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")

plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "final_label")
#########################################################################
# Extraer coordenadas UMAP y etiquetas
umap_df <- as.data.frame(reducedDims(cds)$UMAP)
umap_df$final_label <- colData(cds)$final_label

# Calcular centroides
centroids <- aggregate(cbind(umap_1, umap_2) ~ final_label, data = umap_df, FUN = median)

# Graficar UMAP con etiquetas grandes y en negrita
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = final_label)) +
  geom_point(size = 0.2, alpha = 0.8) +
  geom_text_repel(
    data = centroids,
    aes(x = umap_1, y = umap_2, label = final_label),
    size = 6,
    fontface = "bold"
  ) +
  theme_minimal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "none")
##########################################################################
# Aprender el grafo de trayectorias sin particionar
cds <- learn_graph(cds, use_partition = FALSE)

# Visualización del grafo aprendido
plot_cells(cds, show_trajectory_graph = TRUE, color_cells_by = "cluster")


# Aprender el grafo sin particionar
cds <- learn_graph(cds, use_partition = FALSE)


# Ordenar las células (definir la raíz automáticamente)
cds <- order_cells(cds)

# Visualizar las células coloreadas por pseudotiempo
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_branch_points = FALSE, 
           label_leaves = FALSE, 
           show_trajectory_graph = TRUE)

################################
# ENSEMBL TO GENE_SYMBOL (Opcional pero util)
################################
library(biomaRt)
DefaultAssay(seurat_obj) <- "originalexp"
seurat_obj[["originalexp"]]@counts <- seurat_obj[["originalexp"]]@counts  # Asegura la existencia

cds <- as.cell_data_set(seurat_obj)
library(biomaRt)
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = rownames(cds),
                     mart = mart)
symbol_map <- setNames(annotations$hgnc_symbol, annotations$ensembl_gene_id)
rowData(cds)$gene_short_name <- symbol_map[rownames(cds)]
# 
# ####################################

# # Crear vector con símbolos
symbol_map <- setNames(annotations$hgnc_symbol, annotations$ensembl_gene_id)

# Asignar gene_short_name en Monocle3
rowData(cds)$gene_short_name <- symbol_map[rownames(cds)]

# Asegurarse de que los nombres de los genes estén correctamente asignados
rowData(cds)$gene_name <- rownames(seurat_obj)

# Asignar los nombres cortos de los genes
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
plot_cells(cds,
           genes=c('ENSG00000142227', 'ENSG00000005020', 'ENSG00000011600'),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, 
           min_expr = 3)

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE)
cds_pt_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)


cds_subset <- choose_cells(cds)
cds_subset

cds_subset_pt_res <- graph_test(cds_subset, neighbor_graph = "principal_graph", cores = 30)

# Filtrar resultados significativos
cds_subset_pt_res <- na.omit(cds_subset_pt_res)
cds_subset_pt_res <- cds_subset_pt_res[cds_subset_pt_res$p_value < 0.05 &
                                         cds_subset_pt_res$status == "OK", ]

# Mostrar genes asociados al pseudotiempo
cds_subset_pt_res

cds_subset_pt_res[order(-cds_subset_pt_res$morans_test_statistic), ]

plot_cells(cds,
           genes = c('ENSG00000197728', 'ENSG00000120885', 'ENSG00000168484'),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           reduction_method = "UMAP")

# Subconjunto de genes usando Ensembl IDs
cds_subset_subset <- cds_subset[rowData(cds_subset)$gene_name %in% c(
  'ENSG00000197728', 'ENSG00000120885', 'ENSG00000168484'), ]

# Volver a ordenar las células del subconjunto
cds_subset <- order_cells(cds_subset)
cds_subset_subset <- cds_subset[rowData(cds_subset)$gene_name %in% c(
  'ENSG00000197728', 'ENSG00000120885', 'ENSG00000168484'), ]

plot_genes_in_pseudotime(cds_subset_subset,
                         color_cells_by = "final_label",
                         min_expr = 0.5)

plot_genes_in_pseudotime(cds_subset_subset,
                         color_cells_by = "lung_condition",
                         min_expr = 0.5)

plot_genes_in_pseudotime(cds_subset_subset,
                         color_cells_by = "condition_binary",
                         min_expr = 0.5)