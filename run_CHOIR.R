library(Matrix)
library(Seurat)
library(CHOIR)

# ESTA PARTE QUE ESTA COMENTADA, ES RECOMENDABLE CORRERLA EN UN SCRIPT DE R EN LA TERMINAL Y DEJAR CORRIENDO DE FONDO 
# (ej.con ~80,000 células toma aproximadamente 20 horas en terminar el proceso)

# cat(format(Sys.time(), "%H:%M:%S"), " - Iniciando lectura de archivos...\n")
# 
# # Leer matriz y nombres
# counts <- readMM("counts.mtx")
# rownames(counts) <- readLines("genes.tsv")
# colnames(counts) <- readLines("barcodes.tsv")
# 
# cat(format(Sys.time(), "%H:%M:%S"), " - Construyendo objeto Seurat...\n")
# object <- CreateSeuratObject(counts, min.features = 100, min.cells = 5)
# 
# cat(format(Sys.time(), "%H:%M:%S"), " - Normalizando datos...\n")
# object <- NormalizeData(object)
# 
# cat(format(Sys.time(), "%H:%M:%S"), " - Corriendo buildTree...\n")
# object <- buildTree(object, n_cores = 35)
# 
# cat(format(Sys.time(), "%H:%M:%S"), " - Corriendo pruneTree...\n")
# object <- pruneTree(object, n_cores = 35)
# 
# cat(format(Sys.time(), "%H:%M:%S"), " - Ejecutando runCHOIRumap...\n")
# object <- runCHOIRumap(object, reduction = "P0_reduction")
# 
# cat(format(Sys.time(), "%H:%M:%S"), " - Guardando objeto Seurat completo...\n")
# saveRDS(object, file = "CHOIR_Seurat_object.rds")
# 
# cat(format(Sys.time(), "%H:%M:%S"), " - Proceso completado correctamente.\n")

object <- readRDS("CHOIR_Seurat_object.rds")
Reductions(object)
plotCHOIR(object)
Embeddings(object, "CHOIR_P0_reduction_UMAP")





