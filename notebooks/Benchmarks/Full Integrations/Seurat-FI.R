library(Seurat)

######
## The target batches and the data in this script is for the brain dataset. Change them and the addresses for your dataset.
######

#batches <- c("Pancreas inDrop.csv", "Pancreas CelSeq2.csv", 
 #            "Pancreas CelSeq.csv", "Pancreas Fluidigm C1.csv", "Pancreas SS2.csv")

batches = c("Saunders.csv", "Rosenberg.csv", "Tabula_muris.csv", "Zeisel.csv")

batch.list <- list()

for (i in 1:4) {
  batch <- read.csv(paste("./", batches[i] ,sep=""), row.names = 1)
  batch_matrix <- data.matrix(batch)
  seurat_batch <- CreateSeuratObject(batch_matrix)
  seurat_batch@meta.data$batch = toString(i)
  VariableFeatures(seurat_batch) <- rownames(batch_matrix)
  batch.list[[i]] <- seurat_batch
  gc()
}

names(batch.list) <- 1:4

features <- SelectIntegrationFeatures(object.list = batch.list)
batch.list <- lapply(X = batch.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

t1 <- Sys.time()
s.anchors <- FindIntegrationAnchors(object.list = batch.list, dims = 1:30, reduction = "rpca")
t2 <- Sys.time()

t3 <- Sys.time()
s.integrated <- IntegrateData(anchorset = s.anchors, dims = 1:30)
t4 <- Sys.time()


corrected = as.matrix(s.integrated[["integrated"]]@data)

write.csv(corrected, "./results/Seurat/mouse_brain/corrected.csv")


write(as.numeric(t4 - t1, units = "secs"), file = "./seurat_time_mouse_brain.txt")





