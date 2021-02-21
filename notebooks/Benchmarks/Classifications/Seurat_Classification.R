library(Seurat)

conditions <- c('droplet - 24m', 'droplet - 18m', 'droplet - 21m',
                'droplet - 1m', 'droplet - 30m', 'facs - 18m', 'facs - 24m', 'facs - 21m')



conditions.list <- list()



for (i in 1:8) {
  cond <- read.csv(paste0("./", conditions[i], ".csv"), row.names = 1)
  cond_matrix <- data.matrix(cond)
  seurat_cond <- CreateSeuratObject(cond_matrix)
  seurat_cond@meta.data$batch = toString(i)
  seurat_cond@meta.data$celltype <- as.character(read.csv(paste0("./", conditions[i], "_celltype.csv"))[,1])
  VariableFeatures(seurat_cond) <- rownames(cond_matrix)
  conditions.list[[i]] <- seurat_cond
  gc()
}

names(conditions.list) <- 1:8

features <- SelectIntegrationFeatures(object.list = conditions.list)
conditions.list <- lapply(X = conditions.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


s.anchors <- FindIntegrationAnchors(object.list = conditions.list, dims = 1:30, reduction = "rpca")
s.integrated <- IntegrateData(anchorset = s.anchors, dims = 1:30)


target_conditions <- c("droplet - 3m", "facs - 3m")
target.list <- list()


for (i in 1:2) {
  cond <- read.csv(paste0("./", target_conditions[i], ".csv"), row.names = 1)
  cond_matrix <- data.matrix(cond)
  seurat_cond <- CreateSeuratObject(cond_matrix)
  seurat_cond@meta.data$batch = toString(i)
  seurat_cond@meta.data$celltype <- as.character(read.csv(paste0("./", target_conditions[i], "_celltype.csv"))[,1])
  VariableFeatures(seurat_cond) <- rownames(cond_matrix)
  target.list[[i]] <- seurat_cond
  gc()
}

names(target.list) <- 1:2

query <- target.list[[2]]
s.anchors <- FindTransferAnchors(reference = s.integrated, query = query, dims = 1:30)

predictions <- TransferData(anchorset = s.anchors, refdata = s.integrated$celltype, dims = 1:30)




write.csv(predictions[,c("predicted.id","prediction.score.max")], "facs - 3m_pred.csv")

