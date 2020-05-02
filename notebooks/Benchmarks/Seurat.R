library(Seurat)



batch1 <- read.csv("./mouse_brain_batch1.csv", row.names = 1)
batch2 <- read.csv("./mouse_brain_batch2.csv", row.names = 1)

target1 <- "Tabula_muris"
target2 <- "Zeisel"

for (i in 0:4){
  for (subsample_frac in c("0.1", "0.2", "0.4", "0.6", "0.8", "1.0")){
    
    keep_idx1 <- read.csv(paste("./data/subsample/mouse_brain/", target1, "/", subsample_frac, "/", i, ".csv", sep=""), header = FALSE)[, 1] + 1
    keep_idx2 <- read.csv(paste("./data/subsample/mouse_brain/", target2, "/", subsample_frac, "/", i, ".csv", sep=""), header = FALSE)[, 1] + 1
    
    batch1_matrix <- data.matrix(batch1)[, keep_idx1]
    batch2_matrix <- data.matrix(batch2)[, keep_idx2]
    
    seurat_batch1 = CreateSeuratObject(batch1_matrix)
    seurat_batch1@meta.data$batch = '1'
    VariableFeatures(seurat_batch1) <- rownames(batch1_matrix)
    
    
    seurat_batch2 = CreateSeuratObject(batch2_matrix)
    seurat_batch2@meta.data$batch = '2'
    VariableFeatures(seurat_batch2) <- rownames(batch2_matrix)
    
    
    
    pancreas.list <- list(seurat_batch1, seurat_batch2)
    names(pancreas.list) <- c(1, 2)
    
    k = seurat_batch1@assays$RNA@data@Dim[2]
    if(k < 200){
      pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, dims = 1:30, k.filter = k)
    }
    else{
      pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, dims = 1:30)
    }
    
    pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
    
    
    corrected = as.matrix(pancreas.integrated[["integrated"]]@data)
    
    write.csv(corrected, paste("./results/Seurat/mouse_brain/corrected-", subsample_frac, "-", i, ".csv", sep=""))

  }
}


