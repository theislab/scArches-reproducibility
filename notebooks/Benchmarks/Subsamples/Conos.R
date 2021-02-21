library(Seurat)
library(conos)

######
## The target batches and the data in this script is for the pancreas dataset. Change the targets and the addresses for your dataset.
######

seuratPreprocess <- function(sb) {
  so <- ScaleData(object = sb)
  so <- RunPCA(object = so, npcs = 100)
  
  so <- FindNeighbors(object = so, dims = 1:100)
  so <- FindClusters(object = so, n.iter = 500, n.start = 10)
  
  so <- RunTSNE(object = so, dims = 1:100)
  so <- RunUMAP(object = so, dims = 1:100)
  
  return(so)
}



batch1 <- read.csv("./pancreas_batch1.csv", row.names = 1)
batch2 <- read.csv("./pancreas_batch2.csv", row.names = 1)

target1 <- "Pancreas SS2"
target2 <- "Pancreas CelSeq2"

for (i in 0:4){
  for (subsample_frac in c("0.1", "0.2", "0.4", "0.6", "0.8", "1.0")){

    keep_idx1 <- read.csv(paste("./data/subsample/pancreas/", target1, "/", subsample_frac, "/", i, ".csv", sep=""), header = FALSE)[, 1] + 1
    keep_idx2 <- read.csv(paste("./data/subsample/pancreas/", target2, "/", subsample_frac, "/", i, ".csv", sep=""), header = FALSE)[, 1] + 1
        
    batch1_matrix <- data.matrix(batch1)[, keep_idx1]
    batch2_matrix <- data.matrix(batch2)[, keep_idx2]
        
    seurat_batch1 = CreateSeuratObject(batch1_matrix)
    seurat_batch1@meta.data$batch = '1'
    VariableFeatures(seurat_batch1) <- rownames(batch1_matrix)
        
    
    seurat_batch2 = CreateSeuratObject(batch2_matrix)
    seurat_batch2@meta.data$batch = '2'
    VariableFeatures(seurat_batch2) <- rownames(batch2_matrix)
    
    so1 <- seuratPreprocess(seurat_batch1)
    so2 <- seuratPreprocess(seurat_batch2)
    
    batch.list <- list(so1, so2)
    names(batch.list) <- c(1, 2)
    
    con <- Conos$new(batch.list, n.cores=2)
    
    con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30,
                   n.odgenes=1000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)
    
    con$findCommunities(method=leiden.community, resolution=1)
    
    save.dir <- paste("./results/Conos/pancreas/corrected-", subsample_frac, "-", i, sep="")
    dir.create(save.dir)
    saveConosForScanPy(con, output.path=save.dir, embed=TRUE, pseudo.pca=TRUE, n.dims = 10, verbose=TRUE)

  }
}

