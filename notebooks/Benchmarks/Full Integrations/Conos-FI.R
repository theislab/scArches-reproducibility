library(Seurat)
library(conos)

######
## The batches and the data in this script is for the pancreas dataset. Change them and the addresses for your dataset.
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

batches <- c("Pancreas inDrop.csv", "Pancreas CelSeq2.csv", 
            "Pancreas CelSeq.csv", "Pancreas Fluidigm C1.csv", "Pancreas SS2.csv")

batch.list <- list()

for (i in 1:5) {
  batch <- read.csv(paste("./", batches[i] ,sep=""), row.names = 1)
  batch_matrix <- data.matrix(batch)
  seurat_batch <- CreateSeuratObject(batch_matrix)
  seurat_batch@meta.data$batch = toString(i)
  VariableFeatures(seurat_batch) <- rownames(batch_matrix)
  so <- seuratPreprocess(seurat_batch)
  batch.list[[i]] <- so
}

names(batch.list) <- 1:5

t1 <- Sys.time()
con <- Conos$new(batch.list, n.cores=2)

con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30,
               n.odgenes=1000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)

con$findCommunities(method=leiden.community, resolution=1)

t2 <- Sys.time()

save.dir <- "./results/Conos/pancreas/corrected"
dir.create(save.dir)
saveConosForScanPy(con, output.path=save.dir, embed=TRUE, pseudo.pca=TRUE, n.dims = 10, verbose=TRUE)


write(as.numeric(t2 - t1, units = "secs"), file = paste0("./conos_time_pancreas.txt"))
