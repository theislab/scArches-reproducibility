library(liger)

######
## The target batches and the data in this script is for the pancreas dataset. Change the targets and the addresses for your dataset.
######

batch1 <- read.csv("./pancreas_batch1.csv", row.names = 1)
batch2 <- read.csv("./pancreas_batch2.csv", row.names = 1)

target1 <- "Pancreas SS2"
target2 <- "Pancreas CelSeq2"

for (i in 0:4){
  for (subsample_frac in c("0.1", "0.2", "0.4", "0.6", "0.8", "1.0")){

    keep_idx1 <- read.csv(paste("./data/subsample/pancreas/", target1, "/", subsample_frac, "/", i, ".csv", sep=""), header = FALSE)[, 1] + 1
    keep_idx2 <- read.csv(paste("./data/subsample/pancreas/", target2, "/", subsample_frac, "/", i, ".csv", sep=""), header = FALSE)[, 1] + 1
    
    batch1_df <- batch1[, keep_idx1]
    batch2_df <- batch2[, keep_idx2]
    
    data.list = list(b1=batch1_df, b2=batch2_df)
    data.liger <- createLiger(data.list, remove.missing = F)
    
    data.liger <- normalize(data.liger)
    
    var.genes <- row.names(batch1_df)
    data.liger@var.genes <- var.genes
    
    data.liger <- scaleNotCenter(data.liger)
    
    data.liger <- optimizeALS(data.liger, k=20, thresh = 5e-5, nrep = 5)
    
    
    data.liger <- quantileAlignSNF(data.liger, resolution = 0.4)
    
    
    corrected1 = as.matrix(data.liger@H[["b1"]])
    corrected2 = as.matrix(data.liger@H[["b2"]])
    
    
    write.csv(corrected1, paste("./results/Liger/pancreas/corrected1-", subsample_frac, "-", i, ".csv", sep=""))
    write.csv(corrected2, paste("./results/Liger/pancreas/corrected2-", subsample_frac, "-", i, ".csv", sep=""))

  }
}


