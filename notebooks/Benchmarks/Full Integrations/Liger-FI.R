library(liger)

######
## The target batches and the data in this script is for the brain dataset. Change the targets and the addresses for your dataset.
######

#batches <- c("Pancreas inDrop.csv", "Pancreas CelSeq2.csv", 
 #            "Pancreas CelSeq.csv", "Pancreas Fluidigm C1.csv", "Pancreas SS2.csv")

batches = c("Saunders.csv", "Rosenberg.csv", "Tabula_muris.csv", "Zeisel.csv")

batch.list <- list()

for (i in 1:4) {
  batch <- read.csv(paste("./", batches[i] ,sep=""), row.names = 1)
  batch.list[[paste("b",i,sep = "")]] <- batch
  gc()
}

data.liger <- createLiger(batch.list, remove.missing = F)
  
data.liger@norm.data <- data.liger@raw.data
    
var.genes <- row.names(batch)
data.liger@var.genes <- var.genes

data.liger <- scaleNotCenter(data.liger)

t1 <- Sys.time()

data.liger <- optimizeALS(data.liger, k=20, thresh = 5e-5, nrep = 5)
  
data.liger <- quantileAlignSNF(data.liger, resolution = 0.4)
  
t2 <- Sys.time()

write(as.numeric(t2 - t1, units = "secs"), file = paste0("./liger_time_mouse_brain.txt"))


corrected1 = as.matrix(data.liger@H[["b1"]])
corrected2 = as.matrix(data.liger@H[["b2"]])
corrected3 = as.matrix(data.liger@H[["b3"]])
corrected4 = as.matrix(data.liger@H[["b4"]])


  
write.csv(corrected1, "./results/Liger/mouse_brain/corrected1.csv")
write.csv(corrected2, "./results/Liger/mouse_brain/corrected2.csv")
write.csv(corrected3, "./results/Liger/mouse_brain/corrected3.csv")
write.csv(corrected4, "./results/Liger/mouse_brain/corrected4.csv")

