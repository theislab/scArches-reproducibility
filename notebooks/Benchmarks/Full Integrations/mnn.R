library(batchelor)

data_name <- "pancreas"


batch1 <- data.matrix(read.csv("./Pancreas inDrop.csv", row.names = 1))
batch2 <- data.matrix(read.csv("./Pancreas CelSeq2.csv", row.names = 1))
batch3 <- data.matrix(read.csv("./Pancreas CelSeq.csv", row.names = 1))
batch4 <- data.matrix(read.csv("./Pancreas Fluidigm C1.csv", row.names = 1))
batch5 <- data.matrix(read.csv("./Pancreas SS2.csv", row.names = 1))

t1 <- Sys.time()
mnnCorrected = mnnCorrect(batch1, batch2, batch3, batch4, batch5)
t2 <- Sys.time() 

write(as.numeric(t2 - t1, units = "secs"), file = "./mnn_time_pancreas.txt")

corrected = data.frame(mnnCorrected@assays@data@listData[["corrected"]])

write.csv(corrected, paste("./results/mnnCorrect/", data_name, "/corrected.csv", sep=""))