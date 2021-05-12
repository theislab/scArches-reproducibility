library(harmony)

for (data_name in c("pancreas", "brain")){
  
  x <- read.csv(paste("./data/", data_name, "/", "/X.csv", sep=""), row.names=1)
  batch <- read.csv(paste("./data/", data_name, "/", "/Batch.csv", sep=""), row.names=1)
  
  harmonyCorrected <- HarmonyMatrix(x, batch, 'study', do_pca=F)
  
  harmonyCorrected <- data.frame(harmonyCorrected)
  
  write.csv(harmonyCorrected, paste("./results/Harmony/", data_name, "/harmonyCorrected.csv", sep=""))
  
}