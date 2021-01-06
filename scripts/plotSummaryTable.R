library(tibble)
library(RColorBrewer)
library(dynutils)
library(stringr)
library(Hmisc)
library(plyr)


source("/home/python_scRNA/Munich/MoMo/MoMo_knit_table.R")# You will need to have in the same folder knit_table.R and this plotSingleAtlas.R

# parameters: 
# - 'csv_file_path' would be the path of the csv file 
# - 'type' would be one between "ratio" and "tuning"

plotSummaryTable <- function(csv_file_path, outdir = ".", weight_batch = 0.4, type, plot_time = TRUE){
  
  
  metrics_tab_lab <- read.csv(csv_file_path, sep = ",")
  add_columns <- c("data", "rqr", "method")
  
  
  # rename columns
  col_names <- colnames(metrics_tab_lab)
  col_names <- gsub("\\.", "/", col_names)
  col_names <- gsub("_", " ", col_names)
  col_names <- plyr::mapvalues(col_names, from = c("ASW label", "ASW label/batch", "cell cycle conservation", "hvg overlap",  "graph conn", "KNN"), 
                                               to = c("Cell type ASW", "Batch ASW", "CC conservation", "HVG conservation", "graph connectivity", "kNN accuracy"))
  colnames(metrics_tab_lab) <- col_names
  
  
  
  # metrics names as they are supposed to be ordered
  group_batch <- c("Batch ASW", "PCR batch",  "graph iLISI", "graph connectivity", "kBET", "EBM", "kNN accuracy")
  group_bio <- c("NMI cluster/label", "ARI cluster/label", "Cell type ASW", 
                 "isolated label F1", "isolated label silhouette", "graph cLISI", "CC conservation", "HVG conservation", "trajectory conservation")
  group_time <- c("reference time", "query time")
  # set original values of number of metrics
  n_metrics_batch_original <- sum(group_batch %in% col_names)
  n_metrics_bio_original <- sum(group_bio %in% col_names)
  n_metrics_time_original <- sum(group_time %in% col_names)
  
  
  # order metrics present in the table
  matching.order <- match(c(add_columns, group_batch, group_bio, group_time), col_names)
  col_names.ord <- col_names[matching.order[!is.na(matching.order)]]
  metrics_tab_lab <- metrics_tab_lab[, col_names.ord]
  
  # data scenarios to be saved in file name
  dataset <- unique(as.character(metrics_tab_lab$data))
    
  # info on ratio
  ratio <- metrics_tab_lab$rqr
    
  # time: reference and query scaled together
  time_scaled <- scale_minmax(as.vector(as.matrix(metrics_tab_lab[, group_time])))
  metrics_tab_lab[, group_time[1]] <- time_scaled[1:nrow(metrics_tab_lab)]
  metrics_tab_lab[, group_time[2]] <- time_scaled[(1+nrow(metrics_tab_lab)):length(time_scaled)]
  
  methods <- as.character(metrics_tab_lab$method)
    
    
  ##### Create dataframe of only metrics (excluding time, which is already scaled)
  metrics_tab <- as.data.frame(metrics_tab_lab[, !colnames(metrics_tab_lab) %in% c(add_columns, group_time)])
  metrics_tab[metrics_tab == ""] <- NA
  
  
  ## Remove columns that are full NAs
  na.col <- apply(metrics_tab, 2, function(x) sum(is.na(x)) == nrow(metrics_tab))
  
  # redefine numbers of metrics per group
  if(sum(colnames(metrics_tab)[na.col] %in%  group_batch) > 0){
    n_metrics_batch <- n_metrics_batch_original - sum(colnames(metrics_tab)[na.col] %in%  group_batch)
  } else {
    n_metrics_batch <- n_metrics_batch_original
  }

  if(sum(colnames(metrics_tab)[na.col] %in% group_bio) > 0){
    n_metrics_bio <- n_metrics_bio_original - sum(colnames(metrics_tab)[na.col] %in% group_bio)
  } else{
    n_metrics_bio <- n_metrics_bio_original
  }
    
  metrics_tab <- metrics_tab[, !na.col]
  
  # Scale scores [0,1]
  scaled_metrics_tab <- apply(metrics_tab, 2, function(x) scale_minmax(x))
  
  # calculate average score by group and overall
  score_group_batch <- rowMeans(scaled_metrics_tab[, colnames(scaled_metrics_tab) %in% group_batch], na.rm = T)
  score_group_bio <- rowMeans(scaled_metrics_tab[, colnames(scaled_metrics_tab) %in% group_bio], na.rm = T)
  
  score_all <- (weight_batch*score_group_batch + (1-weight_batch)*score_group_bio)
  
  metrics_tab <- as.data.frame(scaled_metrics_tab)
  metrics_tab <- add_column(metrics_tab, "Method" = methods, .before = 1)
  metrics_tab <- add_column(metrics_tab, "Overall Score" = score_all, .after = "Method")
  metrics_tab <- add_column(metrics_tab, "Batch Correction" = score_group_batch, .after = "Overall Score")
  metrics_tab <- add_column(metrics_tab, "Bio conservation" = score_group_bio, .after = 3+n_metrics_batch)
  
  if(type == "ratio"){
    metrics_tab <- add_column(metrics_tab, "Ratio" = ratio, .after = "Method")
    column_2_width <- 1.5
  } else if(type == "tuning"){
    metrics_tab <- add_column(metrics_tab, "Fine tuning" = ratio, .after = "Method")
    column_2_width <- 4
  }
  
  metrics_tab <- add_column(metrics_tab, "Reference Time" = metrics_tab_lab$`reference time`, .after = ncol(metrics_tab))
  metrics_tab <- add_column(metrics_tab, "Query Time" = metrics_tab_lab$`query time`, .after = ncol(metrics_tab))
  
  # order methods by the overall score
  metrics_tab <- metrics_tab[order(metrics_tab$`Overall Score`,  decreasing = T), ]
  write.csv(metrics_tab, file = paste0(outdir, "/", dataset, "_summary_scores.csv"), quote = F)
  
  
  # Delete rows that are empty
  rowsNA <- which(is.na(metrics_tab$`Overall Score`))
  if(length(rowsNA) >0){
    metrics_tab <- metrics_tab[-rowsNA, ]
  }
  
  
  # Defining column_info, row_info and palettes
  row_info <- data.frame(id = metrics_tab$Method)
  
  column_info <- data.frame(id = colnames(metrics_tab),
                            group = c("Text", "Text", "Score overall", 
                                      rep("Removal of batch effects", (1 + n_metrics_batch)),
                                      rep("Cell type label variance", (1 + n_metrics_bio)),
                                      rep("Time", 2)), 
                            geom = c("text", "text", "bar", "bar", 
                                     rep("circle", n_metrics_batch), "bar", rep("circle", n_metrics_bio), 
                                     rep("bar", n_metrics_time_original)),
                            width = c(7,column_2_width,2,2, rep(1,n_metrics_batch), 2, rep(1,n_metrics_bio), rep(2, n_metrics_time_original)),
                            overlay = F)
  
  # defining colors palette
  palettes <- list("Score overall" = "YlGnBu",
                   "Removal of batch effects" = "BuPu",
                   "Cell type label variance" = "RdPu",
                   "Time" = "YlOrRd")
  
  if(!plot_time){
    column_info <- column_info[1:(nrow(column_info)-2),]
    palettes <- palettes[-which(names(palettes) == "Time")]
  }
  
  
  g <- scIB_knit_table(data = metrics_tab, column_info = column_info, row_info = row_info, palettes = palettes, plot_time)  
  now <- Sys.time()
  ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dataset, "_summary_metrics.pdf"), g, device = cairo_pdf, width = 297, height = 420, units = "mm")
  ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dataset, "_summary_metrics.tiff"), g, device = "tiff", dpi = "retina", width = 297, height = 420, units = "mm")
  ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dataset, "_summary_metrics.png"), g, device = "png", dpi = "retina", width = 210, height = 297, units = "mm")
  
  
}
  
  

