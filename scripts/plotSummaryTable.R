library(tibble)
library(RColorBrewer)
library(dynutils)
library(stringr)
library(Hmisc)
library(plyr)


source("/home/python_scRNA/Munich/MoMo/MoMo_knit_table.R")# You will need to have in the same folder knit_table.R and this plotSingleAtlas.R

# parameters: 
# - 'csv_file_path' would be the path of the csv file 
# - 'type' would be one between "ratio", "best_ratio" and "tuning"



plotSummaryTable <- function(csv_file_path, 
                             outdir = ".", 
                             weight_batch = 0.4, 
                             type, 
                             plot_time = TRUE,
                             tag){
  
  
  metrics_tab_lab <- read.csv(csv_file_path, sep = ",")
  add_columns <- c("data", "rqr", "method")
  
  
  
  # rename columns
  col_names <- colnames(metrics_tab_lab)
  col_names <- gsub("\\.", "/", col_names)
  col_names <- gsub("_", " ", col_names)
  col_names <- plyr::mapvalues(col_names, 
                               from = c("ASW label",
                                        "ASW label/batch", 
                                        "graph conn", 
                                        "KNN", 
                                        "isolated label silhouette", 
                                        "isolated label F1"),
                               to = c("Cell type ASW",
                                      "Batch ASW", 
                                      "Graph connectivity", 
                                      "kNN accuracy", 
                                      "Isolated label silhouette", 
                                      "Isolated label F1"))
  colnames(metrics_tab_lab) <- col_names
  
  
  
  # metrics names as they are supposed to be ordered
  group_batch <- c("Batch ASW", "PCR batch", "Graph connectivity", "EBM")
  group_bio <- c("NMI cluster/label", "ARI cluster/label", "Cell type ASW", 
                 "Isolated label F1", "Isolated label silhouette", "kNN accuracy")
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
  
  # add method name in filename for tuning
  if(type == "tuning"){
    dataset <- paste0(unique(metrics_tab_lab$method), "_", dataset)
  }
  
  # set all zeros in time = NAs
  metrics_tab_lab$`reference time` <- ifelse(metrics_tab_lab$`reference time` == 0, NA, metrics_tab_lab$`reference time`)
  metrics_tab_lab$`query time` <- ifelse(metrics_tab_lab$`query time` == 0, NA, metrics_tab_lab$`query time`)
  
  
  if(type == "ratio"){
    # get rows where ratio == full
    ind.full <- which(metrics_tab_lab$rqr == "full")
    full_df <- metrics_tab_lab[ind.full, ]
    #rename scArches methods
    full_df$method <- as.character(full_df$method)
    ind.sca <- which(full_df$method %in% c("CVAE (MSE)", "trVAE", "CVAE (NB)", "scVI", "scANVI"))
    full_df$method[ind.sca] <- paste0("scArches ", as.character(full_df$method[ind.sca]))
    
    ####---- Ratio benchmark per method
    ratio_df <- metrics_tab_lab[-ind.full, ]
    
    # rename methods
    ratio_df$method <- paste0("scArches ", ratio_df$method)
    
    # refactor ratio
    ratio_df$rqr <- as.numeric(as.character(ratio_df$rqr))
    ratio_df$rqr <- paste(as.numeric(ratio_df$rqr)*100, "%", sep = "")
    
    models <- unique(ratio_df$method)
    rqr_list <- list()
    best_rqr <- data.frame()
    for(m in models){
      model_df <- ratio_df[which(ratio_df$method == m), ]
      metrics_tab_list <- getMetrics_tab(model_df, type, group_time, add_columns, group_batch, n_metrics_batch_original, group_bio,
                                         n_metrics_bio_original, weight_batch, outdir, dataset, tag)
      g <- getTable(metrics_tab_list$metrics_tab, 
                    metrics_tab_list$n_metrics_batch, 
                    metrics_tab_list$n_metrics_bio,
                    n_metrics_time_original, 
                    metrics_tab_list$column_2_width, 
                    plot_time)
      now <- Sys.time()
      ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), m, "_", dataset,"_", tag, "_summary_metrics.pdf"), g, device = cairo_pdf, width = 210, height = 297, units = "mm")
      ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), m, "_", dataset,"_", tag, "_summary_metrics.tiff"), g, device = "tiff", dpi = "retina", width = 210, height = 297, units = "mm")
      ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), m, "_", dataset,"_", tag, "_summary_metrics.png"), g, device = "png", dpi = "retina", width = 210, height = 297, units = "mm")
      rqr_list[[m]]<- g
      
      # Get best RQR
      if(nrow(best_rqr)==0){
        best_rqr <- model_df[model_df$rqr == metrics_tab_list$metrics_tab[1, "Reference fraction"], ]
      } else{
        best_rqr <- rbind(best_rqr, model_df[model_df$rqr == metrics_tab_list$metrics_tab[1, "Reference fraction"], ])
      }
    }
    
    # plot all models together
    all_g <- plot_grid(plotlist = rqr_list, ncol = 1)
    ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dataset,"_", tag, "_summary_metrics.pdf"), all_g, device = cairo_pdf, width = 210, height = 297, units = "mm")
    
    
    # plot all methods (including full) together
    full_best_df <- rbind(best_rqr, full_df)
    full_metrics_tab_list <- getMetrics_tab(full_best_df, type, group_time, add_columns, group_batch, n_metrics_batch_original, group_bio,
                                            n_metrics_bio_original, weight_batch, outdir, dataset, tag)
    full_g <- getTable(full_metrics_tab_list$metrics_tab, 
                       full_metrics_tab_list$n_metrics_batch, 
                       full_metrics_tab_list$n_metrics_bio,
                       n_metrics_time_original, 
                       full_metrics_tab_list$column_2_width, 
                       plot_time)
    
    now <- Sys.time()
    ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dataset,"_",  "_final_summary_metrics.pdf"), full_g, device = cairo_pdf, width = 210, height = 297, units = "mm")
    ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dataset,"_",  "_final_summary_metrics.tiff"), full_g, device = "tiff", dpi = "retina", width = 210, height = 297, units = "mm")
    ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dataset,"_", "_final_summary_metrics.png"), full_g, device = "png", dpi = "retina", width = 210, height = 297, units = "mm")
    
  }
  
  if(type == "tuning"){
    metrics_tab_list <- getMetrics_tab(metrics_tab_lab, type, group_time, add_columns, group_batch, n_metrics_batch_original, group_bio,
                                            n_metrics_bio_original, weight_batch, outdir, dataset, tag)
    g <- getTable(metrics_tab_list$metrics_tab, 
                  metrics_tab_list$n_metrics_batch, 
                  metrics_tab_list$n_metrics_bio,
                  n_metrics_time_original, 
                  metrics_tab_list$column_2_width, 
                  plot_time)
    
    now <- Sys.time()
    ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dataset,"_", tag, "_summary_metrics.pdf"), g, device = cairo_pdf, width = 210, height = 297, units = "mm")
    ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dataset,"_", tag, "_summary_metrics.tiff"), g, device = "tiff", dpi = "retina", width = 210, height = 297, units = "mm")
    ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dataset,"_", tag, "_summary_metrics.png"), g, device = "png", dpi = "retina", width = 210, height = 297, units = "mm")
    
  }
}
  
  
getMetrics_tab <- function(metrics_tab_lab, type, group_time, add_columns, group_batch, n_metrics_batch_original, group_bio,
                           n_metrics_bio_original, weight_batch, outdir, dataset, tag){
  # info on ratio
  ratio <- as.character(metrics_tab_lab$rqr)
  
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
  metrics_tab <- add_column(metrics_tab, "Overall score" = score_all, .after = "Method")
  metrics_tab <- add_column(metrics_tab, "Batch correction" = score_group_batch, .after = "Overall score")
  metrics_tab <- add_column(metrics_tab, "Bio conservation" = score_group_bio, .after = 3+n_metrics_batch)
  
  if(type == "ratio"){
    metrics_tab <- add_column(metrics_tab, "Reference fraction" = ratio, .after = "Method")
    column_2_width <- 1.5
  } else if(type == "tuning"){
    metrics_tab <- add_column(metrics_tab, "Fine tuning" = ratio, .after = "Method")
    column_2_width <- 7
  }
  
  metrics_tab <- add_column(metrics_tab, "Reference time" = metrics_tab_lab$`reference time`, .after = ncol(metrics_tab))
  metrics_tab <- add_column(metrics_tab, "Query time" = metrics_tab_lab$`query time`, .after = ncol(metrics_tab))
  
  if(type == "ratio"){
    # order methods by the overall score
    metrics_tab <- metrics_tab[order(metrics_tab$`Overall score`,  decreasing = T), ]
  } else if(type == "tuning"){
    metrics_tab <- metrics_tab[match(c("query study labels", "input layers", "all weights"), metrics_tab$`Fine tuning`),]
  }
  
  write.csv(metrics_tab, file = paste0(outdir, "/", dataset, "_", tag, "_summary_scores.csv"), quote = F)
  
  
  # Delete rows that are empty
  rowsNA <- which(is.na(metrics_tab$`Overall score`))
  if(length(rowsNA) >0){
    metrics_tab <- metrics_tab[-rowsNA, ]
  }
  
  return(list(metrics_tab=metrics_tab, n_metrics_batch=n_metrics_batch, n_metrics_bio=n_metrics_bio, column_2_width=column_2_width))
}


getTable <- function(metrics_tab, n_metrics_batch, n_metrics_bio,n_metrics_time_original, column_2_width, plot_time){
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
                            width = c(8,column_2_width,2,2, rep(1,n_metrics_batch), 2, rep(1,n_metrics_bio), rep(2, n_metrics_time_original)),
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
  return(g)
}
