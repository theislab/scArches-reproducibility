
library(ggplot2)
library(scales)
# csv_file_path <- "./data/mouse_brain.csv"
# n_cells_file_path <- "./data/cell_numbers_brain.csv"

plot_scatterTime <- function(csv_file_path, n_cells_file_path){
  csv_mat <- read.csv(csv_file_path, sep = ",")
  
  # cell numbers
  n_cells_df <- read.csv(n_cells_file_path, sep = ",")
  
  dataset <- unique(as.character(csv_mat$data))
  
  df <- data.frame(time_min = c(csv_mat$reference_time, csv_mat$query_time)/60,
                   ratio = as.character(csv_mat$rqr),
                   type = c(rep("reference", nrow(csv_mat)), rep("query", nrow(csv_mat))),
                   ratio_type = c(paste0(csv_mat$rqr, "_ref"), 
                                  paste0(csv_mat$rqr, "_que")),
                   method = csv_mat$method,
                   stringsAsFactors = FALSE)
  # refactor ratio
  df$ratio[df$ratio != "full"] <- paste0(as.numeric(df$ratio[df$ratio != "full"])*100, "%")
  
  # remove rows with NA or zero as time
  df <- df[!is.na(df$time_min),]
  df <- df[df$time_min != 0,]
  
  df_merged <- merge(df, n_cells_df, by.x = "ratio_type", by.y = "ratio_ds")
  df_merged$type <- factor(df_merged$type, levels = c("reference", "query"))
  

  
  tiff(paste0("scArches_scatterTime_", dataset, ".tiff"), width = 700, height = 600)
  print(ggplot(df_merged, aes(x=n_cells, y=time_min, shape = ratio, color=method)) +
    geom_point(size = 4) +
    facet_wrap(~type, nrow = 2) +
    theme_bw() +
    labs(shape = "Reference fraction", color = "Method", x = "Number of cells", y = "Time (minutes)") +
    theme(text = element_text(size = 20)) +
    scale_x_continuous(labels = comma)
  )
  dev.off()
}



