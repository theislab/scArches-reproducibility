
library(ggplot2)

csv_file_path <- "./data/rqr_mouse_brain.csv"
csv_mat <- read.csv(csv_file_path, sep = ",")

dataset <- unique(as.character(csv_mat$data))

df <- data.frame(time_sec = c(csv_mat$reference_time, csv_mat$query_time),
                 ratio = csv_mat$rqr,
                 type = c(rep("reference", nrow(csv_mat)), rep("query", nrow(csv_mat))),
                 ratio_type = c(paste0(csv_mat$rqr, "-reference"), 
                                paste0(csv_mat$rqr, "-query")),
                 method = csv_mat$method)

cell_number <- data.frame(ratio_type = c(paste0(unique(csv_mat$rqr), "-reference"), 
                                    paste0(unique(csv_mat$rqr), "-query")))
cell_number$cells = sample(1:20000, nrow(cell_number))

df2<- merge(df, cell_number)


# Change point shapes and colors
ggplot(df2, aes(x=cells, y=time_sec, shape = ratio_type, color=method)) +
  geom_point(size = 3) +
  theme_bw()
