TClu2 <- read.xlsx("E:/Human - 10 clusters/yale-method/0815-new heatmap/ADvsP-T-AD gene - Copy.xlsx", sheet = 2, colNames = T, rowNames = T)


Heatmap(as.matrix(TClu2[ ,1:17]),
               name = "Rash", #### annotate legend
               col = colorRamp2(c(-1, 0, 1), c("navyblue", "white", "red1")), #### set the color scales
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               clustering_distance_columns = "euclidean",
               column_title = "PV_Rash T",
               show_column_names = T,
               cluster_columns = T,
               show_row_names = T,
               row_dend_side = "left", column_dend_side = "top",
               column_title_gp = gpar(fontsize = 12))



TClu2 <- read.xlsx("E:/Human - 10 clusters/yale-method/0815-new heatmap/ADvsP-T-AD gene - Copy.xlsx", sheet = 1, colNames = T, rowNames = T)


Heatmap(as.matrix(TClu2[ ,1:17]),
               name = "Rash", #### annotate legend
               col = colorRamp2(c(-1, 0, 1), c("navyblue", "white", "red1")), #### set the color scales
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               clustering_distance_columns = "euclidean",
               column_title = "AD_Rash T",
               show_column_names = T,
               cluster_columns = T,
               show_row_names = T,
               row_dend_side = "left", column_dend_side = "top",
               column_title_gp = gpar(fontsize = 12))

