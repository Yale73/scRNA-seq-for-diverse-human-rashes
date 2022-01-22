library(openxlsx)
library(ComplexHeatmap)
library(circlize)



F3C <- read.xlsx("F:/4-7-20m-Final-1030/Figure 3/AP_80%_T.xlsx", sheet = 13, colNames = T, rowNames = T)
F3D <- read.xlsx("F:/4-7-20m-Final-1030/Figure 3/AP_80%_T.xlsx", sheet = 14, colNames = T, rowNames = T)

AD <- rownames(F3C)
Pso <-rownames(F3D)


Idents(our) <- our$ID

set1 = list()
for (i in 1:17){
  tryCatch({
    print(i)
    gc() 
    set1[[paste0("C",i)]] <- FindMarkers(our, ident.1 = "AD1", ident.2 ="NML", group.by = "dis", verbose =TRUE,
                                         assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                                         subset.ident = i, logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = AD)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.xlsx(set1, file = "F:/4-7-20m-Final-1030/Figure 3/AD_T_selected.xlsx", rowNames = TRUE, colNames=TRUE)

#############
set2 = list()
for (i in 1:17){
  tryCatch({
    print(i)
    gc() 
    set2[[paste0("C",i)]] <- FindMarkers(our, ident.1 = "Pso", ident.2 ="NML", group.by = "dis", verbose =TRUE,
                                         assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                                         subset.ident = i, logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = Pso)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.xlsx(set2, file = "F:/4-7-20m-Final-1030/Figure 3/Pso_T_selected.xlsx", rowNames = TRUE, colNames=TRUE)



################################### heatmap ##############################
F3C <- read.xlsx("F:/4-7-20m-Final-1030/Figure 3/AD_T_selected.xlsx", sheet = 18, colNames = T, rowNames = T)
Heatmap(as.matrix(F3C),
               name = "Rash", #### annotate legend
               col = colorRamp2(c(-1, 0, 1), c("navyblue", "white", "red1")), #### set the color scales
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               clustering_distance_columns = "canberra",#euclidean
               column_title = "AD-specific Genes",
               show_column_names = T,
               cluster_columns = T,
               show_row_names = T,
               row_dend_side = "left", column_dend_side = "top",
               column_title_gp = gpar(fontsize = 12))




F3D <- read.xlsx("F:/4-7-20m-Final-1030/Figure 3/Pso_T_selected.xlsx", sheet = 18, colNames = T, rowNames = T)
Heatmap(as.matrix(F3D),
               name = "Rash", #### annotate legend
               col = colorRamp2(c(-1, 0, 1), c("navyblue", "white", "red1")), #### set the color scales
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               clustering_distance_columns = "canberra",
               column_title = "Pso-specific Genes",
               show_column_names = T,
               cluster_columns = T,
               show_row_names = T,
               row_dend_side = "left", column_dend_side = "top",
               column_title_gp = gpar(fontsize = 12))

