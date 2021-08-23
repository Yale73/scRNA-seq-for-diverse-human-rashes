############################# get DEGs ####################################
our <- readRDS("E:/Human - 10 clusters/yale-method/000AA-Final objects/our_new_umap_no_218_203.rds")
Idents(our) <- our$donor
our <- RenameIdents(our, `150` = "NML", `154` = "NML", `155` = "NML",`169` = "NML",`195` = "NML", `204` = "NML", `207` = "NML",
                    `170` = "170",`198` = "198", `230` = "230", `231` = "231", `232` = "232", `233` = "233", `236` = "236",
                    `165` = "165", `173` = "173", `194` = "194", `199` = "199", `211` = "211", `222` = "222", `234` = "234", 
                    `235` = "235", `167` = "167", `174` = "174", `192` = "192", `200` = "200", `202` = "202",  `219` = "219", 
                    `141` = "141", `163` = "163", `175` = "175")

our[["donor2"]] <- Idents(object = our)

Idents(our) <- our$ID

donor <- c("167", "174", "192", "200", "202", "219")

my_genes <- read.xlsx("E:/Human - 10 clusters/yale-method/T50-B50-0807/AD-PV-80%-AvsP-p0.05-fc0.25-0809-C2.xlsx", sheetIndex = 1)
Genes <- my_genes$Genes

set1 = list()
for (j in donor){
        tryCatch({
                gc() 
                set1[[paste0(j)]] <- FindMarkers(our, ident.1 = j, ident.2 ="NML", group.by = "donor2", verbose =TRUE, assay = "RNA", 
                                                 slot = "data", subset.ident = "2", test.use = "MAST",
                                                 min.cells.feature = 0, min.cells.group = 0, logfc.threshold = 0, 
                                                 min.pct = 0, min.diff.pct = 0, features = Genes)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.xlsx(set1, file = "E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 6/AN-PN-80pct-merge-AP-C2-TB100.xlsx", rowNames = TRUE, colNames=TRUE)

########################## make figures ###########################
library(xlsx)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

AE1 <- read.xlsx(file = "E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 6/AN-PN-80pct-merge-AP-C2-TB100.xlsx", sheetIndex = 1)
rownames(AE1) <- AE1[,1]# change first column to rownames
AE1[,1] <- NULL# remove the first column
colnames(AE1)#double check colnames

Heatmap(as.matrix(AE1[, 1:7]),
        name = "avg_log2FC", #### annotate legend
        col = colorRamp2(c(-2,-1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2), c("#2D004B","#542788","#8073AC","#B2ABD2","#D8DAEB","#F7F7F7","#FEE0B6","#FDB863","#E08214","#B35806","#7F3B08")), #### set the color scales
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        #cluster_columns = F,
        clustering_distance_columns = "euclidean", #canberra, euclidean
        column_title = "ICR with AD genes",
        show_column_names = T,
        show_row_names = T,
        row_dend_side = "left", column_dend_side = "top",
        column_title_gp = gpar(fontsize = 12))


############################# Pso genes
AE2 <- read.xlsx(file = "E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 6/AN-PN-80pct-merge-AP-C2-TB100.xlsx", sheetIndex = 2)
rownames(AE2) <- AE2[,1]# change first column to rownames
AE2[,1] <- NULL# remove the first column
colnames(AE2)#double check colnames

Heatmap(as.matrix(AE2[, 1:7]),
        name = "avg_log2FC", #### annotate legend
        col = colorRamp2(c(-2,-1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2), c("#7F3B08","#B35806","#E08214","#FDB863","#FEE0B6","#F7F7F7","#D8DAEB","#B2ABD2","#8073AC","#542788","#2D004B")), #### set the color scales
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        #cluster_columns = F,
        clustering_distance_columns = "euclidean",
        column_title = "ICR with Pso genes",
        show_column_names = T,
        show_row_names = T,
        row_dend_side = "left", column_dend_side = "top",
        column_title_gp = gpar(fontsize = 12))

