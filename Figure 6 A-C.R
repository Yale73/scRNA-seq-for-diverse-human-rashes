#Script for scRNA-seq classification of external inflammatory skin disease samples

#1. load packages    
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(metap)
library(MAST)
library(stringr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(rhdf5)
library(scDblFinder)
library(tidyverse)
library(concaveman) 
library(ggforce)
library(monocle3)
library(dplyr)
library(Rmisc)
library(ggrepel)
library(raster)
library(vegan)

our <- readRDS("F:/4-7-20m-Final-1030/our_no_218_203.rds")

genes <- c("TWIST1",	"LGALS1",	"IL32",	"CAPG",	"ITM2C",	"MFHAS1",	"ANXA1",	"SOS1",	"CSGALNACT1",		"LMO4",	"IFITM2",	"S100A10",	"MT-ND5",	"CYSLTR1",	"PLA2G16",	"SYNE2",	"THADA",	"NEAT1",	"IL17RB",	
           "RPL36A",	"ARHGAP21",	"NBAS",	"ACTG1",	"PRKX",	"TGFBR3",	"TIMP1",	"TNFSF10",	"AHNAK",	"MT-ND2",	"ISG15",	"RPL17",	"LONRF2",	"CD99",	"TSHZ2",	"MMP25",	"IFITM1",	"MT-ND1",	"BIRC3",	"FAM102A",	
           "LPCAT2",	"NRIP3",	"CRIP1",	"CLU",	"PLP2",	"ZFP36",	"ZFP36L2",	"TUBA1B",	"GATA3",	"SLC5A3",	"SFXN1",	"FANK1",	"TAGLN2",								
           
           "CXCL13",	"MTRNR2L12",	"CD7",	"MGAT4A",	"FTH1",	"LAYN",	"IL17F",	"KLRB1",	"GNLY",	"CPM",	"CTSH",	"GBP5",	"SOX4",	"CLEC2B",	"GZMB",	"CD2",	"CEBPD",	"ODF2L",	"LAG3",	"LRRN3",	
           "ARHGEF12",	"PTPN13",	"TNFAIP3",	"TRPS1",	"SNX9",	"METRNL",	"BTG1",	"JUN",	"SPOCK2",	"GABARAPL1",	"PMEPA1",	"HIST1H1E",	"RBPJ",	"LINC01871",	"MAP3K4",	"H1FX",	"UBC",	"GALNT1",	"PNRC1",	"GABPB1-AS1",	
           "RPS26",	"MUC20-OT1",	"CHN1",	"NAP1L4",	"PTMS",	"F2R",	"CTLA4",	"DAPK2",	"RAP1B",	"CCR6",	"B3GALT2",	"YPEL2",	"FYN",	"PPDPF",	"SLA2",	"CBLB",	"ADGRG1",	"SARAF")

ad_genes <- c("TWIST1",	"LGALS1",	"IL32",	"CAPG",	"ITM2C",	"MFHAS1",	"ANXA1",	"SOS1",	"CSGALNACT1",		"LMO4",	"IFITM2",	"S100A10",	"MT-ND5",	"CYSLTR1",	"PLA2G16",	"SYNE2",	"THADA",	"NEAT1",	"IL17RB",	
              "RPL36A",	"ARHGAP21",	"NBAS",	"ACTG1",	"PRKX",	"TGFBR3",	"TIMP1",	"TNFSF10",	"AHNAK",	"MT-ND2",	"ISG15",	"RPL17",	"LONRF2",	"CD99",	"TSHZ2",	"MMP25",	"IFITM1",	"MT-ND1",	"BIRC3",	"FAM102A",	
              "LPCAT2",	"NRIP3",	"CRIP1",	"CLU",	"PLP2",	"ZFP36",	"ZFP36L2",	"TUBA1B",	"GATA3",	"SLC5A3",	"SFXN1",	"FANK1",	"TAGLN2")
pv_genes <- c("CXCL13",	"MTRNR2L12",	"CD7",	"MGAT4A",	"FTH1",	"LAYN",	"IL17F",	"KLRB1",	"GNLY",	"CPM",	"CTSH",	"GBP5",	"SOX4",	"CLEC2B",	"GZMB",	"CD2",	"CEBPD",	"ODF2L",	"LAG3",	"LRRN3",	
              "ARHGEF12",	"PTPN13",	"TNFAIP3",	"TRPS1",	"SNX9",	"METRNL",	"BTG1",	"JUN",	"SPOCK2",	"GABARAPL1",	"PMEPA1",	"HIST1H1E",	"RBPJ",	"LINC01871",	"MAP3K4",	"H1FX",	"UBC",	"GALNT1",	"PNRC1",	"GABPB1-AS1",	
              "RPS26",	"MUC20-OT1",	"CHN1",	"NAP1L4",	"PTMS",	"F2R",	"CTLA4",	"DAPK2",	"RAP1B",	"CCR6",	"B3GALT2",	"YPEL2",	"FYN",	"PPDPF",	"SLA2",	"CBLB",	"ADGRG1",	"SARAF")

# RashX DEGs
RashX <-  c("167", "174", "192", "200", "202", "219")

set1 = list()
for (j in RashX){
  tryCatch({
    gc() 
    set1[[paste0(j)]] <- FindMarkers(our, ident.1 = j, ident.2 = c("150", "154", "155", "169", "195", "204", "207"), 
                                     group.by = "donor2", verbose =TRUE, assay = "RNA", 
                                     slot = "data", subset.ident = "2", test.use = "MAST",
                                     min.cells.feature = 0, min.cells.group = 0, logfc.threshold = 0, 
                                     min.pct = 0, min.diff.pct = 0, features = genes)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})        
}


# ref DEGs
marker <- FindMarkers(our, ident.1 = "AD1", ident.2 ="Pso", group.by = "dis", verbose = TRUE, assay = "RNA", 
                      slot = "data", test.use = "MAST", min.cells.feature = 0, subset.ident = "2",
                      min.cells.group = 0, logfc.threshold = 0,  min.pct = 0, min.diff.pct = 0, features = genes)

# Add a column to split AD and PV


marker$group <- ifelse(rownames(marker) %in% ad_genes, "AD", "PV")

set1[["ADvsPso"]] <- marker


library(plyr)

for(i in 1:length(set1)){
  colnames(set1[[i]]) <- paste0( names(set1)[i], "_", colnames(set1[[i]]) )
  set1[[i]]$ROWNAMES  <- rownames(set1[[i]])
}

data <- join_all(set1, by="ROWNAMES", type="full" )
rownames(data) <- data$ROWNAMES
data$ROWNAMES <- NULL

RashX1 <- subset(data, ADvsPso_group=="AD")
RashX1 <- RashX1[, which(str_detect(colnames(RashX1), "_log2FC"))]
RashX2 <- subset(data, ADvsPso_group=="PV")
RashX2 <- RashX2[, which(str_detect(colnames(RashX2), "_log2FC"))]
RashX2$ADvsPso_avg_log2FC = RashX2$ADvsPso_avg_log2FC*(-1)

# make heatmap figures
#AD-specific genes heatmap
colnames(RashX1) <- c("E", "C", "B", "D", "F", "A", "ADvsPV")
Heatmap(as.matrix(RashX1),
        name = "avg_log2FC", #### annotate legend
        col = colorRamp2(c(-2,-1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2), c("#2D004B","#542788","#8073AC","#B2ABD2","#D8DAEB","#F7F7F7","#FEE0B6","#FDB863","#E08214","#B35806","#7F3B08")), #### set the color scales
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        #cluster_columns = F,
        clustering_distance_columns = "canberra", #canberra, euclidean
        column_title = "AD-specific genes",
        show_column_names = T,
        show_row_names = T,
        row_dend_side = "left", column_dend_side = "top",
        column_title_gp = gpar(fontsize = 12))
write.xlsx(RashX1, "F:/4-7-20m-Final-1030/Figure 6/double check-1226/ADvsPV_matrix_Final.xlsx", rowNames=T, colNames=T)
#Pso-specific genes heatmap
colnames(RashX2) <- c("E", "C", "B", "D", "F", "A", "PVvsAD")
Heatmap(as.matrix(RashX2),
        name = "avg_log2FC", #### annotate legend
        col = colorRamp2(c(-2,-1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2), c("#7F3B08","#B35806","#E08214","#FDB863","#FEE0B6","#F7F7F7","#D8DAEB","#B2ABD2","#8073AC","#542788","#2D004B")), #### set the color scales
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        #cluster_columns = F,
        clustering_distance_columns = "canberra",
        column_title = "Pso-specific genes",
        show_column_names = T,
        show_row_names = T,
        row_dend_side = "left", column_dend_side = "top",
        column_title_gp = gpar(fontsize = 12))

write.xlsx(RashX2, "F:/4-7-20m-Final-1030/Figure 6/double check-1226/PVvsAD_matrix_Final.xlsx", rowNames=T, colNames=T)

#8. Hyperdimensionality plot and statistical analysis for external query dataset DEGs for  predicted.id = 2 (or Trm1) cells
#########
Idents(our) <- our$donor
our <- RenameIdents(our, `150` = "NML1", `154` = "NML2", `155` = "NML3",`169` = "NML4",`195` = "NML5", `204` = "NML6", 
                    `207` = "NML7", `170` = "AD1",`198` = "AD2", `230` = "AD3", `231` = "AD4", `232` = "AD5", `233` = "AD6",
                    `236` = "AD7",  `165` = "PV1", `173` = "PV2", `194` = "PV3", `199` = "PV4", `211` = "PV5", 
                    `222` = "PV6", `234` = "PV7", `235` = "PV8",  `167` = "E", `174` = "C", `192` = "B", `200` = "D",
                    `202` = "F", `219` = "A")


our[["STATUS"]] <- Idents(object = our)


#Seurat object --> CDS
exp_mat <- our@assays[["RNA"]]@data #pull out NORMALIZED counts from Seurat object
cell_metadat <- our@meta.data #pull out cell meta-data from Seurat object
gene_annot = data.frame(our@assays[["RNA"]]@counts@Dimnames[[1]])#pull out gene names from Seurat object
names(gene_annot) = "gene_short_name"
row.names(gene_annot) = gene_annot$gene_short_name #row.names of gene_metadata must be equal to row.names of expression_data

data.frame(names(cell_metadat))

human_human_big_cds <- new_cell_data_set(exp_mat,
                                         cell_metadata = cell_metadat,
                                         gene_metadata = gene_annot)

#create 1 large object from all the cds's created above

#Subset to patients of interest
pts = c("AD1", "AD2", "AD3", "AD4", "AD5", "AD6", "AD7", "PV1", "PV2","PV3","PV4","PV5","PV6","PV7","PV8",
        "F", "A", "C", "E", "D", "B")

human_human_big_cds = human_human_big_cds[,colData(human_human_big_cds)$STATUS %in% pts & colData(human_human_big_cds)$ID %in% "2"] #subset to all normal patients and a diseased patient
colData(human_human_big_cds)$dis_updated = colData(human_human_big_cds)$donor
unique(colData(human_human_big_cds)$dis_updated)

colData(human_human_big_cds)$dis_updated = gsub("170","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("198","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("230","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("231","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("232","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("233","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("236","AD",colData(human_human_big_cds)$dis_updated)

colData(human_human_big_cds)$dis_updated = gsub("165","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("173","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("194","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("199","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("211","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("222","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("234","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("235","PV",colData(human_human_big_cds)$dis_updated)

colData(human_human_big_cds)$dis_updated = gsub("167","CIR",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("174","CIR",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("192","CIR",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("200","CIR",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("202","CIR",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("219","CIR",colData(human_human_big_cds)$dis_updated)

unique(colData(human_human_big_cds)$dis_updated)
#--------------------------------------------------------------------#

#--------------------------------------------------------------------#
#AD gene data_matrix
ad_mat = t(human_human_big_cds@assays@data$counts[ad_genes,])
ad_mat[1:10,1:10]
ad_int = data.frame(colData(human_human_big_cds)$dis_updated, colData(human_human_big_cds)$donor, rowSums(ad_mat), colData(human_human_big_cds)$STATUS)
names(ad_int) = c("dis","sample","gene_sig", "status")
ad_dfc = summarySE(ad_int, measurevar='gene_sig', groupvars=c("dis","status"))
head(ad_dfc)
names(ad_dfc)[c(4,6)] = c("ad_gene_sig","ad_se")

#--------------------------------------------------------------------#
#PV gene data_matrix
pv_mat = t(human_human_big_cds@assays@data$counts[pv_genes,])
pv_mat[1:10,1:10]
pv_int = data.frame(colData(human_human_big_cds)$dis_updated, colData(human_human_big_cds)$donor, rowSums(pv_mat), colData(human_human_big_cds)$STATUS)
names(pv_int) = c("dis","sample","gene_sig", "status")
pv_dfc = summarySE(pv_int, measurevar='gene_sig', groupvars=c("dis","status"))
head(pv_dfc)
names(pv_dfc)[c(4,6)] = c("pv_gene_sig","pv_se")

#--------------------------------------------------------------------#
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
re_int = cbind(ad_dfc[,c(1,2,4,6)], pv_dfc[,c(4,6)]) #combine the gene signatures
re_int
write.xlsx(re_int, "F:/4-7-20m-Final-1030/Figure 6/double check-1226/hyperdimensioanl_matrix_Final.xlsx", rowNames=T, colNames=T)

ggplot(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig)) +
  geom_point(alpha=1,aes(color=dis), size=2) +
  geom_text_repel(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig,label=status),box.padding = 0.4) + #add labels
  geom_errorbarh(aes(xmax = ad_gene_sig + ad_se, xmin = ad_gene_sig - ad_se, color=dis)) +
  geom_errorbar(aes(ymax = pv_gene_sig + pv_se, ymin = pv_gene_sig - pv_se, color=dis)) +
  #stat_ellipse(aes(color=dis)) +
  theme_classic() +
  ggtitle("CIR") + 
  xlab("AD-specific genes") +
  ylab("PV-specific genes") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkblue",
                                        "PV" = "darkred",
                                        "CIR" = "darkorange")) +
  scale_fill_manual(name="",values = c("AD" = "darkblue",
                                       "PV" = "darkred",
                                       "CIR" = "darkorange")) +
  geom_mark_hull(aes(fill=dis, label=dis), concavity=5) +
  theme(axis.text=element_text(size=15,color='black'),
        axis.title=element_text(size=15,color='black'),
        plot.title=element_text(size=15,color='black'))


# Statistic analysis
#################################vegdist
#pointDistance

all_coords_mat = as.matrix(data.frame(re_int$ad_gene_sig,re_int$pv_gene_sig))
row.names(all_coords_mat) = re_int$sstatus

dist_mat2 =vegdist(all_coords_mat, method = "canberra", diag = FALSE, upper = TRUE, p = 2)
dist_mat2 <- as.matrix(dist_mat2)
rownames(dist_mat2) = re_int$status
colnames(dist_mat2) = re_int$status

ind_samps = subset(re_int, dis=="CIR")$status
ind_samps <-as.character(ind_samps)

statistic = 
  #Loop through the indeterminate samples and get averages 
  
  all_res = do.call(rbind, lapply(1:length(ind_samps), function(i){
    #i=1
    dist_df = data.frame(re_int$dis,re_int$status, dist_mat2[ind_samps[i],]) #get distance of a sample from the matrix and combine with metadata
    names(dist_df) = c("dis","status","dist")
    
    ad_dist_df = subset(dist_df, dis=="AD") #remove other indeterminate samples
    pv_dist_df = subset(dist_df, dis=="PV") #remove other indeterminate samples
    
    #calculate means to see which sided test to use
    ad_mean = mean(ad_dist_df$dist)
    pv_mean = mean(pv_dist_df$dist)
    
    #now test based on mean direction, find which is least
    if(ad_mean > pv_mean) {
      wiltest = wilcox.test(ad_dist_df$dist, pv_dist_df$dist, alternative="greater")
      wilrest = data.frame("PV",ad_mean,pv_mean,ad_mean-pv_mean,wiltest$statistic,wiltest$p.value)
    } else {wiltest = wilcox.test(ad_dist_df$dist, pv_dist_df$dist, alternative="less")
    wilrest = data.frame("AD",ad_mean,pv_mean,ad_mean-pv_mean,wiltest$statistic,wiltest$p.value)
    }
    
    #comine with data and return
    res = data.frame(ind_samps[i],wilrest)
    names(res) = c("sample","proximity","AD_dist_mean","PV_dist_mean","AD_PV","W","p")
    res
  }))


#### visualize the results directly
statistic

###save the statistic results as xlsx format
write.xlsx(statistic, "F:/4-7-20m-Final-1030/Figure 6/double check-1226/statistic_matrix_Final.xlsx", rowNames=T, colNames=T) 
