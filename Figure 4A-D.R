################## Figure 4C
Tour <- subset(our, idents = c(1:17))
Idents(Tour)<- Tour$dis 

Tour <- subset(Tour, idents = c("AD1", "Pso"))
Idents(Tour) <- Tour$ID
Tour$ID <- factor(x=Tour$ID, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17"))


DimPlot(Tour, reduction = "umap", raster=FALSE, cols = c("grey", "navyblue", "grey", "deepskyblue2", "grey", "orange", "grey",
                                           "grey", "seagreen", "grey", "grey", "grey", "grey", "grey", "grey", 
                                           "grey", "grey"))


DimPlot(Tour, reduction = "umap", raster=FALSE, cols = c("grey", "navyblue", "grey", "grey", "grey", "grey", "grey", "grey",
                                           "grey", "yellowgreen", "royalblue2", "grey", "grey", "grey",
                                           "grey", "darkgreen", "grey"))



################ Figure 4AB
Pso_signature_selected <- list(c("ADGRG1",	"ARHGEF12",	"CALR",	"CBLB",	"CCL5",	"CCSER2",	"CD2",	"CEBPD",	"CHN1",	"CLEC2B",	"CPM",	"CTLA4",	"CTSH",	"CXCL13",
                                 "DAPK2",	"EEF1B2",	"FTH1",	"GABPB1-AS1",	"GALNT1",	"GBP5",	"GNLY",	"GZMB",	"H1FX",	"HIST1H1E",	"IL17F",	"IL7R",	"JAML",	"JUN",
                                 "JUND",	"KIAA0319L",	"KLRB1",	"LAG3",	"LAYN",	"LRRN3",	"MAGEH1",	"MAP3K4",	"METRNL",	"MGAT4A",	"MTRNR2L12",	"MUC20-OT1",	"NAP1L4",	"NSA2",
                                 "ODF2L",	"PIK3R1",	"PNRC1",	"PPDPF",	"PTMS",	"PTPN13",	"RAP1B",	"RBPJ",	"SLA2",	"SOX4",	"SPOCK2",	"SRSF7",	"STK17A",	"TBC1D4",
                                 "TCF25",	"TNFAIP3",	"TNFRSF18",	"TNFRSF4",	"TRPS1",	"UBC",	"YPEL2",	"ZEB2"))

AD_signature_selected <- list(c("ACTG1",	"AHNAK",	"ALOX5AP",	"ANXA1",	"ARHGAP21",	"BDP1",	"BIRC3",	"BTG2",	"CAPG",	"CD99",	"CEP350",	"CLNK",	"CLU",	"CSGALNACT1",
                                "CTSW",	"CYSLTR1",	"GOLGA8A",	"HLA-DQA1",	"IFITM2",	"IL17RB",	"IL32",	"ISG15",	"ISG20",	"ITM2A",	"ITM2C",	"KIAA1551",	"KMT2A",	"LGALS1",
                                "LMNA",	"LMO4",	"LPCAT2",	"MFHAS1",	"MMP25",	"NBAS",	"NEAT1",	"PASK",	"PKM",	"PLA2G16",	"PLP2",	"PRKX",	"RARRES3",	"RCBTB2",
                                "RGCC",	"S100A10",	"S100A4",	"SESN3",	"SLC5A3",	"SLFN5",	"SMDT1",	"SOS1",	"SYNE2",	"TAGLN2",	"TGFBR3",	"THADA",	"TIMP1",	"TNFRSF18",
                                "TNFRSF4",	"TNFSF10",	"TSHZ2",	"TTN",	"TWIST1",	"TXNIP"))

Tour <- AddModuleScore(Tour , features = Pso_signature_selected, name = "Pso_signature_selected", search = T)
P1 <- FeaturePlot(Tour , features = "Pso_signature_selected1", order=T, min.cutoff = 0.3, cols = c("black", "orange", "red"), split.by = "dis")


Tour  <- AddModuleScore(Tour , features = AD_signature_selected, name = "AD_signature_selected", search = T)
P2 <- FeaturePlot(Tour, features = "AD_signature_selected1", order=T, min.cutoff = 0.4, cols = c("black", "orange", "red"), split.by = "dis")

plot_grid(P1, P2, ncol=1)

FeaturePlot(Tour, features = "AD_signature_selected1", order=T, min.cutoff = 0.4, cols = c("black", "orange", "red"))

########################## 4c
science <- readRDS("F:/4-7-20m-Final-1030/science.rds")

Tcluster <- subset(science, idents ="2")
Idents(Tcluster) <- Tcluster$Status

Tcluster <- subset(Tcluster, idents=c("Eczema", "Psoriasis"))
Idents(Tcluster) <- Tcluster$ID

Tcluster <- AddModuleScore(Tcluster, features = Pso_signature_selected, name = "Pso_signature_selected", search = T)
  P1 <- FeaturePlot(Tcluster, features = "Pso_signature_selected1", order=T, pt.size = 0.5, min.cutoff = 0.05, cols = c("black", "orange", "red"), split.by = "dis")

Tcluster <- AddModuleScore(Tcluster, features = AD_signature_selected, name = "AD_signature_selected", search = T)
P2 <- FeaturePlot(Tcluster, features = "AD_signature_selected1", order=T, pt.size = 0.5, min.cutoff =0.2, cols = c("black", "orange", "red"), split.by = "dis")

plot_grid(P2, P1, ncol=1)


FeaturePlot(Tcluster, features = "AD_signature_selected1", order=T, pt.size = 0.5, min.cutoff = 0.11, cols = c("black", "orange", "red"))


###############################4D
library(concaveman)
library(ggforce)
library(Seurat)
library(monocle3)
library(dplyr)
library(Rmisc)
library(ggrepel)
#Read in this object

humandat = Tcluster
#Seurat object --> CDS
exp_mat <- humandat@assays[["RNA"]]@data #pull out NORMALIZED counts from Seurat object
cell_metadat <- humandat@meta.data #pull out cell meta-data from Seurat object
gene_annot = data.frame(humandat@assays[["RNA"]]@counts@Dimnames[[1]])#pull out gene names from Seurat object
names(gene_annot) = "gene_short_name"
row.names(gene_annot) = gene_annot$gene_short_name #row.names of gene_metadata must be equal to row.names of expression_data

data.frame(names(cell_metadat))

human_human_big_cds <- new_cell_data_set(exp_mat,
                                         cell_metadata = cell_metadat,
                                         gene_metadata = gene_annot)

#create 1 giant ass cds object from all the cds's created above
rm(humandat,exp_mat)

#Subset to patients of interest
pts = c("E2","E3","E4","E1","P1","P2","P3")
human_human_big_cds = human_human_big_cds[,colData(human_human_big_cds)$donor3 %in% pts] #subset to all normal patients and a diseased patient
colData(human_human_big_cds)$dis_updated = colData(human_human_big_cds)$donor3
colData(human_human_big_cds)$dis_updated = gsub("E2","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("E3","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("E4","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("E1","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("P1","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("P2","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("P3","PV",colData(human_human_big_cds)$dis_updated)

unique(colData(human_human_big_cds)$dis_updated)

#--------------------------------------------------------------------#
#Get the genes of interest
ad_genes <- c("TWIST1",	"LGALS1",	"IL32",	"CAPG",	"ITM2C",	"MFHAS1",	"ANXA1",	"SOS1",	"CSGALNACT1",		"LMO4",	"IFITM2",	"S100A10",	"MT-ND5",	"CYSLTR1",	"PLA2G16",	"SYNE2",	"THADA",	"NEAT1",	"IL17RB",	
              "RPL36A",	"ARHGAP21",	"NBAS",	"ACTG1",	"PRKX",	"TGFBR3",	"TIMP1",	"TNFSF10",	"AHNAK",	"MT-ND2",	"ISG15",	"RPL17",	"LONRF2",	"CD99",	"TSHZ2",	"MMP25",	"IFITM1",	"MT-ND1",	"BIRC3",	"FAM102A",	
              "LPCAT2",	"NRIP3",	"CRIP1",	"CLU",	"PLP2",	"ZFP36",	"ZFP36L2",	"TUBA1B",	"GATA3",	"SLC5A3",	"SFXN1",	"FANK1",	"TAGLN2")
pv_genes <- c("CXCL13",	"MTRNR2L12",	"CD7",	"MGAT4A",	"FTH1",	"LAYN",	"IL17F",	"KLRB1",	"GNLY",	"CPM",	"CTSH",	"GBP5",	"SOX4",	"CLEC2B",	"GZMB",	"CD2",	"CEBPD",	"ODF2L",	"LAG3",	"LRRN3",	
              "ARHGEF12",	"PTPN13",	"TNFAIP3",	"TRPS1",	"SNX9",	"METRNL",	"BTG1",	"JUN",	"SPOCK2",	"GABARAPL1",	"PMEPA1",	"HIST1H1E",	"RBPJ",	"LINC01871",	"MAP3K4",	"H1FX",	"UBC",	"GALNT1",	"PNRC1",	"GABPB1-AS1",	
              "RPS26",	"MUC20-OT1",	"CHN1",	"NAP1L4",	"PTMS",	"F2R",	"CTLA4",	"DAPK2",	"RAP1B",	"CCR6",	"B3GALT2",	"YPEL2",	"FYN",	"PPDPF",	"SLA2",	"CBLB",	"ADGRG1",	"SARAF")
human_human_big_cds@assays@data$counts[1:10,1:10]

#--------------------------------------------------------------------#
#Now subset to cluster 2 cells
human_human_big_cds = human_human_big_cds[,colData(human_human_big_cds)$ID %in% "2"] #subset to all normal patients and a diseased patient

#--------------------------------------------------------------------#
#AD gene data_matrix
ad_mat = t(human_human_big_cds@assays@data$counts[ad_genes,])
ad_mat[1:10,1:10]
ad_int = data.frame(colData(human_human_big_cds)$dis_updated, colData(human_human_big_cds)$donor3, rowSums(ad_mat) )
names(ad_int) = c("dis","sample","gene_sig")
names(ad_int) = c("dis","sample","gene_sig")
ad_dfc = summarySE(ad_int, measurevar='gene_sig', groupvars=c("dis","sample"))
head(ad_dfc)
names(ad_dfc)[c(4,6)] = c("ad_gene_sig","ad_se")

#--------------------------------------------------------------------#
#PV gene data_matrix
pv_mat = t(human_human_big_cds@assays@data$counts[pv_genes,])
pv_mat[1:10,1:10]
pv_int = data.frame(colData(human_human_big_cds)$dis_updated, colData(human_human_big_cds)$donor3, rowSums(pv_mat) )
names(pv_int) = c("dis","sample","gene_sig")
names(pv_int) = c("dis","sample","gene_sig")
pv_dfc = summarySE(pv_int, measurevar='gene_sig', groupvars=c("dis","sample"))
head(pv_dfc)
names(pv_dfc)[c(4,6)] = c("pv_gene_sig","pv_se")

#--------------------------------------------------------------------#
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
re_int = cbind(ad_dfc[,c(1,2,4,6)], pv_dfc[,c(4,6)]) #combine the gene signatures

ggplot(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig)) +
  geom_point(alpha=1,aes(color=dis), size=2) +
  geom_text_repel(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig,label=sample),box.padding = 0.4) + #add labels
  geom_errorbarh(aes(xmax = ad_gene_sig + ad_se, xmin = ad_gene_sig - ad_se, color=dis)) +
  geom_errorbar(aes(ymax = pv_gene_sig + pv_se, ymin = pv_gene_sig - pv_se, color=dis)) +
  #stat_ellipse(aes(color=dis)) +
  theme_classic() +
  ggtitle("Hyperdimensional proximity of Trm1 cells") + 
  xlab("AD-specific genes") +
  ylab("PV-specific genes") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkblue",
                                        "PV" = "darkred")) +
  scale_fill_manual(name="",values = c("AD" = "darkblue",
                                       "PV" = "darkred")) +
  geom_mark_hull(aes(fill=dis, label=dis), concavity=5) +
  theme(axis.text=element_text(size=15,color='black'),
        axis.title=element_text(size=15,color='black'),
        plot.title=element_text(size=15,color='black'))

write.xlsx(re_int, "F:/4-7-20m-Final-1030/Figure 4/re_int.xlsx")
