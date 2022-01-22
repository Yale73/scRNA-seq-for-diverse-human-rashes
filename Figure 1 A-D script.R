our <- readRDS("E:/Human - 10 clusters/yale-method/000AA-Final objects/our_new_umap_no_218_203.rds")

############################Figure 1A############################
Idents(our) <- our$ID
our <- RenameIdents(our, `1` = "1 Tcm", `2` = "2 Trm1", `3` = "3 eTreg1", `4` = "4 Trm2", `5` = "5 CTLex", `6` = "6 CTLem",
                    `7` = "7 Tet", `8` = "8 Tmm1", `9` = "9 ILC/NK", `10` = "10 NK", `11` = "11 Trm3", `12` = "12 cmTreg", 
                    `13` = "13 Tmm2", `14` = "14 eTreg2", `15` = "15 Tmm3", `16` = "16 ILC2",  `17` = "17 Tn", 
                    `18` = "18 moDC1", `19` = "19 LC1", `20` = "20 InfMono", `21` = "21 Mac1", `22` = "22 Mono",
                    `23` = "23 LC2", `24` = "24 moDC2", `25` = "25 moDC3", `26` = "26 DC1",  `27` = "27 DC2",
                    `28` = "28 LC3", `29` = "29 Mac2", `30` = "30 Plasma", `31` = "31 B", `32` = "32 Mac3", 
                    `33` = "33 DC3", `34` = "34 migDC", `35` = "35 Mac4", `36` = "36 Mast", `37` = "37 Trm-c",
                    `38` = "38 Treg-c", `39` = "39 Mast-c", `40` = "40 ILC/NK-c", `41` = "41 CTL-c")

our[["Ident1"]] <- Idents(object = our)

Idents(our) <- our$ID
our <- RenameIdents(our, `1` = "Tcm", `2` = "Trm1", `3` = "eTreg1", `4` = "Trm2", `5` = "CTLex", `6` = "CTLem",
                    `7` = "Tet", `8` = "Tmm1", `9` = "ILC/NK", `10` = "NK", `11` = "Trm3", `12` = "cmTreg", 
                    `13` = "Tmm2", `14` = "eTreg2", `15` = "Tmm3", `16` = "ILC2",  `17` = "Tn", 
                    `18` = "moDC1", `19` = "LC1", `20` = "InfMono", `21` = "Mac1", `22` = "Mono",
                    `23` = "LC2", `24` = "moDC2", `25` = "moDC3", `26` = "DC1",  `27` = "DC2",
                    `28` = "LC3", `29` = "Mac2", `30` = "Plasma", `31` = "B", `32` = "Mac3", 
                    `33` = "DC3", `34` = "migDC", `35` = "Mac4", `36` = "Mast", `37` = "Trm-c",
                    `38` = "Treg-c", `39` = "Mast-c", `40` = "ILC/NK-c", `41` = "CTL-c")

our[["Ident2"]] <- Idents(object = our)

our$Ident2 <- factor(our$Ident2, levels = c("Tcm", "Trm1", "eTreg1", "Trm2", "CTLex", "CTLem", "Tet",
                                            "Tmm1", "ILC/NK", "NK", "Trm3", "cmTreg", "Tmm2", "eTreg2",
                                            "Tmm3", "ILC2", "Tn", "moDC1", "LC1", "InfMono", "Mac1",
                                            "Mono", "LC2", "moDC2", "moDC3", "DC1", "DC2", "LC3",
                                            "Mac2", "Plasma", "B", "Mac3", "DC3", "migDC", "Mac4",
                                            "Mast", "Trm-c", "Treg-c", "Mast-c", "ILC/NK-c", "CTL-c"))
library(RColorBrewer)
DimPlot(our, reduction="umap", label=F, pt.size = .1, label.size=4.5, group.by = "Ident1", repel = T)+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(41))+
  theme(text=element_text(family="Arial",size=15)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=15),
        axis.title.y=element_text(colour='black', size=15),
        axis.text=element_text(colour='black',size=15),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=15),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=22,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))


our$ID <- factor(x = our$ID, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                                        "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25",
                                        "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", 
                                        "38", "39", "40", "41"))

DimPlot(our, reduction="umap", label=T, pt.size = .1, label.size=4.5, group.by = "ID")+NoLegend()+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(41))+
  theme(text=element_text(family="Arial",size=15)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=15),
        axis.title.y=element_text(colour='black', size=15),
        axis.text=element_text(colour='black',size=15),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=15),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=22,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))


############################Figure 1B############################
Idents(our) <- our$donor
our <- RenameIdents(our, `150` = "NML", `154` = "NML", `155` = "NML",`169` = "NML",`195` = "NML", `204` = "NML", `207` = "NML",
                    `170` = "AD",`198` = "AD", `230` = "AD", `231` = "AD", `232` = "AD", `233` = "AD", `236` = "AD", 
                    `165` = "Pso", `173` = "Pso", `194` = "Pso", `199` = "Pso", `211` = "Pso", `222` = "Pso", `234` = "Pso", 
                    `235` = "Pso",  `167` = "Miscellaneous", `174` = "Miscellaneous", `192` = "Miscellaneous", `200` = "Miscellaneous", `202` = "Miscellaneous", 
                    `219` = "Miscellaneous", `141` = "Miscellaneous", `163` = "Miscellaneous", `175` = "Miscellaneous")


our[["Rash"]] <- Idents(object = our)
Idents(our) <- our$predicted.celltype

our$Rash <- factor(x = our$Rash, levels = c("NML", "AD", "Pso", "Miscellaneous"))

DimPlot(our, reduction="umap", label=F, pt.size = .1,  label.size=6, group.by = "Rash", split.by = "Rash", ncol = 2)+NoLegend()+
    scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(4))+
  theme(text=element_text(family="Arial",size=15)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=15,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))



############################Figure 1C###############################
##markers.to.plot <- c("CD3D", "CCR7",	"SELL",	"TCF7",	"CD69",	"ITGAE",	"CXCR6",	"CD4",	"IL2RA",	"TIGIT",	"FOXP3",	"CTLA4",	"CD8A",	"CD8B",	"GZMA",	"GZMH",	"PRF1",	"CD44",
##                     "ID2",	"SNHG12", "ZFAS1",	"RORA",	"BATF",	"KLRB1",	"NCAM1",	"KLRD1",	"GNLY",	"FCGR3A",	"RGS1",	"PPP1R15A",	"TNFRSF18",	"CD96",	"TSC22D3",	"DUSP4",	"GZMK",	"PRDM1", "GATA3",	"PTGDR2",
##                     "IL7R", "HLA-DRA",	"HLA-DPB1",	"HLA-DRB1",	"LAMP3",	"BIRC3",	"FSCN1",	"CCL22",	"TYROBP",	"IDO1",	"CD74",	"CSF2RA",	"MYO1G",	"CD207",	"ACOT7",	"CD83",	"PTGS2",	"IL1RN",	"G0S2",	"FCGR2A",
##                     "LYZ",	"CXCL8", "CXCL2",	"TNFAIP8",	"CD68",	"C1QB",	"C1QC",	"CD163",	"CD14",	"S100A8",	"S100A9",	"CD1A",	"CLEC10A",	"FCER1A",	"CD1E",	"FCER1G",	"F13A1",	"IGKC", "JCHAIN",
##                     "IGHG1",	"CD79A",	"MS4A1", "KLF4",	"NR4A1",	"NR4A2",	"CLEC9A",	"CEBPB",	"TPSB2",	"TPSAB1",	"CPA3",	"KIT",	"GATA2", "MKI67", "TOP2A")

markers.to.plot <- c("CD3D", "CCR7","SELL",	"KLF2", "CD69",	"ITGAE",	"CXCR6",	"CD4", "TIGIT",	"FOXP3",	"IL2RA", "CTLA4", "CD8A",	"CD8B",	"GZMB", "PDCD1", "LAG3", "KLRB1", "PRF1", "KLRD1",	"GNLY",	
                     "TNFRSF18",	"PRDM1", "BATF", "TRAT1", "RORA", "GATA3",	"PTGDR2",
                     "IL7R", "HLA-DRA",	"HLA-DRB1","CD83","IDO1",	"CD207", "EPCAM",
                     "CD68",	"C1QB",	"C1QC",	"CD163", "CLEC9A", "CD1C", "CLEC10A",	"THBD", "XCR1", "SIRPA",	"F13A1", "IGKC", "JCHAIN",
                     "CD79A",	"MS4A1", "NR4A1",	"NR4A2", "KLF4", "CEBPB","LYZ", "MS4A7", "SERPINA1", "CD14", "S100A9", "IL23A", "TPSB2", "TPSAB1", "MKI67", "TOP2A")

Idents(our) <- our$Ident2
library(viridis)
library(scales)
DotPlot(our, features = markers.to.plot, cols=c("#5F4B8BFF", "#ED2B33FF"), assay = "RNA", col.min = 0.3, col.max = 0.8, dot.min=0.12, dot.scale = 1,
        cluster.idents=F)+
  scale_size(range = c(0, 5))+ 
  scale_size_area(max_size = 5)+ 
  #scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8, family="TT Times New Roman"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 8.5, family="TT Times New Roman"),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9)) +
  scale_color_gradientn(colours = viridis::magma(20), limits = c(0,1), oob = scales::squish, name = 'log2 (count + 1)')


############################Figure 1D############################
library(openxlsx)
ADT <- read.xlsx("F:/Human - 10 clusters/yale-method/0624-ADT cutoff/1.5SD+Mean-ADT process/ADT_cutoff_new_integration.xlsx", sheet = 1, colNames = T, rowNames = T)


breaksList = as.numeric(seq(0, 1, by = .1))
library(ComplexHeatmap)
### heatmap
Heatmap(ADT, col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(breaksList)),
        name="Exp%",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        clustering_distance_columns = "euclidean",
        cluster_rows = T,
        show_column_names = T,
        show_row_names = T,
        row_dend_side = "left", column_dend_side = "top",
        column_title_gp = gpar(fontsize = 12))

