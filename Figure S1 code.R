library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)

DataMix <- readRDS("F:/4-7-20m-Final-1030/Figure_S1AB_object.rds")
our <- subset(DataMix, source=="Our")
rm(DataMix)
DataMix <- our
rm(our)
################### Figure S1A
DimPlot(DataMix, reduction="umap", label=F, pt.size = .1,  label.size=4.5)+
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(26))+
  theme(text=element_text(family="Arial",size=15)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour  = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=15),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=22,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))


################### Figure S1B
p1 <- FeaturePlot(DataMix, reduction="umap", features = "CD3E", pt.size=1,
            min.cutoff =0, max.cutoff = 5)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family= "italic",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family= "italic", size=12),
        legend.key=element_blank())

p2 <- FeaturePlot(DataMix, reduction="umap", features = "KLRB1", pt.size=1,
            min.cutoff =0, max.cutoff = 5)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())

p3 <- FeaturePlot(DataMix, reduction="umap", features = "MS4A1", pt.size=1,
                  min.cutoff =0, max.cutoff = 5)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())

p4 <- FeaturePlot(DataMix, reduction="umap", features = "HLA-DRA", pt.size=1,
            min.cutoff =1, max.cutoff = 5)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())

p5 <- FeaturePlot(DataMix, reduction="umap", features = "TPSAB1", pt.size=1,
                  min.cutoff =0, max.cutoff = 5)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())

p6 <- FeaturePlot(DataMix, reduction="umap", features = "KRT1", pt.size=1,
                  min.cutoff =0, max.cutoff = 5)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())

p7 <- FeaturePlot(DataMix, reduction="umap", features = "COL6A1", pt.size=1,
                  min.cutoff =0, max.cutoff = 5)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())

p8 <- FeaturePlot(DataMix, reduction="umap", features = "MKI67", pt.size=1,
                  min.cutoff =0, max.cutoff = 5)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())


DataMix[["mouse"]] <- PercentageFeatureSet(DataMix, pattern = "^mm10---")


p9 <- FeaturePlot(DataMix, reduction="umap", features = "mouse", pt.size=1,
                  min.cutoff =50, order=T)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())



plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=3)


##########################Figure S1D
our <- readRDS("F:/4-7-20m-Final-1030/our_with_218_203.rds")

Idents(our) <- our$donor
our <- subset(our, idents=c("218", "203"), invert=T)

Idents(our) <- our$donor
our <- RenameIdents(our, `150` = "NML1", `154` = "NML2", `155` = "NML3",`169` = "NML4",`195` = "NML5", `204` = "NML6", `207` = "NML7",
                    `170` = "AD1",`198` = "AD2", `230` = "AD3", `231` = "AD4", `232` = "AD5", `233` = "AD6", `236` = "AD7", 
                    `165` = "Pso1", `173` = "Pso2", `194` = "Pso3", `199` = "Pso4", `211` = "Pso5", `222` = "Pso6", `234` = "Pso7", 
                    `235` = "Pso8",  `167` = "ICR1", `174` = "ICR2", `192` = "ICR3", `200` = "ICR4", `202` = "ICR5", `219` = "ICR7", 
                    `141` = "LP1", `163` = "LP2", `175` = "BP1")

our[["STATUS"]] <- Idents(object = our)
Idents(our) <- our$ID

DimPlot(our, reduction="umap", label=F, pt.size = .1,  label.size=6, group.by = "STATUS")+
  scale_color_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(33))+
  theme(text=element_text(family="Arial",size=15)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=10),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))


saveRDS(our, "F:/4-7-20m-Final-1030/our_no_218_203.rds")

################################## Figure S1E ##################################
DataMix <- DataMix <- readRDS("F:/4-7-20m-Final-1030/three.rds")

Idents(DataMix) <- DataMix$source
DataMix <- RenameIdents(object = DataMix, `Our` = "Manuscript Data", `S` = "Reynolds et al.", `SI` = "Bangert et al.")
DataMix[["source2"]] <- Idents(object = DataMix)

Idents(DataMix) <- DataMix$ID

DimPlot(DataMix, reduction="umap", label=F, pt.size = .1, label.size=4.5, split.by =  "source2")+NoLegend()+
  scale_color_manual(values =  colorRampPalette(brewer.pal(12, "Paired"))(41))+
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
        legend.text=element_text(family="Arial", size=15),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=22,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))

################################## Figure S1F ##################################
DimPlot(DataMix, reduction="umap", label=T, pt.size = .1, label.size=8, group.by = "final_clustering", split.by =  "source")+NoLegend()+
  scale_color_manual(values =  colorRampPalette(brewer.pal(12, "Paired"))(41))+
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
        legend.text=element_text(family="Arial", size=15),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=22,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))
