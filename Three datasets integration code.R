################################ our Skin data
Skin <- readRDS("E:/Human - 10 clusters/yale-method/0218-three object/Skin.rds")
ScienceImmunol <- readRDS("E:/Human - 10 clusters/yale-method/SciImmunol-data analysis/new integration-0304/Science_Immunol.rds")
Science <- readRDS("E:/Human - 10 clusters/yale-method/0218-three object/Science.rds")

############## remove doublets
Idents(Skin) <- Skin$scDblFinder.class
Skin <- subset(Skin, idents = "singlet")

Idents(ScienceImmunol) <- ScienceImmunol$scDblFinder.class
ScienceImmunol <- subset(ScienceImmunol, idents = "singlet")
###############
DataMix <- merge(x = Skin, y = c(Science, ScienceImmunol), project = "DataMix")

rm(Skin, Science, ScienceImmunol)

##### run pca for harmony
DefaultAssay(DataMix) <- 'RNA'
DataMix <- NormalizeData(DataMix, verbose = T, normalization.method = "LogNormalize") %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% RunPCA()

DefaultAssay(DataMix) <- 'ADT'
VariableFeatures(DataMix) <- rownames(DataMix[["ADT"]])
DataMix <- NormalizeData(DataMix, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

#### visualize before harmony
DefaultAssay(DataMix) <- 'RNA'
options(repr.plot.height = 12, repr.plot.width = 30)
p1 <- DimPlot(object = DataMix, reduction = "pca", pt.size = .1, group.by = "donor")+NoLegend()
p2 <- VlnPlot(object = DataMix, features = "PC_1", group.by = "donor", pt.size = .1)+NoLegend()
plot_grid(p1,p2)+NoLegend()

###########run harmony
library(harmony)
options(repr.plot.height = 10, repr.plot.width = 15)
DataMix <- DataMix %>% 
  RunHarmony(c("donor", "chem", "source", "version"), plot_convergence = TRUE, max.iter.harmony = 25)


options(repr.plot.height = 12, repr.plot.width = 30)
p1 <- DimPlot(object = DataMix, reduction = "harmony", pt.size = .1, group.by = "donor")+NoLegend()
p2 <- VlnPlot(object = DataMix, features = "harmony_1", group.by = "donor", pt.size = .1)+NoLegend()
plot_grid(p1,p2)+NoLegend()

options(repr.plot.height = 6, repr.plot.width = 15)
ElbowPlot(DataMix, ndims = 50)

##### check dims for further dim choose
DefaultAssay(DataMix) <- 'RNA'
DimHeatmap(DataMix, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(DataMix, dims = 16:30, cells = 500, balanced = TRUE)
DimHeatmap(DataMix, dims = 31:45, cells = 500, balanced = TRUE)

########Decide the PC numbers
# Determine percent of variation associated with each PC
pct <- DataMix[["harmony"]]@stdev / sum(DataMix[["harmony"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
###42

####The first metric returns PC42 as the PC matching these 
##requirements. Let's check the second metric, which identifies
##the PC where the percent change in variation between 
##consecutive PCs is less than 0.1%:
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
#31

##This second metric returns PC31. Usually, we would choose the minimum of 
###these two metrics as the PCs covering the majority of the variation in 
##the data.
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
#31

##Based on these metrics, for the clustering of cells in Seurat we 
##will use the first 31 PCs to generate the clusters. We 
##can plot the elbow plot again and overlay the information 
##determined using our metrics:
# Create a dataframe with values
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


####downstream analysis
DataMix <- DataMix %>% 
  RunUMAP(reduction = "harmony", dims = 1:31) %>% 
  RunTSNE(reduction = "harmony", dims = 1:31) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:31) 

#####FindCluster for 3 object
DataMix <- DataMix %>%
  FindClusters(resolution = 0.1, algorithm = 1)%>%
  FindClusters(resolution = 0.2, algorithm = 1)%>%
  FindClusters(resolution = 0.3, algorithm = 1)%>%
  FindClusters(resolution = 0.4, algorithm = 1)%>%
  FindClusters(resolution = 0.5, algorithm = 1)%>%
  FindClusters(resolution = 0.6, algorithm = 1)%>%
  FindClusters(resolution = 0.7, algorithm = 1)%>%
  FindClusters(resolution = 0.8, algorithm = 1)%>%
  FindClusters(resolution = 0.9, algorithm = 1)%>%
  FindClusters(resolution = 1, algorithm = 1)%>%
  FindClusters(resolution = 1.1, algorithm = 1)%>%
  FindClusters(resolution = 1.2, algorithm = 1)%>%
  FindClusters(resolution = 1.3, algorithm = 1)%>%
  FindClusters(resolution = 1.4, algorithm = 1)%>%
  FindClusters(resolution = 1.5, algorithm = 1)%>%
  FindClusters(resolution = 1.6, algorithm = 1)%>%
  FindClusters(resolution = 1.7, algorithm = 1)%>%
  FindClusters(resolution = 1.8, algorithm = 1)%>%
  FindClusters(resolution = 1.9, algorithm = 1)%>%
  FindClusters(resolution = 2.0, algorithm = 1)

library(clustree)
clustree(DataMix)

DataMix <- DataMix %>%
  FindClusters(resolution = 0.2, algorithm = 1)

DimPlot(DataMix, reduction="umap", label=T, pt.size = .5, 
        label.size=5)+
  NoLegend()+
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(43))+
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

saveRDS(DataMix, "E:/Human - 10 clusters/yale-method/0304-new three datasets/DataMix_Res01.rds")

########run cluster markers
marker_i <- FindAllMarkers(DataMix_Res01, verbose = TRUE, assay = "RNA", slot = "data",only.pos = T, logfc.threshold = 0.25)
write.xlsx(marker_i, file = "E:/Human - 10 clusters/yale-method/0304-new three datasets/FindAllMarker-13Clusters.xlsx")

