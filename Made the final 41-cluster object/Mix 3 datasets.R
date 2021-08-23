################################ our Skin data
Skin <- readRDS("E:/Human - 10 clusters/yale-method/0208-2021-Allsamples/Skin_All_Dis.rds")

Idents(Skin) <- Skin$donor
Skin <- RenameIdents(object = Skin, `150` = "Our", `154` = "Our", `155` = "Our",`169` = "Our",`195` = "Our", `204` = "Our", `207` = "Our", `170` = "Our",`198` = "Our", `230` = "Our",`231` = "Our",
                     `232` = "Our", `233` = "Our", `236` = "Our", `165` = "Our", `173` = "Our", `194` = "Our", `199` = "Our", `211` = "Our", `222` = "Our", `234` = "Our", `235` = "Our", 
                     `167` = "Our", `174` = "Our", `192` = "Our", `200` = "Our", `202` = "Our", `218` = "Our", `219` = "Our", `141` = "Our", `163` = "Our", `175` = "Our", `203` = "Our")
Skin[["source"]] <- Idents(object = Skin)

###add cellranger version
Idents(Skin) <- Skin$donor
Skin <- RenameIdents(object = Skin, `150` = "v3.0.2", `154` = "v3.0.2", `155` = "v3.0.2",`169` = "v3.0.2",`195` = "v3.0.2", `204` = "v3.0.2", `207` = "v3.0.2", `170` = "v3.0.2",`198` = "v3.0.2",
                     `230` = "v3.0.2",`231` = "v3.0.2", `232` = "v3.0.2", `233` = "v3.0.2", `236` = "v3.0.2", `165` = "v3.0.2", `173` = "v3.0.2", `194` = "v3.0.2", `199` = "v3.0.2", `211` = "v3.0.2",
                     `222` = "v3.0.2", `234` = "v3.0.2", `235` = "v3.0.2", `167` = "v3.0.2", `174` = "v3.0.2", `192` = "v3.0.2",  `200` = "v3.0.2", `202` = "v3.0.2", `218` = "v3.0.2", `219` = "v3.0.2",
                     `141` = "v3.0.2", `163` = "v3.0.2", `175` = "v3.0.2", `203` = "v3.0.2")
Skin[["version"]] <- Idents(object = Skin)
Idents(Skin) <- Skin$seurat_clusters

saveRDS(Skin, "E:/Human - 10 clusters/yale-method/02182021-three object/Skin.rds")
################################ Science Immunology data

ScienceImmunol <- readRDS("E:/Human - 10 clusters/yale-method/SciImmunol-data analysis/SciImmunol_harmony.rds")
###Change Status to dis
Idents(ScienceImmunol) <- ScienceImmunol$status
ScienceImmunol <- RenameIdents(object = ScienceImmunol, `HC` = "N", `AD` = "AD", `D16w` = "D16w",`D1y` = "D1y")
ScienceImmunol[["dis"]] <- Idents(object = ScienceImmunol)

###Change chemistry to chem
Idents(ScienceImmunol) <- ScienceImmunol$donor
ScienceImmunol <- RenameIdents(object = ScienceImmunol, `HC1` = "V2", `HC2` = "V2", `HC3` = "V2",`HC4` = "V2", `HC5` = "V3", `AD1` = "V2", `AD2` = "V2", `AD3` = "V2", `AD4` = "V3",
                               `AD13` = "V3", `AD10` = "V3", `AD11` = "V3", `AD12` = "V3", `AD17` = "V3", `AD18` = "V3", `AD14` = "V3", `AD15` = "V3", `AD16` = "V3", `AD19` = "V3")
ScienceImmunol[["chem"]] <- Idents(object = ScienceImmunol)


Idents(ScienceImmunol) <- ScienceImmunol$donor
ScienceImmunol <- RenameIdents(object = ScienceImmunol, `HC1` = "SI", `HC2` = "SI", `HC3` = "SI",`HC4` = "SI", `HC5` = "SI", `AD1` = "SI", `AD2` = "SI", `AD3` = "SI", `AD4` = "SI",
                               `AD13` = "SI", `AD10` = "SI", `AD11` = "SI", `AD12` = "SI", `AD17` = "SI", `AD18` = "SI", `AD14` = "SI", `AD15` = "SI", `AD16` = "SI", `AD19` = "SI")

ScienceImmunol[["source"]] <- Idents(object = ScienceImmunol)


###add cellranger version
Idents(ScienceImmunol) <- ScienceImmunol$donor
ScienceImmunol <- RenameIdents(object = ScienceImmunol, `HC1` = "v3.0.2", `HC2` = "v3.0.2", `HC3` = "v3.0.2",`HC4` = "v3.0.2", `HC5` = "v3.0.2", `AD1` = "v3.0.2", `AD2` = "v3.0.2", 
                               `AD3` = "v3.0.2", `AD4` = "v3.0.2", `AD13` = "v3.0.2", `AD10` = "v3.0.2", `AD11` = "v3.0.2", `AD12` = "v3.0.2", `AD17` = "v3.0.2", `AD18` = "v3.0.2",
                               `AD14` = "v3.0.2", `AD15` = "v3.0.2", `AD16` = "v3.0.2", `AD19` = "v3.0.2")

ScienceImmunol[["version"]] <- Idents(object = ScienceImmunol)

Idents(ScienceImmunol) <- ScienceImmunol$seurat_clusters

saveRDS(ScienceImmunol, "E:/Human - 10 clusters/yale-method/02182021-three object/ScienceImmunol.rds")
###############
Idents(Science) <- Science$donor
Science <- subset(Science, downsample=5000)


skin.list <- SplitObject(Skin, split.by = "donor")
rm(Skin)
skin.list2 <- SplitObject(ScienceImmunol, split.by = "donor")
rm(ScienceImmunol)
skin.list3 <- SplitObject(Science, split.by = "donor")
rm(Science)


DataMix <- merge(x = skin.list$`150`, y = c(skin.list$`154`, skin.list$`155`, skin.list$`169`, skin.list$`195`, skin.list$`204`, skin.list$`207`, skin.list$`170`, skin.list$`198`,
                                            skin.list$`230`, skin.list$`231`, skin.list$`232`, skin.list$`233`, skin.list$`236`, skin.list$`165`, skin.list$`173`, skin.list$`194`,
                                            skin.list$`199`, skin.list$`211`, skin.list$`222`, skin.list$`234`, skin.list$`235`, skin.list$`167`, skin.list$`174`, skin.list$`192`,
                                            skin.list$`200`, skin.list$`202`, skin.list$`218`, skin.list$`219`, skin.list$`141`, skin.list$`163`, skin.list$`175`, skin.list$`203`,
                                            skin.list2$`HC1`, skin.list2$`HC2`, skin.list2$`HC3`, skin.list2$`HC4`, skin.list2$`HC5`, skin.list2$`AD1`, skin.list2$`AD2`, skin.list2$`AD3`,
                                            skin.list2$`AD4`, skin.list2$`AD13`, skin.list2$`AD10`, skin.list2$`AD11`, skin.list2$`AD12`,skin.list2$`AD17`, skin.list2$`AD18`, skin.list2$`AD14`,
                                            skin.list2$`AD15`,skin.list2$`AD16`, skin.list2$`AD19`, skin.list3$`S1`, skin.list3$`S2`, skin.list3$`S3`, skin.list3$`S4`, skin.list3$`S5`,
                                            skin.list3$`E1`, skin.list3$`E2`, skin.list3$`E3`, skin.list3$`E4`, skin.list3$`P1`, skin.list3$`P2`, skin.list3$`P3`), project = "DataMix")


rm(skin.list, skin.list2, skin.list3)
saveRDS(DataMix, "E:/Human - 10 clusters/yale-method/0208-2021-Allsamples/DataMix_only_Merge.rds")

##### run pca for harmony
DefaultAssay(DataMix) <- 'RNA'
DataMix <- NormalizeData(DataMix, verbose = T, normalization.method = "LogNormalize") %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData() %>% RunPCA()

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
  RunHarmony(c("donor", "chem", "source", "version"), plot_convergence = TRUE, max.iter.harmony = 20)


options(repr.plot.height = 12, repr.plot.width = 30)
p1 <- DimPlot(object = DataMix, reduction = "harmony", pt.size = .1, group.by = "donor")+NoLegend()
p2 <- VlnPlot(object = DataMix, features = "harmony_1", group.by = "donor", pt.size = .1)+NoLegend()
plot_grid(p1,p2)+NoLegend()

options(repr.plot.height = 6, repr.plot.width = 15)
ElbowPlot(DataMix, ndims = 50)

##### check dims for further dim choose
DefaultAssay(DataMix) <- 'RNA'
DimHeatmap(DataMix, dims = 1:15, cells = 500, balanced = TRUE, ncol=5)
DimHeatmap(DataMix, dims = 15:30, cells = 500, balanced = TRUE, ncol=5)
DimHeatmap(DataMix, dims = 31:45, cells = 500, balanced = TRUE, ncol=5)

########Decide the PC numbers
# Determine percent of variation associated with each PC
pct <- DataMix[["harmony"]]@stdev / sum(DataMix[["harmony"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1

####The first metric returns PC43 as the PC matching these 
##requirements. Let's check the second metric, which identifies
##the PC where the percent change in variation between 
##consecutive PCs is less than 0.1%:
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
#46
##This second metric returns PC12. Usually, we would choose 
##the minimum of these two metrics as the PCs covering the 
##majority of the variation in the data.
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
#12

##Based on these metrics, for the clustering of cells in Seurat we 
##will use the first fourteen PCs to generate the clusters. We 
##can plot the elbow plot again and overlay the information 
##determined using our metrics:
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


####downstream analysis
DataMix <- DataMix %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) 

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
  FindClusters(resolution = 1.2, algorithm = 1)

DimPlot(DataMix, reduction="umap", label=F, pt.size = .5, 
        label.size=5, split.by="source")+
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

saveRDS(DataMix, "E:/Human - 10 clusters/yale-method/02182021-three object/DataMix_Res04.rds")

########run cluster markers
marker_i <- FindAllMarkers(DataMix, verbose = TRUE, assay = "RNA", slot = "data",only.pos = T, logfc.threshold = 0.25)
write.xlsx(marker_i, file = "E:/Human - 10 clusters/yale-method/02182021-three object//FindAllMarker-43Clusters.xlsx")

