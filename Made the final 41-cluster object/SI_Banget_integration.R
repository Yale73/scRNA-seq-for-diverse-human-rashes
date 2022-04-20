library('Seurat')
library(dplyr)
library(Matrix)
library(cowplot)
library(data.table)
library(magrittr)
library(tidyverse)
library(RCurl)
library("readxl")
library(purrr)
library(pheatmap)
library(MAST)
library(RColorBrewer)
library(grid)
library(gtable)
library(scds)
library(scater)
library(rsvd)
library(Rtsne)
library(cowplot)
library(SingleCellExperiment)
library(scDblFinder)
library(ggplot2)

####### Normal samples
HC1.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/HC1")
HC2.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/HC2")
HC3.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/HC3")
HC4.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/HC4")
HC5.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/HC5")

####AD untreated
AD1.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD1")
AD2.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD2")
AD3.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD3")
AD4.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD4")
AD13.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD13")

####AD16W
AD10.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD10")
AD11.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD11")
AD12.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD12")
AD17.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD17")
AD18.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD18")

####AD1Y
AD14.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD14")
AD15.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD15")
AD16.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD16")
AD19.data <- Read10X(data.dir = "E:/Human - 10 clusters/yale-method/SciImmunol-data matrix/AD19")

### normal 
HC1 <- CreateSeuratObject(counts = HC1.data, project = "HC1")
HC2 <- CreateSeuratObject(counts = HC2.data, project = "HC2")
HC3 <- CreateSeuratObject(counts = HC3.data, project = "HC3")
HC4 <- CreateSeuratObject(counts = HC4.data, project = "HC4")
HC5 <- CreateSeuratObject(counts = HC5.data, project = "HC5")

####AD untreated
AD1 <- CreateSeuratObject(counts = AD1.data, project = "AD1")
AD2 <- CreateSeuratObject(counts = AD2.data, project = "AD2")
AD3 <- CreateSeuratObject(counts = AD3.data, project = "AD3")
AD4 <- CreateSeuratObject(counts = AD4.data, project = "AD4")
AD13 <- CreateSeuratObject(counts = AD13.data, project = "AD13")

####AD16W
AD10 <- CreateSeuratObject(counts = AD10.data, project = "AD10")
AD11 <- CreateSeuratObject(counts = AD11.data, project = "AD11")
AD12 <- CreateSeuratObject(counts = AD12.data, project = "AD12")
AD17 <- CreateSeuratObject(counts = AD17.data, project = "AD17")
AD18 <- CreateSeuratObject(counts = AD18.data, project = "AD18")

####AD16W
AD14 <- CreateSeuratObject(counts = AD14.data, project = "AD14")
AD15 <- CreateSeuratObject(counts = AD15.data, project = "AD15")
AD16 <- CreateSeuratObject(counts = AD16.data, project = "AD16")
AD19 <- CreateSeuratObject(counts = AD19.data, project = "AD19")

### normal
table(Idents(HC1))
table(Idents(HC2))
table(Idents(HC3))
table(Idents(HC4))
table(Idents(HC5))

### AD
table(Idents(AD1))
table(Idents(AD2))
table(Idents(AD3))
table(Idents(AD4))
table(Idents(AD13))

####AD16W
table(Idents(AD10))
table(Idents(AD11))
table(Idents(AD12))
table(Idents(AD17))
table(Idents(AD18))

####AD1y
table(Idents(AD14))
table(Idents(AD15))
table(Idents(AD16))
table(Idents(AD19))


### NML samples
HC1.sce <- as.SingleCellExperiment(HC1)
HC2.sce <- as.SingleCellExperiment(HC2)
HC3.sce <- as.SingleCellExperiment(HC3)
HC4.sce <- as.SingleCellExperiment(HC4)
HC5.sce <- as.SingleCellExperiment(HC5)

####AD untreated
AD1.sce <- as.SingleCellExperiment(AD1)
AD2.sce <- as.SingleCellExperiment(AD2)
AD3.sce <- as.SingleCellExperiment(AD3)
AD4.sce <- as.SingleCellExperiment(AD4)
AD13.sce <- as.SingleCellExperiment(AD13)

####AD16w
AD10.sce <- as.SingleCellExperiment(AD10)
AD11.sce <- as.SingleCellExperiment(AD11)
AD12.sce <- as.SingleCellExperiment(AD12)
AD17.sce <- as.SingleCellExperiment(AD17)
AD18.sce <- as.SingleCellExperiment(AD18)

####AD1y
AD14.sce <- as.SingleCellExperiment(AD14)
AD15.sce <- as.SingleCellExperiment(AD15)
AD16.sce <- as.SingleCellExperiment(AD16)
AD19.sce <- as.SingleCellExperiment(AD19)


##############label doublets
### NML samples
HC1.sce <- scDblFinder(HC1.sce)
HC2.sce <- scDblFinder(HC2.sce)
HC3.sce <- scDblFinder(HC3.sce)
HC4.sce <- scDblFinder(HC4.sce)
HC5.sce <- scDblFinder(HC5.sce)

####AD untreated
AD1.sce <- scDblFinder(AD1.sce)
AD2.sce <- scDblFinder(AD2.sce)
AD3.sce <- scDblFinder(AD3.sce)
AD4.sce <- scDblFinder(AD4.sce)
AD13.sce <- scDblFinder(AD13.sce)

####AD16w
AD10.sce <- scDblFinder(AD10.sce)
AD11.sce <- scDblFinder(AD11.sce)
AD12.sce <- scDblFinder(AD12.sce)
AD17.sce <- scDblFinder(AD17.sce)
AD18.sce <- scDblFinder(AD18.sce)

####AD1y
AD14.sce <- scDblFinder(AD14.sce)
AD15.sce <- scDblFinder(AD15.sce)
AD16.sce <- scDblFinder(AD16.sce)
AD19.sce <- scDblFinder(AD19.sce)


##############check doublets percentage
### NML samples
table(call=HC1.sce$scDblFinder.class)
table(call=HC2.sce$scDblFinder.class)
table(call=HC3.sce$scDblFinder.class)
table(call=HC4.sce$scDblFinder.class)
table(call=HC5.sce$scDblFinder.class)

####AD untreated
table(call=AD1.sce$scDblFinder.class)
table(call=AD2.sce$scDblFinder.class)
table(call=AD3.sce$scDblFinder.class)
table(call=AD4.sce$scDblFinder.class)
table(call=AD13.sce$scDblFinder.class)

####AD16w
table(call=AD10.sce$scDblFinder.class)
table(call=AD11.sce$scDblFinder.class)
table(call=AD12.sce$scDblFinder.class)
table(call=AD17.sce$scDblFinder.class)
table(call=AD18.sce$scDblFinder.class)

####AD1y
table(call=AD14.sce$scDblFinder.class)
table(call=AD15.sce$scDblFinder.class)
table(call=AD16.sce$scDblFinder.class)
table(call=AD19.sce$scDblFinder.class)


###############
HC1 <- as.Seurat(HC1.sce, project = "HC1")
HC2 <- as.Seurat(HC1.sce, project = "HC2")
HC3 <- as.Seurat(HC1.sce, project = "HC3")
HC4 <- as.Seurat(HC1.sce, project = "HC4")
HC5 <- as.Seurat(HC1.sce, project = "HC5")

AD1 <- as.Seurat(AD1.sce, project = "AD1")
AD2 <- as.Seurat(AD2.sce, project = "AD2")
AD3 <- as.Seurat(AD3.sce, project = "AD3")
AD4 <- as.Seurat(AD4.sce, project = "AD4")
AD13 <- as.Seurat(AD13.sce, project = "AD13")

AD10 <- as.Seurat(AD10.sce, project = "AD10")
AD11 <- as.Seurat(AD11.sce, project = "AD11")
AD12 <- as.Seurat(AD12.sce, project = "AD12")
AD17 <- as.Seurat(AD17.sce, project = "AD17")
AD18 <- as.Seurat(AD18.sce, project = "AD18")

AD14 <- as.Seurat(AD14.sce, project = "AD14")
AD15 <- as.Seurat(AD15.sce, project = "AD15")
AD16 <- as.Seurat(AD16.sce, project = "AD16")
AD19 <- as.Seurat(AD19.sce, project = "AD19")

### normal
HC1[["percent.mt"]] <- PercentageFeatureSet(HC1, pattern = "^MT-")
HC2[["percent.mt"]] <- PercentageFeatureSet(HC2, pattern = "^MT-")
HC3[["percent.mt"]] <- PercentageFeatureSet(HC3, pattern = "^MT-")
HC4[["percent.mt"]] <- PercentageFeatureSet(HC4, pattern = "^MT-")
HC5[["percent.mt"]] <- PercentageFeatureSet(HC5, pattern = "^MT-")


### AD
AD1[["percent.mt"]] <- PercentageFeatureSet(AD1, pattern = "^MT-")
AD2[["percent.mt"]] <- PercentageFeatureSet(AD2, pattern = "^MT-")
AD3[["percent.mt"]] <- PercentageFeatureSet(AD3, pattern = "^MT-")
AD4[["percent.mt"]] <- PercentageFeatureSet(AD4, pattern = "^MT-")
AD13[["percent.mt"]] <- PercentageFeatureSet(AD13, pattern = "^MT-")

####AD16W
AD10[["percent.mt"]] <- PercentageFeatureSet(AD10, pattern = "^MT-")
AD11[["percent.mt"]] <- PercentageFeatureSet(AD11, pattern = "^MT-")
AD12[["percent.mt"]] <- PercentageFeatureSet(AD12, pattern = "^MT-")
AD17[["percent.mt"]] <- PercentageFeatureSet(AD17, pattern = "^MT-")
AD18[["percent.mt"]] <- PercentageFeatureSet(AD18, pattern = "^MT-")

####AD16W
AD14[["percent.mt"]] <- PercentageFeatureSet(AD14, pattern = "^MT-")
AD15[["percent.mt"]] <- PercentageFeatureSet(AD15, pattern = "^MT-")
AD16[["percent.mt"]] <- PercentageFeatureSet(AD16, pattern = "^MT-")
AD19[["percent.mt"]] <- PercentageFeatureSet(AD19, pattern = "^MT-")

### normal
HC1 <- subset(HC1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
HC2 <- subset(HC2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
HC3 <- subset(HC3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
HC4 <- subset(HC4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
HC5 <- subset(HC5, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)

### AD
AD1 <- subset(AD1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
AD2 <- subset(AD2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
AD3 <- subset(AD3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
AD4 <- subset(AD4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
AD13 <- subset(AD13, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)

####AD16W
AD10 <- subset(AD10, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
AD11 <- subset(AD11, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
AD12 <- subset(AD12, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
AD17 <- subset(AD17, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
AD18 <- subset(AD18, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)

####AD16W
AD14 <- subset(AD14, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
AD15 <- subset(AD15, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
AD16 <- subset(AD16, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
AD19 <- subset(AD19, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)

##############
### normal
table(Idents(HC1))
table(Idents(HC2))
table(Idents(HC3))
table(Idents(HC4))
table(Idents(HC5))

### AD
table(Idents(AD1))
table(Idents(AD2))
table(Idents(AD3))
table(Idents(AD4))
table(Idents(AD13))

####AD16W
table(Idents(AD10))
table(Idents(AD11))
table(Idents(AD12))
table(Idents(AD17))
table(Idents(AD18))

####AD16W
table(Idents(AD14))
table(Idents(AD15))
table(Idents(AD16))
table(Idents(AD19))

### individual sample
HC1$donor <- "HC1"
HC2$donor <- "HC2"
HC3$donor <- "HC3"
HC4$donor <- "HC4"
HC5$donor <- "HC5"

### AD
AD1$donor <- "AD1"
AD2$donor <- "AD2"
AD3$donor <- "AD3"
AD4$donor <- "AD4"
AD13$donor <- "AD13"

####AD16W
AD10$donor <- "AD10"
AD11$donor <- "AD11"
AD12$donor <- "AD12"
AD17$donor <- "AD17"
AD18$donor <- "AD18"

HC1$status <- "HC"
HC2$status <- "HC"
HC3$status <- "HC"
HC4$status <- "HC"
HC5$status <- "HC"

### AD
AD1$status <- "AD"
AD2$status <- "AD"
AD3$status <- "AD"
AD4$status <- "AD"
AD13$status <- "AD"

####AD16W
AD10$status <- "D16w"
AD11$status <- "D16w"
AD12$status <- "D16w"
AD17$status <- "D16w"
AD18$status <- "D16w"

####AD16W
AD14$status <- "D1y"
AD15$status <- "D1y"
AD16$status <- "D1y"
AD19$status <- "D1y"

### chemistry
HC1$chem <- "V2"
HC2$chem <- "V2"
HC3$chem <- "V2"
HC4$chem <- "V2"
HC5$chem <- "V3"

### AD
AD1$chem <- "V2"
AD2$chem <- "V2"
AD3$chem <- "V2"
AD4$chem <- "V3"
AD13$chem <- "V3"

####AD16W
AD10$chem <- "V3"
AD11$chem <- "V3"
AD12$chem <- "V3"
AD17$chem <- "V3"
AD18$chem <- "V3"

####AD16W
AD14$chem <- "V3"
AD15$chem <- "V3"
AD16$chem <- "V3"
AD19$chem <- "V3"


### source
HC1$source <- "SI"
HC2$source <- "SI"
HC3$source <- "SI"
HC4$source <- "SI"
HC5$source <- "SI"

### AD
AD1$source <- "SI"
AD2$source <- "SI"
AD3$source <- "SI"
AD4$source <- "SI"
AD13$source <- "SI"

####AD16W
AD10$source <- "SI"
AD11$source <- "SI"
AD12$source <- "SI"
AD17$source <- "SI"
AD18$source <- "SI"

####AD16W
AD14$source <- "SI"
AD15$source <- "SI"
AD16$source <- "SI"
AD19$source <- "SI"

### verv3.0.2on
HC1$version <- "v3.0.2"
HC2$version <- "v3.0.2"
HC3$version <- "v3.0.2"
HC4$version <- "v3.0.2"
HC5$version <- "v3.0.2"

### AD
AD1$version <- "v3.0.2"
AD2$version <- "v3.0.2"
AD3$version <- "v3.0.2"
AD4$version <- "v3.0.2"
AD13$version <- "v3.0.2"

####AD16W
AD10$version <- "v3.0.2"
AD11$version <- "v3.0.2"
AD12$version <- "v3.0.2"
AD17$version <- "v3.0.2"
AD18$version <- "v3.0.2"

####AD16W
AD14$version <- "v3.0.2"
AD15$version <- "v3.0.2"
AD16$version <- "v3.0.2"
AD19$version <- "v3.0.2"


### verv3.0.2on
HC1$dis <- "HC"
HC2$dis <- "HC"
HC3$dis <- "HC"
HC4$dis <- "HC"
HC5$dis <- "HC"

### AD
AD1$dis <- "AD2"
AD2$dis <- "AD2"
AD3$dis <- "AD2"
AD4$dis <- "AD2"
AD13$dis <- "AD2"

####AD16W
AD10$dis <- "D16w"
AD11$dis <- "D16w"
AD12$dis <- "D16w"
AD17$dis <- "D16w"
AD18$dis <- "D16w"

####AD16W
AD14$dis <- "D1y"
AD15$dis <- "D1y"
AD16$dis <- "D1y"
AD19$dis <- "D1y"


Merge <- merge(x = HC1, y = c(HC2, HC3, HC4, HC5, AD1, AD2, AD3, AD4, AD13, AD10, AD11, AD12, AD17, AD18, AD14, AD15, AD16, AD19))


Merge <- Merge %>% 
  NormalizeData(verbose = T)%>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Merge <- CellCycleScoring(Merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Merge <- Merge %>% 
  ScaleData(verbose = T)%>% 
  RunPCA(pc.genes = Merge@var.genes, npcs = 50, verbose = T, nfeatures.print = 15)

library(harmony)
options(repr.plot.height = 6, repr.plot.width = 10)
Merge <- Merge %>% 
  RunHarmony("donor", plot_convergence = TRUE, max.iter.harmony = 20)

ElbowPlot(Merge, ndims = 50)
options(repr.plot.height =20, repr.plot.width = 15)
DimHeatmap(Merge, dims = 1:15, cells = 500)
DimHeatmap(Merge, dims = 16:30, cells = 500)


########Decide the PC numbers
# Determine percent of variation associated with each PC
pct <- Merge[["pca"]]@stdev / sum(Merge[["pca"]]@stdev) * 100

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


Merge <- Merge %>% 
  RunUMAP(reduction = "harmony", dims = 1:16) %>% 
  RunTSNE(reduction = "harmony", dims = 1:16) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:16)%>% 
  FindClusters(resolution = 0.4, algorithm = 1)

saveRDS(Merge, "E:/Human - 10 clusters/yale-method/SciImmunol-data analysis/new integration-0304/Science_Immunol.rds")
