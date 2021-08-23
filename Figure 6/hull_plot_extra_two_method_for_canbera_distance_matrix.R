library(concaveman)
library(ggforce)
library(Seurat)
library(monocle3)
library(dplyr)
library(Rmisc)
library(ggrepel)
#Read in this object

humandat = readRDS("E:/Human - 10 clusters/yale-method/000AA-Final objects/our_new_umap_no_218_203.rds")

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
pts = c("170","198","230","231","232","233","236","165","173","194","199","211","222","234","235","167","174","192","200","202","219")
human_human_big_cds = human_human_big_cds[,colData(human_human_big_cds)$donor3 %in% pts] #subset to all normal patients and a diseased patient
colData(human_human_big_cds)$dis_updated = colData(human_human_big_cds)$donor3
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

colData(human_human_big_cds)$dis_updated = gsub("167","Indeterminate",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("174","Indeterminate",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("192","Indeterminate",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("200","Indeterminate",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("202","Indeterminate",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("219","Indeterminate",colData(human_human_big_cds)$dis_updated)
unique(colData(human_human_big_cds)$dis_updated)

#--------------------------------------------------------------------#
#Get the genes of interest
gene_sigs = data.frame(ad_genes = c("SLC5A3",	"SFXN1",	"GATA3",	"TUBA1B",	"ZFP36L2",	"ZFP36",	"PLP2",	"CLU",	"CRIP1",	"NRIP3",	"LPCAT2",	"BIRC3",
                                     "FAM102A",	"IFITM1",	"MT-ND1",	"TSHZ2",	"LONRF2",	"CD99",	"RPL17",	"ISG15",	"MT-ND2",	"AHNAK",	"TNFSF10",	"TIMP1",
                                     "TGFBR3",	"PRKX",	"ACTG1",	"NBAS",	"ARHGAP21",	"RPL36A",	"NEAT1",	"IL17RB",	"THADA",	"SYNE2",	"PLA2G16",	"CYSLTR1",
                                     "MT-ND5",	"S100A10",	"IFITM2",	"LMO4",	"XIST",	"CSGALNACT1",	"SOS1",	"ANXA1",	"MFHAS1",	"ITM2C",	"CAPG",	"IL32",
                                     "LGALS1",	"TWIST1"), 
                       pv_genes = c("CXCL13",	"MTRNR2L12",	"CD7",	"MGAT4A",	"FTH1",	"LAYN",	"IL17F",	"KLRB1",	"GNLY",	"CPM",	"CTSH",	"GBP5",
                                    "CLEC2B",	"SOX4",	"GZMB",	"CD2",	"CEBPD",	"ODF2L",	"LAG3",	"LRRN3",	"ARHGEF12",	"PTPN13",	"TNFAIP3",	"TRPS1",
                                    "SNX9",	"METRNL",	"BTG1",	"JUN",	"TMIGD2",	"SPOCK2",	"GABARAPL1",	"PMEPA1",	"HIST1H1E",	"RBPJ",	"LINC01871",	"MAP3K4",
                                    "H1FX",	"UBC",	"GALNT1",	"PNRC1",	"MUC20-OT1",	"RPS26",	"GABPB1-AS1",	"CHN1",	"NAP1L4",	"PTMS",	"F2R",	"DAPK2",
                                    "CTLA4",	"CCR6"))									
                                     

human_human_big_cds@assays@data$counts[1:10,1:10]

#--------------------------------------------------------------------#
#Now subset to cluster 2 cells
#human_human_big_cds = human_human_big_cds[,colData(human_human_big_cds)$ID %in% "2"] #subset to all normal patients and a diseased patient

#--------------------------------------------------------------------#
#AD gene data_matrix
ad_mat = t(human_human_big_cds@assays@data$counts[gene_sigs$ad_genes,])
ad_mat[1:10,1:10]
ad_int = data.frame(colData(human_human_big_cds)$dis_updated, colData(human_human_big_cds)$donor3, rowSums(ad_mat) )
names(ad_int) = c("dis","sample","gene_sig")
names(ad_int) = c("dis","sample","gene_sig")
ad_dfc = summarySE(ad_int, measurevar='gene_sig', groupvars=c("dis","sample"))
head(ad_dfc)
names(ad_dfc)[c(4,6)] = c("ad_gene_sig","ad_se")

#--------------------------------------------------------------------#
#PV gene data_matrix
pv_mat = t(human_human_big_cds@assays@data$counts[gene_sigs$pv_genes,])
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
  ggtitle("CIR") + 
  xlab("AD-specific genes") +
  ylab("PV-specific genes") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkblue",
                                        "PV" = "darkred",
                                        "Indeterminate" = "darkorange")) +
  scale_fill_manual(name="",values = c("AD" = "darkblue",
                                       "PV" = "darkred",
                                       "Indeterminate" = "darkorange")) +
  geom_mark_hull(aes(fill=dis, label=dis), concavity=5) +
  theme(axis.text=element_text(size=15,color='black'),
        axis.title=element_text(size=15,color='black'),
        plot.title=element_text(size=15,color='black'))


colData(human_human_big_cds)$dis_updated = colData(human_human_big_cds)$donor3
colData(human_human_big_cds)$dis_updated = gsub("170","",colData(human_human_big_cds)$dis_updated)


#----------------------------------------------#
#Test Now calculate distance between each indeterminate sample, and the PV/AD samples

p1 = as.matrix(subset(re_int))
library(raster)
#pointDistance

all_coords_mat = as.matrix(data.frame(re_int$ad_gene_sig,re_int$pv_gene_sig))
row.names(all_coords_mat) = re_int$sample
dist_mat = pointDistance(all_coords_mat, allpairs=T, lonlat=F) #calculate all-vs-all distance

rownames(dist_mat) = re_int$sample
colnames(dist_mat) = re_int$sample

ind_samps = subset(re_int, dis=="Indeterminate")$sample
pv_sample = 
#Loop through the inderminate samples and get those sweet averages betch

all_res = do.call(rbind, lapply(1:length(ind_samps), function(i){
  #i=1
  dist_df = data.frame(re_int$dis,re_int$sample, dist_mat[ind_samps[i],]) #get distance of a sample from the matrix and combine with metadata
  names(dist_df) = c("dis","sample","dist")
  
  ad_dist_df = subset(dist_df, dis=="AD") #remove other indeterminate samples
  pv_dist_df = subset(dist_df, dis=="PV") #remove other indeterminate samples
  
  #calculate means to see which sided test i wanna do
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



############################################rlof ###########################################
library(Rlof)
dist_mat2 =distmc(all_coords_mat, method = "canberra", diag = FALSE, upper = TRUE, p = 2)
dist_mat2 <- as.matrix(dist_mat2)
rownames(dist_mat2) = re_int$sample
colnames(dist_mat2) = re_int$sample

ind_samps = subset(re_int, dis=="Indeterminate")$sample
pv_sample = 
  #Loop through the inderminate samples and get those sweet averages betch
  
  all_res = do.call(rbind, lapply(1:length(ind_samps), function(i){
    #i=1
    dist_df = data.frame(re_int$dis,re_int$sample, dist_mat2[ind_samps[i],]) #get distance of a sample from the matrix and combine with metadata
    names(dist_df) = c("dis","sample","dist")
    
    ad_dist_df = subset(dist_df, dis=="AD") #remove other indeterminate samples
    pv_dist_df = subset(dist_df, dis=="PV") #remove other indeterminate samples
    
    #calculate means to see which sided test i wanna do
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


#################################vegdist
p1 = as.matrix(subset(re_int))
library(raster)
#pointDistance

all_coords_mat = as.matrix(data.frame(re_int$ad_gene_sig,re_int$pv_gene_sig))
row.names(all_coords_mat) = re_int$sample

library(vegan)
dist_mat2 =vegdist(all_coords_mat, method = "canberra", diag = FALSE, upper = TRUE, p = 2)
dist_mat2 <- as.matrix(dist_mat2)
rownames(dist_mat2) = re_int$sample
colnames(dist_mat2) = re_int$sample

ind_samps = subset(re_int, dis=="Indeterminate")$sample
pv_sample = 
  #Loop through the inderminate samples and get those sweet averages betch
  
  all_res = do.call(rbind, lapply(1:length(ind_samps), function(i){
    #i=1
    dist_df = data.frame(re_int$dis,re_int$sample, dist_mat2[ind_samps[i],]) #get distance of a sample from the matrix and combine with metadata
    names(dist_df) = c("dis","sample","dist")
    
    ad_dist_df = subset(dist_df, dis=="AD") #remove other indeterminate samples
    pv_dist_df = subset(dist_df, dis=="PV") #remove other indeterminate samples
    
    #calculate means to see which sided test i wanna do
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

write.xlsx(pv_sample, "E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 6/vegdist-canberra-all-cluster -together.xlsx", colNames=T, rowNames=T)
