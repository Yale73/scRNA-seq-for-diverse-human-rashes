Trm <- readRDS("E:/Human - 10 clusters/yale-method/000AA-Final objects/our_Trm_no_218_203.rds")
DimPlot(Trm)
#################### Network
table(Idents(Trm))# check cell numbers
# Trm 
# 13976
my_genes <- read.xlsx("E:/Human - 10 clusters/yale-method/T50-B50-0807/AD-PV-80%-AvsP-p0.05-fc0.25-0809-C2.xlsx", sheetIndex = 1)
Genes <- my_genes$Genes
a <- as.data.frame(Trm[["RNA"]]@data[Genes, ])#????????????????????????data.frame
a1 <- as.data.frame(t(a))#????????????????????????data.matrix, ???????????????
a1[a1 == 0] <- NA#??????0???????????????"NA",????????????

not_all_na <- function(x) any(!is.na(x))#????????????NA????????????
library(tidyverse)
b <- a1 %>% select_if(not_all_na)#??????NA????????????
b[is.na(b)] = 0#???NA?????????0,??????????????????

CorrelationMatrix <- cor(b, method = "spearman")#???????????????data.matrix

groups <- list(PV = 1:50,
               AD = 51:100)

library(qgraph)
qgraph(CorrelationMatrix, graph = "cor", layout = "spring", vsize=4, 
       directed = FALSE, posCol = "brown4", negCol = "blue4", 
       edge.labels = F, 
       labels=colnames(CorrelationMatrix), 
       label.prop=1,
       minimum="sig", sampleSize=13976,
       groups = groups,
       color = c("#FF3366", "#99FF00"),
       border.color = "white",border.width = 1,
       #edge.color = "darkgray", 
       edge.width = 0.8,
       curve = 0.2, curveAll = T,
       label.color = "black")

################################# STRING
######### STRING
library(STRINGdb)
string_db <- STRINGdb$new(version="11", species=9606,score_threshold=200, input_directory="")

my_genes <- read.xlsx("E:/Human - 10 clusters/yale-method/0807-T50-B50/AD-PV-80%-AvsP-p0.05-fc0.25-0809-C2.xlsx", sheetIndex = 1)

example1_mapped <- string_db$map(my_genes, "Genes", removeUnmappedRows = FALSE)
#Warning:  we couldn't map to STRING 4% of your identifiers

hits <- example1_mapped$STRING_id[1:200]
string_db$plot_network( hits)

#### label Pso genes
write.xlsx(example1_mapped, "E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 5/example1_mapped_0809.xlsx") ### add order numbericly and then other for color.

example1_mapped <- read.xlsx("E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 5/example1_mapped_0809.xlsx", sheetIndex = 1)

example1_mapped_logFC <- string_db$add_diff_exp_color(subset(example1_mapped, order <103), logFcColStr="other")

unique(example1_mapped_logFC$color)
example1_mapped_logFC$color2[example1_mapped_logFC$color == "#FF0000FF"] <- '#FF3366'
example1_mapped_logFC$color2[example1_mapped_logFC$color == "#00FF00FF"] <- '#99FF00'


payload_id <- string_db$post_payload(example1_mapped_logFC$STRING_id, 
                                     colors = example1_mapped_logFC$color2)

string_db$plot_network(hits, payload_id = payload_id)

