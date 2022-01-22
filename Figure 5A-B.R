our <- readRDS("F:/4-7-20m-Final-1030/our_no_218_203.rds")
DimPlot(our)
Trm <- subset(our, idents="2")
DimPlot(Trm)
#################### Network
table(Idents(Trm))# check cell numbers
# Trm 
# 13961

library(xlsx)
my_genes <- read.xlsx("F:/4-7-20m-Final-1030/Figure 5/AP_80%_p0.001_fc043.xlsx", sheetIndex = 2)
Genes <- my_genes$NA.
a <- as.data.frame(Trm[["RNA"]]@data[Genes, ])#????????????????????????data.frame
a1 <- as.data.frame(t(a))#????????????????????????data.matrix, ???????????????
a1[a1 == 0] <- NA#??????0???????????????"NA",????????????

not_all_na <- function(x) any(!is.na(x))#????????????NA????????????
library(tidyverse)
b <- a1 %>% select_if(not_all_na)#??????NA????????????
b[is.na(b)] = 0#???NA?????????0,??????????????????

CorrelationMatrix <- cor(b, method = "spearman")#???????????????data.matrix
write.xlsx(CorrelationMatrix, "F:/4-7-20m-Final-1030/Figure 5/correlation_matrix.xlsx")
groups <- list(AD = 1:52,
               PV = 53:110)

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

options(download.file.method="libcurl")

example1_mapped <- string_db$map(my_genes, "NA.", removeUnmappedRows = TRUE)

hits <- example1_mapped$STRING_id[1:200]
string_db$plot_network( hits)

#### label Pso genes
write.xlsx(example1_mapped, "F:/4-7-20m-Final-1030/Figure 5/example1_mapped.xlsx")
example1_mapped <- read.xlsx("F:/4-7-20m-Final-1030/Figure 5/example1_mapped.xlsx", sheetIndex = 1)
example1_mapped_logFC <- string_db$add_diff_exp_color(subset(example1_mapped, order <111), logFcColStr="other")

unique(example1_mapped_logFC$color)
example1_mapped_logFC$color2[example1_mapped_logFC$color == "#FF0000FF"] <- '#FF3366'
example1_mapped_logFC$color2[example1_mapped_logFC$color == "#00FF00FF"] <- '#99FF00'


payload_id <- string_db$post_payload(example1_mapped_logFC$STRING_id, 
                                     colors = example1_mapped_logFC$color2)

string_db$plot_network(hits, payload_id = payload_id)

