############################## ADvsN, PvsN, AvsP  ####################
Idents(our) <- our$ID

set1 = list()
for (i in 1:41){
  tryCatch({
    print(i)
    gc() 
    set1[[paste0("C",i)]] <- FindMarkers(our, ident.1 = "AD1", ident.2 ="NML", group.by = "dis", verbose =TRUE,
                                         assay = "RNA", slot = "data", subset.ident = i, logfc.threshold = 0.2, 
                                         test.use = "MAST")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.xlsx(set1, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AN_1031.xlsx", rowNames = TRUE, colNames=TRUE)


set2 = list()
for (i in 1:41){
  tryCatch({
    print(i)
    gc() 
    set2[[paste0("C",i)]] <- FindMarkers(our, ident.1 = "Pso", ident.2 ="NML", group.by = "dis", verbose =TRUE,
                                         assay = "RNA", slot = "data", subset.ident = i,  logfc.threshold = 0.2,
                                         test.use = "MAST")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.xlsx(set2, file = "F:/4-7-20m-Final-1030/Figure 3/double check/PN_1031.xlsx", rowNames = TRUE, colNames=TRUE)

###############################
###########################function
library(readxl)    
library(stringr)
library(gridExtra)
library(openxlsx)

## Task 1

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


## Task 2  

PV_N = read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/PN_1031.xlsx")
AD_N = read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/AN_1031.xlsx")

set1 = list()
for(i in names(PV_N)){
  # select DEG in each contrast that has p.value<0.05 and abs(lfc) > 0.25
  sub_PVN = PV_N[[i]]
  sub_PVN = sub_PVN[sub_PVN$p_val_adj<0.05 & abs(sub_PVN$avg_log2FC)>0.25, ]
  set1[[paste0(i)]] = list("PN" = sub_PVN)
}


write.xlsx(set1, file = "F:/4-7-20m-Final-1030/Figure 3/double check/PN-p0.05-fc025.xlsx", colNames=T, rowNames=T)

detach("package:openxlsx", unload = TRUE)
library(xlsx)
for(i in names(set1)){
  len <- as.data.frame(length(rownames(set1[[i]]$PN)))
  colnames(len) <- "number"
  tab  <- paste0(i)
  write.xlsx(len, file = "F:/4-7-20m-Final-1030/Figure 3/double check/PN-p0.05-fc025-number.xlsx", sheetName = tab, append = TRUE)
}

detach("package:xlsx", unload = TRUE)
library(openxlsx)

set2 = list()
for(i in names(AD_N)){
  # select DEG in each contrast that has p.value<0.05 and abs(lfc) > 0.25
  sub_AD = AD_N[[i]]
  sub_AD = sub_AD[sub_AD$p_val_adj<0.05 & abs(sub_AD$avg_log2FC)>0.25,]
  set2[[paste0(i)]] = list("AN" = sub_AD)
}

write.xlsx(set2, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AN-p0.05-fc025.xlsx",colNames=T, rowNames=T)


detach("package:openxlsx", unload = TRUE)
library(xlsx)
for(i in names(set2)){
  len <- as.data.frame(length(rownames(set2[[i]]$AN)))
  colnames(len) <- "number"
  tab  <- paste0(i)
  write.xlsx(len, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AN-p0.05-fc025-number.xlsx", sheetName = tab, append = TRUE)
}

################################# AD-PV-no 80%##########################
set3=list()
for(i in names(set1)){
  set3[[i]] <- c(set1[[i]][["PN"]][["...1"]], set2[[i]][["AN"]][["...1"]])
}

names(set3)
names(set3) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
                 "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34",
                 "35", "36", "37", "38", "39")

set6 = list()
for (i in names(set3)){
  tryCatch({
    gc() 
    set6[[paste0("C",i)]] <- FindMarkers(our, ident.1 = "AD1", ident.2 = "Pso", group.by = "dis", verbose =TRUE, assay = "RNA", 
                                         slot = "data", subset.ident = i, test.use = "MAST", min.cells.feature = 0, min.pct = 0, 
                                         min.diff.pct = 0, min.cells.group = 0, features = set3[[i]][!duplicated(set3[[i]])])
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.xlsx(set6, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AP_1031.xlsx", rowNames = TRUE, colNames=TRUE)



set1 = list()
for(i in names(set6)){
  # select DEG in each contrast that has p.value<0.05 and abs(lfc) > 0.25
  sub_AP = set6[[i]]
  sub_AP = sub_AP[sub_AP$p_val_adj<0.05 & abs(sub_AP$avg_log2FC)>0.25,]
  set1[[paste0(i)]] = sub_AP
}

#################### remove matrix with empty row, when openxlsx doesn't work ##############
set1 <- set1[sapply(set1, nrow)>0]

#####openxlsx
write.xlsx(set1, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AP-p0.05-fc025-openxlsx.xlsx", rowNames = TRUE, colNames=TRUE)


detach("package:openxlsx", unload = TRUE)
library(xlsx)

for(i in names(set1)){
  len <- as.data.frame(length(rownames(set1[[i]]$AP)))
  colnames(len) <- "number"
  tab  <- paste0(i)
  write.xlsx(len, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AP-p0.05-fc025-number.xlsx", sheetName = tab, append = TRUE)
}


##################################### individula 80% #####################################################
################################# load data for AD and Pso sample DEGs ########################################
Idents(our) <- our$donor
our <- RenameIdents(our, `150` = "NML", `154` = "NML", `155` = "NML",`169` = "NML",`195` = "NML", `204` = "NML", `207` = "NML",
                    `170` = "170",`198` = "198", `230` = "230", `231` = "231", `232` = "232", `233` = "233", `236` = "236",
                    `165` = "165", `173` = "173", `194` = "194", `199` = "199", `211` = "211", `222` = "222", `234` = "234", 
                    `235` = "235", `167` = "167", `174` = "174", `192` = "192", `200` = "200", `202` = "202", `219` = "219",
                    `141` = "141", `163` = "163", `175` = "175")

our[["donor2"]] <- Idents(object = our)

Idents(our) <- our$ID

############################ AD
set2 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/AN-p0.05-fc025.xlsx")
donor <- c("170", "198", "230", "231", "232", "233", "236")
names(set2)
names(set2) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
                 "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", 
                 "27", "28", "29", "30", "31", "32", "33", "35", "36", "37", "38", "39")

for (j in donor){
  set1 = list()
  for (i in names(set2)){
    tryCatch({
      print(j)
      print(i)
      gc()
      set1[[paste0("C", i)]] <- FindMarkers(our, ident.1 = j, ident.2 ="NML", group.by = "donor2", verbose =TRUE, assay = "RNA", 
                              slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                              subset.ident = i, logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = set2[[i]][["...1"]])
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
    write.xlsx(set1, file = paste0("F:/4-7-20m-Final-1030/Figure 3/double check/", j, "-", "AN-p0.05-fc0.25", ".xlsx"), colNames=T, rowNames=T)
}
  
  
############################## Pso
donor <- c("165", "173", "194", "199", "211", "222", "234", "235")

set1 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/PN-p0.05-fc025.xlsx")

names(set1)
names(set1) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
                 "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", 
                 "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39")

for (j in donor){
  set2 = list()
  for (i in names(set1)){
    tryCatch({

      print(j)
      print(i)
      gc()
      set2[[paste0("C", i)]] <- FindMarkers(our, ident.1 = j, ident.2 ="NML", group.by = "donor2", verbose =TRUE,
                                            assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                                            subset.ident = i, logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = set1[[i]][["...1"]])
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  write.xlsx(set1, file = paste0("F:/4-7-20m-Final-1030/Figure 3/double check/", j, "-", "PN-p0.05-fc0.25", ".xlsx"), colNames=T, rowNames=T)
}


#####################expressed by 80%
library(readxl)    
library(stringr)
library(gridExtra)
library(openxlsx)


read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

###############################################p0.05-fc0.25##############################################################
##################################### AD 80% #####################################################
S170 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/170-AN-p0.05-fc0.25.xlsx")
S198 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/198-AN-p0.05-fc0.25.xlsx")
S230 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/230-AN-p0.05-fc0.25.xlsx")
S231 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/231-AN-p0.05-fc0.25.xlsx")
S232 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/232-AN-p0.05-fc0.25.xlsx")
S233 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/233-AN-p0.05-fc0.25.xlsx")
S236 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/236-AN-p0.05-fc0.25.xlsx")

##### names(S170), choose the longest one ############
combine_data_AD=list()
for(i in names(S230)){
  combine_data_AD[[i]] = list(S170[[i]], S198[[i]], S230[[i]], S231[[i]], S232[[i]], S233[[i]], S236[[i]])
  names(combine_data_AD[[i]]) = c("S170", "S198", "S230", "S231", "S232", "S233", "S236")
  combine_data_AD[[i]] = combine_data_AD[[i]][!unlist(lapply(combine_data_AD[[i]], is.null))]
}


library(reshape2)
library(plyr)

AD_results = list()
for(i in names(combine_data_AD)){
  multi_full <- ldply(combine_data_AD[[i]])
  colnames(multi_full)[1] = "sample"
  colnames(multi_full)[2] = "Gene"
  multi_wide = dcast(melt(multi_full, id.vars=c("Gene", "sample")), Gene~variable+sample)
  AD_results[[i]] = multi_wide
}

AD_final = list()
for(i in names(AD_results)){
  Data <- AD_results[[i]]
  
  colnames(Data)[1] = "Gene"
  colnames(Data)[which(str_detect(colnames(Data), "p_val_adj"))] = paste0("p.adj_", colnames(Data)[which(str_detect(colnames(Data), "p"))-1])
  
  # determine significant p values
  qvals = Data[, which(str_detect(colnames(Data), "p.adj_"))]
  qvals_sig = ifelse(qvals<0.05, "S", "NS")
  colnames(qvals_sig) = str_replace(colnames(qvals_sig), "p.adj_", "")
  
  # determine if same directions
  combine_sig = (qvals_sig[,-1] == "S")
  
  genes = Data$Gene[apply(combine_sig, 1, function(x) sum(x) >= length(x)*0.8)]
  
  AD_final[[i]] = genes
}

write.xlsx(AD_final, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AN-80%-p0.05-fc0.25.xlsx", 
           colNames=TRUE, rowNames=TRUE)


for(i in names(AD_final)){
  len <- as.data.frame(length(AD_final[[i]]))
  colnames(len) <- "number"
  tab  <- paste0(i)
  write.xlsx(len, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AN-individual-numbers-p0.05-fc0.25-80%.xlsx", sheetName = tab, append = TRUE)
}

##################################### Psoriasis 80% #####################################################
S165 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/165-PN-p0.05-fc0.25.xlsx")
S173 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/173-PN-p0.05-fc0.25.xlsx")
S194 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/194-PN-p0.05-fc0.25.xlsx")
S199 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/199-PN-p0.05-fc0.25.xlsx")
S211 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/211-PN-p0.05-fc0.25.xlsx")
S222 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/222-PN-p0.05-fc0.25.xlsx")
S234 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/234-PN-p0.05-fc0.25.xlsx")
S235 <- read_excel_allsheets("F:/4-7-20m-Final-1030/Figure 3/double check/235-PN-p0.05-fc0.25.xlsx")

##### names(S199), choose the longest one ############
combine_data_PV=list()
for(i in names(S222)){
  combine_data_PV[[i]] = list(S165[[i]], S173[[i]], S194[[i]], S199[[i]], S211[[i]], S222[[i]], S234[[i]], S235[[i]])
  names(combine_data_PV[[i]]) = c("S165", "S173", "S194", "S199", "S211", "S222", "S234", "S235")
  combine_data_PV[[i]] = combine_data_PV[[i]][!unlist(lapply(combine_data_PV[[i]], is.null))]
}


library(reshape2)

PV_results = list()
for(i in names(combine_data_PV)){
  multi_full <- ldply(combine_data_PV[[i]])
  colnames(multi_full)[1] = "sample"
  colnames(multi_full)[2] = "Gene"
  multi_wide = dcast(melt(multi_full, id.vars=c("Gene", "sample")), Gene~variable+sample)
  PV_results[[i]] = multi_wide
}

PV_final = list()
for(i in names(PV_results)){
  Data <- PV_results[[i]]
  
  colnames(Data)[1] = "Gene"
  colnames(Data)[which(str_detect(colnames(Data), "p_val_adj"))] = paste0("p.adj_", colnames(Data)[which(str_detect(colnames(Data), "p"))-1])
  
  # determine significant p values
  qvals = Data[, which(str_detect(colnames(Data), "p.adj_"))]
  qvals_sig = ifelse(qvals<0.05, "S", "NS")
  colnames(qvals_sig) = str_replace(colnames(qvals_sig), "p.adj_", "")
  
  # determine if same directions
  combine_sig = (qvals_sig[,-1] == "S")
  
  genes = Data$Gene[apply(combine_sig, 1, function(x) sum(x) >= length(x)*0.8)]
  genes = 
    PV_final[[i]] = genes
}



write.xlsx(PV_final, file = "F:/4-7-20m-Final-1030/Figure 3/double check/PN-80%-p0.05-fc0.25.xlsx",
           colNames=TRUE, rowNames=TRUE)

for(i in names(PV_final)){
  len <- as.data.frame(length(PV_final[[i]]))
  colnames(len) <- "number"
  tab  <- paste0(i)
  write.xlsx(len, file = "F:/4-7-20m-Final-1030/Figure 3/double check/PN-individual-numbers-p0.05-fc0.25-80%.xlsx", sheetName = tab, append = TRUE)
}


################################# AD-PV-after 80%##########################
set3=list()
for(i in names(AD_final)){
  set3[[i]] <- c(AD_final[[i]], PV_final[[i]])
}

names(set3)
names(set3) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 
                 "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
                 "30", "31", "32", "33", "36", "37", "38")

set6 = list()
for (i in names(set3)){
  tryCatch({
    gc() 
    set6[[paste0("C",i)]] <- FindMarkers(our, ident.1 = "AD1", ident.2 = "Pso", group.by = "dis", verbose =TRUE, assay = "RNA", 
                                         slot = "data", min.cells.feature = 0, min.cells.group = 0, 
                                         subset.ident = i, logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, 
                                         test.use = "MAST", features = set3[[i]][!duplicated(set3[[i]])])
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.xlsx(set6, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AP_80%.xlsx", rowNames = TRUE, colNames=TRUE)



set1 = list()
for(i in names(set6)){
  # select DEG in each contrast that has p.value<0.05 and abs(lfc) > 0.25
  sub_AP = set6[[i]]
  sub_AP = sub_AP[sub_AP$p_val_adj<0.05 & abs(sub_AP$avg_log2FC)>0.25,]
  set1[[paste0(i)]] = sub_AP
}

#################### remove matrix with empty row, when openxlsx doesn't work ##############
set1 <- set1[sapply(set1, nrow)>0]
#####xlsx
for (i in names(set1)) {
  marker_i <- set1[i]
  tab <- paste(i)
  write.xlsx(marker_i, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AP_80%_p0.05_fc025.xlsx", sheetName = tab, append = TRUE)
}

#####openxlsx
write.xlsx(set1, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AP-p0.05-fc025-openxlsx.xlsx", rowNames = TRUE, colNames=TRUE)


detach("package:openxlsx", unload = TRUE)
library(xlsx)
for(i in names(set1)){
  len <- as.data.frame(length(rownames(set1[[i]])))
  colnames(len) <- "number"
  tab  <- paste0(i)
  write.xlsx(len, file = "F:/4-7-20m-Final-1030/Figure 3/double check/AP_80%_p0.05_fc025-number.xlsx", sheetName = tab, append = TRUE)
}
