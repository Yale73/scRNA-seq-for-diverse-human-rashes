# library
library(ggplot2)
library(cowplot)

##################################### barplot ###############################
AD <- read.xlsx("E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 3/DEG numbers for Barplot.xlsx", sheetIndex = 1)
PV <- read.xlsx("E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 3/DEG numbers for Barplot.xlsx", sheetIndex = 2)
AP <- read.xlsx("E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 3/DEG numbers for Barplot.xlsx", sheetIndex = 3)

AD_pct <- read.xlsx("E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 3/DEG numbers for Barplot.xlsx", sheetIndex = 4)
PV_pct <- read.xlsx("E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 3/DEG numbers for Barplot.xlsx", sheetIndex = 5)
AP_pct <- read.xlsx("E:/Human - 10 clusters/yale-method/000AA-Final figures/Figure 3/DEG numbers for Barplot.xlsx", sheetIndex = 6)

# Grouped
AD$Cluster <- factor(AD$Cluster, levels=c("Tcm",	"Trm1",	"eTreg1",	"Trm2",	"CTLex",	"CTLem",	"Tet",	"Tmm1",	"ILC/NK",	"NK",	"Trm3",
                                          "cmTreg",	"Tmm2",	"eTreg2",	"Tmm3",	"ILC2",	"Tn",	"moDC1",	"LC1",	"InfMono",	"Mac1",	"Mono",
                                          "LC2",	"moDC2",	"moDC3",	"DC1",	"DC2",	"LC3",	"Mac2",	"Plasma",	"B",	"Mac3",	"DC3",
                                          "migDC",	"Mac4",	"Mast",	"Trm-c",	"Treg-c",	"Mast-c",	"ILC/NK-c",	"CTL-c"))
PV$Cluster <- factor(PV$Cluster, levels=c("Tcm",	"Trm1",	"eTreg1",	"Trm2",	"CTLex",	"CTLem",	"Tet",	"Tmm1",	"ILC/NK",	"NK",	"Trm3",
                                          "cmTreg",	"Tmm2",	"eTreg2",	"Tmm3",	"ILC2",	"Tn",	"moDC1",	"LC1",	"InfMono",	"Mac1",	"Mono",
                                          "LC2",	"moDC2",	"moDC3",	"DC1",	"DC2",	"LC3",	"Mac2",	"Plasma",	"B",	"Mac3",	"DC3",
                                          "migDC",	"Mac4",	"Mast",	"Trm-c",	"Treg-c",	"Mast-c",	"ILC/NK-c",	"CTL-c"))
AP$Cluster <- factor(AP$Cluster, levels=c("Tcm",	"Trm1",	"eTreg1",	"Trm2",	"CTLex",	"CTLem",	"Tet",	"Tmm1",	"ILC/NK",	"NK",	"Trm3",
                                          "cmTreg",	"Tmm2",	"eTreg2",	"Tmm3",	"ILC2",	"Tn",	"moDC1",	"LC1",	"InfMono",	"Mac1",	"Mono",
                                          "LC2",	"moDC2",	"moDC3",	"DC1",	"DC2",	"LC3",	"Mac2",	"Plasma",	"B",	"Mac3",	"DC3",
                                          "migDC",	"Mac4",	"Mast",	"Trm-c",	"Treg-c",	"Mast-c",	"ILC/NK-c",	"CTL-c"))
AD_pct$Cluster <- factor(AD_pct$Cluster, levels=c("Tcm",	"Trm1",	"eTreg1",	"Trm2",	"CTLex",	"CTLem",	"Tet",	"Tmm1",	"ILC/NK",	"NK",	"Trm3",
                                          "cmTreg",	"Tmm2",	"eTreg2",	"Tmm3",	"ILC2",	"Tn",	"moDC1",	"LC1",	"InfMono",	"Mac1",	"Mono",
                                          "LC2",	"moDC2",	"moDC3",	"DC1",	"DC2",	"LC3",	"Mac2",	"Plasma",	"B",	"Mac3",	"DC3",
                                          "migDC",	"Mac4",	"Mast",	"Trm-c",	"Treg-c",	"Mast-c",	"ILC/NK-c",	"CTL-c"))
PV_pct$Cluster <- factor(PV_pct$Cluster, levels=c("Tcm",	"Trm1",	"eTreg1",	"Trm2",	"CTLex",	"CTLem",	"Tet",	"Tmm1",	"ILC/NK",	"NK",	"Trm3",
                                          "cmTreg",	"Tmm2",	"eTreg2",	"Tmm3",	"ILC2",	"Tn",	"moDC1",	"LC1",	"InfMono",	"Mac1",	"Mono",
                                          "LC2",	"moDC2",	"moDC3",	"DC1",	"DC2",	"LC3",	"Mac2",	"Plasma",	"B",	"Mac3",	"DC3",
                                          "migDC",	"Mac4",	"Mast",	"Trm-c",	"Treg-c",	"Mast-c",	"ILC/NK-c",	"CTL-c"))
AP_pct$Cluster <- factor(AP_pct$Cluster, levels=c("Tcm",	"Trm1",	"eTreg1",	"Trm2",	"CTLex",	"CTLem",	"Tet",	"Tmm1",	"ILC/NK",	"NK",	"Trm3",
                                          "cmTreg",	"Tmm2",	"eTreg2",	"Tmm3",	"ILC2",	"Tn",	"moDC1",	"LC1",	"InfMono",	"Mac1",	"Mono",
                                          "LC2",	"moDC2",	"moDC3",	"DC1",	"DC2",	"LC3",	"Mac2",	"Plasma",	"B",	"Mac3",	"DC3",
                                          "migDC",	"Mac4",	"Mast",	"Trm-c",	"Treg-c",	"Mast-c",	"ILC/NK-c",	"CTL-c"))

p1 <- ggplot(AD, aes(fill=Comparison, y=number, x=Cluster)) + 
  geom_bar(position="dodge", stat="identity", fill="palevioletred2")+ 
  scale_x_discrete(labels = unique(AD$Cluster)) + labs(x = "")+
  theme(#axis.title.x=element_blank(), 
        axis.text.x=element_text(angle = 90),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("AD-total")

p2 <- ggplot(PV, aes(fill=Comparison, y=number, x=Cluster)) + 
  geom_bar(position="dodge", stat="identity", fill="dodgerblue")+ 
  scale_x_discrete(labels = unique(PV$Cluster)) + labs(x = "")+
  theme(#axis.title.x=element_blank(), 
    axis.text.x=element_text(angle = 90),
    axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())+
  ggtitle("PV-total")

p3 <- ggplot(AP, aes(fill=Comparison, y=number, x=Cluster)) + 
  geom_bar(position="dodge", stat="identity", fill="mediumorchid3")+ 
  scale_x_discrete(labels = unique(AP$Cluster)) + labs(x = "")+
  theme(#axis.title.x=element_blank(), 
    axis.text.x=element_text(angle = 90),
    axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())+
  ggtitle("AP-total")


p4 <- ggplot(AD_pct, aes(fill=Comparison, y=number, x=Cluster)) + 
  geom_bar(position="dodge", stat="identity", fill="palevioletred2")+ 
  scale_x_discrete(labels = unique(AD_pct$Cluster)) + labs(x = "")+
  theme(#axis.title.x=element_blank(), 
    axis.text.x=element_text(angle = 90),
    axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())+
  ggtitle("AD-80%")


p5 <- ggplot(PV_pct, aes(fill=Comparison, y=number, x=Cluster)) + 
  geom_bar(position="dodge", stat="identity", fill="dodgerblue")+ 
  scale_x_discrete(labels = unique(PV_pct$Cluster)) + labs(x = "")+
  theme(#axis.title.x=element_blank(), 
    axis.text.x=element_text(angle = 90),
    axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())+
  ggtitle("PV-80%")

p6 <- ggplot(AP_pct, aes(fill=Comparison, y=number, x=Cluster)) + 
  geom_bar(position="dodge", stat="identity", fill="dodgerblue")+ 
  scale_x_discrete(labels = unique(AP_pct$Cluster)) + labs(x = "")+
  theme(#axis.title.x=element_blank(), 
    axis.text.x=element_text(angle = 90),
    axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())+
  ggtitle("AP-80%")


plot_grid(p1, p2, p3, p4, p5, p6,ncol=3)