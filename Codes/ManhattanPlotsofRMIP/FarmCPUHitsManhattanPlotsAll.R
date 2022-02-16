#ManhattanPlost using ggplot
setwd("~/Documents/FCPUManhattanPlot")
rm(list = ls())
library(tidyverse)
library(data.table)
colors <- c(Agronomic='#E41A1C', CellularOrBiochemical='#377EB8', Disease='#4DAF4A', FloweringTime='#808080', Inflorescence='#984EA3', Root='#FF7F00', SeedComposition='#F781BF', Vegetative='#A65628')

map <- fread("WiDiv1051.geno.map", data.table = T)
res <- fread("Phenotype_162_MaizeFarmCPUResampSummaryCombined.csv")

gwas.dat <- merge(map[,1:3], res[,c(1, 9, 11)], by="SNP", all.x = T)
gwas.dat <- as.data.frame(gwas.dat)

nCHR <- length(unique(gwas.dat$CHROM))
gwas.dat$BPcum <- NA
s <- 0
nbp <- c()
for (i in sort(unique(gwas.dat$CHROM))){
  nbp[i] <- max(gwas.dat[gwas.dat$CHROM == i,]$POS)
  gwas.dat[gwas.dat$CHROM == i,"BPcum"] <- gwas.dat[gwas.dat$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}

axis.set <- gwas.dat %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

#F1 <- 
F1 <- ggplot(gwas.dat, aes(x = BPcum, y = Freq/100, color = PhenoGroup)) +
  geom_point(alpha = 0.75, size=4) +
  geom_hline(yintercept = c(0.05, 0.1), color = c("grey40", "darkred"), linetype = "dashed") + 
  geom_vline(xintercept = c(0, axis.set$maxBP), alpha=0.1) + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  labs(x = NULL,  y = "RMIP") + 
  scale_color_manual(values = colors) + # coment this out if want continious colors
  theme_classic(base_size = 20, base_family = "NimbusSan") +
  theme( legend.position = "top", #top or none
         panel.border = element_blank(),
         panel.grid.minor=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         axis.text.x = element_text(size = 20, vjust = 0.5))
ggsave(F1, file="FarmCPUHitsManhattanPlotAll.png", width = 18, height = 8, dpi = 300)
