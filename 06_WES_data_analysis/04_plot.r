library(GenVisR)
library(ggplot2)

samples <- c(...)
plot_theme_text <- theme(strip.text.x = element_text(size = 15), strip.text.y = element_blank(),
                         strip.background.x = element_rect(colour="#8A8A8A", fill="white"),
                         axis.title.x = element_blank(), axis.title.y = element_blank(),
                         panel.spacing.x = unit(0, "lines"))

csize <- read.table("~/pack/chrom_sizes/hg38.chrom.sizes")
csize <- csize[c(1:7, 9:21, 23:24),]
colnames(csize) <- c("chromosome", "end")
csize$start <- 0

pdf("./hmm_segment.pdf", width=14, height=2)
for(sample in samples){
    cns <- read.table(paste0("./output/",sample,".mappings.sorted.cns"), header=1)
    
    cns$segmean <- cns$log2
    cns <- subset(cns, !chromosome %in% c("chrX","chrY"))
    
    cns$sample <- 1
    ylim <- c(min(cns$segmean), max(cns$segmean))
    q <- as.vector(quantile(cns$segmean))
    iqr <- q[4]-q[2]
    ylim <- c(q[2] - 1.5*iqr, q[4]+1.5*iqr)
    
    plot_theme <- scale_fill_gradient2(low="#2937BB", mid="white", high="#A81909", 
                                       midpoint=0, limit=ylim, name='LogCN')
    p <- cnSpec(cns, y=csize, genome = "hg38", CNscale='relative', plotLayer=plot_theme, plot_title=sample) + plot_theme_text
    plot(p)
}
dev.off()
