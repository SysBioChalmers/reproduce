#!/usr/bin/env Rscript
# Figures from the manuscript
# Benjamin Sanchez


## @knitr figure1
setEPS()
postscript("../results/figures/figure1.eps", height = 3.5, width = 8.5)
par(mfrow = c(1,3), mar = c(4,4,1,1), cex = 0.9)
plotCumulativeDistrib(ups2FC,'UPS2')
text(0.05, 0.9, 'A', pos=4, cex=1.5)
plotCumulativeDistrib(rpFC,'Ribo')
text(0.05, 0.9, 'B', pos=4, cex=1.5)
plotCumulativeDistrib(techFC,'Tech')
text(0.05, 0.9, 'C', pos=4, cex=1.5)
dev.off()


## @knitr figure2
setEPS()
postscript("../results/figures/figure2.eps", height = 3, width = 8.5)
par(mfrow = c(1,3), mar = c(4,4,1,1), cex = 1)
plotVariability(abundanceTPAN[,-1],c('.R1.1','.R2.1','.R3.1'),'',
                bquote('log'['10'] ~ '(abundance)'),
                bquote('log'['10'] ~ '(abundance)'))
text(-3.5, 3, 'A', pos=4, cex=1.5)
plotVariability(abundanceTPAN[,-1],c('_batch1','_batch2','_batch3'),'',
                bquote('log'['10'] ~ '(abundance)'),
                bquote('log'['10'] ~ '(abundance)'))
text(-3.5, 3, 'B', pos=4, cex=1.5)
plotPCA(abundanceTPAN[,-1],'',TRUE)
text(-6, 3.5, 'C', pos=4, cex=1.5)
dev.off()

