#!/usr/bin/env Rscript
# Figures from the manuscript
# Benjamin Sanchez


## @knitr figure1
setEPS()
postscript("../results/figures/figure1.eps", height = 3.5, width = 8.5)
par(mfrow = c(1,3), mar = c(4,4,1,1), cex = 0.9)
plotCumulativeDistrib(ups2FC,'')
plotCumulativeDistrib(rpFC,'')
plotCumulativeDistrib(techFC,'')
dev.off()


## @knitr figure2
setEPS()
postscript("../results/figures/figure2.eps", height = 3, width = 8.5)
par(mfrow = c(1,3), mar = c(4,4,1,1), cex = 1)
plotVariability(abundanceMSR[,-1],c('.R1.1','.R2.1','.R3.1'),'',
                bquote('log'['10'] ~ '(abundance)'),
                bquote('log'['10'] ~ '(abundance)'))
plotVariability(abundanceMSR[,-1],c('_batch1','_batch2','_batch3'),'',
                bquote('log'['10'] ~ '(abundance)'),
                bquote('log'['10'] ~ '(abundance)'))
plotPCA(abundanceMSR[,-1],'',TRUE)
dev.off()

