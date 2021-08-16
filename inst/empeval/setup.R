
library("reshape")
library("lattice")
library("latticeExtra")
library("tactile")
library("memisc")
library("multcomp")
library("plyr")
library("glmmTMB")
library("dplyr")
trellis.par.set(list(plot.symbol = list(col="black",pch=18, cex=0.75),
                     box.rectangle = list(col=1),
                     box.umbrella = list(lty=1, col=1),
                     strip.background = list(col = "white")))
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)
