#!/usr/bin/env Rscript

### height of the "real" bar in the histograms. Higher number for higher num matrices
YHEIGHT = 25
## width of the histogram bins. For higher num matrices, reduce to increase resolution of histogram
BINWIDTH=10
library(optparse)

option_list = list(
  make_option(c("-a", "--analysis_id"),        type = "character", default = "no_id",  help = "")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)
aid = opt$analysis_id

print("...loading libraries")
suppressMessages(suppressWarnings(require(ggplot2, quietly=TRUE)))
suppressMessages(suppressWarnings(require(data.table, quietly=TRUE)))

make_plots = function(RFILE, suffix) {
  
  if (file.exists(RFILE) && file.info(RFILE)$size > 0) {
    
    ## get the distribution of background values
    print(paste("...reading results CSV file:", RFILE))
    rt <- as.data.table(read.delim(RFILE, sep="\t", header=FALSE))
    colnames(rt) <- c("reg1", "reg2","Overlap", "No_Overlap", "ID")
    suppressWarnings(rt[, EXP := paste(reg1, reg2, sep="--")])
    setkey(rt, EXP)
    
    ## get the original data values. dum0 holds the red line
    dum0 <- data.table(EXP=rep(rt[ID==0, EXP], each=2), x=rep(rt[ID==0, Overlap], each=2),y=rep(c(0,YHEIGHT), sum(rt$ID==0)))
    setkey(dum0, EXP)
    
    ## get the odds ratios
    ##Rand is mean of overlaps
    dt2 <- unique(dum0)[setkey(setnames(rt[ID != 0, mean(Overlap), by=EXP], "V1", "Rand"),EXP)]
    dt2[, odds := ifelse(Rand > 0, x / Rand, 0), by=EXP]
    setkey(dt2, odds)
    
    rt$EXP = factor(rt$EXP, levels=unique(dt2$EXP))
    dum0$EXP = factor(dum0$EXP, levels=unique(dt2$EXP))
    
    ## fit normal distribution
    ## dd is for plotting normal dist, dt2 is for odds ratio
    NUM_MATS = (length(unique(rt$ID))-1)*BINWIDTH
    dd=rt[ID > 0, as.list(MASS::fitdistr(Overlap, "normal")$estimate), by=EXP]
    dd.error=rt[ID > 0, as.list(MASS::fitdistr(Overlap, "normal")$sd), by=EXP] ## standard errors on estimates
    setnames(dd.error, c("mean", "sd"), c("mean.error","sd.error"))
    
    setkey(dd, EXP)
    setkey(dd.error, EXP)    
    setkey(dt2, EXP)
    dt2 = dt2[dd]
    setkey(dt2, EXP)
    dt2 = dt2[dd.error]
    dd <- dd[dd.error]
    setkey(dd, EXP)

    ## dd is for the blue normal fit line . the 3* sd is just for graphing purposes (how much to show)
    dd2=dd[, list(x=seq(mean - 3*sd, mean+3*sd, by=BINWIDTH), y=dnorm(x=seq(mean - 3*sd, mean+3*sd, by=BINWIDTH), mean=mean, sd=sd)*NUM_MATS), by=EXP]
    ##dd2=dd[, list(x=seq( (mean - 3 * mean.error) - 3 * (sd + 3 * sd.error), (mean + 3 * mean.error) + 3 * (sd - 3 * sd.error), by=BINWIDTH), y=dnorm(x=seq( (mean - 3 * mean.error) - 3 * (sd + 3 * sd.error), (mean + 3 * mean.error) + 3 * (sd + 3 * sd.error), by=BINWIDTH), mean=mean, sd=sd)*NUM_MATS), by=EXP]
    setkey(dd2, EXP)
    dd = dd[dd2]

    ## get the pvalues of the real overlaps using the fitted normal
    rt2 <- unique(rt[ID > 0])
    setkey(rt2, EXP)
    setkey(dt2, EXP)
    dt2[, pval := pnorm(x, mean , sd)]
    dt2[, hi_odds := x / ( (mean - 1.96 * mean.error) - 1.96 * (sd + sd.error * 1.96) )]
    dt2[, low_odds := x / ( (mean + 1.96 * mean.error) + 1.96 * (sd - sd.error * 1.96) )]
    dt2[, hi_odds_nonconservative := x / ( (mean) - 1 * (sd) )]
    dt2[, low_odds_nonconservative := x / ( (mean) + 1 * (sd) )]
    dt2$sig = "no"
    dt2$sig[dt2$hi_odds < 1] = 'depleted'
    dt2$sig[dt2$low_odds > 1] = 'enriched'
    ##dt2$sig = factor(dt2$sig, levels=c("no","depleted","enriched"))

    ## make the odds ratio plot
    setkey(dt2, odds)
    dt2$EXP = factor(dt2$EXP, levels=unique(dt2$EXP))
    write.table(dt2, row.names=FALSE, col.names=TRUE, quote=FALSE, file=paste0(aid,".odds.", suffix, ".csv"))
    
    g.odds <- ggplot(data=dt2, aes(x=EXP, y=odds, color=sig)) + geom_point() + geom_errorbar(aes(ymin=low_odds, ymax=hi_odds)) +
      theme(axis.text.x = element_text(angle=90), legend.position='none') + xlab("") + ylab("Odds over NULL") + coord_flip() +
        scale_color_manual(values=c("no"="black", "depleted"="darkred", "enriched"="darkgreen"))
    
    g.res <- ggplot() + geom_histogram(data=rt[ID > 0], aes(x=Overlap), binwidth=10) +
      theme_bw() + xlab("Overlapping events") + ylab("Count") + facet_wrap(~ EXP, scale='free', nrow=10) + geom_line(data=dum0, aes(x=x, y=y), color="red") +
        geom_line(data=dd, aes(x=x,y=y), color='blue')

    if (grepl("intra", suffix))
      pdf(paste0(aid,".results.", suffix, ".pdf"), height=8, width=6)
    else
      pdf(paste0(aid,".results.", suffix, ".pdf"), height=20, width=20)      
    print(g.res)
    dev.off()

    pdf(paste0(aid, ".odds.", suffix, ".pdf"), height=10, width=8)
    print(g.odds)
    dev.off()
    
  }
  
}

## make the results histogram
RFILE <- paste0(aid, ".results.interbin.csv")
make_plots(RFILE, "inter")
RFILE2 <- paste0(aid, ".results.intrabin.csv")
make_plots(RFILE2, "intra")
