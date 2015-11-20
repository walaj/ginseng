#!/usr/bin/env Rscript

## check for ImageMagick
if (!grepl("imagemagick", Sys.getenv("PATH"))) 
  stop("Need to have ImageMagick installed for use with R animation package. Broad: use ImageMagick")

library(optparse)

option_list = list(
    make_option(c("-i", "--input"),        type = "character", default = "animation.csv",  help = "Animation data file from swap (default: animation.csv"),
  make_option(c("-r", "--results"),        type = "character", default = "results.csv",  help = "Animation data file from swap (default: animation.csv"),
    make_option(c("-o", "--output"),       type = "character", default = "swap.gif",  help = "Output GIF to generate"),
  make_option(c("-t", "--interval"),     type = "numeric", default = 0.2,  help = "GIF interval (in seconds)"),
    make_option(c("-w", "--width"),        type = "numeric", default = 2000,  help = "Animation width (default 2000)"),
    make_option(c("-d", "--downsample"),   type = "numeric", default = 1.0,  help = "Downsample events by to this fraction of original. Default 1"),
    make_option(c("-f", "--firstandlast"), type="logical", default=FALSE, help="GIF data and final swap only"),
    make_option(c("-a", "--bedA"),         type="character", default=NULL, help="BED file for test regions A"),
    make_option(c("-b", "--bedB")        , type="character", default=NULL, help="BED file for test regions B"),
    make_option(c("-s", "--scale")        , type="numeric", default=1e6, help="Amount to round coordinates to for plotting. Default 1e6"),
    make_option(c("-k", "--height"),       type = "numeric", default = 1000,  help = "Animation height (default 1000)")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

opt$output = paste(normalizePath(opt$output))

if (!file.exists(opt$input)) {
  print(print_help(parseobj))
  stop(paste("Animation file does not exist", opt$input))
}

if (!file.exists(opt$results)) {
  print(print_help(parseobj))
  stop(paste("Results file does not exist", opt$results))
}

if (!is.null(opt$bedA) && !file.exists(opt$bedA)) {
  print(print_help(parseobj))
  stop(paste("BED file A does not exist", opt$bedA))
}

if (!is.null(opt$bedB) && !file.exists(opt$bedB)) {
  print(print_help(parseobj))
  stop(paste("BED file B does not exist", opt$bedB))
}

SCALE=opt$scale

print(opt)

print("...loading libraries")
suppressMessages(suppressWarnings(require(ggplot2, quietly=TRUE)))
suppressMessages(suppressWarnings(require(animation, quietly=TRUE)))
suppressMessages(suppressWarnings(require(gridExtra, quietly=TRUE)))
suppressMessages(suppressWarnings(require(data.table, quietly=TRUE)))

## store genome lenghts
hg19_len <- c(249250621,243199373,198022430,191154276,180915260,
              171115067,159138663,146364022,141213431,135534747,
              135006516,133851895,115169878,107349540,
              102531392,90354753,81195210,78077248,59128983,
              63025520,48129895,51304566,155270560,59373566) 

hg19_clen <- c(0, cumsum(hg19_len))
names(hg19_clen) <- c(seq(22), 'X', 'Y', 'M')

fancy_scientific <- function(l) {
  l <- format(l)
  l <- paste("10^", l)
  parse(text=l)
       # turn in to character string in scientific notation
       #l <- format(l, scientific = TRUE)
       #     # quote the part before the exponent to keep all the digits
       #     l <- gsub("^(.*)e", "'\\1'e", l)
       #     # turn the 'e+' into plotmath format
       #     l <- gsub("e", "%*%10^", l)
       #     # return this as an expression
       #     parse(text=l)
     } 

## get the csv files
#f <- dir()
#f <- f[grepl("^anim.*.csv", f)]
#if (length(f) == 0)
#  stop(paste("No animation csv files found in directory:", opt$input))

## get the temperatures
#t <- gsub("^anim[0-9]+_T_(.*?).csv$", '\\1', f)

## get the steps
#s <- as.numeric(gsub("^anim([0-9]+).*", "\\1", f))
#f <- f[order(s)]
#t <- t[order(s)]

## set the non-scrambled
print("...reading CSV file")
#bt <- read.delim(opt$input, sep=",", header=FALSE)
bt <- fread(opt$input)
setnames(bt, c("V1","V2","V3","V4","V5","V6"), c("chr1","pos1","chr2", "pos2", "count", "step")) ##, "T", "accepted","shared")
bt[, type := ifelse(step==1, "Data", "Swap")]
bt[, d := ifelse(chr1==chr2, abs(pos1-pos2),-1)]
data.ix = bt$type == "Data"
#if (any(bt$d > 0))
#  becdf <- ecdf(bt$d[data.ix & bt$d > 0])

mx = max(bt$step)

## read in the histogram data
print('...reading histogram file')
ht <- read.delim(gsub("animation.csv", "animation.histogram.csv", opt$input), sep=",", header=FALSE)
ht$step = rep(unique(bt$step), each=abs(diff(which(ht$V1==0)[1:2]))) ## add the step information
ht <- ht[ht$V1 != 250e6, ]# don't plot inter-chrom events

## read in the small histogram data
print('...reading small histogram file')
ht_s <- read.delim(gsub("animation.csv", "animation.histogram.small.csv", opt$input), sep=",", header=FALSE)
ht_s$step = rep(unique(bt$step), each=abs(diff(which(ht_s$V1==0)[1:2]))) ## add the step information
ht_s <- ht_s[ht_s$V1 != 250e6, ] # don't plot inter-chrom events

## read in the BED files
#print('...reading BEDs')
#tt <- data.frame(xmin=0,xmax=0,ymin=0,ymax=0)
#if (!is.null(opt$bedA) && !is.null(opt$bedB)) {
#  print('...reading bed file A')
#  tt_a <- read.delim(opt$bedA, sep='\t', header=FALSE, comment.char='#')
#  tt_b <- read.delim(opt$bedB, sep='\t', header=FALSE, comment.char='#')
#  tt_a$V1 <- gsub("chr", "", tt_a$V1)
#  tt_b$V1 <- gsub("chr", "", tt_b$V1)
#  tt_a <- tt_a[tt_a$V1 %in% names(hg19_clen), ]
#  tt_b <- tt_b[tt_a$V1 %in% names(hg19_clen), ]
#  print(paste("Read in", nrow(tt_a), "regions from BED A and", nrow(tt_b), "regions from BED B"))
#  tt <- data.frame(xmin=round((hg19_clen[tt_a$V1] + tt_a$V2)/SCALE),
#                   xmax=round((hg19_clen[tt_a$V1] + tt_a$V3)/SCALE),
#                   ymin=round((hg19_clen[tt_b$V1] + tt_b$V2)/SCALE),
#                   ymax=round((hg19_clen[tt_b$V1] + tt_b$V3)/SCALE))
#}

## convert to absolute dist
bt[, full_pos1 := hg19_clen[chr1+1] + pos1] ## +1 for 1 indexed vecs in R
bt[, full_pos2 := hg19_clen[chr2+1] + pos2] ## +1 for 1 indexed vecs in R

## round to nearest Mb for plotting
bt[, dsam1 := round(full_pos1/SCALE)]
bt[, dsam2 := round(full_pos2/SCALE)]

## downsample
#ne <- sum(bt$type == "Data")
#sss <- rep(FALSE, ne)
#sss[sample(ne, ceiling(opt$downsample*ne))] <- TRUE
#sssa <- rep(sss, length(unique(bt$step)))
#bt <- bt[sssa,]

## set the bounds for chr box drawing
chrdf <- data.frame(xmin=hg19_clen[1:24], xmax=hg19_clen[2:25])
chrdf$ymax = chrdf$xmax
chrdf$ymin = chrdf$xmin
chrdf <- ceiling(chrdf/SCALE)

max.d = max(bt$d[bt$d >=0])

## set the accepted data frame
#df.a <- as.data.frame(bt[!duplicated(step)])
#df.a <- df.a[, c("accepted", "step")]
#df.a$accepted = df.a$accepted/max(df.a$step)

## set the unmoved data frame
##df.u<- as.data.frame(bt[!duplicated(step)])
##df.u <- df.u[, c("unmoved", "step")]
##df.u$unmoved = df.a$accepted/max(df.u$unmoved)

## set the shared data frame
#df.s <- as.data.frame(bt[!duplicated(step)])
#df.s <- df.s[, c("shared","step")]
#df.s$shared = df.s$shared / sum(data.ix)

## make the results histogram
print("...reading results CSV file")
if (file.info(opt$results)$size > 0) {

  rt <- read.delim(opt$results, sep=",", header=FALSE)
  colnames(rt) <- c("reg1", "reg2","Overlap", "No_Overlap", "ID")
  rt$EXP = paste(rt$reg1, rt$reg2, sep="--")

  print(max(rt$ID))

  size = 1
  #id.to.keep <- c("FRAG--FRAG") #,"GENE--GENE","FRAG--GENE","GENE_SHUF--GENE_SHUF","FRAG--GENE_SHUF") #, "GENE--PROM_SHUF","GENE_SHUF--PROM_SHUF")
  id.to.keep <- unique(rt$EXP)
  rt <- rt[rt$EXP %in% id.to.keep,]
  dt.rt <- data.table(rt)
  zscore = (rt$Overlap[1] - mean(rt$Overlap[-1]))/sd(rt$Overlap)
  pvalue = max(sum(rt$Overlap[-1] > rt$Overlap[1]) / (nrow(rt)-1), 1/(nrow(rt)-1))

  #levels(rt$EXP <- factor(rt$EXP))
  #levels(rt$EXP) <- c("Fragile-Fragile", "Fragile-Gene", "Fragile-Gene Shuffle", "Gene-Gene", "Gene Shuffle-Gene Shuffle")
  dum0 <- data.frame(EXP=rep(rt$EXP[rt$ID==-1], each=2), x=rep(rt$Overlap[rt$ID==-1], each=2),y=rep(c(0,300), sum(rt$ID==-1)))
  g.res <- ggplot() + geom_histogram(data=rt[rt$ID >= 0,c(3,4,5,6)], aes(x=Overlap)) +
    ##ggtitle(paste("Data Overlap:", rt$Overlap[1], "Z-score:", zscore, "Matrices:", nrow(rt)-1, "P-val:", pvalue)) +
      theme_bw() + xlab("Overlapping events") + ylab("Count") + facet_wrap(~ EXP, scale='free', nrow=10) + geom_line(data=dum0, aes(x=x, y=y), color="blue")
      #theme(text=element_text(size+10))
  pdf("result.pdf", height=20, width=20)
  print(g.res)
  dev.off()
  pdf("Swap_result.pdf", height=3, width=9)
  print(g.res)
  dev.off()

}

MIN = 0
MAX = hg19_clen[length(hg19_clen)]/SCALE
SIZE = 1
TEXT_SIZE = 12

afunc <- function(bt, steps, data.ix) {
  #sapply(seq_along(steps), function(x) {
    
    mx.u <-max(unique(steps))
    if (steps) { #[x] == 1 && x==1) {
      print(paste("Working on DATA"))
    } else if (steps == mx.u) { ## && x==which(steps==mx.u)[1]) {
      print(paste("Working on FINAL"))
    } else if (steps > 1 && steps < mx.u) {
      print(paste("Working on step:", steps, "of", mx))
    }

    tix <- bt$step == steps
#    tix <- bt$step == steps
    ix = tix | data.ix
    temp <- paste("Temperature:", bt[ix]$T[sum(ix)])

    ## do the histogram plot
    g.dhist <- ggplot() + geom_rect(data=ht[ht$step==steps, ], aes(xmin=log(V1+1,10), xmax=log(V2+1,10), ymin=0, ymax=V3), fill=NA, color='darkred') +
      geom_rect(data=ht[ht$step==1,], aes(xmin=log(V1+1,10), xmax=log(V2+1,10), ymin=0, ymax=V3), fill=NA, color='darkgreen') +
      geom_rect(data=ht_s[ht_s$step==steps,], aes(xmin=log(V1+1,10), xmax=log(V2+1,10), ymin=0, ymax=V3), fill="red", color="red", alpha=0.2) + 
      geom_rect(data=ht_s[ht_s$step==1,], aes(xmin=log(V1+1,10), xmax=log(V2+1,10), ymin=0, ymax=V3), fill="green", color="green", alpha=0.2) +
        theme_bw() + xlab("Span") + ylab("Count") + coord_cartesian(xlim=c(2,8.5), ylim=c(0, max(ht$V3[ht$step==1])*1.1)) + scale_x_continuous(labels=fancy_scientific, breaks=1:8)
    
    ## set the orig plot
    g.orig <- ggplot() + geom_point(data=bt[data.ix,], aes(x=dsam1, y=dsam2), size=SIZE) + theme_bw() + ggtitle(paste("Original data:", sum(data.ix), "events")) +
      xlab("") + ylab("") + geom_rect(data=chrdf, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill='green', alpha=0.10, color='black') 
        #geom_rect(data=tt, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="purple", alpha=0.35, color='purple') + ## plot the bed regions
        #theme(legend.position="none") + coord_cartesian(xlim=c(MIN,MAX), ylim=c(MIN,MAX)) 
    
    ## set the swap plot
    g.swap <- ggplot() + geom_point(data=bt[tix,], aes(x=dsam1, y=dsam2), size=SIZE) + theme_bw() + ggtitle(paste("Swapped, at Temp:", temp, "Step:", steps)) +
      xlab("") + ylab("") + geom_rect(data=chrdf, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill='green', alpha=0.15, color='black') 
        #geom_rect(data=tt, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="purple", alpha=0.35, color='purple') + ## plot the bed regions
        #geom_rect(data=tt, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="yellow", alpha=0.3, color='black') + 
        #theme(legend.position="none") + coord_cartesian(xlim=c(MIN,MAX), ylim=c(MIN,MAX)) 

    ## set the accepted plot
    #df.a$type = "other";
    #df.a$type[df.a$step == steps] = "this"
    #g.accept <- ggplot() + geom_point(data=df.a, aes(x=step, y=accepted, color=type), size=3) + theme_bw() + 
    #  xlab("Steps") + ylab("Step Acceptance %") + theme(legend.position="none") + scale_colour_manual(values=c("black", "green")) + ylim(c(0,max(df.a$accepted)*1.1)) +
    #    theme(text = element_text(size=TEXT_SIZE))

    ## set the shared plot
    #df.s$type = "other";
    #df.s$type[df.s$step == steps] = "this"
    #g.shared <- ggplot() + geom_point(data=df.s, aes(x=step, y=shared, color=type), size=3) + theme_bw() + 
    #  xlab("Steps") + ylab("Shared %") + theme(legend.position="none") + scale_colour_manual(values=c("black", "green")) + ylim(c(0,max(df.s$shared)*1.1)) +
    #    theme(text = element_text(size=TEXT_SIZE)) + coord_cartesian(ylim=c(0,0.3))
    
    ## get inter/inter
    inter = sum(bt$d[tix]<0)
    intra = sum(bt$d[tix]>=0)
    ratio = inter / (intra+0.0001)

    ## plot the non-binned size histogram
    if (any(bt$d > 0)) {
      g.ohist <- ggplot(bt[tix & d > 0]) + geom_histogram(aes(x=log10(d)), fill='gray80', color='black') +
        scale_x_continuous(name="Size (bp)", breaks=0:8, label=parse(text=paste("10", 0:8, sep="^"))) + theme_bw() +
          ylab("Event count") + ggtitle("Original Histogram")
      
      g.thist <- ggplot(bt[data.ix & d > 0]) + geom_histogram(aes(x=log10(d)), fill='gray80', color='black') +
        scale_x_continuous(name="Size (bp)", breaks=0:8, label=parse(text=paste("10", 0:8, sep="^"))) + theme_bw() +
          ylab("Event count") + ggtitle("Swapped Histogram")
    }
      
    #secdf <- ecdf(bt$d[tix & bt$d >= 0])
    #sq <- c(0,1,10,100,200,500,1000,2000,5000,10000,20000,seq(from=30000, to=max.d, by=1e5))
    #oecdf <- ecdf(bt$d[bt$type=="Data" & bt$d >= 0])
    ##iiix <- c(1,sort(sample(length(sq)*2, 10000)))
    #df <- data.frame(x=c(sq,sq), y=c(secdf(sq),oecdf(sq)),type=c(rep("Swap 1M", length(sq)), rep("Swap 5M", length(sq)),rep("Data", length(sq))))
    #gs <- ggplot() + geom_line(data=df, aes(x=x, y=y, color=type), size=2) + theme_bw() + xlab("Distance") + ylab("CDF")  +
      #ggtitle(paste("Intra-", intra, "Inter-", inter, "Ratio: ", ratio)) + theme(legend.title=element_blank()) + ylim(c(0,1)) + xlim(c(0,175e6)) +
    #          theme(text = element_text(size=TEXT_SIZE)) + scale_colour_manual(values=c("darkgreen", "darkred"), name='') 

    if (any(bt$d > 0))
      suppressWarnings(grid.arrange(arrangeGrob(g.orig, g.swap, ncol=2), arrangeGrob(g.ohist, g.thist, ncol=2), g.dhist, nrow=3, ncol=1, heights=c(1/2,1/4,1/4)))
    else
      suppressWarnings(grid.arrange(arrangeGrob(g.orig, g.swap, ncol=2), g.dhist, nrow=2, ncol=1, heights=c(1/2,1/2)))
  }

## make the animation
r <- floor(1/opt$interval)
steps <- c(rep(1,r),unique(bt$step), rep(unique(bt$step)[length(unique(bt$step))], r))
if (opt$firstandlast)
    steps = c(steps[1], steps[(length(steps))])
#ppdf(print(afunc(bt, steps[1], data.ix)))
pdf("swap.pdf")
print(afunc(bt, steps[length(steps)], data.ix))
dev.off()
#saveGIF(afunc(bt, steps, data.ix), movie.name=opt$output, interval=opt$interval, ani.width=opt$width, ani.height=opt$height, ani.dev='pdf')

#####
if (FALSE)
  {
    system.time(bt <- read.delim(gzfile("/xchip/gistic/Jeremiah/Projects/Significance/Sanger2500/Runs/150427/all.matrices.csv.gz"), sep="\t", header=FALSE)    )
    dt <- data.table(bt)
    setnames(dt, c("V1","V2","V3","V4","V5"), c("r_chr","r_pos","c_chr","c_pos","id"))
    setkey(dt, id, r_chr, r_pos)
  }

