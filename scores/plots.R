library(jsonlite)
library(boot)

# Basic setup
source("../centres.R")   #defines "centres"
source("scorespec.R")    #defines "spec"
set.seed(123)

# User configuration
field <- "t"
level <- 500
lead <- 240
scores <- c('rmse-ref', 'rmse', 'bias')

# Utility functions
meanfun <- function(data, i){
    return(mean(data[i]))
}
tcol <- function(col, trans=0.8){
    rgb.val <- col2rgb(col)
    return(rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max=255, alpha=(1-trans) * 255))
}
compute.stats <- function(d.fcst, d.ref, is.ref, lscore, fcsts){
    ci.type <- "basic"
    fmean <- c()
    fmean.low <- c()
    fmean.high <- c()
    fcsts <- as.integer(names(d.fcst[[lscore]]))
    for (fcst in fcsts){
        sfcst <- as.character(fcst)
        d <- d.fcst[[lscore]][[sfcst]]
        if (is.ref){d <- d - d.ref[[lscore]][[sfcst]]}            
        draws <- boot(d, statistic=meanfun, R=1000)
        ci <- boot.ci(draws, 0.95, ci.type)[[ci.type]]
        fmean <- c(fmean, draws$t0)
        fmean.low <- c(fmean.low, ci[4])
        fmean.high <- c(fmean.high, ci[5])
    }
    return(list('mean'=fmean, 'low'=fmean.low, 'high'=fmean.high, 'fcsts'=fcsts))
}
plot.line <- function(fmean, col, lwd, lty=1){
    polygon(c(fmean$fcsts, rev(fmean$fcsts)), c(fmean$low, rev(fmean$high)), col=tcol(col), border=NA)
    lines(fmean$fcsts, fmean$mean, col=col, lwd=lwd, lty=lty)
}

# Main plotting function
genPlot <- function(score, stream, mtype){

    # Plot properties
    charsize <- 1.4
    lwd <- 3
    
    # Prepare output file
    ofile <- paste(paste(score, stream, mtype, sep='_'), 'pdf', sep='.')
    pdf(ofile, paper='special', width=8, height=5, bg='white')
    default <- par(mar=c(5,5,4,2))

    # Generate plot background
    units <- spec[[field]][[score]]$units
    ylab <- spec[[field]][[score]]$name
    if (!is.na(units)){
        ylab <- paste(ylab, ' (', spec[[field]][[score]]$units, ')', sep='')
    }
    plot(c(0, lead), c(0, 0), type='n', xaxt='n', yaxt='n', xaxs='i',
         xlab='Forecast Hour', ylab=ylab, cex.lab=charsize, ylim=spec[[field]][[score]]$range)
    axis(1, at=seq(0, lead, by=24), cex.axis=charsize)
    axis(2, cex.axis=charsize, las=1)
    hline <- spec[[field]][[score]]$hline
    if (!is.na(hline)) abline(h=hline, lwd=2, lty=2)
    leg.centre <- c()
    leg.col <- c()

    # Retrieve reference data
    is.ref <- grepl('-ref', score)
    lscore <- sub('-ref', '', score)
    fname <- paste(paste(field, level, 'oic', 'ecmf', 'pm', '00', sep='_'), 'json', sep='.')
    d.ref <- fromJSON(fname)

    # Plot reference data
    if (! is.ref && mtype != 'pm'){
        fmean <- compute.stats(d.ref, d.ref, FALSE, lscore)
        plot.line(fmean, 'grey', 1, lty=5)
    }
    
    # Process data by centre
    for (c in seq(1, length(centre))){
        
        # Retrieve input data
        lstream <- stream
        if (centre[[c]]$id == 'ecmf'){lstream <- 'oic'}
        fname <- paste(paste(field, level, lstream, centre[[c]]$id, mtype, '00', sep='_'), 'json', sep='.')
        if (! file.exists(fname)){ next }
        d.fcst <- fromJSON(fname)
        
        # Boostrap for confidence intervals
        fmean <- compute.stats(d.fcst, d.ref, is.ref, lscore)
        
        # Add lines
        plot.line(fmean, centre[[c]]$col, lwd)

        # Set legend values
        leg.centre <- c(leg.centre, centre[[c]]$name)
        leg.col <- c(leg.col, centre[[c]]$col)

    }
    
    # Finalize plot
    if (length(leg.centre) > 0){
        legend(x=spec[[field]][[score]]$leg, inset=0.01, legend=as.matrix(leg.centre), col=leg.col, lwd=lwd,
               cex=1+(charsize-1)/3, ncol=spec[[field]][[score]]$ncol)
    }
    par(default)
    dev.off()
}

# Main calculations

genPlot('rmse', 'oic', 'ai')
genPlot('bias', 'oic', 'ai')
stop()

for (stream in c('oic', 'sic')){
    for (mtype in c('ai', 'hy', 'pm')){
        for (score in scores){
            genPlot(score, stream, mtype)
        }
    }
}
