library(jsonlite)
library(boot)

# Basic setup
source("../centres.R")   #defines "centres"
source("specspec.R")     #defines "spec"
set.seed(123)

# User configuration
field <- "ke"
level <- 250
lead <- 240
nboot <- 1000

# Utility functions
meanfun <- function(data, i){
    return(mean(data[i]))
}
tcol <- function(col, trans=0.8){
    rgb.val <- col2rgb(col)
    return(rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max=255, alpha=(1-trans) * 255))
}
myfloor <- function(x){return(max(x, 1e-10))}
refslope <- function(wn, slope, pt){return((pt[2] / pt[1]**slope) * wn**slope)}
towl <- function(wn){return(2*pi*6371/wn)}
get_axp <- function(x){return(10^c(ceiling(x[1]), floor(x[2])))}
get_fname <- function(centre_id, stream, mtype){
    lstream <- stream
    if (centre_id == 'ecmf'){lstream <- 'oic'}
    fnames <- Sys.glob(paste(paste('spec', field, level, lstream, centre_id, mtype, '??', lead, sep='_'), 'json', sep='.'))
    if (length(fnames) < 1){return('No_File')}
    return(fnames[1])
}
compute.stats <- function(d.fcst, d.ref, is.ratio, lscore){
    ci.type <- "basic"
    small <- 1e-10
    fmean <- c()
    fmean.low <- c()
    fmean.high <- c()
    wnlist <- c()
    wns <- seq(1, ncol(d.fcst[[lscore]]))
    for (wn in wns){
        d <- d.fcst[[lscore]][,wn]
        if (is.ratio){d <- sapply(sqrt(d), myfloor) / sapply(sqrt(d.ref[[lscore]][,wn]), myfloor)}
        draws <- boot(d, statistic=meanfun, R=nboot)
        fmean <- c(fmean, draws$t0)
        wnlist <- c(wnlist, d.fcst[['deg']][wn])
        if (length(unique(d)) == 1){
            fmean.low <- c(fmean.low, draws$t0)
            fmean.high <- c(fmean.high, draws$t0)
        } else {
            ci <- boot.ci(draws, 0.95, ci.type)[[ci.type]]
            fmean.low <- c(fmean.low, ci[4])
            fmean.high <- c(fmean.high, ci[5])
        }
    }
    return(list('mean'=fmean, 'low'=fmean.low, 'high'=fmean.high, 'wns'=wnlist))
}
plot.line <- function(fmean, col, lwd, lty=1){
    polygon(c(fmean$wns, rev(fmean$wns)), c(fmean$low, rev(fmean$high)), col=tcol(col), border=NA)
    lines(fmean$wns, fmean$mean, col=col, lwd=lwd, lty=lty)
}

# Main plotting function
genPlot <- function(score, stream, mtype){

    # Plot properties
    charsize <- 1.4
    lwd <- 3
    
    # Convert wavenumber to wavelength
    xlims <- towl(spec[[field]][[score]]$xrange)
    usr.i <- log10(xlims)
    wl.ticks <- axTicks(3, usr=usr.i, axp=c(get_axp(usr.i), n=3), log=TRUE, nintLog=5)
    
    # Prepare output file
    ofile <- paste(paste(score, stream, mtype, as.character(lead), sep='_'), 'pdf', sep='.')
    pdf(ofile, paper='special', width=8, height=5, bg='white')
    default <- par(mar=c(5,7,4,2))

    # Generate plot background
    plot(c(1, 2), c(1, 2), type='n', xaxt='n', yaxt='n', xaxs='i', log=spec[[field]][[score]]$log,
         xlab='Global Spherical Wavenumber', ylab=NA, cex.lab=charsize, xlim=spec[[field]][[score]]$xrange,
         ylim=spec[[field]][[score]]$yrange)
    axis(1, cex.axis=charsize)
    axis(2, cex.axis=charsize, las=1)
    axis(3, at=towl(wl.ticks), labels=wl.ticks)
    mtext('km', side=3, at=max(spec[[field]][[score]]$xrange - 65), line=1)
    mtext(spec[[field]][[score]]$name, 2, line=5, cex=charsize)
    hline <- spec[[field]][[score]]$hline
    if (!is.na(hline)) abline(h=hline, lwd=2, lty=2)
    leg.centre <- c()
    leg.col <- c()
    
    # Retrieve reference data
    is.ratio <- grepl('-ratio', score)
    lscore <- sub('-ratio', '', score)
    d.ref <- fromJSON(get_fname('ecmf', 'oic', 'pm'))

    # Add filter cutoffs if needed
    if (mtype == 'hy' && is.ratio){
        for (c in seq(1, length(centre))){
            if (!is.na(centre[[c]]$cutoff) && file.exists(get_fname(centre[[c]]$id, stream, mtype))){
                abline(v=centre[[c]]$cutoff, col=centre[[c]]$col, lty=3, lwd=2)
            }
        }
    }
    
    # Generate reference lines
    if (! is.ratio){
        wnl <- seq(20, 100)
        kel <- refslope(wnl, -3, c(min(wnl),20))
        lines(wnl, kel)
        text(50, 5, labels='-3', cex=1.2)
        wnl <- seq(120,600)
        kel <- refslope(wnl, -5/3, c(min(wnl),0.1))
        text(270, 0.1, labels='-5/3', cex=1.2)
        lines(wnl, kel)
    }
    
    # Plot reference data
    if (! is.ratio && mtype != 'pm'){
        fmean <- compute.stats(d.ref, d.ref, FALSE, lscore)
        plot.line(fmean, 'grey', 3, lty=5)
    }
    
    # Process data by centre
    for (c in seq(1, length(centre))){

        # Retrieve input data
        fname <- get_fname(centre[[c]]$id, stream, mtype)
        if (! file.exists(fname)){ next }
        d.fcst <- fromJSON(fname)
        
        # Boostrap for confidence intervals
        fmean <- compute.stats(d.fcst, d.ref, is.ratio, lscore)

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

#genPlot('en-ratio', 'oic', 'hy')
#stop()


for (stream in c('oic', 'sic')){
    for (mtype in c('ai', 'hy', 'pm')){
        for (score in c('en', 'en-ratio', 'enrot', 'endiv')){
            genPlot(score, stream, mtype)
        }
    }
}
