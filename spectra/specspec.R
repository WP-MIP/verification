set.specs <- function(name, yrange, xrange=c(1,700), leg='bottomleft', ncol=2, log='xy', hline=NA){
    return(list(name=name, yrange=yrange, xrange=xrange, leg=leg, ncol=ncol, log=log, hline=hline))
}

spec <- list()

v <- 'ke'
spec[[v]] <- list()
spec[[v]][['en']] <- set.specs(expression(paste("Kinetic Energy (",m^2, s^-2, ")", sep='')), c(10e-7,100))
spec[[v]][['en-ratio']] <- set.specs("Kinetic Energy Amplitude Ratio", c(0,2), log='x', leg='topleft', hline=1)
spec[[v]][['enrot']] <- set.specs(expression(paste("Rotational Kinetic Energy (",m^2, s^-2, ")", sep='')), c(10e-7,100))
spec[[v]][['endiv']] <- set.specs(expression(paste("Divergence Kinetic Energy (",m^2, s^-2, ")", sep='')), c(10e-7,100))
