set.specs <- function(name, units, yrange, xrange=c(1,700), leg='bottomleft', ncol=2, log=NA){
    return(list(name=name, units=units, yrange=yrange, xrange=xrange, leg=leg, ncol=ncol, log=log))
}

spec <- list()

v <- 'ke'
spec[[v]] <- list()
spec[[v]][['en']] <- set.specs("Kinetic Energy", "m^2 s^-1", c(10e-7,100), log='xy')
spec[[v]][['en-ratio']] <- set.specs("Kinetic Energy Ratio", NA, c(0,2))
spec[[v]][['enrot']] <- set.specs("Rotational Kinetic Energy", "m^2 s^-1", c(10e-7,100), log='xy')
spec[[v]][['endiv']] <- set.specs("Divergent Kinetic Energy", "m^2 s^-1", c(10e-7,100), log='xy')
