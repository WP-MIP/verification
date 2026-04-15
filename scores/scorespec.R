set.specs <- function(name, units, range, mult=1, hline=NA, leg='topleft', ncol=2){
    return(list(name=name, units=units, range=range, mult=mult, hline=hline, leg=leg, ncol=ncol))
}

spec <- list()

v <- 't'
spec[[v]] <- list()
spec[[v]][['rmse']] <- set.specs("RMSE", "K", c(0,4.5), leg='topleft', ncol=2)
spec[[v]][['rmse-ref']] <- set.specs("RMSE", "K", c(-0.7,0.7), hline=0, leg='bottomleft', ncol=3)
spec[[v]][['bias']] <- set.specs("Bias", "K", c(-0.6,0.6), hline=0, ncol=3)
spec[[v]][['acor']] <- set.specs("Anomaly Correlation", NA, c(0,1), hline=1)

v <- 'gh'
spec[[v]] <- list()
spec[[v]][['rmse']] <- set.specs("RMSE", "dam", c(0,12), mult=0.1)
spec[[v]][['rmse-ref']] <- set.specs("RMSE", "dam", c(-2,2), mult=0.1, hline=0, leg='bottomleft')
spec[[v]][['bias']] <- set.specs("Bias", "dam", c(-3,3), mult=0.1, hline=0)
spec[[v]][['acor']] <- set.specs("Anomaly Correlation", NA, c(0,1), hline=1)


