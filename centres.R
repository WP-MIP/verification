set_centre <- function(id, name, col, cutoff=NA){
    return(list(id=id, name=name, col=col, cutoff=cutoff))
}

centre <- list()
centre[[length(centre)+1]] <- set_centre('babj', 'CMA', '#6d5d01')
centre[[length(centre)+1]] <- set_centre('sbsj', 'CPTEC/INPE', '#e69b27')
centre[[length(centre)+1]] <- set_centre('cwao', 'ECCC', '#f62323', cutoff=16)
centre[[length(centre)+1]] <- set_centre('ecmf', 'ECMWF', '#000098', cutoff=21)
centre[[length(centre)+1]] <- set_centre('edzw', 'DWD', '#28be9c')
centre[[length(centre)+1]] <- set_centre('gfdl', 'GFDL', '#006f00')
centre[[length(centre)+1]] <- set_centre('rjtd', 'JMA', '#856feb')
centre[[length(centre)+1]] <- set_centre('rksl', 'KMA/KIAPS', '#ac02d6')
centre[[length(centre)+1]] <- set_centre('kwbc', 'NOAA', '#67b320')
centre[[length(centre)+1]] <- set_centre('rums', 'RAS', '#c5bf46')
centre[[length(centre)+1]] <- set_centre('egrr', 'UKMO', '#307391', cutoff=21)
