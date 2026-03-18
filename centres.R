set_centre <- function(id, name, col, cutoff=NA){
    return(list(id=id, name=name, col=col, cutoff=cutoff))
}

centre <- list()
centre[[length(centre)+1]] <- set_centre('babj', 'CMA', 'orange')
centre[[length(centre)+1]] <- set_centre('sbsj', 'CPTEC/INPE', 'brown')
centre[[length(centre)+1]] <- set_centre('cwao', 'ECCC', 'red', cutoff=16)
centre[[length(centre)+1]] <- set_centre('ecmf', 'ECMWF', 'black', cutoff=21)
centre[[length(centre)+1]] <- set_centre('edzw', 'DWD', 'cyan')
centre[[length(centre)+1]] <- set_centre('gfdl', 'GFDL', 'darkgreen')
centre[[length(centre)+1]] <- set_centre('rjtd', 'JMA', 'purple')
centre[[length(centre)+1]] <- set_centre('rksl', 'KMA/KIAPS', 'magenta')
centre[[length(centre)+1]] <- set_centre('kwbc', 'NOAA', 'green')
centre[[length(centre)+1]] <- set_centre('rums', 'RAS', 'pink')
centre[[length(centre)+1]] <- set_centre('egrr', 'UKMO', 'blue')
