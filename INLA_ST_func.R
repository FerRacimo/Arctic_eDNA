
library("tidyr")
library("INLA")
library("splancs")
library("rgeos")
library("ggplot2")
library("reshape2")
library("dplyr")
library("gridExtra")
library("deldir")
library("viridis")
library(rnaturalearth)
library(ggspatial)
source("spde-book-functions.R")

# Find nearest value in an array
FindNearest <- function(array, value){
    dists <- abs((array - value))
    idx <- which(dists == min(dists))[1]
    nearest <- array[idx]
    return(nearest)
}

MinMaxNorm <- function(x,a=0,b=1){
    numerator <- (x - min(x,na.rm=TRUE))*(b-a)
    denominator <- max(x,na.rm=TRUE) - min(x,na.rm=TRUE)
    final <- a + numerator/denominator
    return(final)
}


LonLatToAzi <- function(test){
    colnames(test) <- c("Lon","Lat")
    test <- data.frame(test)
    coordinates(test) <- ~Lon+Lat
    proj4string(test) <- CRS("+init=epsg:4326") 
    CRS.new <- CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
    test2 <- data.frame(spTransform(test, CRS.new))
    test2 <- cbind(test2[,1],test2[,2])
    return(test2)
}


AziToLonLat <- function(m1){
    colnames(m1) <- c("X","Y")
    coordinates(m1) <- ~X+Y
    proj4string(m1) <- CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
CRS.new <- CRS("+init=epsg:4326")
    m2 <- data.frame(spTransform(m1, CRS.new))
    m2 <- cbind(m2[,1],m2[,2])
    colnames(m2) <- c("Lon","Lat")
    m2 <- data.frame(m2)
    return(m2)
}


RunInlaBin <- function(var=NA,sitetab=NA,xytab=NA,mesh.s=NA,mesh.t=NA,Ntrials=1,namescov=c(),normcov=FALSE){

    print(var)
    sitetab <- data.frame(sitetab)
    #Z <- c(sitetab[,var])
    n <- dim(sitetab)[1]

    spde <- inla.spde2.pcmatern(mesh = mesh.s, 
        prior.range = c(3000, 0.5), # P( range < 3000) = 0.5
        prior.sigma = c(10, 0.01)) # P( sigma > 10) = 0.01

    iset <- inla.spde.make.index('i', n.spde = spde$n.spde, n.group = k)

    A <- inla.spde.make.A(mesh = mesh.s, loc = xytab, group = sitetab$Age, group.mesh = mesh.t) 

    #namescov <- c("annualtemp","annualprec","tempseaso","precseaso","human")
    tablist <- as.list(sitetab)
    covlist <- tablist[namescov]

   # Normalize covariates to (-1,1) if (0,1); Stanrdize them otherwise  
    if(normcov == TRUE){
    	covlist <- lapply(covlist,function(x){
       		if( all(unique(x) %in% c(0,1,NA)) ){ return(MinMaxNorm(x,-1,1))
        	} else{ return(c(base::scale(x))) }
    	})
    }

    sdat <- inla.stack( 
    	 data = tablist[var], 
    	 A = list(A,1), 
    	 effects = list( iset, c( list(b0 = rep(1, n)), covlist ) ), 
    	 tag = "stdata")

    # Prior for temporal auto-regressive parameter
    h.spec <- list(theta = list(prior = 'pccor1', param = c(0, 0.9)))

    # Precision prior for default data model - NOT USED
    #prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))

    # Gaussian process model formula
    formula <- as.formula( paste(var, " ~ 0 + b0 + ",paste(namescov, collapse=" + ")," + f(i, model = spde, group = i.group,control.group = list(model = 'ar1', hyper = h.spec))",sep=""))

    # Binomial data model
    res <- inla(formula,
    	family="binomial", Ntrials=1, data = inla.stack.data(sdat), 
    	control.family = list(),
	control.predictor = list(link=1,compute = TRUE,A = inla.stack.A(sdat)),
	control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))

    return(res)
}

# Obtain projected mesh
ProjMesh <- function(res,xytab,mesh.s,mesh.t,spde.s,PlotBox=NULL){

    if(is.null(PlotBox)){
	limx = range(xytab[, 1])
    	limy = range(xytab[, 2])
    } else {
      	limx = range(PlotBox[,1])
	limy = range(PlotBox[,2])
    }      	 
    
    r0 <- diff(limx) / diff(limy)
    prj <- inla.mesh.projector(mesh.s, xlim = limx,
        ylim = limy, dims = c(100 * r0, 100))
    #in.pr <- inout(prj$lattice$loc, xytab)
    m.prj <- lapply(1:mesh.t$n, function(j) {
        idx <- 1:spde.s$n.spde + (j - 1) * spde.s$n.spde
        #print(idx)
	#valid <- which(sign(res$summary.fixed[,"0.025quant"]) == sign(res$summary.fixed[,"0.975quant"]))
	#field <- sum(res$summary.fixed[valid,"mean"]) + res$summary.ran$i$mean[idx]
	field <- res$summary.fixed[1,"mean"] + res$summary.ran$i$mean[idx]
        r <- inla.mesh.project(prj,field = field)
        #r[!in.pr] <- NA
        prob <- exp(r)/(1 + exp(r))
        return(prob)
    })
    return(list(m.prj,prj))
}






PlotPolar <- function(m1,title=NA,minplot=NA,maxplot=NA){

# custom theme without axes and annotations
theme_polar <- function(){
   list(
      theme_bw(base_size=10),
      theme(
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text = element_blank(),
         axis.ticks = element_blank()),
      labs(x='',y=''))
}

dt <- rnaturalearth::ne_coastline()
clip_boundary <- sp::SpatialPolygons(
  list(sp::Polygons(
    list(sp::Polygon(
      data.frame(lon = c(-180, 180, 180, -180), lat = c(60, 60, 90, 90)))), ID = 1)
  ), proj4string = sp::CRS(sp::proj4string(dt)))

arctic <- raster::crop(dt, clip_boundary)
arctic <- sp::spTransform(arctic, sp::CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
arctic_sf <- sf::st_as_sf(arctic)

ggplot(data=m1,aes(x = Lon, y = Lat, z = Prob,fill=Prob)) + 
  #geom_point(data = m1, aes(x = Lon, y = Lat, fill = Prob,colour=Prob)) +
  geom_tile(aes(fill = Prob)) +
  ggspatial::layer_spatial(data = arctic_sf,colour="white") +
  ggtitle(title) +
  scale_fill_viridis(limits = c(minplot,maxplot)) +
  theme_polar()
}


