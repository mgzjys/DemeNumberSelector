#' This function allows you to calculate the maximum edge length for each triangle side
#' Edgelength()
#' @param Gene_diff genetic distance file (.diffs)
#' @param Sample_location genetic distance file (.diffs)
#' @param Automode parameter to determine whether to use the defualt settings of density clustering, TRUE or FALSE
#' @keywords Edgelength
#' @export
#' @examples
#' Edgelength()
Edgelength = function(Gene_diff, Sample_location,Automode) {
  # Loading genetic distance data
  Gene_diff <- read.table(Gene_diff, header = FALSE)
  Sample_location<-read.table(Sample_location,header = FALSE)
  
  # Loading reqiured library
  library(densityClust)
  library(rgeos)
  library(spatstat.utils)
  library(geosphere)
  
  # density clustering based on genetic distance data
  Gene_dist <- as.dist(Gene_diff)
  Gene_Clust <- densityClust(Gene_dist, gaussian = TRUE)
  Gene_Clust <-
    findClusters(
      Gene_Clust,
      rho = max(Gene_Clust$rho) / 2 ,
      delta = max(Gene_Clust$delta) / 2 ,
      plot = TRUE
    )
  rehalo <- clusters(Gene_Clust, as.list = FALSE, halo.rm = TRUE)
  
  # add reference line
  abline(v = max(Gene_Clust$rho) / 2)
  abline(h = max(Gene_Clust$delta) / 2)
  
  # plotting the clustering result
  plotMDS(Gene_Clust)
  
  # calculate the spatail mean for each cluster
  Nclusters <-
    length(Gene_Clust$clusters[!duplicated(Gene_Clust$clusters)])
  ClusterID_Coord <- cbind(rehalo, Sample_location)
  colnames(ClusterID_Coord) <- c("clusterID", "lon", "lat")
  clustercenters <- list()
  
 
  
  
  for (i in 1:Nclusters)
  {
    cluster.member <- ClusterID_Coord[ClusterID_Coord$clusterID == i,]
    cluster.center <-
      cbind(mean(cluster.member$lon), mean(cluster.member$lat))
    clustercenters <- rbind(clustercenters, cluster.center)
  }
  clustercenters<-matrix(as.numeric(clustercenters),ncol=2)
  clustercenters<-SpatialPoints(cbind(clustercenters[,1], clustercenters[,2]),CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  # Calculated the maximum edege length for each triangle side
  max_edge_length <-
    min(distm(clustercenters)[distm(clustercenters) > 0])
  return(max_edge_length)
  #max_edge_length
  
}


#' This function calculate the deme number for EEMS based on the optimized edge length calculated from function Edgelength()
#'
#' @param outlinefile outline file
#' @param max_edge_length optimized edge length calculated from function Edgelength()
#' @keywords GridGenerator
#' @export
#' @examples
#' DemeNumber()
DemeNumber=function(outlinefile, max_edge_length) {
  edge_eems_triangle <- function(x, outlinefile) {
    # Loading boundary data and sample location data
    outline <- read.table(outlinefile, header = F)
    colnames(outline) <- c("Lon", "Lat")
    
    # creating ourline polygon
    library(sp)
    library(rgeos)
    Outline_Polygon <-
      SpatialPolygons(list(Polygons(list(
        Polygon(cbind(outline$Lon, outline$Lat))
      ), "x")))
    outline_area <- gArea(Outline_Polygon)
    
    #outline_area <- areaPolygon(Outline_Polygon)
    
    
    
    
    # calculate the edge length based on triangulation method from eems
    x_span <- max(outline$Lon) - min(outline$Lon)
    y_span <- max(outline$Lat) - min(outline$Lat)
  #  
    
   # x_span <- distm(c(max(outline$Lon),max(outline$Lat)),c(min(outline$Lon),max(outline$Lat)))
   # y_span <- distm(c(max(outline$Lon),max(outline$Lat)),c(max(outline$Lon),min(outline$Lat)))
    
    
    xDemes <- as.integer(sqrt(x * x_span * x_span / outline_area))
    yDemes <- as.integer(sqrt(x * y_span * y_span / outline_area))
    
    
    scalex <- x_span / (xDemes - 0.5)
    scaley <- y_span / (yDemes - 1.0)
    
    edge1 <- distm(cbind(max(outline$Lon),min(abs(outline$Lat))+yDemes*scaley),cbind(max(outline$Lon)-(scalex/2),min(abs(outline$Lat))+(yDemes-1)*scaley))
    edge2 <- distm(cbind(max(outline$Lon),min(abs(outline$Lat))+yDemes*scaley),cbind(max(outline$Lon)-(scalex/2),min(abs(outline$Lat))+yDemes*scaley))
    
    edge1 <-as.numeric(edge1[1])
    edge2 <-as.numeric(edge2[1])
    
    #edge1 <- sqrt((scalex / 2) * (scalex / 2) + scaley * scaley)
    #edge2 <- scalex
    
    if (edge1 >= edge2)
    {
      return(edge1)
    }
    else
    {
      return(edge2)
    }
    
  }
  # calculate deme number
  deme_number <-ceiling(
    uniroot((
      function (x)
        edge_eems_triangle(x, outlinefile) - max_edge_length
    ),
    lower = 0,
    upper = 1000,
    tol = 1e-10
    )[1]$root)
  
  # plotting the result
  
  edge_eems_plot<- Vectorize(edge_eems_triangle);
  suppressWarnings(curve(
    edge_eems_plot(x,outlinefile),
    deme_number - 10,
    deme_number + 10,
    main = "deme number"
  ))
  abline(h = max_edge_length)
  abline(v = deme_number)
  points(
    deme_number,
    edge_eems_triangle(deme_number, outlinefile),
    pch = 18,
    col = "red",
    cex = 1.8
  )
  return(deme_number)
}
