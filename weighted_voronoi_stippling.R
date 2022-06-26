##Ressource: http://dahtah.github.io/imager/stippling.html

library(dplyr)
library(imager)

##Compute Voronoi diagram of point set xy with width w and height h.
voronoi <- function(xy,w,h)
{
  v <- imfill(w,h)
  ind <- round(xy) %>% index.coord(v,.)
  v[ind] <- seq_along(ind)
  d <- distance_transform(v>0,1)
  watershed(v,-d,fill_lines=FALSE)
}

##Compute Voronoi diagram of point set xy, and return center of mass of each cell (with density given by image dens)
cvt <- function(xy,dens)
{
  voronoi(xy,width(dens),height(dens)) %>% as.data.frame %>%
    mutate(vim=c(dens)) %>%
    group_by(value) %>%
    dplyr::summarise(x=weighted.mean(x,w=vim),y=weighted.mean(y,w=vim)) %>%
    select(x,y) %>%
    filter(x %inr% c(1,width(dens)),y %inr% c(1,height(dens)))
}

##Initializing sample Points for a given Image and loops the cvt calculations.
##nDots: Amount of Sample Points, gamma: Density of Input Values can be changed, nSteps: Number of Voronoi Steps
stipple <- function(im,nDots=1e3,gamma=1,nSteps=5)
{
  dens <- (1-im)^gamma
  xy <- sample(nPix(im),nDots,replace=TRUE,prob=dens) %>% coord.index(im,.) %>% select(x,y)
  
  for (i in 1:nSteps)
  {
    xy <- cvt(xy,dens)

  }
  xy
}

##Use of Algorithm in the following lines:

##Load Image and safe it grayscaled
im <- load.image("Input/Zebra.jpg")
im <- grayscale(im) 

##Stippling Algorithm with the use of weighted centroid coronoi diagrams
out <- stipple(im,nDots=50000,gamma=2,nSteps=2)

##Save image as pdf
pdf("Output/Zebra.pdf") 
plot(out,ylim=c(height(im),1),main=paste("Zebra mit 50.000 Punkten und 2 Voronoi Durchläufen"), xlim=c(1, width(im)), cex=.1,pch=19,axes=FALSE,xlab="",ylab="")
dev.off() 
