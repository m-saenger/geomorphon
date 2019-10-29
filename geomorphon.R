##' ...............................................................................
##'  Landform Classification after Jasiewicz and Stepinkski (2013) 
##' ...............................................................................
##'  M. Sänger  10/2019
##' ...............................................................................
##' 
##' Raster data to be on equally spaced grid (e.g. km)
##' Output very sensitive to flatness threshold
##' Function not optmised for performance - run time very long for large data sets


## ------------------------------------- Definitions --------------------------------------

geomorph.def <- data.frame(
  num_lf = 1:10,
  id_lf = c("PK", "RI", "SH", "SP", "SL", "FS", "FL", "HL", "VL", "PT"),
  name_en = c("Peak", "Ridge", "Shoulder", "Spur", "Slope", "Footslope", "Flat", "Hollow", "Valley", "Pit"),
  colour = c("magenta", "red", "orange", "yellow", "grey40",  "grey70", "grey90", "skyblue1", "dodgerblue", "royalblue3"),
  stringsAsFactors = F
)

geomorph.lut <- data.frame(
  V0 = c("FL", "FL", "FL", "SH", "SH", "RI", "RI", "RI", "PK"),
  V1 = c("FL", "FL", "SH", "SH", "SH", "RI", "RI", "RI", NA),
  V2 = c("FL", "FS", "SL", "SL", "SP", "SP", "RI", NA, NA),
  V3 = c("FS", "FS", "SL", "SL", "SL", "SP", NA, NA, NA),
  V4 = c("FS", "FS", "HL", "SL", "SL", NA, NA, NA, NA),
  V5 = c("VL", "VL", "HL", "HL", NA, NA, NA, NA, NA),
  V6 = c("VL", "VL", "VL", NA, NA, NA, NA, NA, NA),
  V7 = c("VL", "VL", NA, NA, NA, NA, NA, NA, NA),
  V8 = c("PT", NA, NA, NA, NA, NA, NA, NA, NA)
)
geomorph.lut <- as.matrix(geomorph.lut)
geomorph.lut.num <- matrix(match(geomorph.lut, geomorph.def$id_lf), nrow = nrow(geomorph.lut)) 

## ------------------------------------- Function --------------------------------------

geomorph <- function(x, flatness.thresh = NA, res = NA, ncell = NA, geomorph.lut.num, verbose = F){
  #' @description Note, that no performance optimisation has been done to this function, yet.
  #' @author M. Sänger 2018
  #' @source Jasiewicz, Stepinkski 2013
  #' @param x vector from raster::focal function
  #' @param flatness.thresh Flatness threshold in degrees
  #' @param res resolution, same unit as values in r
  #' @param ncell total number of cells in raster (used for progress bar only)
  #' @param geomorph.lut.num look-up table to derive landform class from ternary patterns

  # Breaks for flatness threshold
  brks <- c(-Inf, -flatness.thresh, flatness.thresh, Inf)
  brks.ind <- c(-1, 0, 1)
  
  # Progress bar (if verbose = TRUE)
  if(verbose){
    counter <<- counter + 1
    if(counter %in% round(seq(0, ncell, length.out = 20))) 
      cat("=", round(counter/ncell*100)); if(counter == ncell) cat("\n")
  } 
  
  # Create matrix from incoming vector x
  size = sqrt(length(x))
  m <- matrix(x, nrow = size)
  
  # Distance from central point to edge (number of cells)
  mid <- ceiling(size/2)
  
  # Matrix of all vectors from the central point to the octants
  oct <- rbind(
    ne = cbind(mid:size, mid:size),
    e = cbind(mid:size, mid),
    se = cbind(mid:size, mid:1),
    s = cbind(mid, mid:1),
    sw = cbind(mid:1, mid:1),
    w = cbind(mid:1, mid),
    nw = cbind(mid:1, mid:size),
    n = cbind(mid, mid:size)
  )
  
  # Coordinates and cell distance (sqrt(2) for diagonals)
  oct.vector <- m[oct]
  cell.scaling <- rep(c(sqrt(2), 1), 4) # Horizontal cell distance in all 8 directions
  cell.size <- res * cell.scaling
  
  # Matrix octants vs. cell values
  m1 <- matrix(oct.vector, nrow = 8, byrow = T)
  
  # z diff from central point
  m.diff <-  m1[, -1] - m1[, 1]
  
  # Calculate slope angle and transform to degrees
  m.slope <- atan(m.diff/(cell.size * 1:ncol(m.diff)))
  m.angle <- m.slope * 180/pi
  
  # Calculate zenith and nadir angles for each octant
  nadir <- 90 + apply(m.angle, 1, min, na.rm = T)
  zenith <- 90 - apply(m.angle, 1, max, na.rm = T)
  
  # Derive ternary pattern
  ternary.pattern <- brks.ind[findInterval(nadir - zenith, brks)]
  
  plus.ind <- length(which(ternary.pattern == 1))
  neg.ind <- length(which(ternary.pattern == -1))
  
  # Look up ternarity pattern and assign landform class
  geomorph.lut.num[neg.ind + 1, plus.ind + 1]  
}

## ------------------------------------- Example --------------------------------------

library(raster)
library(tidyverse)
library(reshape2)

# Definition
focal.window.size <- 7
flatness.thresh <- 1
counter <- 1
verbose <- FALSE

# Data
data("volcano")
dat <- volcano/1e2 # Scale z axis to get meaningful slope values
r <- raster(t(dat), xmn=0, xmx=nrow(dat), ymn=0, ymx=ncol(dat))

# Focal function
focal.function <- function(x){
  geomorph(x,  flatness.thresh, res(r), ncell(r), geomorph.lut.num, verbose = verbose)
}
# Focal matrix
focal.matrix <- matrix(1, nrow = focal.window.size, ncol = focal.window.size)

# Apply focal function
r.volcano.lf <- raster::focal(r, fun = focal.function, w = focal.matrix, pad = T, padValue = NA)

# Melt raster to data frame
df.lf <- r.volcano.lf %>% 
  flip(direction = 2) %>% 
  as.data.frame(xy = T)

ggplot(df.lf, aes(x, y, fill = factor(layer, geomorph.def$num_lf))) +
  geom_raster() +
  scale_fill_manual("Landform", values = geomorph.def$colour, labels = geomorph.def$name, drop = F) +
  coord_cartesian(expand = F) +
  theme_bw(16) +
  theme(plot.background = element_rect(colour = "black", size = .3)) +
  labs(x = "x", y = "y", title = "Volcano Landform Classes")


