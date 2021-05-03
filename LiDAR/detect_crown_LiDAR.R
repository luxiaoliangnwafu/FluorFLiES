# Attach the 'ForestTools' and 'raster' libraries
library(ForestTools)
library(raster)
library(stringr)

# functions ---------------------------------------------------------------
detect.tree <- function(fname, label){
  r <- raster(fname)
  # change extent
  extent(r) <- extent(1, 30, 1, 30)
  chm <- raster(nrow = 30, ncol = 30)
  extent(chm) <- extent(r)
  chm <- resample(r, chm, method = 'bilinear') # resample output
  
  plot(chm)
  
  
  # Remove plot margins (optional)
  #par(mar = rep(0.5, 4))
  
  # Plot CHM (extra optional arguments remove labels and tick marks from the plot)
  # plot(chm, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
  
  lin <- function(x){x * 0.05 + 0.6}
  ttops <- vwf(CHM = chm, winFun = lin, minHeight = 2)
  
  # Create polygon crown map
  crownsPoly <- mcws(treetops = ttops, CHM = chm, format = "polygons", minHeight = 1.5, verbose = FALSE)
  crownsPoly[["crownDiameter"]] <- sqrt(crownsPoly[["crownArea"]]/ pi) * 2
  
  # Plot tree area
  fname = paste0(label, "_CHM", ".png")
  png(fname)
  plot(chm, xlab = "", ylab = "", xaxt='n', yaxt = 'n', main = fname)
  plot(crownsPoly, border = "blue", lwd = 0.5, add = TRUE)
  plot(ttops, col = "blue", pch = 20, cex = 0.5, add = TRUE)
  dev.off()
  
  # Create crown map
  crowns <- mcws(treetops = ttops, CHM = chm, minHeight = 1.5, verbose = FALSE)
  fname = paste0(label, "_CHM_SHAPE", ".png")
  png(fname)
  plot(crowns, col = sample(rainbow(50), length(unique(crowns[])), replace = TRUE), 
       legend = FALSE, xlab = "", ylab = "", xaxt='n', yaxt = 'n', main = fname)
  dev.off()
  
  
  # And then convert it (back) to a data.frame
  df.xy <- coordinates(ttops)
  df.xy <- as.data.frame(df.xy)
  df.tree <- data.frame(df.xy, height = crownsPoly[["height"]], radius = crownsPoly[["crownDiameter"]], area = crownsPoly[["crownArea"]])
  colnames(df.tree) <- c("x", "y", "z", "ca", "crad")
  head(df.tree)
  write.csv(df.tree, paste0(label, "_pre_treeData.csv"), row.names = F)
}


# work --------------------------------------------------------------------
file <- "US-HA1_p_01.tif"
detect.tree(file, "US_HA1")
