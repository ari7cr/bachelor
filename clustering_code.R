# iCCL clustering loop (ehemals multiple dimensions and resolutions)

#funktion 1: clustering loop über 2 for schleifen.
##Die erste für Dimensions von x bis y und die zweite für resolutions 0.25, 0.5, 0.75 und 1


iCCL <- function(SeuratObject, min.dim, max.dim){
  wd <- getwd()
  dir.create("iCCL_results")
  dim.list = min.dim:max.dim
  resol = c(0.25, 0.50, 0.75, 1)
  
  for(i in dim.list){
    dimens <- Seurat::FindNeighbors(SeuratObject, dims = 1:i)
    for(r in resol){
      clust <- Seurat::FindClusters(dimens, resolution = r)
      clust <- Seurat::RunUMAP(clust, dims = 1:i)
      nam <- paste("Dim_", i, "_Res_", r, sep = "")
      print(nam)
      assign(nam, clust)
      pl <- Seurat::DimPlot(clust, label = 1, repel = 1) + ggplot2::ggtitle(nam) # +NoLegend()
      mypath <- file.path(wd, projectname, paste(nam, ".png", sep = ""))
      png(file=mypath)
      print(pl)
      dev.off()
    }
  }
  
  
}

# Funktion 2: Dimensionalität (ElbowPlot) abschätzen. Ich habe zeitlich leider keine Clustering Scores einbauen können aber
# aufgrund deines inputs, dass das ablesen im elbowplot immer ein bisschen geschätzt ist, habe ich die Lösung hier gefunden 
# die die Varianz über das Seurat Objekt ermittelt und einen Punkt empfiehlt, den man dann mit obrigem Command kontrollieren kann
# Die Methode ist aus einem Harvard scRNAseq Kurs https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html

predictdimension <- function(SeuratObject){
  
  
  # Determine percent of variation associated with each PC
  pct <- SeuratObject[["pca"]]@stdev / sum(SeuratObject[["pca"]]@stdev) * 100
  
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  #co1
  
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  # last point where change of % of variation is more than 0.1%.
  #co2
  
  # proceed with Minimum of the two calculation
  pcs <- min(co1, co2)
  #pcs
  
  # Create a dataframe with values
  plot_df <- data.frame(pct = pct,
                        cumu = cumu,
                        rank = 1:length(pct))
  
  # Elbow plot to visualize
  prediction <-   ggplot2::ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
    geom_text() +
    geom_vline(xintercept = 90, color = "grey") +
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw() + ggplot2::ggtitle("prediction: ", pcs)
  plot(prediction)
  
  
}



