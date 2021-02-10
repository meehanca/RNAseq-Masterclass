###################################
#Modified plotPCA
###################################
#Allows us to specify which PCs to plot and how many genes to take into account


plotPCAEx = function(object, PCx = 1, PCy = 2, cond="condition", ntop=500, labels = TRUE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(cond %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  cond.df <- as.data.frame(colData(object)[, cond, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(cond) > 1) {
    factor(apply( cond.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[cond]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PCa=pca$x[,PCx], PCb=pca$x[,PCy], group=group, cond.df, name=colnames(object))
  
  pc1 <- ggplot(data=d, aes_string(x="PCa", y="PCb", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC",PCx,": ",round(percentVar[PCx] * 100),"% variance")) +
    ylab(paste0("PC",PCy,": ",round(percentVar[PCy] * 100),"% variance")) +
    coord_fixed() #+
  #geom_text(show.legend = F)
  
  pc2 <- pc1 + geom_point()
  pc3 <- pc2  +     theme(legend.position="none")
  
  #  Finally add the labels, using ggrepel so that they dont write over each other or the points  
  if (labels)
  {
    library("ggrepel")
    pc3 + geom_label_repel(aes(label = name),
                           color = "gray20",
                           data = d,
                           force = 10)
  }
  else{
    pc3
  }
}

###################################
#Distance Matrix
###################################
 
#Overview of sample distances
plotDistanceMatrix <- function(vst, colname = name){
  names <- sampleTable[,colname]
  sampleDists <- dist(t(assay(vst)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- names
  colnames(sampleDistMatrix) <- names
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  #pdf("Figures/Distance_matrix.pdf", height = 9, width = 15,onefile=FALSE)
  pheatmap <- pheatmap(sampleDistMatrix,
                       clustering_distance_rows=sampleDists,
                       clustering_distance_cols=sampleDists,
                       col=colors, fontsize = 13)
}

###################################
#XY plots
###################################

plot_xy_scatter<-function(dds,x='one',y='two'){
  dds<-estimateSizeFactors(dds)
  rep_values<- counts(dds, normalized=TRUE)[,c(x,y)]
  # Take logs
  vals<-log2(rep_values[,c(x,y)] + 0.5)
  plot(vals,pch=16, cex=0.4,xlab=paste('sample',x),ylab=paste('sample',y))
  grid(col = "darkgray", lty = "solid",lwd = par("lwd"), equilogs = TRUE)
  abline(coef=c(0,1))
  return(vals)
}

###################################
#Filter differentially expressed genes
###################################

filter_degs <- function(res, padj = 0.05, logFC=0){
  res2 = res[!(is.na(res$padj)),]
  res2 = res2[res2$padj < padj,]
  res2 = res2[abs(res2$log2FoldChange) > logFC,]
  summary(res)
  return(res2)
}

###################################
#Heatmap Clusters
###################################

plot_tree <- function(pheatmap_Data){
  plot(pheatmap_Data$tree_row)
  max <- round(range(pheatmap_Data$tree_row$height))
  #change seq() to make different intersections on the tree
  for(i in seq(0.5,max[2],1)){
    abline(a = i, b = 0, col = "blue", lty = i+1 )
  }
}
produce_clusters <- function(pheatmap_Data, value){
  #grab gene names + cluster information from heatmap, melt + transform into dataframe
  
  melted.pheatmap_Data <- melt(cutree(pheatmap_Data$tree_row, h = value))
  melted.pheatmap_Data <- as.data.frame(melted.pheatmap_Data)
  #get counts for gene names + scale counts
  
  interesting_counts <- counts[rownames(melted.pheatmap_Data),]
  interesting_counts = t(scale(t(as.matrix(interesting_counts)),center = TRUE, scale = TRUE))
  interesting_counts <- as.data.frame(interesting_counts)
  #move columns of count table into same order produced by pheatmap clustering
  interesting_counts <- interesting_counts[,pheatmap_Data$tree_col$order]
  #assign cluster ids to count table
  interesting_counts$cluster <- melted.pheatmap_Data$value
  interesting_counts$ids <- rownames(interesting_counts)
  melt_counts <<- melt(interesting_counts, id.vars = c("cluster", "ids"))
  #plot that shit
  c <- ggplot(melt_counts,aes(alpha = 0.5))
  
  c <- c + geom_line(aes(x = variable, y = value,
                         group = ids)) + theme_bw() + facet_grid(cluster~., scales = "fixed") +
    scale_color_manual(values = c("lightblue","goldenrod1", "orange", "red")) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = -90,
                                     hjust = 0),
          axis.title = element_blank(),
          legend.position = "top",
          axis.ticks.y = element_blank()) 
  c
}
cluster_to_heatmap <- function( cluster, rownames = F, scale = "row"){
  specific_genes <<- unique((melt_counts[melt_counts$cluster==cluster,])[,2])
  pheatmap(log2((counts[specific_genes,])+1), scale = "row", main = paste("Cluster", cluster), show_rownames = rownames, border_color = NA)
}