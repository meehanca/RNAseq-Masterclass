###################################
#Load libraries
###################################
#Check whether we have installed BioCManager, and install if not
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Vector of packages we will use today
packages <- c("DESeq2","pheatmap","reshape2",
              "VennDiagram","ggplot2","ggrepel",
              "grid","gridExtra","data.table",
              "RColorBrewer","biomaRt","tidyverse",
              "EnhancedVolcano")

#You can individually install packages you may not have with this command
BiocManager::install("DESeq2")

#Alternatively, this command will check whether you have previously installed these packages and if not install in one go
BiocManager::install(setdiff(packages, rownames(installed.packages())))

#Load packages all in one go
lapply(packages, require, character.only=TRUE)
              

###################################
# set up directories 
###################################

#Directory containing counts
htseqDir <- ("./ALL_COUNTS")

#Load file containing some useful functions
source(file='Useful_functions.R')

#load our sample_table file with metadata
sampleTable <- read.table("sample_table.txt",header=TRUE) 

#have a look at the metadata
View(sampleTable) 

#Load counts
dds<-DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = htseqDir,design = ~condition) 

#Initialise DESeq 
dds<- DESeq(dds) 

###################################
# PCA and transformation
###################################

# vst transform for PCA 
vst <- vst(dds) 


#plot (and write) PCA
dir.create("Plots")
#pdf("./Plots/PCA.pdf", height = 5, width = 9,onefile=FALSE)
o <- plotPCAEx(vst, cond ="condition",PCx=1,PCy=2, label = T)
o +      theme(legend.position="right") + theme(text = element_text(size = 12))
#dev.off()

#pdf("Plots/Distance_matrix.pdf", height = 9, width = 15,onefile=FALSE)
plotDistanceMatrix(vst, colname = "sample")
#dev.off()


###################################
# scatter
###################################

#change plot window to 1 row + 3 columns of plots
par(mfrow=c(1,3)) 

#Plot scatterplots of gene expression
plot_xy_scatter(dds, 'col_0_1','col_0_2')
plot_xy_scatter(dds, 'col_0_1','col_0_3')
plot_xy_scatter(dds, 'col_0_2','col_0_3')

#See if you can work out how to plot all replicate scatterplots without duplicates in one plot
#(Maybe in your own time)

# change plot window back to default
par(mfrow=c(1,1))

###################################
# DEG
###################################
#Set comparisons
level <- "condition"
control <- "col_0" #col_O is our wild-type condition
condition <- "met1" #met1 is one of our mutant lines

#Calculates differential expression 
res = results(dds, contrast=c(level,condition,control),alpha=0.1) 
summary(res)
head(res)

#Filters differentially expressed genes
met1_DEGs <- res
met1_filtered_DEGs <- filter_degs(res, padj = 0.1, logFC = 0)

#We should produce sets of differentially expressed genes for the other comparisons too


#See if you can produce separate objects containing the upregulated or downregulated genes

###################################
# DEG visualisation
###################################

#Volcano plot
#type  ?EnhancedVolcano  into your console to view how you can change your output
EnhancedVolcano(met1_DEGs,
                lab = rownames(met1_DEGs),
                pLabellingCutoff = NULL,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-5, 20),
                ylim = c(0,10),
                pCutoff = 0.01)

#MA plot
DESeq2::plotMA(met1_DEGs, ylim=c(-10,15))


#Heatmap
counts <- counts(dds)
met1_deg_counts <- counts[rownames(met1_filtered_DEGs),]
pheatmap(log10(met1_deg_counts+1), show_rownames = F,
         scale = "row", border_color = NA,
         fontsize_col = 19)

#Which plot (or plots) do you think best describe the data?
#Try to produce and save plots for each comparison

###################################
# mart
###################################
#Connect to the plants biomart database
ensembl = useEnsembl(biomart="plants_mart",host= "plants.ensembl.org")

#List organisms available
listDatasets(ensembl)

#Load the arabidopsis dataset
ensembl = useEnsembl(biomart="plants_mart", dataset="athaliana_eg_gene", host= "plants.ensembl.org")

#Have a quick at what information is available
head((listAttributes(ensembl)),30)

#Lets grab the gene ids and the gene descriptions
genedesc <- getBM(attributes=c('ensembl_gene_id','description'), filters ='ensembl_gene_id'
                  ,values =rownames(met1_filtered_DEGs) , mart =ensembl)

#We should probably join this information up with our DEG statistics
met1_filtered_DEGs$ensembl_gene_id <- rownames(met1_filtered_DEGs)
met1_DEGS_descriptions_stats <- met1_filtered_DEGs[,c(7,1:6)]
met1_DEGS_descriptions_stats <- left_join(as.data.frame(met1_DEGS_descriptions_stats), genedesc)
View(met1_DEGS_descriptions_stats)

#Lets write these met1 statistics
dir.create("GeneTables")
write.table(met1_DEGS_descriptions_stats, file = "GeneTables/met1_DEGS_descriptions_stats.txt",quote = F, row.names = F, col.names = T, sep = '\t')

#We'll have to produce tables for each comparison


###################################
# Single expression plot
###################################

#Boring plot for plotting expression of one single gene
d <- plotCounts(dds, 
                gene="AT5G13930", 
                intgroup="condition",
                returnData=TRUE)
x = ggplot(d, aes(x=condition, y=log2(count), size = 10, group = condition)) +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.line = element_line(colour = "black")) +
  geom_point() +
  #geom_line(size = 1.5, linetype = "dotted") +
  ylab("log2(Counts)") +
  ggtitle("AT5G13930") +
  geom_label_repel(aes(label = row.names(d))) +
  guides(size=FALSE) +
  theme(legend.text = element_text( size = 14),axis.text.x = element_text(size = 12)
        ,legend.title = element_text( size = 16),axis.title=element_text(size=14,face="bold")) 
print(x)



###################################
# Heatmap Clusters
###################################
#Bit of extra stuff down here if you want
#You can explore the hierarchical clustering of differentially expressed genes calculated during the pheatmap function
pheatmap_Data <- pheatmap(log10(met1_deg_counts+1), show_rownames = F,
                          scale = "row", border_color = NA,
                          fontsize_col = 19)
plot_tree(pheatmap_Data)
produce_clusters(pheatmap_Data, 6)
pheatmap_Data <- cluster_to_heatmap(3)

