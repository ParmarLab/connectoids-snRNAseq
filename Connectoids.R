#' Reconstruction of the human nigrostriatal pathway in vitro reveals target-dependent dopamine neuron maturation
#' @author Edoardo Sozzi
#' @date 2026-09-15
# Connectoid analysis

# Packages ---------------------------------------------------------------------
{
  library(tidyverse)
  library(Seurat)
  library(writexl)
  library(harmony)
  library(msigdbr)
  library(msigdb)
  library(pheatmap)
}

# Colors -----------------------------------------------------------------------
{
  Undefined.c <- c("#cccccc")
  
  #Group
  d60_vMB_vMB.c <- c("#e5acab")
  d60_vMB_STR.c <- c("#acdae2")
  d90_vMB_vMB.c <- c("#cc4949")
  d90_vMB_STR.c <- c("#49bacc")
  
  Group_timepoint.c <- c(d60_vMB_vMB.c,d90_vMB_vMB.c,d60_vMB_STR.c,d90_vMB_STR.c)
  
  #Celltype
  FP1.c <- c("#e3abcd")
  FP2.c <- c("#db80b3")
  FP3.c <- c("#c74e9b")
  DA_Neurons.c <- c("#49bacc")
  Cholinergic.c <- c("#5078ba")
  Celltype.c <- c(FP1.c,FP2.c,FP3.c,DA_Neurons.c,Cholinergic.c)
  
  #Gradient
  Max_gradient.c <- c("#5b2466")
  colfunc_white <- colorRampPalette(c("white",Max_gradient.c))
  Gradient_white.c<-colfunc_white(100)
}
#-------------------------------------------------------------------------------



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1 Preprocessing and Quality control -------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("//Users/edoardo/OneDrive - Lund University/_Research/PhD_ParmarLab/Parmar Lab/Projects/Connectoids_project/snRNAseq/Connectoids/Raw_data/")

#Loading Data
directories<-c("./Connectoid_THCre_VM_FlexGFP_H9_VM_VMside_d60_1/",
               "./Connectoid_THCre_VM_FlexGFP_H9_VM_VMside_d60_2/",
               "./Connectoid_THCre_VM_FlexGFP_H9_VM_VM_side_d90_1/filtered_feature_bc_matrix/",
               "./Connectoid_THCre_VM_FlexGFP_H9_STR_VMside_d60_1/",
               "./Connectoid_THCre_VM_FlexGFP_H9_STR_VMside_d60_2/",
               "./Connectoid_THCre_VM_FlexGFP_H9_STR_VM_side_d90_1/filtered_feature_bc_matrix/")

# Define orig.ident names
names(directories)<-c("Con-vMB-vMB-d60-1","Con-vMB-vMB-d60-2","Con-vMB-vMB-d90-1",
                      "Con-vMB-STR-d60-1","Con-vMB-STR-d60-2","Con-vMB-STR-d90-1")

# Load the dataset
Connectoid.data <- Read10X(data.dir = directories) 

# Create the Seurat object, normalize and run PCA
Connectoid <- CreateSeuratObject(Connectoid.data, 
                                 project = "ES002", 
                                 min.cells = 0, 
                                 min.features=200) %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = Connectoid.data@var.genes, npcs = 50, verbose = FALSE)

# Reorder idents
Connectoid$orig.ident <- factor(Connectoid$orig.ident, 
                                levels = c("Con-vMB-vMB-d60-1","Con-vMB-vMB-d60-2","Con-vMB-vMB-d90-1",
                                           "Con-vMB-STR-d60-1","Con-vMB-STR-d60-2","Con-vMB-STR-d90-1"))

# Quality control
Connectoid[["percent.mt"]] <- PercentageFeatureSet(Connectoid,pattern = "^.*?MT-")
Connectoid[["percent.rb"]] <- PercentageFeatureSet(Connectoid,pattern = "^.*RP[SL]")
Connectoid <- subset(Connectoid, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 1)

# Group
Connectoid$Group <- dplyr::case_match(Connectoid$orig.ident,
                                      "Con-vMB-vMB-d60-1" ~ "Con_vMB_vMB",
                                      "Con-vMB-vMB-d60-2" ~ "Con_vMB_vMB",
                                      "Con-vMB-vMB-d90-1" ~ "Con_vMB_vMB",
                                      "Con-vMB-STR-d60-1" ~ "Con_vMB_STR",
                                      "Con-vMB-STR-d60-2" ~ "Con_vMB_STR",
                                      "Con-vMB-STR-d90-1" ~ "Con_vMB_STR")

# Group_timepoint
Connectoid$Group_timepoint <- dplyr::case_match(Connectoid$orig.ident,
                                                "Con-vMB-vMB-d60-1" ~ "d60_vMB_vMB",
                                                "Con-vMB-vMB-d60-2" ~ "d60_vMB_vMB",
                                                "Con-vMB-vMB-d90-1" ~ "d90_vMB_vMB",
                                                "Con-vMB-STR-d60-1" ~ "d60_vMB_STR",
                                                "Con-vMB-STR-d60-2" ~ "d60_vMB_STR",
                                                "Con-vMB-STR-d90-1" ~ "d90_vMB_STR")

# Reorder Group
Connectoid$Group <- factor(Connectoid$Group, levels = c("Con_vMB_vMB","Con_vMB_STR"))

# Reorder Group_timepoint
Connectoid$Group_timepoint <- factor(Connectoid$Group_timepoint, levels = c("d60_vMB_vMB","d90_vMB_vMB","d60_vMB_STR","d90_vMB_STR"))
Connectoid$Group_timepoint_2 <- factor(Connectoid$Group_timepoint, levels = c("d60_vMB_vMB","d60_vMB_STR","d90_vMB_vMB","d90_vMB_STR"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Genetic demultiplexing based on demuxlet
Con_vMB_vMB_d60_1_genetic_demultiplexing <- read.csv("./Connectoid_THCre_VM_FlexGFP_H9_VM_VMside_d60_1.best", sep="\t")
Con_vMB_vMB_d60_2_genetic_demultiplexing <- read.csv("./Connectoid_THCre_VM_FlexGFP_H9_VM_VMside_d60_2.best", sep="\t")
Con_vMB_vMB_d90_1_genetic_demultiplexing <- read.csv("./Connectoid_THCre_VM_FlexGFP_H9_VM_VM_side_d90_1.best", sep="\t")
Con_vMB_STR_d60_1_genetic_demultiplexing <- read.csv("./Connectoid_THCre_VM_FlexGFP_H9_STR_VMside_d60_1.best", sep="\t")
Con_vMB_STR_d60_2_genetic_demultiplexing <- read.csv("./Connectoid_THCre_VM_FlexGFP_H9_STR_VMside_d60_2.best", sep="\t")
Con_vMB_STR_d90_1_genetic_demultiplexing <- read.csv("./Connectoid_THCre_VM_FlexGFP_H9_STR_VM_side_d90_1.best", sep="\t")

Con_vMB_vMB_d60_1_genetic_demultiplexing$sample <- "Con-vMB-vMB-d60-1"
Con_vMB_vMB_d60_2_genetic_demultiplexing$sample <- "Con-vMB-vMB-d60-2"
Con_vMB_vMB_d90_1_genetic_demultiplexing$sample <- "Con-vMB-vMB-d90-1"
Con_vMB_STR_d60_1_genetic_demultiplexing$sample <- "Con-vMB-STR-d60-1"
Con_vMB_STR_d60_2_genetic_demultiplexing$sample <- "Con-vMB-STR-d60-2"
Con_vMB_STR_d90_1_genetic_demultiplexing$sample <- "Con-vMB-STR-d90-1"

# Combine all of the genetic demultiplexing data
Genetic_demultiplexing <- rbind(Con_vMB_vMB_d60_1_genetic_demultiplexing,
                                Con_vMB_vMB_d60_2_genetic_demultiplexing,
                                Con_vMB_vMB_d90_1_genetic_demultiplexing,
                                Con_vMB_STR_d60_1_genetic_demultiplexing,
                                Con_vMB_STR_d60_2_genetic_demultiplexing,
                                Con_vMB_STR_d90_1_genetic_demultiplexing)

# HS1001 and H9 were mislabeled, here is correct
Genetic_demultiplexing$Genome <- ifelse(Genetic_demultiplexing$SNG.BEST.GUESS=="WH-3759-HS1001", "H9",
                                        ifelse(Genetic_demultiplexing$SNG.BEST.GUESS=="WH-3759-RC17", "RC17",
                                               ifelse(Genetic_demultiplexing$SNG.BEST.GUESS=="WH-3759-H9", "HS1001","Undefined")))

Genetic_demultiplexing$Genome_next_guess <- ifelse(Genetic_demultiplexing$SNG.NEXT.GUESS=="WH-3759-HS1001", "H9",
                                                   ifelse(Genetic_demultiplexing$SNG.NEXT.GUESS=="WH-3759-RC17", "RC17",
                                                          ifelse(Genetic_demultiplexing$SNG.NEXT.GUESS=="WH-3759-H9", "HS1001","Undefined")))

# Combine Sample and cell barcode
Genetic_demultiplexing <- Genetic_demultiplexing %>% unite(cellbarcode, c(sample, BARCODE), sep = "_", remove = FALSE, na.rm = T)
Genetic_demultiplexing_filtered <- Genetic_demultiplexing[,c(2,23,24)]
rownames(Genetic_demultiplexing_filtered) <- Genetic_demultiplexing_filtered$cellbarcode
Connectoid <- AddMetaData(Connectoid,Genetic_demultiplexing_filtered) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2 Clustering and annotation ---------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vMB_side <- subset(Connectoid, subset=Genome==c("RC17"))

# Run Harmony integration 
vMB_side.harmony <- RunHarmony(vMB_side, group.by.vars = "orig.ident",lambda = 1)

# UMAP embeddings - selected resolution = 0.14
vMB_side.harmony <- vMB_side.harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.14) %>% 
  identity()

# Annotation
Idents(vMB_side.harmony) <- "seurat_clusters"
vMB_side.harmony <- RenameIdents(vMB_side.harmony, 
                                 `0` = "FP_Late",
                                 `1` = "DA_neurons",
                                 `2` = "DA_neurons",
                                 `3` = "FP_Early",
                                 `4` = "Chol_neurons",
                                 `5` = "FP_Cycling")

vMB_side.harmony$Celltype <- Idents(vMB_side.harmony)
vMB_side.harmony$Celltype <- factor(vMB_side.harmony$Celltype, 
                                    levels = c("FP_Cycling",
                                               "FP_Early",
                                               "FP_Late",
                                               "DA_neurons",
                                               "Chol_neurons"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots
# Selected Markers 
Markers.M <- c("SOX2","ZEB2","TCF7L2","GLI3","CENPK","ETV5","NTN1","CORIN","DCN",
               "SHH","FOXA2","LMX1B","VIM","SLC2A3","SLC3A2","DLK1","SLC5A3","JUN",
               "NRXN1","MYT1L","SYN3","DCX","SNAP25","GAP43","NLGN1","PBX3","MEIS2",
               "EBF3","SLC5A7","ISL1","EPHA6","CHAT","MET","RET")

# Dotplot - Fig. EV4E
DotPlot(vMB_side.harmony, features = Markers.M, group.by = "Celltype") + 
  scale_colour_gradient(low = "white", high = Max_gradient.c) + 
  scale_x_discrete(limits=rev) + coord_flip() + theme_classic()
ggsave("DotPlot_vMB_side_Celltype.pdf", width = 1100, height = 2000, units="px")

# Feature plots - Fig. EV4F & Fig. 4F
Feature_plot_markers <- c("TH","LMX1B","MYT1L","FOXA2",
                          "EN1","KCNJ6","NES","SLC5A7",
                          "NR4A1","GRIA4","HES1","ISL1")
for (i in Feature_plot_markers){
  FeaturePlot(vMB_side.harmony, features = i, order=T, cols = c(Undefined.c,Max_gradient.c), raster = F) 
  ggsave(paste0("Feature_vMB_side_",i,".pdf"), width = 900, height=750, units = "px")
}

# UMAP scatter plot by Celltype - Fig. 4E
DimPlot(vMB_side.harmony, label = F, group.by = "Celltype",cols = Celltype.c, shuffle = T, raster=F) + NoLegend()
ggsave("UMAP_vMB_side_celltypes.pdf", width = 1500, height = 1500, units="px")

# Proportions - Fig. 4E
n_cells <- FetchData(vMB_side.harmony, vars = c("Group", "Celltype")) %>%
  dplyr::count(Group,Celltype) %>%
  tidyr::spread(Celltype, n)
n_cells_melt <- reshape2::melt(n_cells)
Group_cluster_plot <- ggplot(n_cells_melt, aes(x=Group, y=value, fill=variable)) + 
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=c(Celltype.c)) + NoLegend()
ggsave("Barchart_vMB_side_Celltype_by_Group.pdf", width=500, height=1000, units = "px")

# UMAP scatter plot by Group - Fig. EV4G
DimPlot(vMB_side.harmony, label = F, group.by = "Group_timepoint",cols = Group_timepoint.c, shuffle=T, raster=F) + NoLegend()
ggsave("UMAP_vMB_side_Group_timepoint.pdf", width = 1500, height = 1500, units="px")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3 DA neuron analysis ----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Subset DA neurons
DA_neurons <- subset(vMB_side.harmony, subset=Celltype==c("DA_neurons"))
expr <- GetAssayData(DA_neurons, layer = "counts")

# Define GeneSet
geneSet <- "GOBP_MIDBRAIN_DEVELOPMENT" 
Genes.list <- escape::getGeneSets(library = "C5",gene.sets =geneSet)
escape_mat <- escape::escape.matrix(expr,gene.sets = Genes.list,
                            method = "ssGSEA",
                            groups = 1000,
                            min.size = 5)

DA_neurons <- Seurat::AddMetaData(DA_neurons,escape_mat)

# GOBP_MIDBRAIN_DEVELOPMENT
# GOBP_NEURON_DEVELOPMENT
# GOBP_NEUROGENESIS
# GOBP_DEVELOPMENTAL_MATURATION
# GOBP_NEURON_MATURATION
# GOBP_NEUROTRANSMITTER_SECRETION
# GOBP_RESPONSE_TO_DOPAMINE

# Heatmap - Fig. 4G
metadata <- DA_neurons@meta.data
Heatmap_GSEA_Group <- metadata %>% group_by(Group_timepoint_2) %>% 
  summarise(GOBP_MIDBRAIN_DEVELOPMENT=mean(`GOBP-MIDBRAIN-DEVELOPMENT`), 
            GOBP_NEURON_DEVELOPMENT=mean(`GOBP-NEURON-DEVELOPMENT`), 
            GOBP_NEUROGENESIS=mean(`GOBP-NEUROGENESIS`),   
            GOBP_DEVELOPMENTAL_MATURATION=mean(`GOBP-DEVELOPMENTAL-MATURATION`), 
            GOBP_NEURON_MATURATION=mean(`GOBP-NEURON-MATURATION`), 
            GOBP_NEUROTRANSMITTER_SECRETION=mean(`GOBP-NEUROTRANSMITTER-SECRETION`),
            GOBP_RESPONSE_TO_DOPAMINE=mean(`GOBP-RESPONSE-TO-DOPAMINE`))
Heatmap_GSEA_Group <- t(Heatmap_GSEA_Group[,-1])
Heatmap_maturation_by_Group <- pheatmap(Heatmap_GSEA_Group,scale = "row", color = Gradient_white.c, cluster_rows = F, cluster_cols = F)

# Save RDS file ----------------------------------------------------------------
saveRDS(object = DA_neurons, "/Users/edoardo/OneDrive - Lund University/_Research/PhD_ParmarLab/Parmar Lab/Projects/Connectoids_project/snRNAseq/Connectoids/Processed_data/DA_neurons.rds")     
#-------------------------------------------------------------------------------

#sessionInfo()
#R version 4.6.1 (2026-06-24)
#Platform: aarch64-apple-darwin23
#Running under: macOS Tahoe 26.6.2

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#time zone: Europe/Zurich
#tzcode source: internal

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] pheatmap_1.0.13      GSEABase_1.74.0      graph_1.90.0         annotate_1.90.0      XML_3.99-0.23        AnnotationDbi_1.74.0
#[7] IRanges_2.46.0       S4Vectors_0.50.1     Biobase_2.72.0       BiocGenerics_0.58.1  generics_0.1.4       future_1.75.0       
#[13] msigdb_1.20.0        msigdbr_26.1.1       harmony_2.0.5        Rcpp_1.1.2           writexl_2.0.0        Seurat_5.5.1        
#[19] SeuratObject_5.4.0   sp_2.2-3             lubridate_1.9.5      forcats_1.0.1        stringr_1.6.0        dplyr_1.2.1         
#[25] purrr_1.2.2          readr_2.2.0          tidyr_1.3.2          tibble_3.3.1         ggplot2_4.0.3        tidyverse_2.0.0     