#' Reconstruction of the human nigrostriatal pathway in vitro reveals target-dependent dopamine neuron maturation
#' @author Edoardo Sozzi
#' @date 2026-09-15
# Assembloid-Connectoid analysis

# Packages ---------------------------------------------------------------------
{
  library(tidyverse)
  library(Seurat)
  library(writexl)
  library(harmony)
  library(EnhancedVolcano)
  library(org.Hs.eg.db)
  library(pheatmap)
}

# Colors ------------------------------------------------------------------------
{
  Undefined.c <- c("#cccccc")
  
  #Group
  Assembloid_d60.c <- c("#aed9b3")
  Connectoid_d60.c <- c("#acdae2")
  Assembloid_d90.c <- c("#5ebb60")
  Connectoid_d90.c <- c("#49bacc")
  
  Group.c <- c(Assembloid_d90.c,Connectoid_d90.c)
  Group_2.c <- c(Assembloid_d60.c,Assembloid_d90.c,Connectoid_d60.c,Connectoid_d90.c)
  
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
setwd("//Users/edoardo/OneDrive - Lund University/_Research/PhD_ParmarLab/Parmar Lab/Projects/Connectoids_project/snRNAseq/")

# Loading Data
directories<-c("./Assembloids/Raw_data/Assembloid_THCre_VM_FlexGFP_H9_STR_d60_1/",
               "./Assembloids/Raw_data/Assembloid_THCre_VM_FlexGFP_H9_STR_d60_2/",
               "./Assembloids/Raw_data/Assembloid_THCre_VM_FlexGFP_H9_STR_d90_1/filtered_feature_bc_matrix/",
               "./Connectoids/Raw_data/Connectoid_THCre_VM_FlexGFP_H9_STR_VMside_d60_1/",
               "./Connectoids/Raw_data/Connectoid_THCre_VM_FlexGFP_H9_STR_VMside_d60_2/",
               "./Connectoids/Raw_data/Connectoid_THCre_VM_FlexGFP_H9_STR_VM_side_d90_1/filtered_feature_bc_matrix/")

# Define orig.ident names
names(directories)<-c("Ass-vMB-STR-d60-1","Ass-vMB-STR-d60-2","Ass-vMB-STR-d90-1",
                      "Con-vMB-STR-d60-1","Con-vMB-STR-d60-2","Con-vMB-STR-d90-1")

# Load the dataset
vMB_STR.data <- Read10X(data.dir = directories) 

# Create the Seurat object, normalize and run PCA
vMB_STR <- CreateSeuratObject(vMB_STR.data, 
                              project = "ES002", 
                              min.cells = 0, 
                              min.features=200) %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = vMB_STR.data@var.genes, npcs = 50, verbose = FALSE)

# Reorder Orig.ident
vMB_STR$orig.ident <- factor(vMB_STR$orig.ident, 
                             levels = c("Ass-vMB-STR-d60-1","Ass-vMB-STR-d60-2","Ass-vMB-STR-d90-1",
                                        "Con-vMB-STR-d60-1","Con-vMB-STR-d60-2","Con-vMB-STR-d90-1"))

# Quality control
vMB_STR[["percent.mt"]] <- PercentageFeatureSet(vMB_STR,pattern = "^.*?MT-")
vMB_STR[["percent.rb"]] <- PercentageFeatureSet(vMB_STR,pattern = "^.*RP[SL]")
vMB_STR <- subset(vMB_STR, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 1)

# Group
vMB_STR$Group <- dplyr::case_match(vMB_STR$orig.ident,
                                   "Ass-vMB-STR-d60-1" ~ "Assembloid_d60",
                                   "Ass-vMB-STR-d60-2" ~ "Assembloid_d60",
                                   "Ass-vMB-STR-d90-1" ~ "Assembloid_d90",
                                   "Con-vMB-STR-d60-1" ~ "Connectoid_d60",
                                   "Con-vMB-STR-d60-2" ~ "Connectoid_d60",
                                   "Con-vMB-STR-d90-1" ~ "Connectoid_d90")

# Group_type
vMB_STR$Group_type <- dplyr::case_match(vMB_STR$orig.ident,
                                        "Ass-vMB-STR-d60-1" ~ "Assembloid",
                                        "Ass-vMB-STR-d60-2" ~ "Assembloid",
                                        "Ass-vMB-STR-d90-1" ~ "Assembloid",
                                        "Con-vMB-STR-d60-1" ~ "Connectoid",
                                        "Con-vMB-STR-d60-2" ~ "Connectoid",
                                        "Con-vMB-STR-d90-1" ~ "Connectoid")

# Reorder Group
vMB_STR$Group <- factor(vMB_STR$Group, levels = c("Assembloid_d60","Assembloid_d90","Connectoid_d60","Connectoid_d90"))

# Reorder Group_type
vMB_STR$Group_type <- factor(vMB_STR$Group_type, levels = c("Assembloid","Connectoid"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Genetic demultiplexing based on demuxlet
Ass_vMB_STR_d60_1_genetic_demultiplexing <- read.csv("./Assembloids/Raw_data/Assembloid_THCre_VM_FlexGFP_H9_STR_d60_1.best", sep="\t")
Ass_vMB_STR_d60_2_genetic_demultiplexing <- read.csv("./Assembloids/Raw_data/Assembloid_THCre_VM_FlexGFP_H9_STR_d60_2.best", sep="\t")
Ass_vMB_STR_d90_1_genetic_demultiplexing <- read.csv("./Assembloids/Raw_data/Assembloid_THCre_VM_FlexGFP_H9_STR_d90_1.best", sep="\t")
Con_vMB_STR_d60_1_genetic_demultiplexing <- read.csv("./Connectoids/Raw_data/Connectoid_THCre_VM_FlexGFP_H9_STR_VMside_d60_1.best", sep="\t")
Con_vMB_STR_d60_2_genetic_demultiplexing <- read.csv("./Connectoids/Raw_data/Connectoid_THCre_VM_FlexGFP_H9_STR_VMside_d60_2.best", sep="\t")
Con_vMB_STR_d90_1_genetic_demultiplexing <- read.csv("./Connectoids/Raw_data/Connectoid_THCre_VM_FlexGFP_H9_STR_VM_side_d90_1.best", sep="\t")

# Add "sample" column name, should match "orig.ident" name
Ass_vMB_STR_d60_1_genetic_demultiplexing$sample <- "Ass-vMB-STR-d60-1"
Ass_vMB_STR_d60_2_genetic_demultiplexing$sample <- "Ass-vMB-STR-d60-2"
Ass_vMB_STR_d90_1_genetic_demultiplexing$sample <- "Ass-vMB-STR-d90-1"
Con_vMB_STR_d60_1_genetic_demultiplexing$sample <- "Con-vMB-STR-d60-1"
Con_vMB_STR_d60_2_genetic_demultiplexing$sample <- "Con-vMB-STR-d60-2"
Con_vMB_STR_d90_1_genetic_demultiplexing$sample <- "Con-vMB-STR-d90-1"

# Combine all of the genetic demultiplexing data
Genetic_demultiplexing <- rbind(Ass_vMB_STR_d60_1_genetic_demultiplexing,
                                Ass_vMB_STR_d60_2_genetic_demultiplexing,
                                Ass_vMB_STR_d90_1_genetic_demultiplexing,
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
vMB_STR <- AddMetaData(vMB_STR,Genetic_demultiplexing_filtered) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2 Clustering and annotation ---------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vMB_side <- subset(vMB_STR, subset=Genome==c("RC17"))

# Run Harmony integration 
vMB_side.harmony <- RunHarmony(vMB_side, group.by.vars = "orig.ident",lambda = 1)

# UMAP embeddings - selected resolution = 0.12
vMB_side.harmony <- vMB_side.harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.12) %>% 
  identity()

# Annotation
Idents(vMB_side.harmony) <- "seurat_clusters"
vMB_side.harmony <- RenameIdents(vMB_side.harmony, 
                                 `0` = "DA_neurons",
                                 `1` = "FP_Late",
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

# UMAP scatter plot by Celltype - Fig. 4J
DimPlot(vMB_side.harmony, label = F, group.by = "Celltype",cols = Celltype.c, shuffle = T, raster=F) + NoLegend()
ggsave("UMAP_vMB_side_celltypes.pdf", width = 1500, height = 1500, units="px")
n_cells <- FetchData(vMB_side.harmony, vars = c("Group_type", "Celltype")) %>%
  dplyr::count(Group_type,Celltype) %>%
  tidyr::spread(Celltype, n)
n_cells_melt <- reshape2::melt(n_cells)
Group_cluster_plot <- ggplot(n_cells_melt, aes(x=Group_type, y=value, fill=variable)) + 
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=c(Celltype.c)) + NoLegend()
ggsave("Barchart_vMB_side_Celltype_by_Group.pdf", width=700, height=1000, units = "px")

# UMAP scatter plot by Group - Fig. EV4I
DimPlot(vMB_side.harmony, label = F, group.by = "Group",cols = Group_2.c, shuffle=F, raster=F) + NoLegend()
ggsave("UMAP_vMB_side_Group.pdf", width = 1500, height = 1500, units="px")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3 DA neuron analysis ----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DA_neurons <- subset(vMB_side.harmony, subset=Celltype==c("DA_neurons"))

# Volcano Plot - Fig. EV4J
labels <- c("TUBA1A","BEX2","SNAP25","IFT20","VIM","KCNA4","FOS","VEGFA","NFIA",
            "HES4","BCL11A","RERGL","PLCZ1","ZIC1")
Idents(DA_neurons) <- "Group_type"
Differences_DA_neurons_Group <- FindMarkers(DA_neurons,ident.1 = "Assembloid", ident.2 = "Connectoid")
Volcano_DA_neurons_by_Group.c <- ifelse(Differences_DA_neurons_Group$avg_log2FC < 0, Connectoid_d90.c, 
                                        ifelse(Differences_DA_neurons_Group$avg_log2FC > 0, Assembloid_d90.c,
                                               Undefined.c))
Volcano_DA_neurons_by_Group.c[is.na(Volcano_DA_neurons_by_Group.c)] <- Undefined.c
names(Volcano_DA_neurons_by_Group.c)[Volcano_DA_neurons_by_Group.c == Connectoid_d90.c] <- 'Connectoid'
names(Volcano_DA_neurons_by_Group.c)[Volcano_DA_neurons_by_Group.c == Undefined.c] <- "not significant"
names(Volcano_DA_neurons_by_Group.c)[Volcano_DA_neurons_by_Group.c == Assembloid_d90.c] <- 'Assembloid'

labels_Assembloid <- Differences_DA_neurons_Group %>% filter(avg_log2FC > 1)
labels_Assembloid <- rownames(labels_Assembloid)
labels_Connectoid <- Differences_DA_neurons_Group %>% filter(avg_log2FC < -1)
labels_Connectoid <- rownames(labels_Connectoid)
labels_unspecified <- Differences_DA_neurons_Group %>% filter(avg_log2FC > -1 & avg_log2FC < 1)
labels_unspecified <- rownames(labels_unspecified)

EnhancedVolcano(Differences_DA_neurons_Group,
                lab = rownames(Differences_DA_neurons_Group),
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = labels,
                drawConnectors = TRUE,
                colCustom = Volcano_DA_neurons_by_Group.c)
ggsave("Volcano_DA_neurons_Assembloids_vs_Connectoid.pdf",width=3000, height = 3000, units = "px")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GO - Fig. EV4K
# Convert symbol 
Differences_DA_neurons_Group$Entrez <- mapIds(org.Hs.eg.db,
                                              keys=rownames(Differences_DA_neurons_Group),
                                              column="ENTREZID",
                                              keytype="SYMBOL",
                                              multiVals="first")

# Defining up and down regulated genes
Up_genes <- Differences_DA_neurons_Group[Differences_DA_neurons_Group$avg_log2FC > 1 & Differences_DA_neurons_Group$p_val_adj < 0.05, 6] 
Dn_genes <- Differences_DA_neurons_Group[Differences_DA_neurons_Group$avg_log2FC < -1 & Differences_DA_neurons_Group$p_val_adj < 0.05, 6]
Backgroud_genes <- Differences_DA_neurons_Group[,6]

#Gene ontology
egobp_up <- clusterProfiler::enrichGO(
  gene     = Up_genes,
  universe = Backgroud_genes,
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, 
  maxGSSize = 500,
  readable = TRUE)
barplot(egobp_up,showCategory = 10,order = TRUE) +
  scale_fill_gradient(low = Max_gradient.c,high = Undefined.c,name = "Adjusted p-value") 
ggsave("GO_Barplot_Up_Connectoids_vs_Assembloids_LogFC_1_padj_0.05.pdf",width = 2000, height = 2800, units = "px")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig. 4K
expr <- GetAssayData(DA_neurons, layer = "counts")
DA_neurons$Group_2 <- factor(DA_neurons$Group, levels = c("Assembloid_d60","Connectoid_d60","Assembloid_d90","Connectoid_d90"))

# Define GeneSet
geneSet <- "GOBP_DOPAMINE_RECEPTOR_SIGNALING_PATHWAY" 
Genes.list <- escape::getGeneSets(library = "C5",gene.sets =geneSet)
escape_mat <- escape::escape.matrix(expr,gene.sets = Genes.list,
                                    method = "ssGSEA",
                                    groups = 1000,
                                    min.size = 5)

DA_neurons <- Seurat::AddMetaData(DA_neurons,escape_mat)

# GOCC_SYNAPSE
# GOBP_AXON_DEVELOPMENT
# GOBP_DOPAMINE_RECEPTOR_SIGNALING_PATHWAY
# GOBP_NEURON_PROJECTION_EXTENSION
# GOBP_SYNAPSE_ORGANIZATION
# GOBP_SYNAPTIC_SIGNALING
# GOBP_POSTSYNAPSE_ORGANIZATION

# Heatmap - Fig. 4K
metadata <- DA_neurons@meta.data
Heatmap_GSEA_Group <- metadata %>% group_by(Group_2) %>% 
  summarise(GOCC_SYNAPSE=mean(`GOCC-SYNAPSE`),
            GOBP_AXON_DEVELOPMENT=mean(`GOBP-AXON-DEVELOPMENT`),
            GOBP_DOPAMINE_RECEPTOR_SIGNALING_PATHWAY=mean(`GOBP-DOPAMINE-RECEPTOR-SIGNALING-PATHWAY`),
            GOBP_NEURON_PROJECTION_EXTENSION=mean(`GOBP-NEURON-PROJECTION-EXTENSION`),
            GOBP_SYNAPSE_ORGANIZATION=mean(`GOBP-SYNAPSE-ORGANIZATION`),
            GOBP_POSTSYNAPSE_ORGANIZATION=mean(`GOBP-POSTSYNAPSE-ORGANIZATION`),
            GOBP_SYNAPTIC_SIGNALING=mean(`GOBP-SYNAPTIC-SIGNALING`),
            GOBP_PRESYNAPSE_ORGANIZATION=mean(`GOBP-PRESYNAPSE-ORGANIZATION`))
Heatmap_GSEA_Group <- t(Heatmap_GSEA_Group[,-1])
Heatmap_maturation_by_Group <- pheatmap(Heatmap_GSEA_Group,scale = "row", color = Gradient_white.c, cluster_rows = F, cluster_cols = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DotPlot - Fig. 4L
DA_neurons$Group_2 <- factor(DA_neurons$Group, levels = c("Assembloid_d60","Connectoid_d60","Assembloid_d90","Connectoid_d90"))

Markers.M <- c("MAP1B","GRIA2","KIF1A","SNAP25","STMN2","ITGB1","TUBB3","TUBB2A",
               "SCN2A","SEMA3A","GRID2","EPHA5",
               "NCAM1","NPAS3","DCX","ZBTB20","SOX5","NTN1","ZEB2","GLI2")
DotPlot(DA_neurons, features = Markers.M, group.by = "Group") + 
  scale_colour_gradient(low = "white", high = Max_gradient.c) + 
  scale_x_discrete(limits=rev) + coord_flip() + theme_classic()
ggsave("DotPlot_DA_neurons_Group.pdf", width = 1290, height = 1500, units="px")

# Save RDS file ----------------------------------------------------------------
saveRDS(object = DA_neurons, "/Users/edoardo/OneDrive - Lund University/_Research/PhD_ParmarLab/Parmar Lab/Projects/Connectoids_project/snRNAseq/Assembloids_Connectoids_combined_vMB_STR/Processed_data/DA_neurons.rds")     
#-------------------------------------------------------------------------------

sessionInfo()
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
#  [1] org.Hs.eg.db_3.23.1    EnhancedVolcano_1.30.0 ggrepel_0.9.8          pheatmap_1.0.13       
#[5] GSEABase_1.74.0        graph_1.90.0           annotate_1.90.0        XML_3.99-0.23         
#[9] AnnotationDbi_1.74.0   IRanges_2.46.0         S4Vectors_0.50.1       Biobase_2.72.0        
#[13] BiocGenerics_0.58.1    generics_0.1.4         future_1.75.0          msigdb_1.20.0         
#[17] msigdbr_26.1.1         harmony_2.0.5          Rcpp_1.1.2             writexl_2.0.0         
#[21] Seurat_5.5.1           SeuratObject_5.4.0     sp_2.2-3               lubridate_1.9.5       
#[25] forcats_1.0.1          stringr_1.6.0          dplyr_1.2.1            purrr_1.2.2           
#[29] readr_2.2.0            tidyr_1.3.2            tibble_3.3.1           ggplot2_4.0.3         
#[33] tidyverse_2.0.0       

