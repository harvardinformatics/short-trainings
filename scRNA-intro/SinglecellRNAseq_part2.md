---
title: "Harvard Informatics scRNA-seq workshop: part 2"
author: "Adam Freedman and Tim Sackton"
date: "Fall, 2025"
output:
  html_document:
    keep_md: true
  pdf_document: default
---

# Introduction to single-cell RNAseq: part2: moving beyond clustering

**Welcome to part 2 of the workshop on single-cell RNA-seq!**

### Recap of part 1

Last week we covered: 

* potential benefits of scRNAseq over bulk RNAseq analysis
* current sequencing technologies (all available at the Bauer Core)
* artefacts requiring data cleaning
* base-level QC metrics: stuff reported from instrument-associated software
* scRNAseq data structures
  * sparse matrices, cell barcode and feature tables,zero-inflation
* loading data with the R package Seurat
* data processing,filtering, and clustering

Specifically with respect to data pre-processing and clustering using a
publicly available human PBMC data set as an example, we: 

* created Seurat objects for our raw and unfiltered count matrices
* in Seurat, normalized the filtered counts with the SCTransform method
* did a preliminary clustering of cells in the filtered data
* removed ambient RNA contamination with *SoupX*, taking as input
  * raw counts
  * filtered counts
  * provisional clusters from the filtered Seurat object
* removed doublets with scDblFinder
* did threshold-based removal of low UMI count, low feature (gene) count, and high mtDNA% cells

### Today, we will cover: 

* using the PBMC dataset,demonstrate how to identify marker genes for cell clusters
* demonstrate one approach for annotating cell types in the PBMC data set using a reference database
* multi-sample analysis using mouse brain data
  * simple merging versus batch effect correction with sample integration
    * metrics for determining whether integration is necessary
* differential expression analysis within cell clusters across conditions

## Load/install libraries and set working directory

These will have to be installed in order to run this R markdown file.



``` r
installed_packages <- rownames(installed.packages())
for (pkg in c("Seurat", "tidyverse", "patchwork", "sctransform","SoupX",
                "cowplot","metap")) {
  if (!pkg %in% installed_packages) {
    install.packages(pkg, quiet = TRUE)
  }
  
  library(pkg, character.only = TRUE)
}

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

for (pkg in c("glmGamPoi", "scDblFinder", "scater", "SingleR","celldex",
                "multtest")) {
  if (!pkg %in% installed_packages) {
    BiocManager::install(pkg)
  }
  
  library(pkg, character.only = TRUE)
}
```


``` r
data_dir=getwd()
setwd(data_dir)
```

## Marker gene discovery

For much of the history of scRNA-seq, a standard part of data analysis workflows has been the identification of genes that are canonical markers for cell clusters. Lists of marker genes are a necessary prerequisite for any attempt at manual annotation of clusters, e.g. if you have prior information that a set of genes are highly up-regulated in a particular cell type--information obtained from, for example, bulk RNA-seq on sorted cells of that type--one can provisionally infer that a cluster of cells within similar up-regulation of that set of genes belong to that cell type. Detecting marker genes for clusters involves performing tests of differential expression between sets of cells of interest. Detecting marker genes for a particular cluster involves treating cells as replicates and testing for differential expression between cells in the cluster of interest and all other cells. To demonstrate this, we will load up the filtered *pbmc* Seurat object we saved in part 1 of this workshop:


``` r
pbmc_filtered_seurat <- readRDS("data/pbmc/pbmc4k_seurat_obj_filt_singlets_posthocfilt.rds")
filtered_umap_plot<-DimPlot(pbmc_filtered_seurat, label = TRUE) + ggtitle("Filtered")
filtered_umap_plot
```

![](SinglecellRNAseq_part2_v2_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

### Find all marker genes for cluster 8

In this case, marker genes are identified by looking for DE between
cluster 8 and the aggregation of all other clusters. With Seurat, the
standard syntax is to specify an *identity*, such that if only one is
supplied, testing is done for cluster specific genes.


``` r
cluster8.markers <- FindMarkers(pbmc_filtered_seurat, ident.1 = 8)
head(cluster8.markers, n = 5)
```

```
##                p_val avg_log2FC pct.1 pct.2     p_val_adj
## KLRF1   0.000000e+00   6.651095 0.920 0.022  0.000000e+00
## GZMB    0.000000e+00   4.010149 0.879 0.037  0.000000e+00
## CLIC3  4.434973e-294   5.245734 0.816 0.047 6.613876e-290
## GNLY   5.098425e-289   5.458362 0.971 0.087 7.603282e-285
## FGFBP2 1.716290e-283   5.185766 0.805 0.047 2.559504e-279
```

The default DE test method is the Wilcoxon Rank Sum test, although several other methods are available. It is worth looking at the help for the `FindMarkers()` function to determine whether other changes to default settings are appropriate for your data. *pct.1* and *pct.2* specify the fraction of cells expressing the gene in cluster 8 (i.e. *ident.1*) and all other cells, respectively, and log-fold change is with cluster 8 in the numerator, i.e. positive numbers indicate up-regulation in cluster 8. *p_val_adj* is the p-value adjusted for multiple comparisons using the Bonferroni method. Clearly, the top 5 genes are expressed in the majority of cluster 8 cells, and very few cells from all other clusters.

The above implementation of `FindMarkers()` includes genes that are both up and down-regulated in cluster 8, i.e.


``` r
downreg_cluster8_markers <- subset(cluster8.markers,avg_log2FC<0)
head(downreg_cluster8_markers,n=5)
```

```
##               p_val avg_log2FC pct.1 pct.2    p_val_adj
## RPL39  5.300508e-93  -1.354651 1.000     1 7.904648e-89
## RPL18A 2.092823e-84  -1.254406 0.994     1 3.121027e-80
## RPL34  1.065071e-82  -1.254684 1.000     1 1.588340e-78
## RPL13  1.950295e-78  -1.119134 1.000     1 2.908475e-74
## RPL21  6.803242e-77  -1.095427 1.000     1 1.014568e-72
```

One can also provide the option to only look for genes up-regulated in the cluster of interest, by supplying the `only.pos = TRUE` keyword argument in `FindMarkers()`. If one wants to simultaneously find makers for all clusters, you can instead use the `FindAllMarkers()` function. **NOTE**: this will take a lot of time to run.

Finally, one can supply more than one identity to focus on genes that distinguish different clusters (or groups of clusters). Looking at our UMAP plot, we might be interested in knowing which genes have expression patterns that distinguish clusters 2 from 4:


``` r
cluster2Vs4.markers <- FindMarkers(pbmc_filtered_seurat, ident.1 = 2, ident.2=4)
head(cluster2Vs4.markers, n = 5)
```

```
##                      p_val avg_log2FC pct.1 pct.2    p_val_adj
## CD8B          4.946079e-89   5.099766 0.800 0.057 7.376088e-85
## CD8A          1.218662e-57   4.191128 0.617 0.048 1.817390e-53
## RP11-291B21.2 1.029645e-47   5.104648 0.508 0.021 1.535510e-43
## CTSW          6.018199e-24   1.796136 0.527 0.174 8.974941e-20
## CD40LG        1.260292e-18  -2.870340 0.032 0.246 1.879474e-14
```

or genes whose expression distinguish clusters 1 and 12 from 9:


``` r
cluster1.12Vs9.markers <- FindMarkers(pbmc_filtered_seurat, ident.1 = c(1,12), ident.2=9)
head(cluster2Vs4.markers, n = 5)
```

```
##                      p_val avg_log2FC pct.1 pct.2    p_val_adj
## CD8B          4.946079e-89   5.099766 0.800 0.057 7.376088e-85
## CD8A          1.218662e-57   4.191128 0.617 0.048 1.817390e-53
## RP11-291B21.2 1.029645e-47   5.104648 0.508 0.021 1.535510e-43
## CTSW          6.018199e-24   1.796136 0.527 0.174 8.974941e-20
## CD40LG        1.260292e-18  -2.870340 0.032 0.246 1.879474e-14
```

**REMEMBER**: what you supply to *ident.1* will be the numerator in log-fold change calculations, and will be the up-regulated cluster(s) when LFC > 0 (and the down-regulated cluster when LFC < 0).

## Celltype annotation: beyond marker genes

Manual annotation of cell clusters is a time-consuming process, that also requires a great deal of expertise. To the extent that both requiring the leveraging of considerable expertise to annotate cell types, manual annotation is not much different from creating a cell atlas! Increasingly however, annotation of cell clusters observed in scRNA-seq data can be carried out with computational tools. These tools typically take as input your expression data (in R either a Seurat or SingleCellExperiment object), and the specification of a reference. The reference data set contains expression data for individual cells and cell type labels. In other words, your data are the query which is aligned to the reference such that cell type labels can be transferred from reference cells to the query cells that share a similar expression pattern. Of course, for this approach to work, there must be consistency between the gene symbols in the query and the reference databases.

### Annotating the PBMC dataset

We will demonstrate one method `SingleR`, for annotating cells in our filtered PBMC dataset (which we have already loaded as an *Rds* objects at the start of today's workshop.

#### Load the reference data set

The MonacoImmuneData is a comprised of normalized expression data from 114 bulk RNA-seq experiments on sorted immune cell populations, relevant for annotating immune cells such as PBMCs. We load it as follows:


``` r
reference <- celldex::MonacoImmuneData()
```

#### Bookkeeping

We found there can be issues with running `celldex` if the cache directory doesn't exist. If loading the referenceSo, we attempt to clear the cache, and in the console type `yes` if we are informed it doesn't exist and do we want to create it:


``` r
library(ExperimentHub)
eh <- ExperimentHub()
removeCache(eh)
```

And if you want more information, use help!


``` r
help(MonacoImmuneData)
```

#### Convert our pbmc data to an SCE object


``` r
sce <- as.SingleCellExperiment(pbmc_filtered_seurat)
```

#### Run SingleR


``` r
pred <- SingleR(test = sce, ref = reference, labels = reference$label.main)
```

#### Add predicted cell type labels to our data


``` r
pbmc_filtered_seurat$SingleR.labels <- pred$labels
table(pbmc_filtered_seurat$SingleR.labels)
```

```
## 
##         B cells    CD4+ T cells    CD8+ T cells Dendritic cells       Monocytes 
##             576            1077             597             160             536 
##        NK cells     Progenitors         T cells 
##             178              10             478
```

The labels can now be accessed via: 
* `pbmc_filtered_seurat$SingleR.labels` or
* `pbmc_filtered_seurat@meta.data$SingleR.labels`

### a peek at cluster labels and annotations


``` r
pbmc_annots_umap<-DimPlot(pbmc_filtered_seurat, group.by = "SingleR.labels", label = FALSE) +
  ggplot2::ggtitle("Cell annotations")
plot_grid(filtered_umap_plot,pbmc_annots_umap,ncol=2,nrow=1)
```

![](SinglecellRNAseq_part2_v2_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

Note that our clustering reveals finer grained subdivision than is implied simply by the cell type annotations. This can be caused by multiple sources, including: 

* coarser resolution of the cell type reference database
* missing cell types
* not able to distinguish subtypes
  * for example, our annotation does not specify memory and naive CD4+T subtypes
* clustering algorithm over-splitting of functional cell types

### Annotation caveats and future directions

* building a reference database for annotating cell types is a non-trivial challenge
* other "current" tools may not be adequately maintained and hard to implement
  *the Seurat Azimuth method no longer appears to work, likely due to outdated reference file specifications
* "best practice" for cell type annotation and is an area of active development
* we expect that AI-based/transformer methods will likely factor prominently in the future



## Merging, integrating and DE analysis

As we have demonstrated above, marker gene discovery for cell clusters *within* a sample is straightforward and is based upon testing for differential expression between a defined cluster of cells versus all other cells, or between sets of clusters of interest. In a multi-sample scRNA-seq study, where individual samples represent experimental conditions of interest (e.g. treated versus control), differential expression testing is not as straightforward. In this context, one typically is interested in understanding how expression changes *within a cell type* across conditions. The biggest challenge is to correctly assign cells across experimental conditions, e.g. so that cells of cell type A subject to conditions x and y, respectively, are both assigned to the same same cell type. The two methods for doing this are **merging** and **batch integration**.

##### When to merge?

Merging assumes that there are no batch effects other than the effects of the different experimental conditions. This might be the case where samples are taken from a cell line with one sample (or a set of samples) are exposed to a perturbation treatment, and a second sample (or set of samples) is an untreated control, and then scRNA-seq libraries are generated simultaneously for all the samples using the same chemistry, then multiplexed and sequenced. Merging in this case, involved concatenating the count matrices from the samples, while using the cell barcodes track which cells belong to which condition.

##### Whent to integrate?

Often times, a study is conducted in a manner in which batch effects that have nothing to do with the biological variation of interest are unavoidable. This typically is the case when you are merging data you have generated with previously published data (or data you have generated previously) using different experimental protocols, using different library chemistry, a different sequencing instrument, or some combination of all three. Batch integration involves a formal statistical procedure in which the samples are aligned together in low-dimensional space. Batch integration is still an active area of research as there is a well-known tradeoff between optimal removal of batch effects and preservation of biological variation. In other words, it is hard to completely remove batch effects without removing some biological variation--variation which might be relevant for the experiment that has been conducted! A good starting point for understanding the scope of the problem is [Luecken et al. 2022, Nature Methods](https://www.nature.com/articles/s41592-021-01336-8) which compares the performance of integration methods for atlas-scale data sets.

##### Which is best for your data?

There are two complementary approaches for assessing whether integration
is necessary:

1. Construct a UMAP plot of a merged data set, with cells labelled by batch. Batch effects would be indicated by:
  a.  poor mixing of batches across clusters
  b.  clusters that are unique to a batch
2.  Compute summary statistics that describe the degree of mixing between batches

We will demonstrate both approaches below.

### Merging case study: male vs. female MOp brain regions

For our first example, we will use two mouse samples that were used as part of the [Allen BrainAtlas](https://mouse.brain-map.org/static/atlas). Both comprise cells from the primary motor areas (MOp) of the brain, with libraries built from the same 10x v3 chemistry, and both mice have the same genotype.The only difference between them besides sex is that the female sample was harvested at 58 days, while the male sample was harvested at 61 days. We begin with, and test, the assumption that there should be little if any batch effect. To peform this test, we do the following:

1. For each sample separately: 
  a. remove ambient RNA contamination with SoupX
  b. remove doublets with scDblFinder
  c. manually filter out low count and high mtDNA% cells on reasonable thresholds
2. Then merge the samples
  a. examine the UMAP plot
  b. calculate LISI metric to quantify batch effects
  
#### Load data, create Seurat objects and run SoupX


``` r
#### male1
male1_raw <- Seurat::Read10X("data/mousebrain/male_MOp_L8TX_181211_01_G12/raw_feature_bc_matrix")
male1_filtered <- Seurat::Read10X("data/mousebrain/male_MOp_L8TX_181211_01_G12/filtered_feature_bc_matrix") 

male1_seurat<- CreateSeuratObject(counts = male1_filtered)
male1_seurat <- PercentageFeatureSet(male1_seurat, pattern = "^mt-", col.name = "percent.mt")
male1_seurat <- SCTransform(male1_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
male1_seurat <- RunPCA(male1_seurat, verbose = FALSE)
male1_seurat <- RunUMAP(male1_seurat, dims = 1:30)
male1_seurat <- FindNeighbors(male1_seurat, dims = 1:30)
male1_seurat <- FindClusters(male1_seurat)

male1_soup_channel <- SoupX::SoupChannel(tod = male1_raw,toc=male1_filtered,
                is10X = TRUE)
male1_soup_channel$tod<-male1_raw
male1_soup_channel <- SoupX::setClusters(male1_soup_channel, 
                                   clusters = as.factor(Idents(male1_seurat)))
male1_soup_channel <- setDR(male1_soup_channel, 
                DR=Seurat::Embeddings(male1_seurat, "umap"))
male1_soup_channel <- autoEstCont(male1_soup_channel)
```

``` r
male1_corrected_counts <- adjustCounts(male1_soup_channel,roundToInt=TRUE)
male1_seurat_soupx <- CreateSeuratObject(counts = male1_corrected_counts)

male1_seurat <- NULL
male1_corrected_counts <- NULL
male1_soup_channel <- NULL
male1_raw <-NULL
male1_filtered <- NULL

##### male2
male2_raw <- Seurat::Read10X("data/mousebrain/457909/raw_feature_bc_matrix")
male2_filtered <- Seurat::Read10X("data/mousebrain/457909/filtered_feature_bc_matrix") 

male2_seurat<- CreateSeuratObject(counts = male2_filtered)
male2_seurat <- PercentageFeatureSet(male2_seurat, pattern = "mt-", col.name = "percent.mt")
male2_seurat <- SCTransform(male2_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
male2_seurat <- RunPCA(male2_seurat, verbose = FALSE)
male2_seurat <- RunUMAP(male2_seurat, dims = 1:30)
male2_seurat <- FindNeighbors(male2_seurat, dims = 1:30)
male2_seurat <- FindClusters(male2_seurat)

male2_soup_channel <- SoupX::SoupChannel(tod = male2_raw,toc=male2_filtered,
                is10X = TRUE)
male2_soup_channel$tod<-male2_raw
male2_soup_channel <- SoupX::setClusters(male2_soup_channel, 
                                   clusters = as.factor(Idents(male2_seurat)))
male2_soup_channel <- setDR(male2_soup_channel, 
                DR=Seurat::Embeddings(male2_seurat, "umap"))
male2_soup_channel <- autoEstCont(male2_soup_channel)
```

``` r
male2_corrected_counts <- adjustCounts(male2_soup_channel,roundToInt=TRUE)
male2_seurat_soupx <- CreateSeuratObject(counts = male2_corrected_counts)

male2_seurat <- NULL
male2_corrected_counts <- NULL
male2_soup_channel <- NULL
male2_raw <-NULL
male2_filtered <- NULL

##### female1
female1_raw <- Seurat::Read10X("data/mousebrain/female_MOp_L8TX_181211_01_C01/raw_feature_bc_matrix")
female1_filtered <- Seurat::Read10X("data/mousebrain/female_MOp_L8TX_181211_01_C01/filtered_feature_bc_matrix")  
female1_seurat<- CreateSeuratObject(counts = female1_filtered)
female1_seurat <- PercentageFeatureSet(female1_seurat, pattern = "^mt-", col.name = "percent.mt")
female1_seurat <- SCTransform(female1_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
female1_seurat <- RunPCA(female1_seurat, verbose = FALSE)
female1_seurat <- RunUMAP(female1_seurat, dims = 1:30)
female1_seurat <- FindNeighbors(female1_seurat, dims = 1:30)
female1_seurat <- FindClusters(female1_seurat)

female1_soup_channel <- SoupX::SoupChannel(tod = female1_raw,toc=female1_filtered,
                is10X = TRUE)
female1_soup_channel$tod<-female1_raw
female1_soup_channel <- SoupX::setClusters(female1_soup_channel, 
                                   clusters = as.factor(Idents(female1_seurat)))
female1_soup_channel <- setDR(female1_soup_channel, 
                DR=Seurat::Embeddings(female1_seurat, "umap"))
female1_soup_channel <- autoEstCont(female1_soup_channel)
```

``` r
female1_corrected_counts <- adjustCounts(female1_soup_channel,roundToInt=TRUE)
female1_seurat_soupx <- CreateSeuratObject(counts = female1_corrected_counts)

female1_seurat <- NULL
female1_corrected_counts <- NULL
female1_soup_channel <- NULL
female1_raw <-NULL
female1_filtered <- NULL

##### female2
female2_raw <- Seurat::Read10X("data/mousebrain/457911/raw_feature_bc_matrix")
female2_filtered <- Seurat::Read10X("data/mousebrain/457911/filtered_feature_bc_matrix") 

female2_seurat<- CreateSeuratObject(counts = female2_filtered)
female2_seurat <- PercentageFeatureSet(female2_seurat, pattern = "^mt-", col.name = "percent.mt")
female2_seurat <- SCTransform(female2_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
female2_seurat <- RunPCA(female2_seurat, verbose = FALSE)
female2_seurat <- RunUMAP(female2_seurat, dims = 1:30)
female2_seurat <- FindNeighbors(female2_seurat, dims = 1:30)
female2_seurat <- FindClusters(female2_seurat)

female2_soup_channel <- SoupX::SoupChannel(tod = female2_raw,toc=female2_filtered,
                is10X = TRUE)
female2_soup_channel$tod<-female2_raw
female2_soup_channel <- SoupX::setClusters(female2_soup_channel, 
                                   clusters = as.factor(Idents(female2_seurat)))
female2_soup_channel <- setDR(female2_soup_channel, 
                DR=Seurat::Embeddings(female2_seurat, "umap"))
female2_soup_channel <- autoEstCont(female2_soup_channel)
```

``` r
female2_corrected_counts <- adjustCounts(female2_soup_channel,roundToInt=TRUE)
female2_seurat_soupx <- CreateSeuratObject(counts = female2_corrected_counts)

female2_seurat <- NULL
female2_corrected_counts <- NULL
female2_soup_channel <- NULL
female2_raw <-NULL
female2_filtered <- NULL

#### female3
female3_raw <- Seurat::Read10X("data/mousebrain/500199/raw_feature_bc_matrix")
female3_filtered <- Seurat::Read10X("data/mousebrain/500199/filtered_feature_bc_matrix") 

female3_seurat<- CreateSeuratObject(counts = female3_filtered)
female3_seurat <- PercentageFeatureSet(female3_seurat, pattern = "^mt-", col.name = "percent.mt")
female3_seurat <- SCTransform(female3_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
female3_seurat <- RunPCA(female3_seurat, verbose = FALSE)
female3_seurat <- RunUMAP(female3_seurat, dims = 1:30)
female3_seurat <- FindNeighbors(female3_seurat, dims = 1:30)
female3_seurat <- FindClusters(female3_seurat)

female3_soup_channel <- SoupX::SoupChannel(tod = female3_raw,toc=female3_filtered,
                is10X = TRUE)
female3_soup_channel$tod<-female3_raw
female3_soup_channel <- SoupX::setClusters(female3_soup_channel, 
                                   clusters = as.factor(Idents(female3_seurat)))
female3_soup_channel <- setDR(female3_soup_channel, 
                DR=Seurat::Embeddings(female3_seurat, "umap"))
female3_soup_channel <- autoEstCont(female3_soup_channel)
```

``` r
female3_corrected_counts <- adjustCounts(female3_soup_channel,roundToInt=TRUE)
female3_seurat_soupx <- CreateSeuratObject(counts = female3_corrected_counts)

female3_seurat <- NULL
female3_corrected_counts <- NULL
female3_soup_channel <- NULL
female3_raw <-NULL
female3_filtered <- NULL

#### male_cnupal
male_cnupal_raw <- Seurat::Read10X("data/mousebrain/male_CNU-PAL_L8TX_190327_01_E04/raw_feature_bc_matrix")
male_cnupal_filtered <- Seurat::Read10X("data/mousebrain/male_CNU-PAL_L8TX_190327_01_E04/filtered_feature_bc_matrix") 

male_cnupal_seurat<- CreateSeuratObject(counts = male_cnupal_filtered)
male_cnupal_seurat <- PercentageFeatureSet(male_cnupal_seurat, pattern = "^mt-", col.name = "percent.mt")
male_cnupal_seurat <- SCTransform(male_cnupal_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
male_cnupal_seurat <- RunPCA(male_cnupal_seurat, verbose = FALSE)
male_cnupal_seurat <- RunUMAP(male_cnupal_seurat, dims = 1:30)
male_cnupal_seurat <- FindNeighbors(male_cnupal_seurat, dims = 1:30)
male_cnupal_seurat <- FindClusters(male_cnupal_seurat)

male_cnupal_soup_channel <- SoupX::SoupChannel(tod = male_cnupal_raw,toc=male_cnupal_filtered,
                is10X = TRUE)
male_cnupal_soup_channel$tod<-male_cnupal_raw
male_cnupal_soup_channel <- SoupX::setClusters(male_cnupal_soup_channel, 
                                   clusters = as.factor(Idents(male_cnupal_seurat)))
male_cnupal_soup_channel <- setDR(male_cnupal_soup_channel, 
                DR=Seurat::Embeddings(male_cnupal_seurat, "umap"))
male_cnupal_soup_channel <- autoEstCont(male_cnupal_soup_channel)
```

``` r
male_cnupal_corrected_counts <- adjustCounts(male_cnupal_soup_channel,roundToInt=TRUE)
male_cnupal_seurat_soupx <- CreateSeuratObject(counts = male_cnupal_corrected_counts)

male_cnupal_seurat <- NULL
male_cnupal_corrected_counts <- NULL
male_cnupal_soup_channel <- NULL
male_cnupal_raw <-NULL
male_cnupal_filtered <- NULL
rm(list = names(which(sapply(ls(), function(x) is.null(get(x))))))
```

#### Doublet removal with scDblFinder

##### convert objects to SCE


``` r
sce_male1 <- as.SingleCellExperiment(male1_seurat_soupx)
sce_male2 <- as.SingleCellExperiment(male2_seurat_soupx)
sce_female1 <- as.SingleCellExperiment(female1_seurat_soupx)
sce_female2 <- as.SingleCellExperiment(female2_seurat_soupx)
sce_female3 <- as.SingleCellExperiment(female3_seurat_soupx)
sce_malecnupal <- as.SingleCellExperiment(male_cnupal_seurat_soupx)
```

#### Run scDblFinder


``` r
sce_male1 <- scDblFinder(sce_male1)
sce_male2 <- scDblFinder(sce_male2)
sce_female1 <- scDblFinder(sce_female1)
sce_female2 <- scDblFinder(sce_female2)
sce_female3 <- scDblFinder(sce_female3)
sce_malecnupal <- scDblFinder(sce_malecnupal)
```

##### Filter Seurat objects on singlet class


``` r
male1_seurat_soupx$scDblFinder.class <- colData(sce_male1)$scDblFinder.class
male1_seurat_singlets <- subset(male1_seurat_soupx, subset = scDblFinder.class == "singlet")
male1_seurat_soupx <-NULL

male2_seurat_soupx$scDblFinder.class <- colData(sce_male2)$scDblFinder.class
male2_seurat_singlets <- subset(male2_seurat_soupx, subset = scDblFinder.class == "singlet")
male2_seurat_soupx <-NULL

female1_seurat_soupx$scDblFinder.class <- colData(sce_female1)$scDblFinder.class
female1_seurat_singlets <- subset(female1_seurat_soupx, subset = scDblFinder.class == "singlet")
female1_seurat_soupx <-NULL

female2_seurat_soupx$scDblFinder.class <- colData(sce_female2)$scDblFinder.class
female2_seurat_singlets <- subset(female2_seurat_soupx, subset = scDblFinder.class == "singlet")
female2_seurat_soupx <-NULL

female3_seurat_soupx$scDblFinder.class <- colData(sce_female3)$scDblFinder.class
female3_seurat_singlets <- subset(female3_seurat_soupx, subset = scDblFinder.class == "singlet")
female3_seurat_soupx <-NULL

male_cnupal_seurat_soupx$scDblFinder.class <- colData(sce_malecnupal)$scDblFinder.class
male_cnupal_seurat_singlets <- subset(male_cnupal_seurat_soupx, subset = scDblFinder.class == "singlet")
male_cnupal_seurat_soupx <-NULL

sce_male1 <- NULL
sce_male2 <- NULL
sce_female1 <- NULL
sce_female2 <- NULL
sce_female3 <- NULL
sce_male_cnupal <- NULL

rm(list = names(which(sapply(ls(), function(x) is.null(get(x))))))
```

#### add mtDNA counts to object metadata


``` r
male1_seurat_singlets <- PercentageFeatureSet(male1_seurat_singlets, pattern = "^mt-", col.name = "percent.mt")
male2_seurat_singlets <- PercentageFeatureSet(male2_seurat_singlets, pattern = "^mt-", col.name = "percent.mt")
female1_seurat_singlets <- PercentageFeatureSet(female1_seurat_singlets, pattern = "^mt-", col.name = "percent.mt")
female2_seurat_singlets <- PercentageFeatureSet(female2_seurat_singlets, pattern = "^mt-", col.name = "percent.mt")
female3_seurat_singlets <- PercentageFeatureSet(female3_seurat_singlets, pattern = "^mt-", col.name = "percent.mt")
male_cnupal_seurat_singlets <- PercentageFeatureSet(male_cnupal_seurat_singlets, pattern = "^mt-", col.name = "percent.mt")
```

#### post hoc filtering

In part 1 of this workshop, we demonstrated how to use median absolute
deviations (MAD) for threshold-based filtering, but today we will simply
use manually set filters:


``` r
male1_seurat_singlets_posthoc <- subset(male1_seurat_singlets,subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 5)
male2_seurat_singlets_posthoc <- subset(male2_seurat_singlets,subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 5)
female1_seurat_singlets_posthoc <- subset(female1_seurat_singlets,subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 5)
female2_seurat_singlets_posthoc <- subset(female2_seurat_singlets,subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 5)
female3_seurat_singlets_posthoc <- subset(female3_seurat_singlets,subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 5)

male_cnupal_seurat_singlets_posthoc <- subset(male_cnupal_seurat_singlets,subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 5)

male1_seurat_singlets <- NULL
male2_seurat_singlets <- NULL
female1_seurat_singlets <- NULL
female2_seurat_singlets <- NULL
female3_seurat_singlets <- NULL
male_cnupal_seurat_singlets <- NULL

rm(list = names(which(sapply(ls(), function(x) is.null(get(x))))))
```

#### merge data sets

##### add sex,replicate and tissue data to seurat objects


``` r
male1_seurat_singlets_posthoc[["sex"]] <- "male"
male1_seurat_singlets_posthoc[["replicate"]] <- 1
male1_seurat_singlets_posthoc[["tissue"]] <- "mop"

male2_seurat_singlets_posthoc[["sex"]] <- "male"
male2_seurat_singlets_posthoc[["replicate"]] <- 2
male2_seurat_singlets_posthoc[["tissue"]] <- "mop"

female1_seurat_singlets_posthoc[["sex"]] <- "female"
female1_seurat_singlets_posthoc[["replicate"]] <- 1
female1_seurat_singlets_posthoc[["tissue"]] <- "mop"

female2_seurat_singlets_posthoc[["sex"]] <- "female"
female2_seurat_singlets_posthoc[["replicate"]] <- 2
female2_seurat_singlets_posthoc[["tissue"]] <- "mop"

female3_seurat_singlets_posthoc[["sex"]] <- "female"
female3_seurat_singlets_posthoc[["replicate"]] <- 3
female3_seurat_singlets_posthoc[["tissue"]] <- "mop"

male_cnupal_seurat_singlets_posthoc[["sex"]] <- "male"
male_cnupal_seurat_singlets_posthoc[["replicate"]] <- 1
male_cnupal_seurat_singlets_posthoc[["tissue"]] <- "cnupal"
```

##### merge objects


``` r
merged_mop <- merge(male1_seurat_singlets_posthoc, y = c(male2_seurat_singlets_posthoc, female1_seurat_singlets_posthoc, female2_seurat_singlets_posthoc, female3_seurat_singlets_posthoc), add.cell.ids = c("male1", "male2","female1","female2","female3"), project = "merged_mop")
```

If we look at the metadata of the merged object in `@meta.data` we will
now see the replicate-sex combinations appended as cell prefixes, and as
a column variable


``` r
merged_mop[["sex_replicate"]] <- paste(merged_mop$sex,merged_mop$replicate,sep="_")
head(merged_mop@meta.data)
```

```
##                             orig.ident nCount_RNA nFeature_RNA
## male1_AAACCCAAGCTTCATG-1 SeuratProject      34586         6875
## male1_AAACCCAAGTGAGGTC-1 SeuratProject      22711         5816
## male1_AAACCCAGTGAACGGT-1 SeuratProject      61795         8412
## male1_AAACCCAGTGGCATCC-1 SeuratProject      46467         7716
## male1_AAACCCAGTTTGCAGT-1 SeuratProject        867          678
## male1_AAACCCATCTACCTTA-1 SeuratProject      18946         5130
##                          scDblFinder.class percent.mt  sex replicate tissue
## male1_AAACCCAAGCTTCATG-1           singlet  0.8037934 male         1    mop
## male1_AAACCCAAGTGAGGTC-1           singlet  1.8229052 male         1    mop
## male1_AAACCCAGTGAACGGT-1           singlet  2.6150983 male         1    mop
## male1_AAACCCAGTGGCATCC-1           singlet  1.9368584 male         1    mop
## male1_AAACCCAGTTTGCAGT-1           singlet  0.4613610 male         1    mop
## male1_AAACCCATCTACCTTA-1           singlet  1.0609100 male         1    mop
##                          sex_replicate
## male1_AAACCCAAGCTTCATG-1        male_1
## male1_AAACCCAAGTGAGGTC-1        male_1
## male1_AAACCCAGTGAACGGT-1        male_1
## male1_AAACCCAGTGGCATCC-1        male_1
## male1_AAACCCAGTTTGCAGT-1        male_1
## male1_AAACCCATCTACCTTA-1        male_1
```

##### downsample merged object (for workshop only!)


``` r
N<-3000
cell_metadata <- merged_mop@meta.data
cell_metadata$cell_id <- rownames(cell_metadata)

cell_subset <- cell_metadata %>%
  group_by(sex_replicate) %>%
  slice_sample(n = N) %>%
  pull(cell_id)

downsamp_mop_merged <- subset(merged_mop, cells = cell_subset)
```

#### run clustering on downsampled merged dataset


``` r
downsamp_mop_merged <- PercentageFeatureSet(downsamp_mop_merged, pattern = "^mt-", col.name = "percent.mt")
downsamp_mop_merged <- SCTransform(downsamp_mop_merged, vars.to.regress = "percent.mt", verbose = FALSE)
downsamp_mop_merged <- RunPCA(downsamp_mop_merged, verbose = FALSE)
downsamp_mop_merged <- RunUMAP(downsamp_mop_merged, dims = 1:30)
downsamp_mop_merged <- FindNeighbors(downsamp_mop_merged, dims = 1:30)
downsamp_mop_merged <- FindClusters(downsamp_mop_merged)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 15000
## Number of edges: 489048
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9121
## Number of communities: 26
## Elapsed time: 1 seconds
```

#### create UMAP plot labelled by sex


``` r
bysexplot<-DimPlot(downsamp_mop_merged, reduction = "umap", group.by = c("sex"),alpha=0.3)
byclusterplot<-DimPlot(downsamp_mop_merged, reduction = "umap",label=TRUE)
bysexreplicateplot<-DimPlot(downsamp_mop_merged, reduction = "umap",group.by=c("sex_replicate"))
bysexplot
```

![](SinglecellRNAseq_part2_v2_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

``` r
byclusterplot
```

![](SinglecellRNAseq_part2_v2_files/figure-html/unnamed-chunk-27-2.png)<!-- -->

``` r
bysexreplicateplot
```

![](SinglecellRNAseq_part2_v2_files/figure-html/unnamed-chunk-27-3.png)<!-- -->

While the UMAP plot on the left suggests there is overall good mixing between male and female samples, such that batch effects might not be a big problem, it is notable that the large cluster at the center of the plot seems to have a few continuous bands of male cells at its periphery, and at the center of the cluster as well. We have downsampled all samples to the same number of cells, such that the male-specific aggregations at the periphery of clusters cannot simply be due to having sampled more cells from male samples.In fact, we have sampled more female samples than male samples, so one might expect the opposite pattern to hold. That said, the UMAP method tends to inflate the distance between clusters, such that any clusters that are dominated by cells from one batch may superficially inflate the degree of a batch effect. 


#### Differential expression testing

There are two types of differential expression analysis one is
interested in:

1.  tests that look for conserved cell type markers regardless of the experimental condition those cells belong to. This involves using the `FindConservedMarkers()` function in Seurat, which, for the cluster of interest, performs separate differential expression tests for each experimental condition (in our case, sex) and then uses a meta-analysis approach implemented in the R *MetaDE* package for combining p-values. This is analogous to marker gene discovery (as demonstrated above) but combining the results of condition-specific tests.

2. tests that look for DE within clusters across conditions, analogous to what one normally does with bulk RNA-seq data. This involves: 
  a."pseudo-bulking", such that for each replicate/sample, for each defined cluster of cells, for each gene, counts for that gene are summed. This prouduces a single count per gene per replicate. then
  b. DE testing, often using statistical models that are employed with bulk RNA-seq, e.g. `DESeq2`.

##### Look for conserved markers for cluster 0

To use the SCT-transformed counts for downstream differential expression tests, one has to generated counts from the SCT-transformed values, effectively reversing the SCT regression model, so at to populate the `counts` and `data` slots in `SCT` assay, with the former being corrected counts that represent a conversion-to-count of the SCTransformed values. We do this with Seurat's `PrepSCTFindMarkers()` function.


``` r
options(future.globals.maxSize = 8 * 1024^3)
downsamp_mop_merged <- PrepSCTFindMarkers(downsamp_mop_merged, verbose = FALSE)
```

Then, to find conserved makers for a cluster (in our example, cluster 0) such that expression differences between that cluster and all other cells is conserved across sexes, we can do the following


``` r
Idents(downsamp_mop_merged)<- "seurat_clusters"
cluster0_conserved_markers <- FindConservedMarkers(downsamp_mop_merged, ident.1 = 0, grouping.var = "sex", assay = "SCT",verbose = FALSE)
head(cluster0_conserved_markers)
```

```
##         male_p_val male_avg_log2FC male_pct.1 male_pct.2 male_p_val_adj
## Stard8           0        3.236120      0.882      0.165              0
## Prex1            0        2.802730      0.918      0.226              0
## Calb1            0        2.656055      0.912      0.220              0
## Gucy1a1          0        2.854776      0.904      0.237              0
## Tesc             0        2.515893      0.923      0.266              0
## Lamp5            0        2.221590      0.976      0.409              0
##         female_p_val female_avg_log2FC female_pct.1 female_pct.2
## Stard8             0          3.149816        0.859        0.166
## Prex1              0          2.690153        0.908        0.250
## Calb1              0          2.661769        0.908        0.186
## Gucy1a1            0          2.644162        0.873        0.239
## Tesc               0          2.272215        0.912        0.293
## Lamp5              0          2.132868        0.954        0.391
##         female_p_val_adj max_pval minimump_p_val
## Stard8                 0        0              0
## Prex1                  0        0              0
## Calb1                  0        0              0
## Gucy1a1                0        0              0
## Tesc                   0        0              0
## Lamp5                  0        0              0
```

The output of `FindConservedMarkers` includes all the fields you'd expect with `FindMarkers()` and in addition, *minimump_p_vale* which is the combined p-value.

##### Differential expression within clusters between sexes

In our merged, downsampled Seurat object `downsamp_mop_merged`, we have three replicates for female MOp and two for males. The robust approach for testing for DE within clusters is to do *pseudobulking*, where for each sample, counts for each cell within a cluster are aggregated. Then tools such as `DESeq2` can be used in a manner analogous to bulk RNA-seq analysis. Note that, as with the implementation of `FindConservedMarkers()`, we can use the counts back-calculated from the SCTransformed data, rather than the raw, unnormalized counts.


``` r
pseudobulk_counts <- AggregateExpression(downsamp_mop_merged, return.seurat = T,
                                           assays = "SCT",
                                           slot = "counts",
                                           group.by = c("seurat_clusters","sex_replicate"))
```

Prior to DE testing, we need to add a *sex* variable back into the `pseudobulk_counts`, as all that exists is `sex_replicate`.


``` r
sex <- mutate(as_tibble(pseudobulk_counts@meta.data),sex = str_remove(sex_replicate, "-.*")) %>% select(sex)
pseudobulk_counts$sex <- sex$sex
```

Now, we can do DE testing for cluster 0. Note that Seurat's `AggregateExpression` function puts a "g" before the cluster labels, so this has to be added when subsetting the data. Also note, the Seurat default is to have a minimum number of replicates set to 3, and will throw an error if there are fewer replicates, as with the case for our two male samples. There are methodological performance issues with so few replicates, but for the purpose of demonstration, we can still use the `DESeq2` model, by setting `min.cells.group = 2`. Remember too, the default assay is *SCT* so we are using the SCT-converted counts for testing as when we find conserved markers (above). To do DE testing within a cluster, we can subset the data, and then run `FindMarkers`, doing so in a way that contrasts the specified conditions:



``` r
cluster0_pseudobulk <- subset(pseudobulk_counts, seurat_clusters == "g0")
  Idents(cluster0_pseudobulk) <- "sex"
cluster0_de_markers <- FindMarkers(cluster0_pseudobulk, ident.1 = "male", ident.2 = "female", slot = "counts", test.use = "DESeq2", min.cells.group = 2, verbose = F)
head(cluster0_de_markers)
```

```
##                p_val avg_log2FC pct.1 pct.2    p_val_adj
## Xist    6.644944e-23 -3.7768126     0     1 1.474779e-18
## Ddx3y   8.856161e-20  1.7217314     1     0 1.965536e-15
## Uty     1.648811e-19  1.6754842     1     0 3.659370e-15
## Eif2s3y 2.110938e-19  1.6739161     1     0 4.685017e-15
## Tsix    2.479584e-17 -0.1174559     1     1 5.503189e-13
## Kdm5d   1.152393e-15  1.1881395     1     0 2.557620e-11
```

And then we can make a nice volcano plot!


``` r
library(ggrepel)
cluster0_de_markers$gene <- rownames(cluster0_de_markers)
volcano_plot<- ggplot(cluster0_de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
      ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,
      "")), colour = "red", size = 3)
volcano_plot
```

![](SinglecellRNAseq_part2_v2_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

We can use the same `FindMarkers()` function to do differential expression testing across conditions. Remember, the default assay is *SCT* so we are using the SCT-converted counts for testing as when we find conserved markers (above).

An alternative code implementation, and which won't require subsetting for additional cluster tests would be to set a composite variable for cluster+sex, then set that variable with `Idents`, then use the composite variable for specifying the groups for which to test for DE:


``` r
pseudobulk_counts$cluster.sex <- paste(pseudobulk_counts$seurat_clusters, pseudobulk_counts$sex, sep = "_")
Idents(pseudobulk_counts) <- pseudobulk_counts$cluster.sex
cluster0_de_markers_alt <- FindMarkers(pseudobulk_counts, ident.1 = "g0_male", ident.2 = "g0_female", slot = "counts", test.use = "DESeq2", min.cells.group = 2, verbose = F)
```

```
## converting counts to integer mode
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

``` r
head(cluster0_de_markers_alt)
```

```
##                p_val avg_log2FC pct.1 pct.2    p_val_adj
## Xist    6.644944e-23 -3.7768126     0     1 1.474779e-18
## Ddx3y   8.856161e-20  1.7217314     1     0 1.965536e-15
## Uty     1.648811e-19  1.6754842     1     0 3.659370e-15
## Eif2s3y 2.110938e-19  1.6739161     1     0 4.685017e-15
## Tsix    2.479584e-17 -0.1174559     1     1 5.503189e-13
## Kdm5d   1.152393e-15  1.1881395     1     0 2.557620e-11
```

#### Integration case study

To demonstrate how to do integration, we will consider a subset of our samples, namely 1 male MOp brain region sample and another male sample from the CNU-PAL region.


``` r
merged_mop_cnupal <- merge(male1_seurat_singlets_posthoc, y = male_cnupal_seurat_singlets_posthoc, add.cell.ids = c("male_mop", "male_cnupal"), project = "merged_mopcnupal")

N<-3000
cell_metadata_mopcnupal <- merged_mop_cnupal@meta.data
cell_metadata_mopcnupal$cell_id <- rownames(cell_metadata_mopcnupal)

cell_subset <- cell_metadata_mopcnupal %>%
  group_by(tissue) %>%
  slice_sample(n = N) %>%
  pull(cell_id)

downsamp_mop_cnupal_merged <- subset(merged_mop_cnupal, cells = cell_subset)

downsamp_mop_cnupal_merged <- PercentageFeatureSet(downsamp_mop_cnupal_merged, pattern = "^mt-", col.name = "percent.mt")
downsamp_mop_cnupal_merged <- SCTransform(downsamp_mop_cnupal_merged, vars.to.regress = "percent.mt", verbose = FALSE)
downsamp_mop_cnupal_merged <- RunPCA(downsamp_mop_cnupal_merged, verbose = FALSE)
downsamp_mop_cnupal_merged <- RunUMAP(downsamp_mop_cnupal_merged, dims = 1:30)
downsamp_mop_cnupal_merged <- FindNeighbors(downsamp_mop_cnupal_merged, dims = 1:30)
downsamp_mop_cnupal_merged <- FindClusters(downsamp_mop_cnupal_merged)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6000
## Number of edges: 209331
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9179
## Number of communities: 29
## Elapsed time: 0 seconds
```

``` r
bytissueplot<-DimPlot(downsamp_mop_cnupal_merged, reduction = "umap", group.by = c("tissue"),alpha=0.3)
bytissueplot
```

![](SinglecellRNAseq_part2_v2_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

There are rather pronounced batch effects by tissue ... which we might expect given likely tissue-specific cell type composition and expression profiles. One could simply proceed with integration, but it is worth using a quantitative metric to summarize the magnitude of batch effects, such as *LISI*.

##### quantitative assessment of batch effects
After visualization with UMAP, it is advised to use a statistical approach to quantify the batch effect. One useful statistic is the *local inverse Simpson's index* aka *LISI*. This metric is based upon probabilities that two randomly sampled cells come from different batches (in our case, male or female). The score is calculated for each cell and ranges from 1 to the total number of batches, with 1 meaning poor mixing and a value closer to the number of batches indicative of good mixing (and minimal batch effect). We need to use the `reticulate` library to calculate LISI in python, specifically using the `harmonypy` package. 



``` r
library(reticulate)
use_condaenv("r-lisi", conda="/Users/adamfreedman/miniforge3/bin/conda", required = TRUE)
py_discover_config()
```

```
## python:         /Users/adamfreedman/miniforge3/envs/r-lisi/bin/python
## libpython:      /Users/adamfreedman/miniforge3/envs/r-lisi/lib/libpython3.8.dylib
## pythonhome:     /Users/adamfreedman/miniforge3/envs/r-lisi:/Users/adamfreedman/miniforge3/envs/r-lisi
## version:        3.8.20 | packaged by conda-forge | (default, Sep 30 2024, 17:48:42)  [Clang 17.0.6 ]
## numpy:          /Users/adamfreedman/miniforge3/envs/r-lisi/lib/python3.8/site-packages/numpy
## numpy_version:  1.24.4
## 
## NOTE: Python version was forced by use_python() function
```

``` r
py_run_string("import harmonypy")
```


``` r
unintegrated_obj <- downsamp_mop_cnupal_merged
reduction_to_use <- "pca" 
batch_variable <- list("tissue")
embeddings <- Embeddings(unintegrated_obj, reduction = reduction_to_use)
metadata <- unintegrated_obj@meta.data
lisi_module <- import("harmonypy.lisi")
message("Calculating LISI... This may take a moment.")
```

```
## Calculating LISI... This may take a moment.
```

``` r
lisi_scores <- lisi_module$compute_lisi(
  X = embeddings,
  metadata = metadata,
  label_colnames = c(batch_variable),
  perplexity = 30L
)
lisi_col_name <- paste0(batch_variable, "_LISI")
unintegrated_obj[[lisi_col_name]] <- as.numeric(lisi_scores)
message(paste("LISI scores added to metadata column:", lisi_col_name))
```

```
## LISI scores added to metadata column: tissue_LISI
```

``` r
message(paste("Mean LISI score for unintgrated ata is",mean(unintegrated_obj@meta.data$tissue_LISI)))
```

```
## Mean LISI score for unintgrated ata is 1.02247557258345
```
This score is close to 1 indicating a very strong batch effect in the un-integrated data, i.e. on average any cell's neighborhood only contains one batch, i.e. the batch it belongs to. Thus, integration is crucial.


As an example of how to do this, we can wrap integration using the CCA method.  

**NOTE**: not all integration methods work well for all tasks, and so it is worth exploring options if a method produces poor results. While Harmony integration is known to do well for simple integration tasks, a test run using this method on this data set generated a very poor integration, with clear batch effects retained in the UMAP plot.

##### Perform CCA integration


``` r
integrated_downsamp_mop_cnupal <- IntegrateLayers(object = downsamp_mop_cnupal_merged,
    method = "CCAIntegration",
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    normalization.method = "SCT",
    verbose = FALSE)

integrated_downsamp_mop_cnupal <- FindNeighbors(integrated_downsamp_mop_cnupal, reduction = "integrated.cca", dims = 1:30)
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

``` r
integrated_downsamp_mop_cnupal <- FindClusters(integrated_downsamp_mop_cnupal, resolution = 1)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6000
## Number of edges: 187718
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8755
## Number of communities: 28
## Elapsed time: 0 seconds
```

``` r
integrated_downsamp_mop_cnupal <- RunUMAP(integrated_downsamp_mop_cnupal, dims = 1:30, reduction = "integrated.cca")
```

```
## 19:49:03 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## 19:49:03 Read 6000 rows and found 30 numeric columns
```

```
## 19:49:03 Using Annoy for neighbor search, n_neighbors = 30
```

```
## 19:49:03 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 19:49:03 Writing NN index file to temp file /var/folders/y4/5qbzd1h11vdb2lww93vpl_pc0000gq/T//RtmphX7GWS/file17f925b0c782f
## 19:49:03 Searching Annoy index using 1 thread, search_k = 3000
## 19:49:04 Annoy recall = 100%
## 19:49:04 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
## 19:49:05 Initializing from normalized Laplacian + noise (using RSpectra)
## 19:49:05 Commencing optimization for 500 epochs, with 247128 positive edges
## 19:49:05 Using rng type: pcg
## 19:49:09 Optimization finished
```

``` r
cca_integrate_umapplot <- DimPlot(integrated_downsamp_mop_cnupal, reduction = "umap", group.by = c("tissue"))
cca_cluster_integrate_umapplot <- DimPlot(integrated_downsamp_mop_cnupal, reduction = "umap",label=TRUE)
cca_cluster_integrate_umapplot + cca_integrate_umapplot 
```

![](SinglecellRNAseq_part2_v2_files/figure-html/unnamed-chunk-38-1.png)<!-- -->

While some of the batch effects have been corrected, there still appear to be substantial portions of the UMAP plot that are restricted to a particular tissue type. **This might very well be due to underlying biology and differences in cell type composition!**. We can calculate LISI for the integrate data to quantify the remaining degree of batch effect:


``` r
integrated_obj <- integrated_downsamp_mop_cnupal
reduction_to_use <- "integrated.cca" 
batch_variable <- list("tissue")
embeddings <- Embeddings(integrated_obj, reduction = reduction_to_use)
metadata <- integrated_obj@meta.data
lisi_module <- import("harmonypy.lisi")
message("Calculating LISI... This may take a moment.")
```

```
## Calculating LISI... This may take a moment.
```

``` r
lisi_scores <- lisi_module$compute_lisi(
  X = embeddings,
  metadata = metadata,
  label_colnames = c(batch_variable),
  perplexity = 30L
)
lisi_col_name <- paste0(batch_variable, "_LISI")
integrated_obj[[lisi_col_name]] <- as.numeric(lisi_scores)
message(paste("LISI scores added to metadata column:", lisi_col_name))
```

```
## LISI scores added to metadata column: tissue_LISI
```

``` r
message(paste("Mean LISI score for intgrated ata is",mean(integrated_obj@meta.data$tissue_LISI)))
```

```
## Mean LISI score for intgrated ata is 1.13677987175341
```
Indeed, as the UMAP suggested, batch effect correction was not all that successful. **If we were confident that integration had sufficiently corrected the batch effects, and if we had multiple replicates per biological condition, we could proceed with differential expression testing in the same fashion as we did above on the merged data**.

**But**, the poor outcome of integration leads to the following questions:

1. Should we even try to integrate data from distinct brain regions?
2. How much do batch effects have to be corrected to insure that downstream analyses are robust?
  a. Given *n* batches we are trying to integrate, is there a particular LISI score threshold (given a theoretical range of 1 to *n*) below which we should not do further analyses? 

#### The case of unreplicated conditions
While not ideal, there may be cases where one does not have biological replicates for each condition. For example, if integration between the male MOp and CNU-PAL brain regions had been successful, and that was all the preliminary data you had, is there a way to assess differential expression between these brain regions? One can't take a pseudo-bulking approach, because there are no replicates from which to aggregate expression. Instead, one can only use `FindMarkers` in a manner such that cells are treated as replicates--in the same way in which DE testing is done between clusters to identify marker genes for clusters. One could do this:


``` r
options(future.globals.maxSize = 2 * 1024^3)
integrated_downsamp_mop_cnupal <- PrepSCTFindMarkers(integrated_downsamp_mop_cnupal,verbose=TRUE)
```

```
## Found 2 SCT models. Recorrecting SCT counts using minimum median counts: 1377.5
```


``` r
integrated_downsamp_mop_cnupal$cluster_tissue <- paste(integrated_downsamp_mop_cnupal$seurat_clusters, integrated_downsamp_mop_cnupal$tissue, sep = "_")
Idents(integrated_downsamp_mop_cnupal) <- "cluster_tissue"
de_genes_cluster5 <- FindMarkers(
  integrated_downsamp_mop_cnupal,
  ident.1 = "5_mop",  
  ident.2 = "5_cnupal",
  verbose = FALSE
)
head(de_genes_cluster5)
```

```
##                 p_val avg_log2FC pct.1 pct.2    p_val_adj
## Tmem178b 2.030930e-75  -8.737610     0     1 4.546034e-71
## Sema6a   2.030930e-75  -8.737610     0     1 4.546034e-71
## Gad1     2.044068e-75  -8.889613     0     1 4.575441e-71
## Pcp4l1   2.052872e-75  -9.567685     0     1 4.595150e-71
## Lingo2   2.052872e-75  -9.655147     0     1 4.595150e-71
## Erbb4    2.059943e-75 -12.362100     0     1 4.610977e-71
```

**HOWEVER**...this is **NOT** the statistically most robust approach as:

* cells (with noisy counts) are treated as replicates
* lack of independence between cells means p-values are not necessarily accurate
* there is no way to determine how representative the one sample representing male or female of other individuals of that sex

### Conclusions/Future directions
