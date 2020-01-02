# ClusterMole

> See, children, the misguided Mole.  
> He lives down in a deep, dark hole;  
> Sweetness, and light, and good fresh air  
> Are things for which he does not care.  
> He has not even that makeshift  
> Of feeble minds - the social gift.  
> But say not that he has no soul,  
> Lest haply we misjudge the Mole;  
> Nay, if we measure him by men,  
> No doubt he sits in his dark den  
> Instructing others blind as he  
> Exactly how the world should be.  
>
> -- Oliver Herford

Cell type annotation of single-cell RNA sequencing (scRNA-seq) data typically requires a reference dataset, but finding an appropriate one may be challenging.
ClusterMole is an R package that provides a collection of cell type markers for thousands of human and mouse cell populations sourced from a variety of databases.

Install ClusterMole:

```r
BiocManager::install("igordot/clustermole", update = FALSE)
```

Retrieve a table of all cell type markers:

```r
markers_tbl = clustermole_markers()
head(markers_tbl)
```

See a summary of the available cell types:

```r
markers_tbl %>% distinct(celltype, organ, db)
```

Perform cell type enrichment for a given gene expression matrix:

```r
clustermole_enrichment(expr_mat, species)
```
