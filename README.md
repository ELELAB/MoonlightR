Cancer Structural Biology Group, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark 

Cancer Systems Biology, Section of Bioinformatics, Department of Health and Technology, Technical University of Denmark, 2800, Lyngby, Copenhagen


#MoonlightR 


Repository associated to the publication:

Interpreting pathways to discover cancer driver genes with Moonlight.
Colaprico A, Olsen C, Bailey MH, Odom GJ, Terkelsen T, Silva TC, Olsen AV, Cantini L, Zinovyev A, Barillot E, Noushmehr H, Bertoli G, Castiglioni I, Cava C, Bontempi G, Chen XS, Papaleo E.
Nat Commun. 2020 Jan 3;11(1):69. doi: 10.1038/s41467-019-13803-0., PMID: 31900418

contacts for repository: Elena Papaleo, elpap-at-dtu.dk, elenap-at-cancer.dk; Matteo Tiberti: tiberti-at-cancer.dk


# Identify oncogenes and tumor suppressor genes from genomic (gene expression) data.

### Installation from Bioconductor ###
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("MoonlightR")
```

### Installation from GitHub ###
```R
devtools::install_github(repo = "ibsquare/MoonlightR")
```
