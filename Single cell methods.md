Epithelial-Myoepithelial dataset:
Initially, sequencing data was demultiplexed using the bcl-convert tool from llumina and processed using the Cell Ranger (Zheng2017) v.7.1.0 (10X Genomics Pleasanton, CA) cellranger multi function, which was run with the 'Include introns=False' option. The transcriptome reference was mm10-2020-A. Downstream analyses were performed with Seurat v.4.4.0 (Hao2021) with the parameters described below. CreateSeuratObject was run with min.cells=3. Low quality cells were filtered with the following thresholds for gene counts, UMI counts, and percentage mitochondria:  500 < nFeature_RNA < 6000, nCount_RNA < 30000, percent.mt <7. The RNA count data was normalized with the LogNormalize method, while the HTO count data was normalized with the CLR method. HTOs were demultiplexed with the HTODemux function with the parameter positive.quantile = 0.99. RunPCA, RunUMAP, and FindNeighbors were run with the top 28 dimensions, and a chosen resolution of 0.6, and the resulting clusters were manually annotated using cell type specific markers.

[Zheng2017] Zheng GX, Terry JM, Belgrader P, Ryvkin P, Bent ZW, Wilson R, Ziraldo SB, Wheeler TD, McDermott GP, Zhu J, Gregory MT. Massively parallel digital transcriptional profiling of single cells. Nature communications. 2017 Jan 16;8(1):1-2.

[Hao2021] Hao Y, Hao S, Andersen-Nissen E, Mauck III WM, Zheng S, Butler A, Lee MJ, Wilk AJ, 
Darby C, Zager M, Hoffman P. Integrated analysis of multimodal single-cell data. Cell. 2021 Jun 24;184(13):3573-87.

