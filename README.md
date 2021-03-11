# P2_cochlear
It contains the analytical code to analyze P2 scRNA-seq and scATAC-seq multi-omics datasets. 


Abstract
============================

The function of auditory hair cells is to transduce sound to the brain and in mammals these cells reside together with supporting cells in the sensory epithelium of the cochlea, called the organ of Corti. To establish the organ’s delicate function during development and differentiation, spatiotemporal gene expression is strictly controlled by a combination of cell type specific transcription factors representing the regulatory landscape. Investigations on the interplay between transcription factors and chromatin accessibility, were obscured to identify cellular heterogeneity by bulk-sequencing technology . To study the formation of the regulatory landscape in hair cells, we collected single-cell chromatin accessibility profiles accompanied by single-cell RNA data from genetically labeled murine hair cells and supporting cells after birth. Using an integrative approach, we predicted cell type specific activating and repressing functions of developmental transcription factors. Furthermore, by integrating gene expression and chromatin accessibility datasets, we reconstructed gene regulatory networks. Using a comparative approach, 20 hair cell specific activators and repressors, including putative downstream targets genes, were identified. Clustering of target genes resolved groups of related transcription factors and was utilized to infer their developmental functions. Finally, the heterogeneity in the single cell data allowed us to spatially reconstruct transcriptional as well as chromatin accessibility trajectories, indicating that gradual changes in the epigenetic landscape were lagging behind the transcriptional identity of hair cells along the organ’s longitudinal axis. Overall, this study provides a strategy to spatially reconstruct the formation of a lineage specific regulatory landscape using a single-cell multi-omics approach.


Analysis Pipelines
============================

Here we include analysis pipeline code for our manuscript. 
* Seurat pipeline: Seurat pipeline is to analyze P2 scRNA-seq dataset. We identified clusters, annotated cell type for each cluster using marker genes, determined differentially expressed genes across cell types, and identified differentially expressed genes between library IDs (Apical VS Basal cells). We followed the documentations from https://satijalab.org/seurat/
* SnapATAC pipeline: SnapATAC pipeline is to analyze P2 scATAC-seq dataset. We started with Snaptools, which is a python module for pre-processing and working with snap files. We followed the documentation from https://github.com/r3fang/SnapTools. 

```
snaptools snap-pre --input-file=/P2_ATAC_ApexBase_aggr_2872/outs/fragments_sort.bed.gz --output-snap=atac_peak_cell_20000.snap 
--genome-name=mm10 --genome-size=mm10.chrom.sizes --min-mapq=30 --min-flen=50 
--max-flen=1000 --keep-chrm=TRUE --keep-single=FALSE --keep-secondary=False --overwrite=True 
--max-num=20000 --min-cov=100 --verbose=True 
```

```
snaptools snap-add-bmat --snap-file=atac_peak_cell_20000.snap --bin-size-lis 5000 --verbose=True
```
Then we followed the SnapATAC pipeline https://github.com/r3fang/SnapATAC. We identified the scATAC-seq clusters, annotated the clusters based on Jaccard similarity matrix, created gene-based annotation, called peaks for each cluster, and identified differentially accessible regions. 

* LIGER pipeline: LIGER pipeline is to integrate and align scRNA-seq and scATAC-seq datasets. We co-embedded these two datasets onto the same UMAP and validate our Jaccard similarity-based annotation approach. We followed the documentation of LIGER from http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html. 

* ataqv pipeline: ataqv pipeline is to measure and compare ATAC-seq results. We conducted peak calling in SnapATAC for each individual cluster and extracted corresponding bam files to run ataqv. We followed the pipeline from https://github.com/ParkerLab/ataqv. For each cluster, we run `ataqv --peak-file` :

```
ataqv --peak-file /atacqv_results/P2_scATACseq_cluster1_peaks.narrowPeak --name cluster1 --metrics-file /atacqv_results/cluster1.ataqv.json.gz 
--excluded-region-file /atacqv_results/mm10.blacklist.bed.gz --tss-file /atacqv_results/ref_genes_R_package_4_columns.bed 
--tss-extension 2000bp --ignore-read-groups mouse /bam_files/P2_scATACseq_cluster1.bam > /atacqv_results/cluster1.ataqv.out
```

After we run `ataqv --peak-file` for each cluster, we run:

```
mkarv 6_sample cluster1.ataqv.json.gz cluster2.ataqv.json.gz cluster3.ataqv.json.gz cluster4.ataqv.json.gz cluster5.ataqv.json.gz cluster6.ataqv.json.gz
```

* chromVAR pipeline: chromVAR pipeline is to analyze the scATAC-seq dataset. chromVAR takes as inputs aligned fragments from scATAC-seq as well as genomic annotations such as motif positions. Here we used HOCOMOCOv10 motif database. chromVAR computes for each annotation and each cell a bias corrected “deviation” in accessibility from the expected accessibility based on the average of all the cells or samples. Finally, we calculated the z-scores for each TF-motif per cell to represent the TF activity for that cell. We followed the pipeline from https://github.com/GreenleafLab/chromVAR. 

* HINT-ATAC pipeline: HINT-ATAC pipeline is to identify the active transcription factor binding sites using ATAC-seq datasets. HINT-ATAC software is designed for "bulk" datasets. We aggregated the reads from the same cell type and created "pseudo-bulk" data for each cluster to run HINT-ATAC. Here we used HINT-ATAC to compare changes in the activity of transcription factors with differential footprinting between hair cells (cluster 1) and supporting cells (cluster 2) in scATAC-seq data. To make our analysis consistent, we used HOCOMOCOv10 motif database as well for footprinting analysis. We followed the documentations from https://www.regulatory-genomics.org/hint/tutorial/. 

```
python setupGenomicData.py --mm10
```
We created HOCOMOCOv10 motif PWM for further analysis using `python createPwm.py` function:
```
python createPwm.py -f jaspar-2016 -i hocomoco_v10_mouse_motif.txt  -o hocomoco_v10_mouse_motif
```

First, we called footprints `rgt-hint footprinting` for hair cell cluster and supporting cell cluster separately:
```
rgt-hint footprinting --atac-seq --paired-end --organism=mm10 ~/Desktop/P2_analysis/P2_bam_files/P2_scATACseq_cluster1.bam ~/Desktop/P2_analysis/P2_scATACseq_cluster_peaks/P2_scATACseq_cluster1_peaks.narrowPeak
rgt-hint footprinting --atac-seq --paired-end --organism=mm10 ~/Desktop/P2_analysis/P2_bam_files/P2_scATACseq_cluster2.bam ~/Desktop/P2_analysis/P2_scATACseq_cluster_peaks/P2_scATACseq_cluster2_peaks.narrowPeak
```

Secondly, we conducted motif matching `rgt-motifanalysis matching` for footprints for cluster 1 and cluster 2, respectively, using HOCOMOCOv10 motif database:
```
rgt-motifanalysis matching --motif-dbs ~/rgtdata/motifs/hocomoco_v10_mouse_motif/ --organism=mm10 --input-files P2_cluster1_footprints.bed
rgt-motifanalysis matching --motif-dbs ~/rgtdata/motifs/hocomoco_v10_mouse_motif/ --organism=mm10 --input-files P2_cluster2_footprints.bed
```

Finally, we used HINT-ATAC to generate average chromatin accessibility profiles `rgt-hint differential ` around binding sites of particular TF. Additionally, by comparing the cleavage profiles from two scATAC-seq clusters, we can get insights on changes in binding in two cell types.
```
rgt-hint differential --organism=mm10 --bc --nc 2 --mpbs-files=P2_cluster1_hocomoco_v10_match/P2_cluster1_footprints_mpbs.bed,P2_cluster2_hocomoco_v10_match/P2_cluster2_footprints_mpbs.bed 
--reads-files=P2_scATACseq_cluster1.bam,P2_scATACseq_cluster2.bam --conditions=cluster1,cluster2
```

* CellTrails pipeline: CellTrails pipeline is to infer trajectory using TF-by-cell matrix we calculated from chromVAR. Here we concentrated on hair cells which compose of inner hair cells (IHC) and outer hair cells (OHC). IHCs and OHCs are differentiated from same progenitor cells. We inferred scATAC-seq hair cell trajectory as "Y-shape" using CellTrails. We followed the documentation from https://hellerlab.stanford.edu/celltrails/. 

* Slingshot pipeline: Slingshot pipeline is to reconstruct trajectory of hair cells using TF-by-cell matrix we calculated from chromVAR. We concentrated on hair cells which compose of inner hair cells (IHC) and outer hair cells (OHC). IHCs and OHCs are differentiated from same progenitor cells. We inferred scATAC-seq hair cell trajectory as "Y-shape" using Slingshot. The asymmetric developmental result from Slingshot is consistent with our result using CellTrails algorithm. Applying a second trajectory reconstruction algorithm rules out the computational variations. We followed the documentation from https://www.bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html. 

* Monocle pipeline: Monocle pipeline is to reconstruct trajectory of hair cells using TF-by-cell matrix we calculated from chromVAR. The goal of `Monocle pipeline` is the same as `Slingshot pipeline`. We concentrated on hair cells which compose of inner hair cells (IHC) and outer hair cells (OHC). IHCs and OHCs are differentiated from same progenitor cells. We inferred scATAC-seq hair cell trajectory as "Y-shape" using Monocle. The asymmetric developmental result from Monocle is consistent with our result using CellTrails and Slingshot algorithms. Applying a different trajectory reconstruction algorithm rules out the computational variations. We followed the documentation from http://cole-trapnell-lab.github.io/monocle-release/docs/.  

* diffTF pipeline: diffTF pipeline is to classify TFs into activators and repressors based on TF activities. diffTF is designed for "bulk" RNA-seq and ATAC-seq datasets. The RNA-seq and ATAC-seq samples have to match. To apply diffTF on our single-cell multi-omics data, we ranked the cells for scRAN-seq and scATAC-seq based on 1D spatial map, respectively. Then, we separated the cells into 4 partitions based on the rank and aggregated the cells from the same partition. For scRNA-seq, we aggregated the reads for each partition individually. For scATAC-seq, we created 4 "pseudo-bulk" ATAC-seq datasets and called peaks for each "pseudo-bulk" ATAC-seq dataset. We followed diffTF pipeline from https://difftf.readthedocs.io/en/latest/. 

* Cicero pipeline: Cicero pipeline is to predict cis-regulatory interactions in the genome by examining co-accessibility using scATAC-seq data. We run Cicero for hair cell cluster and supporting cell cluster, separately. We followed the Cicero pipeline from https://cole-trapnell-lab.github.io/cicero-release/docs/.  Additionally, we integrated gene regulatory information into Cicero results by running FIMO http://meme-suite.org/doc/fimo.html. Also, we run FIMO for TFs of interest for hair cells and supporting cells, respectively. To make our analysis consistent, we used HOCOMOCOv10 motif database. Take `S100b` gene as an example.  

```
bedtools getfasta -fi /nfs/turbo/umms-joergwal/shuzwang/mm10/refdata-cellranger-atac-mm10-1.1.0/fasta/genome.fa -bed S100b_hc_peaks_up50kb_down50kb.bed -fo S100b_hc_peaks_up50kb_down50kb.fa
bedtools getfasta -fi /nfs/turbo/umms-joergwal/shuzwang/mm10/refdata-cellranger-atac-mm10-1.1.0/fasta/genome.fa -bed S100b_sc_peaks_up50kb_down50kb.bed -fo S100b_sc_peaks_up50kb_down50kb.fa
```

```
fimo --text --parse-genomic-coord HOCOMOCOv10_MOUSE_mono_meme_format.meme S100b_hc_peaks_up50kb_down50kb.fa | gzip > S100b_hc_peaks_up50kb_down50kb_fimo.txt.gz
fimo --text --parse-genomic-coord HOCOMOCOv10_MOUSE_mono_meme_format.meme S100b_sc_peaks_up50kb_down50kb.fa | gzip > S100b_sc_peaks_up50kb_down50kb_fimo.txt.gz
```

Additional Analysis
============================

Besides the existing packages we used, we also generated some novel functions for our manuscript. 

* Jaccard_matrix.R : `Jaccard_matrix.R` script includes the functions to annotate differentially accessible regions to genes and generated Jaccard similarity matrix between differentially expressed genes from scRNA-seq data and annotated differentially accessible regions from scATAC-seq data.

* TF_mode_classification.R : `TF_mode_classification.R` script includes the functions to identify differentially expressed TF genes and identify differentially accessible TF motifs between two populations (HC VS SC). Such functions provide three different statistical tests for users to choose like median value comparison, Student t-test, and Wilcoxon rank sum test to compare two populations. Also, the script includes the function to classify the TF mode into activators and repressors based on relations between gene expression and chromatin accessiblity. 

* GRN.R : `GRN.R` script contains a 3-step pipeline, inspired by SCENIC algorithm https://github.com/aertslab/SCENIC, to reconstruct gene regulatory network (GRN). The first setp is to identify gene co-expression using GENIE3 https://bioconductor.org/packages/release/bioc/vignettes/GENIE3/inst/doc/GENIE3.html. We integrated scATAC-seq putative binding sites using FIMO within a certin window (e.g. 50kb, 100kb) to identify direct downstream target genes of TFs. Once we removed indirect target genes, we determined regulons which represent groups of genes that are regulated as units. The third step is to analyze gene set enrichment in single-cell RNA-seq using AUCell https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html. 

* 1D_spatial_reconstruction.R : `1D_spatial_reconstruction.R` script to reconstruct 1D spatial map using either scRNA-seq or scTAC-seq data. The first step to reconstruct 1D spatial map is to identify centroids of apical and basal cells and to rotate the cells with the apical cells on the top and basal cells on the bottom. The second step is to project the cells along apex to base by identifying the differentially expressed genes/differentially accessibile regions between apical and basal cells. Using the spatial-specific genes as features, we assumed the PC1 would represent the apex-to-base axis and we used the PC1 rank to order the cells. Finally, we made curtain plots (1D PCA plot) to resolve cochlea tonotopic structure. The x-axis for the 1D PCA plot is random jitter for better visualization.

* Utilities.R : `Utilities.R` script contains utility functions. It includes `overlapped_genes` function which calculates the number of overlaps and returns the overlapped genes. `smooth_fragment_length` function visualizes the QC results of scATAC-seq from `ataqv` with smoothing. `atac_single_curtain_plot` function makes a 1D PCA plot using scATAC-seq data. `make_dot_plot` function makes dotplot for differentially expressed genes in scRNA-seq data. `chipseq_annotation_sox2` function annotates ChIP-seq peak file to the nearest genes based on the user-defined window size.


Supplementary Files
============================

* HOCOMOCO_v10_gene_translation_table.txt: `HOCOMOCO_v10_gene_translation_table.txt` is a tab-delimited file inclduing HOCOMOCO ID on the first column and corresponding gene name on the second column.

* mm10.blacklist.bed.gz: `mm10.blacklist.bed.gz` is a zipped tab-delimited file including the defined ENCODE blacklist in mouse (mm10). The removal of the ENCODE blacklist is an essential quality measure when analyzing functional genomics data. `mm10.blacklist.bed.gz` file was downloaded from https://github.com/Boyle-Lab/Blacklist.

* mm10.chrom.sizes.txt: `mm10.chrom.sizes.txt` is a tab-delimited file inclduing the size of chromosomes in mm10. The file was downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/. 






