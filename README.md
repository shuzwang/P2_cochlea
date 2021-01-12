# P2_cochlear
It contains the code to analyze P2 scRNA-seq and scATAC-seq data. 


Abstract
============================

The function of auditory hair cells is to transduce sound to the brain and in mammals these cells reside together with supporting cells in the sensory epithelium of the cochlea, called the organ of Corti. To establish the organ’s delicate function during development and differentiation, spatiotemporal gene expression is strictly controlled by a combination of cell type specific transcription factors representing the regulatory landscape. Investigations on the interplay between transcription factors and chromatin accessibility, were obscured to identify cellular heterogeneity by bulk-sequencing technology . To study the formation of the regulatory landscape in hair cells, we collected single-cell chromatin accessibility profiles accompanied by single-cell RNA data from genetically labeled murine hair cells and supporting cells after birth. Using an integrative approach, we predicted cell type specific activating and repressing functions of developmental transcription factors. Furthermore, by integrating gene expression and chromatin accessibility datasets, we reconstructed gene regulatory networks. Using a comparative approach, 20 hair cell specific activators and repressors, including putative downstream targets genes, were identified. Clustering of target genes resolved groups of related transcription factors and was utilized to infer their developmental functions. Finally, the heterogeneity in the single cell data allowed us to spatially reconstruct transcriptional as well as chromatin accessibility trajectories, indicating that gradual changes in the epigenetic landscape were lagging behind the transcriptional identity of hair cells along the organ’s longitudinal axis. Overall, this study provides a strategy to spatially reconstruct the formation of a lineage specific regulatory landscape using a single-cell multi-omics approach.


Analysis Pipelines
============================

Here we include analysis pipeline code for our manuscript. 
* Seurat pipeline: Seurat pipeline is to analyze P2 scRNA-seq dataset. We identified clusters, annotated cell type for each cluster using marker genes, determined differentially expressed genes across cell types, and identified differentially expressed genes between library IDs (Apical VS Basal cells). We followed the documentations from https://satijalab.org/seurat/
* SnapATAC pipeline: SnapATAC pipeline is to analyze P2 scATAC-seq dataset. We started with Snaptools, which is a python module for pre-processing and working with snap files. We followed the documentation from https://github.com/r3fang/SnapTools. ::
  $ snaptools snap-pre --input-file=/P2_ATAC_ApexBase_aggr_2872/outs/fragments_sort.bed.gz --output-snap=atac_peak_cell_20000.snap --genome-name=mm10 --genome-size=mm10.chrom.sizes --min-mapq=30 --min-flen=50 --max-flen=1000 --keep-chrm=TRUE --keep-single=FALSE --keep-secondary=False --overwrite=True --max-num=20000 --min-cov=100 --verbose=True
  $
