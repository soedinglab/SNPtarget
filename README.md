# SNPtarget (currently under construction)

SNPtarget is a pipeline for predicting enhancer-gene interactions based on
statistical signatures of DNase-I hypersensitivity across more than 100 tissues,
genomic distance and HiC interaction frequency.

## Repository structure

The folder *enh_scripts* contains a python package of the statistical model and
utility scripts.

The *pipelines* directory contains two pipelines for running on an LSF cluster.

* SNPpipeline (prediction of target genes of 1000 Genomes SNPs)
* FullGenomePipeline (prediction of 100bp segmentation of the full human genome)
