# Pathway-Guided Kernel-MLP for Cross-Disease Transcriptomic Analysis

This repository implements a pathway-guided Kernel-MLP framework for the analysis of gene expression data and cross-disease transfer learning. The approach integrates biological pathway information with machine learning to enable interpretable modeling at the pathway level.

## Project Scope

The study focuses on:

* Analysis of microarray gene expression data
* Integration of KEGG pathway knowledge into machine learning models
* Stability analysis through repeated experiments
* Evaluation of cross-disease transfer learning behavior

The repository provides a structured and reproducible implementation of the proposed methodology.

## Computational Pipeline

The workflow consists of the following stages:

### Data Collection

* Public datasets retrieved from GEO
* Multiple cohorts combined for analysis

### Preprocessing

* RMA normalization performed in R
* Probe-to-gene mapping
* Intersection of common genes across datasets
* Batch effect correction using ComBat

### Pathway-Guided Modeling

* Integration of KEGG pathways
* Construction of Gaussian kernels per pathway
* Training of Kernel-MLP model
* Baseline comparisons using Random Forest and SVM

### Stability and Transfer Analysis

* Repeated training runs
* Pathway selection frequency evaluation
* Cross-disease generalization analysis

## Data Processing Notes

### Normalization

Raw data were processed using RMA normalization 

### Batch Effect Correction

Batch effects were corrected using ComBat, with dataset identifiers (GSE IDs) used as batch variables.

### Phenotype Labels

Phenotype labels were manually curated and standardized prior to analysis.
Label consistency is critical and must be ensured before running the pipeline.

Example labels:

* RA
* HE

## Pathway Data

Pathway definitions are based on the KEGG_2021_Human collection.



## Execution Strategy

Due to computational constraints, experiments are designed to run with a single L1 regularization value per execution.

Multiple configurations should be executed separately.



## Methodological Features

* Pathway-level representation instead of gene-level features
* Gaussian kernel-based modeling
* Sparse and interpretable learning via L1 regularization
* Stability analysis through repeated runs
* Integration of biological knowledge into machine learning

## Reproducibility

To reproduce the experiments:

1. Apply identical preprocessing steps
2. Ensure correct label standardization
3. Use consistent random seeds
4. Execute repeated runs

