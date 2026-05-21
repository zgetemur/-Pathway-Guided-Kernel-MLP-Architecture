# Pathway-Guided Kernel-MLP for Cross-Disease Transcriptomic Analysis

A pathway-informed **Kernel-MLP framework** for transcriptomic analysis and cross-disease transfer learning.

This project integrates **prior biological pathway knowledge** with **non-linear machine learning** to improve interpretability, pathway-level feature selection, and generalization across diseases.

---

## Overview

This study presents a biologically guided computational framework designed for:

- High-throughput **microarray gene expression analysis**
- Integration of **KEGG** and **PID** pathway knowledge into neural architectures
- Stability and robustness evaluation through repeated computational experiments
- Investigation of disease-specific molecular relationships using **cross-disease transfer learning**
- Biological specificity validation using an independent solid tumor control

### Transfer Learning Design

**Source Domain:** Rheumatoid Arthritis (RA)  
↓  
**Target Domain:** Acute Myeloid Leukemia (AML)  
↓  
**Specificity Control:** Prostate Cancer (PCa)

---

# Computational Pipeline

## Stage A — Data Collection

### Transcriptomic Data Acquisition

- Public datasets retrieved from **NCBI Gene Expression Omnibus (GEO)**
- Multi-cohort integration:
  - **8 independent RA cohorts**
  - Validation datasets for **AML**
  - Control datasets for **PCa**

---

## Stage B — Preprocessing & Harmonization

### Expression Processing

- Raw **CEL** extraction
- **Robust Multi-array Average (RMA)** normalization
- Transformation into **log2-scale expression values**

### Feature Harmonization

- Probe → Gene annotation
- Cross-platform feature alignment

### Batch Effect Correction

- Empirical Bayes adjustment using **ComBat / neuroCombat**
- Batch labels defined using **GSE identifiers**
- Validation performed via:
  - PCA
  - t-SNE
  - ML-based batch signature verification

---

## Stage C — Pathway-Guided Modeling

### Biological Prior Integration

Expression vectors are decomposed into biologically constrained modules:

| Database | Description | Count |
|----------|-------------|------:|
| KEGG | Global functional pathway maps | 320 |
| PID | Fine-grained signaling networks | 209 |

### Modeling Architecture

1. Pathway-level expression factorization  
2. Construction of pathway-specific **Gaussian (RBF) kernels**  
3. End-to-end **Kernel-MLP training**  
4. Comparative benchmarking against:

- Random Forest
- Support Vector Machine (SVM)

---

## Stage D — Stability & Constrained Transfer Learning

### Stability Analysis

- **80 independent repeated runs**
- Frequency-based pathway selection
- Statistical filtering:

\[
p < 0.05
\]

### Transfer Strategy — *Freeze & Adapt*

Global pathway coefficients remain fixed:

\[
\eta = \text{Frozen}
\]

Pathway-specific weights remain adaptable:

\[
w = \text{Trainable}
\]

This preserves domain-invariant biological structure while adapting to target-specific molecular shifts.

---

# Data Processing Notes

## Normalization

Raw CEL files were summarized using **RMA** to reduce technical variability and generate normalized transcriptomic profiles.

## Batch Effect Mitigation

Batch variables were explicitly modeled using cohort identifiers and corrected via empirical Bayes harmonization.

## Phenotype Labels

| Label | Description |
|-------|-------------|
| RA | Rheumatoid Arthritis (Source Domain) |
| AML | Acute Myeloid Leukemia (Target Domain) |
| PCa | Prostate Cancer (Negative Control) |
| HE | Healthy Controls |

---

# Embedded Biological Knowledge

## KEGG

**Kyoto Encyclopedia of Genes and Genomes**

Captures large-scale functional organization and pathway topology.

**Total pathways:** `320`

---

## PID

**Pathway Interaction Database**

Provides detailed signaling interactions and regulatory cascades.

**Total pathways:** `209`

---

# Key Methodological Contributions

### Pathway-Level Factorization
Reduces dimensional complexity by transforming gene-level inputs into biologically meaningful modules.

### Non-linear Representation Learning
Models complex intra-pathway interactions through Gaussian kernel branches.

### Sparse Biological Selection
Applies non-negative **L1 regularization** to enforce compact pathway selection.

### Catastrophic Forgetting Prevention
Maintains transferable biological structures through constrained adaptation.

### Specificity Validation
Demonstrates that inflammatory → myeloid transfer is biologically meaningful by confirming transfer failure on an unrelated solid tumor domain (**PCa**).

---

# Repository Goals

✔ Reproducible computational pipeline  
✔ Modular architecture  
✔ Interpretable pathway selection  
✔ Robust cross-disease generalization  
✔ Biological validation and specificity analysis

---

## Citation

If you use this repository in your research, please cite the associated study.

```
