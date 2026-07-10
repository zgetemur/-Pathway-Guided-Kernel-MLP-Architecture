# Pathway Guided Kernel MLP for Cross Disease Transcriptomic Analysis

A biologically informed **Kernel-MLP framework** for transcriptomic modeling and cross-disease transfer learning.

This repository provides a modular and reproducible implementation of a pathway-guided machine learning pipeline that integrates **prior biological knowledge** with **non-linear representation learning** to improve interpretability, pathway selection, and disease transferability.

---

## Project Overview

The framework was developed to investigate whether molecular representations learned from inflammatory disease can generalize to hematologic malignancy while preserving biological specificity.

### Study Objectives

- Analyze large-scale **microarray gene expression data**
- Integrate **KEGG** and **PID** pathway information into neural architectures
- Evaluate robustness through repeated computational experiments
- Investigate molecular transfer from **Rheumatoid Arthritis (RA)** to **Acute Myeloid Leukemia (AML)**
- Validate biological specificity using **Prostate Cancer (PCa)** as a negative control

---

# Workflow

## A. Data Collection

Public transcriptomic datasets were obtained from the **NCBI Gene Expression Omnibus (GEO)**.

### Included Cohorts

| Domain | Description |
|--------|-------------|
| RA | 8 independent cohorts (source domain) |
| AML | Target transfer domain |
| PCa | Specificity control domain |
| HE | Healthy controls |

---

## B. Preprocessing & Harmonization

### Expression Processing

- Raw CEL extraction
- Robust Multi-array Average (**RMA**) normalization
- Log2-scale expression summarization

### Feature Alignment

- Probe-to-gene annotation
- Cross-platform harmonization

### Batch Effect Correction

Technical variation across datasets was mitigated using:

- **ComBat**
- **neuroCombat**

Batch labels were defined using GEO series identifiers (GSE).

Validation included:

- PCA visualization
- t-SNE projection
- Machine learning-based batch verification

---

## C. Pathway-Guided Modeling

Gene expression profiles were decomposed into biologically constrained modules before classification.

### Embedded Pathway Knowledge

| Database | Description | Count |
|---------|-------------|------:|
| KEGG | Functional pathway organization | 320 |
| PID | Regulatory and signaling interactions | 209 |

### Modeling Pipeline

```text
Gene Expression
      ↓
Pathway Factorization
      ↓
Pathway-specific RBF Kernels
      ↓
Kernel-MLP
      ↓
Classification
```

### Benchmark Models

Performance was compared against:

- Random Forest
- Support Vector Machine (SVM)

---

## D. Stability Analysis & Transfer Learning

### Stability Assessment

To evaluate reproducibility:

- 80 independent repeated experiments
- Frequency-based pathway ranking
- Binomial statistical filtering (`p < 0.05`)

### Transfer Strategy — Freeze & Adapt

| Component | Strategy |
|----------|----------|
| Global pathway coefficients (η) | Frozen |
| Pathway-specific weights (w) | Trainable |

This mechanism preserves transferable biological structure while adapting to target-specific transcriptomic variation.

---

# Methodological Highlights

### Pathway-Level Factorization

Reduces transcriptomic dimensionality by organizing genes into biologically meaningful modules.

### Non-linear Representation Learning

Captures complex molecular interactions through pathway-specific Gaussian kernel branches.

### Sparse Biological Selection

Applies non-negative L1 regularization for compact and interpretable pathway selection.

### Transfer Without Catastrophic Forgetting

Separates shared biological representations from target-domain adaptation.

### Specificity Validation

Demonstrates that observed transfer behavior reflects biologically meaningful inflammatory–myeloid relationships rather than generic cross-domain effects.

Transfer performance was expected for **AML** but intentionally failed for **PCa**, supporting domain specificity.

---

# Repository Structure

```text
data/
├── raw/
├── processed/

preprocessing/
├── normalization/
├── batch_correction/

model/
├── kernel_mlp/
├── baselines/

experiments/
├── stability/
├── transfer_learning/

results/
```

---

# Reproducibility

This repository is designed to support:

- Modular experimentation
- Transparent biological interpretation
- Reproducible transcriptomic analysis
- Cross-disease generalization studies

---

# Citation

If you use this repository in academic work, please cite the associated publication.
