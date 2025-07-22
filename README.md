# Chocolate Allergy – RNA-seq Simulation

## Introduction

A simulated RNA-seq count matrix for genes involved in allergic response, built to explore immune activation in a chocolate allergy scenario.  
Created as part of a conceptual and ethical exploration in computational biology.

to be reviewed

## Purpose

This simulation was designed to:
- Practice RNA-seq data modeling using Python
- Explore the biological logic of allergic reactions
- Prepare for future DE analysis with real datasets

to be reviewed

## Methods Used

- Poisson distribution for count data generation
- Custom λ values for control vs allergy-exposed samples
- Pandas and NumPy for structuring the data

to be reviewed

## Files Included

- `chocolate_simulation.py` – Main Python script for data generation
- `counts.csv` – Simulated count matrix (optional)

to be reviewed

## Creative Context

Inspired by the question: "What genes wake up when chocolate triggers an immune reaction?"  
Blending scientific rigor with intuitive design, this repo reflects a personal learning journey.

# Data Transformation and Statistical Analyses

A log2 transformation was applied to normalize raw count data and reduce expression scale variance. Samples were grouped as Control and Allergy to reflect contrasting immune environments.
- T-tests were conducted per gene to detect differential expression.
- Log2 fold change values were computed between allergy and control groups.
- Multiple hypothesis testing was corrected using Benjamini-Hochberg (FDR) method.
- Genes with adjusted p-values below 0.05 were selected as statistically significant.
This step ensures reliable identification of genes that respond to the simulated allergic condition.

to be reviewed 

# Visualization

Two visualizations were used to highlight the results:

- Volcano Plot Highlights genes with significant fold changes and low p-values. Genes such as IL4, IL5, IL13 appeared prominently, indicating Th2-related activation.

- Heatmap Reveals expression intensity across samples for significant genes. Allergy samples exhibited distinct upregulation in selected markers.
These plots were generated with Python's matplotlib, seaborn, and statistical libraries.

to be reviewed

# Biological Interpretation

The expression profiles mimic a typical Th2 immune response seen in allergic reactions:
- IL4, IL5, IL13: Promote IgE class switching and eosinophil activation
- CCL17, CCL22: Direct immune cell recruitment to inflamed tissues
- GATA3: Th2 lineage transcription factor, consistent with observed gene patterns
Although based on simulated counts, gene behavior aligns with real biological pathways.

to be reviewed

# Final Insights and Future Description

This RNA-seq simulation offers a scaffold for ethical data exploration and pipeline prototyping.
- Enables DE analysis without experimental or animal-derived data
- Supports scientific storytelling through open and reproducible code
- Prepares the groundwork for future integration with real datasets
Future plans include applying the pipeline to public RNA-seq repositories and extending biological annotation with pathway and ontology analyses.

to be reviewed

## Technologies Used

- `Python` | `pandas` | `numpy`  
- `scipy.stats` | `matplotlib` | `seaborn`  
- `statsmodels` for FDR correction

---
to be reviewed

## Author

Designed and implemented by **Elif**, with a focus on ethical innovation and meaningful scientific impact.
