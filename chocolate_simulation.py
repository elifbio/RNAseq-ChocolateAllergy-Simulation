#Elif_ChocolateAllergy_RNAseq/

#Transcriptomic Signature of Chocolate Allergy: Simulated RNA-Seq

#Import required packages
import numpy as np
import pandas as pd

#Simulated gene expression counts

# Fix random seed for simulation
np.random.seed(42)

# Gene names (personal context can be added if desired)
genes = ['IL4', 'IL13', 'FCER1A', 'HDC', 'CCL2', 'GENEX']


# Gene expression for control and allergy conditions (Poisson distribution)

# Spread lambda values across 6 genes × 3 samples
lam_control = np.tile([50, 60, 30, 20, 40, 100], (3,1)).T  # shape (6,3)
lam_allergy = np.tile([90, 100, 80, 65, 70, 130], (3,1)).T  # shape (6,3)

# Generate count data using Poisson distribution
control = np.random.poisson(lam=lam_control)
allergy = np.random.poisson(lam=lam_allergy)

# Create count matrix
counts = pd.DataFrame(
    data = np.hstack([control, allergy]),
    index = genes,
    columns = ['Ctrl_1', 'Ctrl_2', 'Ctrl_3', 'Allergy_1', 'Allergy_2', 'Allergy_3']
)

counts
print(counts)

genes = ['IL4', 'IL5', 'IL13', 'CCL17', 'CCL22', 'GATA3']
control_lambda = [20, 15, 25, 10, 12, 18]
allergy_lambda = [50, 45, 60, 30, 35, 55]

np.random.seed(42)

control_counts = [np.random.poisson(lam, 3) for lam in control_lambda]
allergy_counts = [np.random.poisson(lam, 3) for lam in allergy_lambda]

df = pd.DataFrame({
    'Gene': genes,
    'Control_1': [c[0] for c in control_counts],
    'Control_2': [c[1] for c in control_counts],
    'Control_3': [c[2] for c in control_counts],
    'Allergy_1': [a[0] for a in allergy_counts],
    'Allergy_2': [a[1] for a in allergy_counts],
    'Allergy_3': [a[2] for a in allergy_counts]
})
df.set_index('Gene', inplace=True)

df.to_csv("counts.csv")

#Data transformation and statistical analyses

#Import required libraries
import pandas as pd
import numpy as np

#Load simulated data from .csv file
df = pd.read_csv("counts.csv", index_col=0)
print(df.head())  # Display the first few rows

#Normalize RNA-seq data using log2 transformation
df_log = np.log2(df + 1)  # +1 to avoid log(0) issues

#Split samples: first 3 columns are control, last 3 are allergy
control_samples = df_log.iloc[:, 0:3]  # Control_1, Control_2, Control_3
allergy_samples = df_log.iloc[:, 3:6]  # Allergy_1, Allergy_2, Allergy_3

#Apply t-test and compute fold change for each gene
from scipy.stats import ttest_ind

p_values = []
fold_changes = []

for gene in df_log.index:
    control = control_samples.loc[gene]
    allergy = allergy_samples.loc[gene]
    
    stat, p = ttest_ind(control, allergy)
    p_values.append(p)
    
    fc = allergy.mean() - control.mean()
    fold_changes.append(fc)
    
#Create results DataFrame
results = pd.DataFrame({
    'Gene': df_log.index,
    'log2FoldChange': fold_changes,
    'pValue': p_values
})

#Transform p-values for visualization
results['-log10(pValue)'] = -np.log10(results['pValue'])
print(results)

#Draw volcano plot – fold change vs p-value
import matplotlib.pyplot as plt

plt.figure(figsize=(8,6))
plt.scatter(results['log2FoldChange'], results['-log10(pValue)'], color='purple', alpha=0.7)

#Add threshold lines: statistical and biological significance
plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--')  # p < 0.05 threshold
plt.axvline(x=1, color='gray', linestyle='--')                # FC > 2
plt.axvline(x=-1, color='gray', linestyle='--')               # FC < 0.5

#Axis labels and title
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot – Chocolate Allergy Simulation')
plt.grid(True)
plt.show()

#Save volcano plot
plt.savefig("volcano_plot.png", dpi=300)

import matplotlib.pyplot as plt

#Volcano plot – labeled version (with significant gene names)
plt.figure(figsize=(8,6))
plt.scatter(results['log2FoldChange'], results['-log10(pValue)'], color='purple', alpha=0.7)

# Add threshold lines
plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--')  # p < 0.05
plt.axvline(x=1, color='gray', linestyle='--')                # FC > 2
plt.axvline(x=-1, color='gray', linestyle='--')               # FC < 0.5

# Label significant genes
for i in range(len(results)):
    fc = results.loc[i, 'log2FoldChange']
    pval = results.loc[i, '-log10(pValue)']
    gene = results.loc[i, 'Gene']
    
    if abs(fc) > 1 and pval > -np.log10(0.05):  # Biologically and statistically significant
        plt.text(fc, pval + 0.1, gene, fontsize=9, ha='center', color='darkred')

#Axis labels and title
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot – Chocolate Allergy Simulation')
plt.grid(True)
plt.tight_layout()
plt.show()

#Save labeled volcano plot
plt.savefig("volcano_plot_labeled.png", dpi=300)

#Why p-value correction is needed?
#A separate t-test is done for each gene → e.g., 1000 genes → 1000 p-values
#Some genes may appear significant just by chance
#Correction controls false discovery rate (FDR)
#Most common method: Benjamini-Hochberg (FDR-BH)

#Multiple testing correction – Benjamini-Hochberg method
from statsmodels.stats.multitest import multipletests

#Apply correction and add new column
results['adj_pValue'] = multipletests(results['pValue'], method='fdr_bh')[1]

#Filter significant genes (adj_p < 0.05)
significant_genes = results[results['adj_pValue'] < 0.05]
print(significant_genes[['Gene', 'log2FoldChange', 'adj_pValue']])

#Import packages for heatmap
import seaborn as sns
import matplotlib.pyplot as plt

#Select log2 expression values for significant genes
selected_genes = results[results['adj_pValue'] < 0.05]['Gene']
heatmap_data = df_log.loc[selected_genes]

#Draw heatmap
sns.heatmap(heatmap_data, cmap='magma', annot=True)
plt.title("Expression Heatmap – Significant Genes")
plt.tight_layout()
plt.show()

#Save heatmap visualization
plt.savefig("heatmap_significant_genes.png", dpi=300)








#Final insights and future description
