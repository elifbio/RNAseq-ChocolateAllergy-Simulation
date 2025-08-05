#Elif_ChocolateAllergy_RNAseq/

#Transcriptomic Signature of Chocolate Allergy: Simulated RNA-Seq

# === 1. Import Required Packages ===
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import subprocess
import sys

# Install gseapy if not available (safe for scripts)
try:
    import gseapy as gp
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "gseapy"])
    import gseapy as gp

# === 2. Simulate RNA-Seq Gene Expression Counts ===
np.random.seed(42)

genes = ['IL4', 'IL5', 'IL13', 'CCL17', 'CCL22', 'GATA3']
control_lambda = [20, 15, 25, 10, 12, 18]
allergy_lambda = [50, 45, 60, 30, 35, 55]

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

# === 3. Log2 Transformation ===
df_log = np.log2(df + 1)  # Avoid log(0)
control_samples = df_log.iloc[:, 0:3]
allergy_samples = df_log.iloc[:, 3:6]

# === 4. Differential Expression Analysis ===
p_values = []
fold_changes = []
for gene in df_log.index:
    ctrl = control_samples.loc[gene]
    alrg = allergy_samples.loc[gene]
    stat, p = ttest_ind(ctrl, alrg)
    p_values.append(p)
    fold_changes.append(alrg.mean() - ctrl.mean())

results = pd.DataFrame({
    'Gene': df_log.index,
    'log2FoldChange': fold_changes,
    'pValue': p_values
})
results['-log10(pValue)'] = -np.log10(results['pValue'])

# === 5. Volcano Plot ===
plt.figure(figsize=(8,6))
plt.scatter(results['log2FoldChange'], results['-log10(pValue)'], color='purple', alpha=0.7)
plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--')
plt.axvline(x=1, color='gray', linestyle='--')
plt.axvline(x=-1, color='gray', linestyle='--')
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot – Chocolate Allergy Simulation')
plt.grid(True)
plt.savefig("volcano_plot.png", dpi=300)
plt.show()

# === 6. Labeled Volcano Plot ===
plt.figure(figsize=(8,6))
plt.scatter(results['log2FoldChange'], results['-log10(pValue)'], color='purple', alpha=0.7)
plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--')
plt.axvline(x=1, color='gray', linestyle='--')
plt.axvline(x=-1, color='gray', linestyle='--')

for i in range(len(results)):
    fc = results.loc[i, 'log2FoldChange']
    pval = results.loc[i, '-log10(pValue)']
    gene = results.loc[i, 'Gene']
    if abs(fc) > 1 and pval > -np.log10(0.05):
        plt.text(fc, pval + 0.1, gene, fontsize=9, ha='center', color='darkred')

plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot – Chocolate Allergy Simulation')
plt.grid(True)
plt.tight_layout()
plt.savefig("volcano_plot_labeled.png", dpi=300)
plt.show()

# === 7. Multiple Testing Correction ===
results['adj_pValue'] = multipletests(results['pValue'], method='fdr_bh')[1]
significant_genes = results[results['adj_pValue'] < 0.05]
significant_genes.to_csv("significant_genes.csv")

# === 8. Heatmap ===
sns.set(style='whitegrid')
heatmap_data = df_log.loc[significant_genes['Gene']]
sns.heatmap(heatmap_data, cmap='magma', annot=True)
plt.title("Expression Heatmap – Significant Genes")
plt.tight_layout()
plt.savefig("heatmap_significant_genes.png", dpi=300)
plt.show()

# === 9. Gene Ontology Enrichment (Enrichr via gseapy) ===
gene_list = significant_genes['Gene'].tolist()

# Save significant genes to a text file for enrichment input
with open("gene_list.txt", "w") as f:
    for gene in gene_list:
        f.write(f"{gene}\n")
        
if gene_list:
    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=['GO_Biological_Process_2023'],
        organism='Human',
        outdir='enrichr_results',
        cutoff=0.05
    )

    top_terms = enr.results.sort_values('Adjusted P-value').head(10)
    plt.figure(figsize=(8, 6))
    plt.barh(top_terms['Term'], -np.log10(top_terms['Adjusted P-value']), color='skyblue')
    plt.xlabel('-log10(Adjusted P-value)')
    plt.title('Top GO Biological Processes')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()

    gp.barplot(enr.results, title='GO Biological Process Enrichment')
    plt.show()

    # === 10. KEGG Pathway Enrichment ===
    enr_kegg = gp.enrichr(
        gene_list=gene_list,
        gene_sets=['KEGG_2021_Human'],
        organism='Human',
        outdir='kegg_results',
        cutoff=0.05
    )

    top_kegg = enr_kegg.results.sort_values('Adjusted P-value').head(10)
    plt.figure(figsize=(8, 6))
    plt.barh(top_kegg['Term'], -np.log10(top_kegg['Adjusted P-value']), color='salmon')
    plt.xlabel('-log10(Adjusted P-value)')
    plt.title('Top KEGG Pathways')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()

    # === 11. Reactome Pathway Enrichment ===
    enr_reactome = gp.enrichr(
        gene_list=gene_list,
        gene_sets=['Reactome_2022'],
        organism='Human',
        outdir='reactome_results',
        cutoff=0.05
    )

    top_reactome = enr_reactome.results.sort_values('Adjusted P-value').head(10)
    plt.figure(figsize=(8, 6))
    plt.barh(top_reactome['Term'], -np.log10(top_reactome['Adjusted P-value']), color='mediumseagreen')
    plt.xlabel('-log10(Adjusted P-value)')
    plt.title('Top Reactome Pathways')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()
else:
    print("No significant genes found (adj_p < 0.05). Enrichment analyses skipped.")
