import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

deseq_results_file = "GSE75070_MCF7_shRUNX1_shNS_RNAseq_log2_foldchange.txt"
deseq_results = pd.read_csv(deseq_results_file, sep="\t")

annotated_peaks_file = "annotated.txt"
annotated_peaks = pd.read_csv(annotated_peaks_file, sep="\t").drop_duplicates(subset=['Gene Name'])

fold_change_cutoff = 1
padj_cutoff = 0.01

upregulated = deseq_results[(deseq_results['log2FoldChange'] > fold_change_cutoff) & (deseq_results['padj'] < padj_cutoff)]
downregulated = deseq_results[(deseq_results['log2FoldChange'] < -fold_change_cutoff) & (deseq_results['padj'] < padj_cutoff)]

upregulated_with_peaks = pd.merge(upregulated, annotated_peaks, left_on='genename', right_on='Gene Name', how='left')
downregulated_with_peaks = pd.merge(downregulated, annotated_peaks, left_on='genename', right_on='Gene Name', how='left')

distance_cutoffs = [5, 20, 100]  

def calculate_bound_genes(df, distance_cutoffs):
    bound_counts = []
    unbound_counts = [] 
    for distance in distance_cutoffs:
        bound_genes = df[(df['Distance to TSS'].notnull()) & (df['Distance to TSS'].abs() <= (distance * 1000))]
        total_unbound_genes = len(df) - len(bound_genes)
        bound_counts.append(len(bound_genes))
        unbound_counts.append(total_unbound_genes)  
    return bound_counts, unbound_counts 

upregulated_bound_counts, upregulated_unbound_counts = calculate_bound_genes(upregulated_with_peaks, distance_cutoffs)
downregulated_bound_counts, downregulated_unbound_counts = calculate_bound_genes(downregulated_with_peaks, distance_cutoffs)
upregulated_percent_bound = [(bound_count / unbound_count * 100) for bound_count, unbound_count in zip(upregulated_bound_counts, upregulated_unbound_counts)]
downregulated_percent_bound = [(bound_count / unbound_count * 100) for bound_count, unbound_count in zip(downregulated_bound_counts, downregulated_unbound_counts)]


fig, ax = plt.subplots()
upregulated_bound_bars = ax.bar(index, upregulated_percent_bound, bar_width, label='Upregulated Bound', color='green', alpha=0.5, edgecolor='black')
upregulated_unbound_bars = ax.bar(index, [100] * len(distance_cutoffs), bar_width, color='gray', alpha=0.5, edgecolor='black')
downregulated_bound_bars = ax.bar(index + bar_width, downregulated_percent_bound, bar_width, label='Downregulated Bound', color='red', alpha=0.5, edgecolor='black')
downregulated_unbound_bars = ax.bar(index + bar_width, [100] * len(distance_cutoffs), bar_width, color='gray', alpha=0.5, edgecolor='black')
for rect, value in zip(upregulated_bound_bars, upregulated_bound_counts):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height/2, f'{int(value)}', ha='center', va='center', color='black')
for rect, value in zip(downregulated_bound_bars, downregulated_bound_counts):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height/2, f'{int(value)}', ha='center', va='center', color='black')
for rect, value in zip(upregulated_unbound_bars, upregulated_unbound_counts):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height/2, f'{int(value)}', ha='center', va='center', color='black')
for rect, value in zip(downregulated_unbound_bars, downregulated_unbound_counts):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height/2, f'{int(value)}', ha='center', va='center', color='black')
ax.set_xlabel('Distance to TSS (kb)')
ax.set_ylabel('Percentage of genes (%)')
ax.set_title('Bound and unbound DE genes')
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(distance_cutoffs)
legend_labels = ['Upregulated bound genes', 'Downregulated bound genes', 'Unbound genes']
legend_colors = ['green', 'red', 'gray']
ax.legend(labels=legend_labels, handles=[plt.Rectangle((0,0),1,1, color=color) for color in legend_colors], loc='upper left', bbox_to_anchor=(1.05, 1))
ax.set_ylim([0, 100])
ax.set_yticks(np.arange(0, 101, 10))

plt.savefig('Project2_barchart.png')
plt.show()
