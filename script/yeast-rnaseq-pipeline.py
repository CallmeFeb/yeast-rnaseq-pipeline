# %%
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import gzip
import gseapy as gp

# %%
table = pd.read_csv('../data/processed/combined_matrix.tsv', sep='\t', index_col='Name' )
table.columns = table.columns.str.replace("_quant", "", regex=False)
table = table.round().astype(int)
Meta = pd.DataFrame({
    'sample': ['ERR458493', 'ERR458494', 'ERR458495', 'ERR458500', 'ERR458501', 'ERR458502'],
    'group': ['Normal','Normal','Normal','Mutated','Mutated','Mutated'],
})
Meta = Meta.set_index("sample")
table_t= table.transpose()
print(all(table_t.index == Meta.index))
#print(table_t.head())
#print(Meta.head())

# %%
dds = DeseqDataSet(
    counts = table_t,
    metadata = Meta,
    design = "~group",
    refit_cooks= True
)

dds.deseq2()

# %%
stat_res = DeseqStats(dds, contrast=['group', 'Mutated', 'Normal'] )
stat_res.summary()

result_df= stat_res.results_df
result_df= result_df.dropna(subset=['padj'])


clean_results =result_df[(result_df['padj']<0.05) & (result_df['log2FoldChange'].abs()>1)]
sorted_results= clean_results.sort_values('padj')
sorted_results.head(20)


# %%
result_df.loc[result_df.index.str.contains("YOR290C")]
table.loc[["YOR290C_mRNA"]]


# %%
log_counts = np.log2(table_t +1 )
scaled = StandardScaler().fit_transform(log_counts)

pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled)

print(pca_result)
print(pca.explained_variance_ratio_)

# %%
pca_dif = pd.DataFrame(pca_result, columns=["PC1", "PC2"], index=table_t.index)
pca_dif['group']= Meta['group']

for group, color in [('Normal', 'blue'), ('Mutated', 'red')]:
    subset= pca_dif[pca_dif['group'] ==group]
    plt.scatter(subset['PC1'], subset['PC2'], color= color, label= group, s=80)

plt.xlabel(f"PC1({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("PCA: wild-type vs SNF2 mutant")
plt.legend()
plt.savefig("../results/pca_plot.png", dpi=300)
plt.close()

# %%
row =[]

with gzip.open('../data/raw/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz','rt') as f:
    for line in f:
        if not line.startswith('>'):
            continue

        line= line.strip()[1:]
        transcript_id = line.split(" ")[0]

        symbol= ""
        if 'gene_symbol' in line:
            symbol = line.split('gene_symbol:')[1].split(" ")[0]

        description = ""
        if 'description' in line:
            description= line.split('description:')[1]

        row.append({
            'transcript_id': transcript_id,
            'symbol': symbol,
            'description': description
        })

annot =pd.DataFrame(row).set_index('transcript_id')
annot.shape
annot.head()


# %%
clean_results.shape
clean_results.head()
print("Up in mutant:  ", (clean_results["log2FoldChange"] > 0).sum())
print("Down in mutant:", (clean_results["log2FoldChange"] < 0).sum())

annotated =sorted_results.join(annot, how='left')
print(annotated['description'].isna().sum())
pd.set_option('display.max_colwidth', 300)
annotated[['log2FoldChange','padj', 'symbol', 'description']]

# %%
annotated[annotated['log2FoldChange'] > 0][['symbol', 'description', 'log2FoldChange', 'padj']].head(15)


# %%
plot_df=result_df.copy()
plot_df['neg_log10_padj']= -np.log10(plot_df['padj'].replace(0, 1e-300))

def local_gene(row):

    if row['padj']<.05 and row['log2FoldChange']> 1:
        return "Up in Mutant"
    elif row['padj']<.05 and row['log2FoldChange']<-1:
        return "Down in Mutant"
    else:
        return "Not Significant"

plot_df['Category']= plot_df.apply(local_gene, axis=1)
print(plot_df['Category'].value_counts())

# %%
colors= {"Up in Mutant":"red", "Down in Mutant":'blue', "Not Significant":'grey'}

for category, color in colors.items():
    subset= plot_df[plot_df['Category']==category]
    plt.scatter(subset['log2FoldChange'],subset['neg_log10_padj'], color=color, s=6, alpha= 0.6, label=f"{category}({len(subset)})")

plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8)
plt.axvline(1, color="black", linestyle="--", linewidth=0.8)
plt.axvline(-1, color="black", linestyle="--", linewidth=0.8)

plt.xlabel('log2 fold change (Mutated vs Normal)')
plt.ylabel('-log10 adjusted p-value')
plt.title("Volcano plot: SNF2 mutant vs wild-type")
plt.legend()
plt.savefig("../results/volcano_plot.png", dpi=300)
plt.show()


# %%
result_df.to_csv("../results/dge_all_genes.csv")
annotated.to_csv("../results/dge_significant_genes_annotated.csv")
print("Saved")

# %%
libraries =gp.get_library_name(organism="yeast")
print(len(libraries))
print([lib for lib in libraries if "GO"in lib])

go_bp=gp.get_library('GO_Biological_Process_2018', organism="yeast")
print(len(go_bp))
first=list(go_bp.keys())[0]
print(first)
print(go_bp[first])

systematic = annot.index.to_series().str.replace("_mRNA","", regex=False)
annot['gene_name']= annot['symbol'].where(annot['symbol'] !="", systematic)
annot.loc[["YHR174W_mRNA", "YPL071C_mRNA"], ["symbol", "gene_name"]]


down_ids = clean_results[clean_results["log2FoldChange"] < 0].index
up_ids = clean_results[clean_results["log2FoldChange"] > 0].index

down_genes = annot.loc[down_ids, "gene_name"].tolist()
up_genes = annot.loc[up_ids, "gene_name"].tolist()
background = annot.loc[result_df.index, "gene_name"].tolist()

library_genes= set()
for go in go_bp.values():
    library_genes.update(go)

found= [g for g in background if g in library_genes]
print(len(library_genes), len(found))

# %% [markdown]
# 

# %%
enr_down = gp.enrichr(
    gene_list= down_genes,
    gene_sets= go_bp,
    background=background,
    outdir=None
)

res_down= enr_down.results.sort_values('Adjusted P-value')
res_down[["Term", "Overlap", "Adjusted P-value", "Genes"]].head(10)
res_down[res_down["Term"].str.contains("glycoly")][["Term", "Overlap", "Adjusted P-value"]]

enr_up= gp.enrichr(
    gene_list= up_genes,
    gene_sets= go_bp,
    background=background,
    outdir=None
)

res_up= enr_up.results.sort_values('Adjusted P-value')
res_up[["Term", "Overlap", "Adjusted P-value", "Genes"]].head(10)

top = res_down.head(10).copy()
top["neg_log10_padj"] = -np.log10(top["Adjusted P-value"])
top["label"] = top["Term"].str.replace(r" \(GO:\d+\)", "", regex=True)

plt.figure(figsize=(8, 5))
plt.barh(top["label"], top["neg_log10_padj"], color="steelblue")
plt.axvline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8)
plt.gca().invert_yaxis()

plt.xlabel("-log10 adjusted p-value")
plt.title("GO Biological Process: genes down in snf2 mutant")
plt.tight_layout()
plt.savefig("results/go_down_barplot.png", dpi=300)
plt.close()

res_down.to_csv("results/go_bp_down.csv", index=False)
res_up.to_csv("results/go_bp_up.csv", index=False)



