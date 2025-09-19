import pandas as pd
from analysis_tasks import munge_sumstats_preprocessing, filter_significant_variants, finemap_region

# Load mock GWAS
df = pd.read_csv("data/gwas_sumstats/mock_gwas.tsv", sep="\t")

# Run Munge (mock)
munged_df, munged_path = munge_sumstats_preprocessing.fn("data/gwas_sumstats/mock_gwas.tsv", "data/out_munge")

# Filter significant SNPs
sig_df, sig_path = filter_significant_variants.fn(munged_df, "data/out_filter")

# Run finemapping (mock)
results = finemap_region.fn(seed=42, sumstats=sig_df, chr_num=22, lead_variant_position=int(sig_df.iloc[0]['BP']))
print(results.head())
