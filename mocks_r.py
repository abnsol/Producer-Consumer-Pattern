import pandas as pd
import numpy as np
import os
from prefect import task

@task
def mock_munge_sumstats(gwas_file_path, output_dir, **kwargs):
    """Mock replacement for MungeSumstats preprocessing."""
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(gwas_file_path, sep="\t")
    out_path = os.path.join(output_dir, "munged_sumstats_processed.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[MOCK MUNGE] Saved {len(df)} SNPs to {out_path}")
    return df, out_path


@task
def mock_finemap_region(seed, sumstats, chr_num, lead_variant_position, **kwargs):
    """Mock replacement for SuSiE fine-mapping."""
    snps = sumstats['SNP'].sample(min(5, len(sumstats)), random_state=seed).tolist()
    df = pd.DataFrame({
        "SNP": snps,
        "PIP": np.random.uniform(0.1, 1, len(snps)),
        "cs": np.random.randint(1, 3, len(snps)),
        "region_id": f"chr{chr_num}:{lead_variant_position}",
        "region_chr": chr_num,
        "region_center": lead_variant_position,
        "converged": True
    })
    print(f"[MOCK FINEMAP] Produced {len(df)} variants for {chr_num}:{lead_variant_position}")
    return df
