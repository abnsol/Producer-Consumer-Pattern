import pandas as pd
import numpy as np
import os
import subprocess


def generate_mock_data(num_snps=5000, num_individuals=100, chrom=22, seed=42):
    """Generates mock GWAS summary stats and PLINK precursor files."""

    np.random.seed(seed)

    # === Create Directories ===
    print("Creating necessary directories...")
    gwas_dir = "data/gwas_sumstats"
    plink_eur_dir = "data/1000Genomes_phase3/plink_format_b37/EUR"
    os.makedirs(gwas_dir, exist_ok=True)
    os.makedirs(plink_eur_dir, exist_ok=True)

    # === 1. Generate Mock GWAS Summary Statistics ===
    print("Generating mock GWAS summary statistics...")
    bp = np.random.randint(16050000, 51304566, size=num_snps)
    bp.sort()
    alleles = ['A', 'T', 'C', 'G']
    a1 = np.random.choice(alleles, size=num_snps)
    a2 = np.random.choice(alleles, size=num_snps)

    # Ensure alleles differ
    for i in range(num_snps):
        while a1[i] == a2[i]:
            a2[i] = np.random.choice(alleles)

    # P-values (force some to be significant)
    p_values = np.random.uniform(0, 1, size=num_snps)
    sig_idx = np.random.choice(num_snps, size=min(50, num_snps), replace=False)
    p_values[sig_idx] = np.random.uniform(1e-10, 5e-9, size=len(sig_idx))

    gwas_df = pd.DataFrame({
        'CHR': chrom,
        'BP': bp,
        'SNP': [f"rs{i}" for i in range(num_snps)],
        'A1': a1,
        'A2': a2,
        'BETA': np.random.normal(0, 0.1, size=num_snps),
        'SE': np.random.uniform(0.01, 0.05, size=num_snps),
        'P': p_values,
        'N': np.random.randint(50000, 150000, size=num_snps),
        'FRQ': np.random.uniform(0.01, 0.99, size=num_snps)
    })

    gwas_output_path = os.path.join(gwas_dir, "mock_gwas.tsv")
    gwas_df.to_csv(gwas_output_path, sep='\t', index=False)
    print(f"Mock GWAS data saved to {gwas_output_path}")

    # === 2. Generate Mock PLINK .map and .ped files ===
    print("\nGenerating mock PLINK text files (.map and .ped)...")
    map_df = gwas_df[['CHR', 'SNP', 'BP']].copy()
    map_df['CM'] = 0
    map_df = map_df[['CHR', 'SNP', 'CM', 'BP']]
    map_path = os.path.join(plink_eur_dir, "mock_plink_data")
    map_df.to_csv(f"{map_path}.map", sep='\t', index=False, header=False)

    ped_data = []
    genotypes = ['1 1', '1 2', '2 2', '0 0']
    for i in range(num_individuals):
        row = [
            f"FAM{i}", f"IND{i}", "0", "0",
            np.random.choice([1, 2]), "-9"
        ]
        snp_genos = np.random.choice(genotypes, size=num_snps,
                                     p=[0.45, 0.10, 0.44, 0.01])
        row.extend(snp_genos)
        ped_data.append(row)

    ped_df = pd.DataFrame(ped_data)
    ped_df.to_csv(f"{map_path}.ped", sep=' ', index=False, header=False)
    print(f"Mock .map and .ped files saved as {map_path}.[map/ped]")

    # === 3. PLINK Conversion Command ===
    print("\n--- ACTION REQUIRED ---")
    output_prefix = os.path.join(plink_eur_dir, "EUR.22.1000Gp3.20130502")
    plink_cmd = f"plink --file {map_path} --make-bed --out {output_prefix} --allow-no-sex"
    print(f"\nRun this command (if PLINK installed):\n{plink_cmd}\n")

    # Optional: try auto-run
    try:
        subprocess.run(plink_cmd, shell=True, check=True)
        print("PLINK binary files generated successfully!")
    except FileNotFoundError:
        print("PLINK not found on PATH — please install it if you want binaries.")
    except subprocess.CalledProcessError as e:
        print("PLINK failed:\n", e.stderr)


if __name__ == "__main__":
    generate_mock_data(num_snps=1000, num_individuals=50)  # smaller default for testing
