import os
import numpy as np
import pandas as pd
from typing import Tuple, Dict, Optional, List
from dataclasses import dataclass
from scipy.stats import norm
import subprocess
from pathlib import Path

@dataclass
class RealisticMockConfig:
    # Cohort and region
    chrom: int = 22
    num_individuals: int = 5000
    num_snps: int = 3000
    region_start_bp: int = 16_050_000
    region_end_bp: int = 51_304_566

    # Genetic architecture
    num_causal: int = 5
    heritability: float = 0.5
    effect_scale: float = 1.0  # multiply causal effects to strengthen signals
    min_maf: float = 0.01
    max_maf: float = 0.5

    # LD model (Markov haplotype copying)
    ld_decay_kb: float = 20.0  # correlation ~ exp(-d_kb / ld_decay_kb)
    hotspot_prob: float = 0.01  # probability of a recombination hotspot between adjacent SNPs
    hotspot_ld_multiplier: float = 0.2  # reduce copying prob in hotspots

    # Outputs
    population: str = "EUR"
    base_dir: str = "data"
    gwas_filename: str = "realistic_mock_gwas.tsv"
    plink_prefix_template: str = "{pop}/{pop}.{chrom}.1000Gp3.20130502"

    # Misc
    seed: int = 42
    clean_text_intermediates: bool = True  # remove .ped/.map after PLINK conversion if successful


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(seed)


def _generate_positions(cfg: RealisticMockConfig, rng: np.random.Generator) -> np.ndarray:
    # Uniformly spread SNPs across the region
    positions = rng.integers(cfg.region_start_bp, cfg.region_end_bp, size=cfg.num_snps)
    positions.sort()
    return positions


def _random_alleles(cfg: RealisticMockConfig, rng: np.random.Generator, m: int) -> Tuple[np.ndarray, np.ndarray]:
    bases = np.array(list("ATCG"))
    a1 = rng.choice(bases, size=m, replace=True)
    a2 = rng.choice(bases, size=m, replace=True)
    # ensure not equal
    same = (a1 == a2)
    while same.any():
        a2[same] = rng.choice(bases, size=same.sum(), replace=True)
        same = (a1 == a2)
    return a1, a2


def _sample_mafs(cfg: RealisticMockConfig, rng: np.random.Generator, m: int) -> np.ndarray:
    # Realistic MAF distribution: Beta(0.5, 0.5) skewed to rare variants, truncated to [min_maf, max_maf]
    maf = rng.beta(0.5, 0.5, size=m)
    maf = np.clip(maf, cfg.min_maf, cfg.max_maf)
    return maf


def _simulate_haplotypes(cfg: RealisticMockConfig, positions: np.ndarray, maf: np.ndarray, n_hap: int, rng: np.random.Generator) -> np.ndarray:
    """
    Simulate haplotypes (0/1) for n_hap chromosomes with distance-based LD.
    Markov copying: with prob rho = exp(-d_kb/ld_decay), copy previous allele; else resample Bernoulli(maf_i).
    Hotspots reduce copying probability between adjacent SNPs.
    """
    m = len(positions)
    H = np.zeros((n_hap, m), dtype=np.int8)

    # Precompute adjacent distances in kb
    d_bp = np.diff(positions, prepend=positions[0])
    d_kb = np.maximum(d_bp / 1000.0, 1e-3)

    # Random hotspots between adjacent SNPs
    hotspot_mask = rng.random(m) < cfg.hotspot_prob

    for h in range(n_hap):
        # First SNP: draw from Bernoulli(maf)
        H[h, 0] = rng.random() < maf[0]
        for i in range(1, m):
            rho = np.exp(-d_kb[i] / cfg.ld_decay_kb)
            if hotspot_mask[i]:
                rho *= cfg.hotspot_ld_multiplier
            if rng.random() < rho:
                H[h, i] = H[h, i-1]
            else:
                H[h, i] = rng.random() < maf[i]
    return H


def _haplotypes_to_genotypes(H: np.ndarray) -> np.ndarray:
    """
    Pair consecutive haplotypes into diploids.
    H shape: (2N, M) of 0/1; returns G shape (N, M) of 0/1/2 counts of allele '1' (A1 placeholder).
    """
    assert H.shape[0] % 2 == 0, "Number of haplotypes must be even"
    N = H.shape[0] // 2
    G = H[0::2, :] + H[1::2, :]
    return G.astype(np.int8)


def _ensure_biallelic_and_maf(G: np.ndarray, maf: np.ndarray, positions: np.ndarray, a1: np.ndarray, a2: np.ndarray,
                              cfg: RealisticMockConfig) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Filter monomorphic or MAF out-of-range SNPs after sampling.
    """
    N = G.shape[0]
    freq = G.mean(axis=0) / 2.0
    keep = (freq >= cfg.min_maf) & (freq <= (1.0 - cfg.min_maf)) & (G.var(axis=0) > 0)
    G2 = G[:, keep]
    maf2 = np.minimum(freq[keep], 1.0 - freq[keep])
    pos2 = positions[keep]
    a1_2 = a1[keep]
    a2_2 = a2[keep]
    return G2, maf2, pos2, a1_2, a2_2


def _choose_causals(cfg: RealisticMockConfig, m: int, rng: np.random.Generator) -> np.ndarray:
    c = min(cfg.num_causal, m)
    return np.sort(rng.choice(m, size=c, replace=False))


def _simulate_trait_and_sumstats(G: np.ndarray,
                                 positions: np.ndarray,
                                 a1: np.ndarray,
                                 a2: np.ndarray,
                                 cfg: RealisticMockConfig,
                                 causal_idx: np.ndarray,
                                 rng: np.random.Generator) -> Tuple[pd.DataFrame, Dict]:
    """
    Simulate quantitative trait with specified heritability and compute marginal GWAS stats.
    Orientation:
      - A1 is effect (counted) allele for computing y and sumstats.
      - If A1 frequency > 0.5, flip A1/A2 and the sign of BETA/Z so that A1 is minor.
    """
    N, M = G.shape
    # Effects for causal SNPs
    beta = rng.normal(loc=0.0, scale=1.0, size=len(causal_idx)) * cfg.effect_scale

    # Genetic value
    gval = (G[:, causal_idx] @ beta).astype(np.float64)
    var_g = np.var(gval)
    if var_g <= 1e-12:
        var_g = 1e-12

    # Residual variance to achieve h2
    h2 = np.clip(cfg.heritability, 1e-6, 0.999999)
    var_e = var_g * (1 - h2) / h2
    e = rng.normal(0.0, np.sqrt(var_e), size=N)

    y = gval + e
    y_mean = y.mean()

    # Precompute for regression
    G_mean = G.mean(axis=0)
    G_var = G.var(axis=0)
    Sxx = ((G - G_mean) ** 2).sum(axis=0)
    # beta_hat = cov(G, y) / var(G)
    cov_Gy = ((G - G_mean) * (y[:, None] - y_mean)).sum(axis=0) / (N - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        beta_hat = np.where(G_var > 0, cov_Gy / G_var, 0.0)

    # Residuals and SE per SNP (simple OLS with intercept)
    alpha_hat = y_mean - beta_hat * G_mean
    y_hat = alpha_hat + G * beta_hat
    resid = y[:, None] - y_hat
    # SSE/(N-2)
    dof = max(N - 2, 1)
    sigma2_i = (resid ** 2).sum(axis=0) / dof
    with np.errstate(divide="ignore", invalid="ignore"):
        se_hat = np.sqrt(sigma2_i / Sxx)
        z = np.where(se_hat > 0, beta_hat / se_hat, 0.0)

    # Orientation: ensure A1 is the minor allele
    frq = G_mean / 2.0  # frequency of counted allele (current A1)
    flip = frq > 0.5
    # Flip alleles where needed
    a1_flipped = a1.copy()
    a2_flipped = a2.copy()
    a1_flipped[flip], a2_flipped[flip] = a2[flip], a1[flip]
    # Adjust stats on flip (count of minor allele = 2 - G)
    frq_adj = np.where(flip, 1.0 - frq, frq)
    beta_hat_adj = np.where(flip, -beta_hat, beta_hat)
    z_adj = np.where(flip, -z, z)
    se_hat_adj = se_hat  # same magnitude under flipping
    # Two-sided p-value (normal approx)
    pvals = 2.0 * norm.sf(np.abs(z_adj))

    # Build DataFrame
    snp_ids = np.array([f"rs{100_000_000 + i}" for i in range(M)], dtype=object)
    df = pd.DataFrame({
        "CHR": cfg.chrom,
        "BP": positions,
        "SNP": snp_ids,
        "A1": a1_flipped,
        "A2": a2_flipped,
        "BETA": beta_hat_adj,
        "SE": se_hat_adj,
        "P": pvals,
        "Z": z_adj,
        "N": N,
        "FRQ": frq_adj,
    })

    meta = {
        "causal_idx": causal_idx.tolist(),
        "beta_causal": beta.tolist(),
        "var_g": float(var_g),
        "var_e": float(var_e),
        "h2": float(h2),
        "effect_scale": float(cfg.effect_scale),
        "num_significant_p5e8": int((pvals < 5e-8).sum()),
        "num_snps": int(M),
        "num_individuals": int(N)
    }
    return df, meta


def _write_ped_map(G: np.ndarray, df: pd.DataFrame, out_prefix: str):
    """
    Write PLINK .ped/.map using df for SNP metadata. A1 is effect allele (minor), A2 other.
    .map: CHR, SNP, CM=0, BP
    .ped: FID IID PID MID SEX PHENOTYPE then pairs of alleles per SNP.
    """
    out_dir = Path(out_prefix).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # .map
    map_df = pd.DataFrame({
        "CHR": df["CHR"].astype(int),
        "SNP": df["SNP"],
        "CM": 0,
        "BP": df["BP"].astype(int),
    })
    map_path = f"{out_prefix}.map"
    map_df.to_csv(map_path, sep="\t", header=False, index=False)

    # .ped
    N, M = G.shape
    # For each SNP, create allele pairs based on genotype counts of A1 (minor)
    # 0 -> A2 A2, 1 -> A1 A2, 2 -> A1 A1
    A1 = df["A1"].to_numpy()
    A2 = df["A2"].to_numpy()

    ped_rows: List[List[str]] = []
    # Minimal family info
    for i in range(N):
        row = [f"FAM{i+1}", f"IND{i+1}", "0", "0", "1", "-9"]
        ped_genos = []
        gi = G[i, :]
        # Build allele pairs
        a1_pair = np.column_stack((A1, A1))
        a2_pair = np.column_stack((A2, A2))
        het_pair = np.column_stack((A1, A2))
        # Map counts to pairs
        pairs = np.where(gi[:, None] == 0, a2_pair,
                 np.where(gi[:, None] == 1, het_pair, a1_pair))
        # Flatten per SNP (two columns per SNP)
        for j in range(pairs.shape[0]):
            ped_genos.extend([pairs[j, 0], pairs[j, 1]])
        row.extend(ped_genos)
        ped_rows.append(row)

    ped_df = pd.DataFrame(ped_rows)
    ped_path = f"{out_prefix}.ped"
    ped_df.to_csv(ped_path, sep=" ", header=False, index=False)
    return ped_path, map_path


def _plink_make_bed(ped_prefix: str, bed_prefix: str) -> bool:
    """
    Convert .ped/.map to .bed/.bim/.fam using PLINK 1.9, if available.
    Returns True if .bed created successfully.
    """
    try:
        subprocess.run([
            "plink",
            "--file", ped_prefix,
            "--make-bed",
            "--out", bed_prefix,
            "--allow-no-sex"
        ], check=True, capture_output=True, text=True)
    except FileNotFoundError:
        print("[MOCK GEN] PLINK 1.9 not found on PATH — skipping binary conversion.")
        return False
    except subprocess.CalledProcessError as e:
        print("[MOCK GEN] PLINK failed:\n", e.stderr)
        return False
    return os.path.exists(f"{bed_prefix}.bed") and os.path.exists(f"{bed_prefix}.bim") and os.path.exists(f"{bed_prefix}.fam")


def generate_realistic_gwas_and_plink(cfg: Optional[RealisticMockConfig] = None) -> Dict[str, str]:
    """
    Generate realistic GWAS summary statistics and a matching small PLINK reference panel
    (ped/map -> bed/bim/fam) for a single chromosome region.

    Returns paths:
      - gwas_sumstats_path
      - plink_prefix (for .bed/.bim/.fam)
    """
    cfg = cfg or RealisticMockConfig()
    rng = _rng(cfg.seed)

    # Prepare output paths
    base_dir = Path(cfg.base_dir)
    gwas_dir = base_dir / "gwas_sumstats"
    plink_dir = base_dir / "1000Genomes_phase3" / "plink_format_b37" / cfg.population
    plink_dir.mkdir(parents=True, exist_ok=True)
    gwas_dir.mkdir(parents=True, exist_ok=True)

    plink_prefix_rel = cfg.plink_prefix_template.format(pop=cfg.population, chrom=cfg.chrom)
    plink_prefix = str(plink_dir / plink_prefix_rel.split("/", 1)[-1] if "/" in plink_prefix_rel else plink_prefix_rel)
    plink_prefix_full = str(plink_dir / Path(plink_prefix_rel).name)

    # 1) Positions, alleles, MAF
    positions = _generate_positions(cfg, rng)
    a1, a2 = _random_alleles(cfg, rng, len(positions))
    maf = _sample_mafs(cfg, rng, len(positions))

    # 2) Haplotypes with LD, then genotypes
    H = _simulate_haplotypes(cfg, positions, maf, n_hap=2 * cfg.num_individuals, rng=rng)
    G = _haplotypes_to_genotypes(H)

    # 3) Filter monomorphic / extreme MAF after sampling
    G, maf_obs, positions, a1, a2 = _ensure_biallelic_and_maf(G, maf, positions, a1, a2, cfg)

    # 4) Choose causal variants and simulate trait + sumstats
    causal_idx = _choose_causals(cfg, G.shape[1], rng)
    sumstats_df, meta = _simulate_trait_and_sumstats(G, positions, a1, a2, cfg, causal_idx, rng)

    # 5) Save GWAS sumstats
    gwas_path = str(gwas_dir / cfg.gwas_filename)
    sumstats_df.to_csv(gwas_path, sep="\t", index=False)
    print(f"[MOCK GEN] GWAS sumstats saved to {gwas_path}  (N={cfg.num_individuals}, M={G.shape[1]})")
    print(f"[MOCK GEN] Significant (p<5e-8): {meta['num_significant_p5e8']} (h2={meta['h2']:.2f}, effect_scale={meta['effect_scale']:.2f})")

    # 6) Write PLINK text and convert to binary (if PLINK installed)
    ped_prefix = str(plink_dir / "mock_plink_data")
    _ = _write_ped_map(G, sumstats_df, ped_prefix)
    ok = _plink_make_bed(ped_prefix, plink_prefix_full)

    if ok and cfg.clean_text_intermediates:
        for ext in (".ped", ".map"):
            try:
                os.remove(f"{ped_prefix}{ext}")
            except OSError:
                pass

    if ok:
        print(f"[MOCK GEN] PLINK binary files generated at prefix: {plink_prefix_full}")
    else:
        print(f"[MOCK GEN] PLINK binaries not created. You can run manually:\n"
              f"  plink --file {ped_prefix} --make-bed --out {plink_prefix_full} --allow-no-sex")

    return {
        "gwas_sumstats_path": gwas_path,
        "plink_prefix": plink_prefix_full
    }


if __name__ == "__main__":
    # Example usage with defaults; adjust as needed
    cfg = RealisticMockConfig(
        chrom=22,
        num_individuals=5000,
        num_snps=3000,
        num_causal=5,
        heritability=0.5,
        effect_scale=1.5,  # bump to see a few genome-wide signals with small N
        seed=42
    )
    paths = generate_realistic_gwas_and_plink(cfg)
    print("[MOCK GEN] Done:", paths)