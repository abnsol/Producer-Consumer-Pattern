# How to run

## Make scripts executable and run the setup
```bash
chmod +x scripts/setup_minimal_real.sh
chmod +x scripts/install_r_packages.R
./scripts/setup_minimal_real.sh
```

## Export environment variables (or configure your Config.from_env to read these)
```sh
export PLINK_DIR="$PWD/data/1000Genomes_phase3/plink_format_b37"
export REF_GENOME="GRCh37" # if your Config uses it
```

## Sanity checks
```sh
plink --version; plink2 --version; gcta64 --version
R -q -e 'sessionInfo()'
source .venv/bin/activate; python -c 'import rpy2; import cyvcf2; print("OK")'
```

## Run a quick end-to-end with the tiny dataset
- Use data/gwas_sumstats/mock_gwas.tsv as the GWAS input file.
- Your PLINK reference prefix created is:
    - $PLINK_DIR/EUR/EUR.22.1000Gp3.20130502.{bed,bim,fam}
- The pipeline will only find chr22; other chromosomes will be skipped with warnings (expected).

## Notes and limitations

- This is the real stack (MungeSumstats, COJO via GCTA, LD via PLINK2, finemapping via susieR) on tiny data. It validates code paths and integration. Results are synthetic and not biologically meaningful.
- Ensure Config.plink_dir pulls from the PLINK_DIR env var (or set it directly in your Config).
- If an R package fails to build (common on fresh systems), re-run setup after checking its error; the script installs the usual system libs required for susieR/MungeSumstats.

Optional: without sudo If your Codespace disallows sudo to /usr/local/bin, replace moves to $HOME/.local/bin and add to PATH:
```sh
mkdir -p $HOME/.local/bin
mv plink plink2 gcta64 $HOME/.local/bin
echo 'export PATH="HOME/.local/bin:PATH"' >> ~/.bashrc && source ~/.bashrc
```

Troubleshooting tips

- rpy2 import errors: ensure R is installed and R_HOME is discoverable by rpy2. pip installing rpy2 after R is usually sufficient.
- COJO finds no signals: that’s okay on synthetic data. The pipeline still exercises the code paths.
- ID mismatches: your munge step already preserves RS_ID; the mock PLINK .bim uses rsIDs, so COJO and PLINK2 should match on RS_ID/SNP. If not, log the IDs used and verify the .bim second column contains rsIDs like rs123.