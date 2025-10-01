#### How to use

- Save this as realistic_mock_data_generator.py at repo root.
- Run it:
```sh
python realistic_mock_data_generator.py
```
- It will produce:
```sh
data/gwas_sumstats/realistic_mock_gwas.tsv
data/1000Genomes_phase3/plink_format_b37/EUR/EUR.22.1000Gp3.20130502.{bed,bim,fam} (if PLINK is installed)
```
- Point Config.plink_dir to data/1000Genomes_phase3/plink_format_b37 and use the GWAS file above as input.

#### Tuning

- effect_scale: increase if you want more genome-wide significant hits with small N.
- num_causal and heritability: adjust polygenicity and signal strength.
- ld_decay_kb and hotspot_*: shape local LD; smaller decay → stronger LD.

This generator ensures:

- rsID-based SNP IDs compatible with the pipelinex’s RS_ID handling.
- Allele alignment: A1 is the minor/effect allele in sumstats and in PLINK (after conversion), making - PLINK2 LD and SuSiE orientation consistent.
- Realistic-looking summary stats that exercise the full “real” pipeline end to end.