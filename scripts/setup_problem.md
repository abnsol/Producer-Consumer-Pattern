### Why is there heavy download ? 

- The thing that blows up disk usage is downloading full reference PLINK panels (e.g., 1000 Genomes across all chromosomes/populations). 
- **Binaries and packages are relatively small.**
- Locally generated tiny PLINK dataset will run the “real” pipeline end-to-end on small tests within limits. It’s functionally real but not biologically meaningful (random synthetic data).

### Estimated sizes 
- These are ballpark ranges for Linux x86_64; actual numbers vary by distro and version.
    - #### PLINK 2 binary
        - Download: 10–20 MB (zip)
        - Installed: 8–15 MB (binary)

    - #### GCTA64 binary
        - Download: 5–20 MB (zip)
        - Installed: 3–10 MB (binary)

    - #### R base + dev toolchain (r-base, r-base-dev, plus runtime deps)
        - Download: 100–300 MB
        - Installed: 500–1500 MB

    - #### rpy2 (pip)
        - Download: 3–10 MB
        - Installed: 10–30 MB

    - #### susieR (R package)
        - Download: 1–5 MB
        - Installed: 5–30 MB (includes Rcpp/Armadillo deps if not already present)

    - #### System dev libs (libcurl4-openssl-dev, libssl-dev, libxml2-dev)
        - Download: 10–50 MB combined
        - Installed: 20–100 MB combined

    - #### local PLINK dataset (from mock_data_generator.py)
        - .bed/.bim/.fam: well under 1–5 MB for 1k–5k SNPs and 50–100 individuals
        - The intermediate .ped/.map text files can be a few–tens of MB; you can delete them after conversion

    - #### Full 1000 Genomes PLINK reference panels
        - Download/Installed: tens of GB (20–60+ GB depending on population coverage and format). This is what we should avoid.

### Notes for your 4 cores / 16 GB RAM Codespace

- Plenty for the small “real” run. SuSiE on small LD matrices and PLINK LD computations will be fast.
- If you keep Optuna enabled (n_trials=10), that adds some CPU time; you can set fewer trials for quick tests.
- If space gets tight after apt installs, you can reclaim some:
```bash 
sudo apt-get clean && sudo rm -rf /var/lib/apt/lists/*
```

### Gotchas

- Results from the synthetic dataset are for pipeline validation, not scientific interpretation.
- Mocking functions can be a second optional way to this set up 
- Ensure rpy2’s version is compatible with the R version you install. If there’s a mismatch, pip may fail to import rpy2.