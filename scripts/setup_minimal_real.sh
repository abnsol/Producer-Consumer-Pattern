#!/usr/bin/env bash
set -euo pipefail

# This script prepares a Codespace for a "real minimal" run:
# - Installs system deps, R, PLINK 1.9, PLINK2, GCTA64
# - Creates a Python venv and installs Python deps (incl. rpy2)
# - Installs R packages (susieR, MungeSumstats) via a companion R script
# - Generates a tiny PLINK dataset for chr22
#
# After this, set PLINK_DIR and run your pipeline on the tiny dataset.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="/usr/local/bin"

echo "[1/8] Update APT and install system dependencies + R"
sudo apt-get update -y
sudo apt-get install -y \
  curl wget unzip ca-certificates \
  build-essential gfortran \
  r-base r-base-dev \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  liblapack-dev libblas-dev

echo "[2/8] Install PLINK 1.9 (needed to convert .ped/.map -> .bed/.bim/.fam)"
# Official builds: https://www.cog-genomics.org/plink/1.9/
pushd /tmp >/dev/null
PLINK19_ZIP="plink_linux_x86_64_20231211.zip"
wget -q https://s3.amazonaws.com/plink1-assets/${PLINK19_ZIP}
unzip -q ${PLINK19_ZIP}
sudo mv plink ${BIN_DIR}/plink
rm -f ${PLINK19_ZIP}
popd >/dev/null
plink --version || { echo "PLINK 1.9 not found on PATH"; exit 1; }

echo "[3/8] Install PLINK 2"
# Official builds: https://www.cog-genomics.org/plink/2.0/
pushd /tmp >/dev/null
PLINK2_ZIP="plink2_linux_x86_64.zip"
wget -q -O ${PLINK2_ZIP} https://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip
unzip -q ${PLINK2_ZIP}
sudo mv plink2 ${BIN_DIR}/plink2
rm -f ${PLINK2_ZIP}
popd >/dev/null
plink2 --version || { echo "PLINK 2 not found on PATH"; exit 1; }

echo "[4/8] Install GCTA64"
# Official builds: https://yanglab.westlake.edu.cn/software/gcta/
pushd /tmp >/dev/null
GCTA_ZIP="gcta_1.94.1beta_linux_kernel_3_x86_64.zip"
wget -q -O ${GCTA_ZIP} https://cnsgenomics.com/software/gcta/${GCTA_ZIP} || \
wget -q -O ${GCTA_ZIP} https://yanglab.westlake.edu.cn/software/gcta/${GCTA_ZIP}
unzip -q ${GCTA_ZIP}
sudo mv gcta64 ${BIN_DIR}/gcta64
rm -f ${GCTA_ZIP}
popd >/dev/null
gcta64 --version || { echo "GCTA64 not found on PATH"; exit 1; }

echo "[5/8] Create Python virtualenv and install Python dependencies"
# If your Codespace already uses a devcontainer venv, you can skip the venv creation.
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip wheel setuptools

# Install Python deps used in the repo
pip install \
  pandas numpy scipy \
  prefect loguru psutil optuna \
  rpy2 cyvcf2

# Optional: install any project-specific extras you need
pip install -r requirements.txt

echo "[6/8] Install R packages (susieR, MungeSumstats, etc.)"
# Uses scripts/install_r_packages.R
Rscript "${ROOT_DIR}/scripts/install_r_packages.R"

echo "[7/8] Generate tiny mock GWAS + PLINK chr22 dataset"
pushd "${ROOT_DIR}" >/dev/null
python - <<'PY'
from mock_data_generator import generate_mock_data
# Smaller numbers keep files tiny; adjust if needed
generate_mock_data(num_snps=1000, num_individuals=50, chrom=22, seed=42)
PY
popd >/dev/null

# The generator already prints a PLINK command; attempt conversion again for safety.
PLINK_EUR_DIR="${ROOT_DIR}/data/1000Genomes_phase3/plink_format_b37/EUR"
MAP_PREFIX="${PLINK_EUR_DIR}/mock_plink_data"
OUT_PREFIX="${PLINK_EUR_DIR}/EUR.22.1000Gp3.20130502"

echo "[7b/8] Converting .ped/.map -> .bed/.bim/.fam with PLINK 1.9"
plink --file "${MAP_PREFIX}" --make-bed --out "${OUT_PREFIX}" --allow-no-sex || true

# Optional: cleanup large text intermediates to save disk space
rm -f "${MAP_PREFIX}.ped" "${MAP_PREFIX}.map" || true

echo "[8/8] Environment hints"
echo
echo "Add to your shell env (bash):"
echo "  export PLINK_DIR='${ROOT_DIR}/data/1000Genomes_phase3/plink_format_b37'"
echo "  export REF_GENOME='GRCh37'   # if your Config uses it"
echo
echo "Verify tools:"
echo "  plink --version && plink2 --version && gcta64 --version"
echo "  R -q -e 'sessionInfo()'"
echo "  python -c 'import rpy2; import cyvcf2; print(\"rpy2/cyvcf2 OK\")'"
echo
echo "Next:"
echo "  - Point Config.plink_dir to \$PLINK_DIR (via env or config)."
echo "  - Use data/gwas_sumstats/mock_gwas.tsv as your GWAS input."
echo "  - Run your pipeline."