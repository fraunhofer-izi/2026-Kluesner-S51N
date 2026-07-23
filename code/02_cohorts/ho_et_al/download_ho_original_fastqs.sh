#!/bin/bash
set -euo pipefail

BASE_DIR="/mnt/ribolution/user_worktmp/michael.rade/work/2026-Kluesner-FredHutch/ho_et_al/fastq/"
mkdir -p "$BASE_DIR"

submit_download_job () {
    local SRR="$1"
    local URL="$2"

    local OUT_DIR="${BASE_DIR}/${SRR}"
    mkdir -p "$OUT_DIR"

    local BASENAME
    BASENAME="$(basename "$URL")"

    # Remove trailing .1 so the local file ends in .fastq.gz
    local OUT_NAME="${BASENAME%.1}"

    local JOB_SCRIPT="${OUT_DIR}/download_${OUT_NAME}.sh"

    cat > "$JOB_SCRIPT" <<EOT
#!/bin/bash
#SBATCH --job-name=dl_${SRR}
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=${OUT_DIR}/dl_${OUT_NAME}.%j.out
#SBATCH --error=${OUT_DIR}/dl_${OUT_NAME}.%j.err

set -euo pipefail

OUT_FILE="${OUT_DIR}/${OUT_NAME}"
URL="${URL}"

echo "Downloading:"
echo "\$URL"
echo "To:"
echo "\$OUT_FILE"

if [[ -s "\$OUT_FILE" ]]; then
    echo "Output already exists and is non-empty, skipping:"
    ls -lh "\$OUT_FILE"
    exit 0
fi

curl -fL --retry 5 --retry-delay 10 -o "\${OUT_FILE}.partial" "\${URL}"

mv "\${OUT_FILE}.partial" "\${OUT_FILE}"

ls -lh "\${OUT_FILE}"
EOT

    sbatch "$JOB_SCRIPT"
}

submit_download_job "SRR35658576" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658576/HTS_SO188_08_492_PBMC_S8_I1_001.fastq.gz.1"
submit_download_job "SRR35658576" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658576/HTS_SO188_08_492_PBMC_S8_I2_001.fastq.gz.1"
submit_download_job "SRR35658576" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658576/HTS_SO188_08_492_PBMC_S8_R1_001.fastq.gz.1"
submit_download_job "SRR35658576" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658576/HTS_SO188_08_492_PBMC_S8_R2_001.fastq.gz.1"

submit_download_job "SRR35658577" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658577/HTS_SO188_07_470_PBMC_S7_I1_001.fastq.gz.1"
submit_download_job "SRR35658577" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658577/HTS_SO188_07_470_PBMC_S7_I2_001.fastq.gz.1"
submit_download_job "SRR35658577" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658577/HTS_SO188_07_470_PBMC_S7_R1_001.fastq.gz.1"
submit_download_job "SRR35658577" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658577/HTS_SO188_07_470_PBMC_S7_R2_001.fastq.gz.1"

submit_download_job "SRR35658578" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658578/HTS_SO188_06_439_PBMC_S6_I1_001.fastq.gz.1"
submit_download_job "SRR35658578" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658578/HTS_SO188_06_439_PBMC_S6_I2_001.fastq.gz.1"
submit_download_job "SRR35658578" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658578/HTS_SO188_06_439_PBMC_S6_R1_001.fastq.gz.1"
submit_download_job "SRR35658578" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658578/HTS_SO188_06_439_PBMC_S6_R2_001.fastq.gz.1"

submit_download_job "SRR35658579" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658579/HTS_SO188_04_492_CSF_S4_I1_001.fastq.gz.1"
submit_download_job "SRR35658579" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658579/HTS_SO188_04_492_CSF_S4_I2_001.fastq.gz.1"
submit_download_job "SRR35658579" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658579/HTS_SO188_04_492_CSF_S4_R1_001.fastq.gz.1"
submit_download_job "SRR35658579" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658579/HTS_SO188_04_492_CSF_S4_R2_001.fastq.gz.1"

submit_download_job "SRR35658581" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658581/HTS_SO188_05_251_PBMC_S5_I1_001.fastq.gz.1"
submit_download_job "SRR35658581" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658581/HTS_SO188_05_251_PBMC_S5_I2_001.fastq.gz.1"
submit_download_job "SRR35658581" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658581/HTS_SO188_05_251_PBMC_S5_R1_001.fastq.gz.1"
submit_download_job "SRR35658581" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658581/HTS_SO188_05_251_PBMC_S5_R2_001.fastq.gz.1"

submit_download_job "SRR35658582" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658582/HTS_SO188_02_439_CSF_S2_I1_001.fastq.gz.1"
submit_download_job "SRR35658582" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658582/HTS_SO188_02_439_CSF_S2_I2_001.fastq.gz.1"
submit_download_job "SRR35658582" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658582/HTS_SO188_02_439_CSF_S2_R1_001.fastq.gz.1"
submit_download_job "SRR35658582" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658582/HTS_SO188_02_439_CSF_S2_R2_001.fastq.gz.1"

submit_download_job "SRR35658583" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658583/HTS_SO188_01_251_CSF_S1_I1_001.fastq.gz.1"
submit_download_job "SRR35658583" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658583/HTS_SO188_01_251_CSF_S1_I2_001.fastq.gz.1"
submit_download_job "SRR35658583" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658583/HTS_SO188_01_251_CSF_S1_R1_001.fastq.gz.1"
submit_download_job "SRR35658583" "https://sra-pub-src-1.s3.amazonaws.com/SRR35658583/HTS_SO188_01_251_CSF_S1_R2_001.fastq.gz.1"

