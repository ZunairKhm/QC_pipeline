#!/bin/bash --login
#SBATCH --job-name=run_bracken
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4


# The project ID which this job should run under:
#SBATCH --account=pawsey1172

module load  singularity/4.1.0-slurm
# 1. Define your directories and parameters
# Update these paths to match your actual directories
REPORTS_DIR="results/kraken2"  # Directory containing Kraken2 report files
DB_DIR="../databases/kraken2_database"
OUTPUT_DIR="results/bracken"

READ_LEN=150 # Must match the kmer length built into your Kraken2 database
TAX_LVL="S"  # S = Species, G = Genus

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# 2. Define the container
# This pulls the official BioContainers image for Bracken 2.9
BRACKEN_CONTAINER="docker://quay.io/biocontainers/bracken:2.9--py311h2a4ad6c_1"

echo "Starting Bracken processing..."

# 3. Loop through all Kraken2 reports
# Adjust the '*.report' extension below if your files are named differently (e.g., *.kraken2.report.txt)
for report in "${REPORTS_DIR}"/*.report.txt; do
    
    # Extract the base filename to keep output names consistent
    base_name=$(basename "${report}" .report.txt)
    
    echo "Processing: ${base_name}"
    
    # Execute Bracken inside the Apptainer container
    singularity exec "${BRACKEN_CONTAINER}" bracken \
        -d "${DB_DIR}" \
        -i "${report}" \
        -o "${OUTPUT_DIR}/${base_name}.bracken" \
        -w "${OUTPUT_DIR}/${base_name}.bracken.report.txt" \
        -r ${READ_LEN} \
        -l ${TAX_LVL}
        
done

echo "All reports processed successfully!"