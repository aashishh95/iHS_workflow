#!/bin/bash

# Define the resource requirements here using #SBATCH
#SBATCH -c 10  # Requesting 10 CPUs
#SBATCH -t 24:00:00  # Max wallTime for the job
#SBATCH --mem=128G  # Requesting 128G memory
#SBATCH -o job.%J.out  # Output file
#SBATCH -e job.%J.err  # Error file

# Check if conda is initialized, if not, initialize it
if ! command -v conda &> /dev/null
then
    echo "Conda not found. Initializing conda..."
    conda init
fi

# Activate the Conda environment only if it's not already activated
if [[ "$CONDA_DEFAULT_ENV" != "selscan" ]]; then
    echo "Activating conda environment: selscan"
    conda activate selscan
else
    echo "Conda environment 'selscan' is already activated."
fi

# Check if required modules are loaded; load them if not
module_list=$(module list 2>&1)

if [[ ! "$module_list" =~ "gencore/2" ]]; then
    echo "Loading module: gencore/2"
    module load all gencore/2
else
    echo "Module 'gencore/2' is already loaded."
fi

if [[ ! "$module_list" =~ "beagle" ]]; then
    echo "Loading module: beagle"
    module load beagle
else
    echo "Module 'beagle' is already loaded."
fi

if [[ ! "$module_list" =~ "bcftools" ]]; then
    echo "Loading module: bcftools"
    module load bcftools
else
    echo "Module 'bcftools' is already loaded."
fi

# Run the Python script for iHS workflow
python iHS_workflow.py \
    --plink_file ./sample_plink \
    --sample_file ./sample_population_iHS \
    --ancestral_allele_file ./ancestral_allele_hg38 \
    --beagle_jar ./beagle_jar/beagle.28Jun21.220.jar \
    --genetic_map_dir ./genetic_map_phasing \
    --genetic_map_dir_shapeit ./genetic_map_shapeIT \
    --genome_version 38

# Remove the intermediate files if they exist
if ls ./sample_plink_* 1> /dev/null 2>&1; then
    rm ./sample_plink_*
fi
if ls ./sample_population_* 1> /dev/null 2>&1; then
    rm ./sample_population_iHS_*
fi
if ls ./chr* 1> /dev/null 2>&1; then
    rm ./chr*
fi
