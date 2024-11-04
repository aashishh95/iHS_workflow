iHS Workflow
Welcome to the iHS Workflow for calculating iHS on genomic data! This workflow is tailored for population genetics research and includes steps for phasing, ancestral allele polarization, and selection scan analysis. Below are detailed instructions for setup and usage, including information on file formats and key components.

Getting Started
Step 1: Set Up the Environment
Before submitting a job to SLURM, activate your Conda environment with the necessary modules and software pre-installed:

Modules: beagle, bcftools
Software: selscan, plink, plink2
Ensure these dependencies are properly configured within your environment to enable smooth execution of the workflow.

Step 2: Check Sample File Format
Review the format of the sample file located in samples.txt. This file is critical for correctly processing and analyzing your data, so verify its structure and content before starting.

Step 3: Understand .map File Formats
This workflow utilizes two different genetic map files:

Genetic Map (.map): Used by Beagle for phasing.
Genetic Map (combined.txt): Used to update the genetic map of the .bim file.
Each map file serves a distinct purpose, so please inspect both files to understand their formats and intended uses.

Step 4: Prepare Ancestral Allele Files
To recode for ancestral alleles, this workflow uses a Python script (recode.py) adapted from the Obed project. You can generate ancestral allele files for both hg18 and hg38 references, based on your specific data needs.


Step 5: Check bash_iHS_workflow.sh file. This shows how to pass argument to run the iHS_wrokflow.py script.


Additional Resources
Feel free to reach out for support or questions regarding this workflow:

Contact: vet.aashish2020@gmail.com
