import subprocess
import os
import glob

def run_command(command):
    """Helper function to run shell commands"""
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True)

def split_chromosomes(plink_file):
    """Step 1: Split chromosomes, create VCF files, and zip them"""
    for i in range(1, 23):
        vcf_out = f"{plink_file}_chr{i}"
        run_command(f"plink2 --bfile {plink_file} --chr {i} --recode vcf --out {vcf_out}")
        run_command(f"gzip -f {vcf_out}.vcf")  # Use -f to force overwrite

def polarize_ancestral_alleles(plink_file, ancestral_allele_file, genome_version):
    """Step 2: Polarize ancestral alleles using a Python script"""
    for i in range(1, 23):
        # Determine the file name based on the genome version
        if genome_version == 19:
            ancestral_file = f"{ancestral_allele_file}/chr{i}_hg19_AA.txt"
        elif genome_version == 38:
            ancestral_file = f"{ancestral_allele_file}/chr{i}_hg38_AA.txt"
        else:
            raise ValueError("Invalid genome version. Use 19 for hg19 or 38 for hg38.")
        
        run_command(f"python recodeAA.py {ancestral_file} {plink_file}_chr{i}.vcf.gz")

def phase_genotypes(plink_file, beagle_jar, genetic_map_dir):
    """Step 3: Phase genotypes using Beagle"""
    for i in range(1, 23):
        genetic_map = f"{genetic_map_dir}/plink.chr{i}.GRCh37.map"  # Use genetic map directory
        output = f"{plink_file}_chr{i}_phased_beagle"
        run_command(f"java -jar {beagle_jar} gt={plink_file}_chr{i}_recodedAA_recodedAA.vcf.gz map={genetic_map} out={output} chrom={i}")
        run_command(f"bcftools index {output}.vcf.gz")

def split_population(sample_file, plink_file):
    """Step 4: Split VCF by population based on the sample.txt file"""
    sample_file_with_ext = f"{sample_file}.txt"  # Add .txt extension
    for i in range(1, 23):
        run_command(f"bcftools view --samples-file {sample_file_with_ext} {plink_file}_chr{i}_phased_beagle.vcf.gz -Oz -o {sample_file}_chr{i}_phased.vcf.gz")

def create_map_files(sample_file, plink_file, genetic_map_dir_shapeit):
    """Step 5: Create MAP files and sync with genetic map"""
    sample_file_with_ext = f"{sample_file}.txt"  # Add .txt extension
    for i in range(1, 23):
        # Extract info for BIM
        run_command(f"bcftools view {sample_file}_chr{i}_phased.vcf.gz -H | awk '{{print $1, $3, \"0\", $2, $4, $5}}' > {sample_file}_chr{i}_phased.bim")
        # Sync BIM with genetic map
        genetic_map = f"{genetic_map_dir_shapeit}/genetic_map_chr{i}_combined_b37.txt"  # Use genetic map directory
        run_command(f"plink --bim {sample_file}_chr{i}_phased.bim --cm-map {genetic_map} {i} --make-just-bim --out chr{i}_gm")
        # Create MAP file
        run_command(f"awk '{{print $1, $2, $3, $4}}' chr{i}_gm.bim > chr{i}_gm.map")

def run_selscan(sample_file, plink_file, output_dir):
    """Step 6: Run selscan for iHS"""
    sample_file_with_ext = f"{sample_file}.txt"  # Add .txt extension
    for i in range(1, 23):
        run_command(f"selscan --ihs --vcf {sample_file}_chr{i}_phased.vcf.gz --map chr{i}_gm.map --max-extend 1000000 --max-gap 500000 --gap-scale 50000 --cutoff 0.05 --alt --out {output_dir}/{sample_file}_chr{i}_iHS")
        run_command(f"norm --ihs --bins 100 --files {output_dir}/{sample_file}_chr{i}_iHS.ihs.alt.out")
        run_command(f"awk 'BEGIN{{OFS=\",\"; FS=\"\\t\"}} {{for(i=1; i<=NF; i++) {{printf \"\\\"%s\\\"\", $i; if(i<NF) printf \",\"}}; printf \"\\n\"}}' {output_dir}/{sample_file}_chr{i}_iHS.ihs.alt.out.100bins.norm > {output_dir}/{sample_file}_chr{i}_iHS.csv")

def concatenate_csvs(sample_file, output_dir):
    """Step 7: Concatenate all CSV files"""
    run_command(f"cat {output_dir}/{sample_file}_chr{{1..22}}_iHS.csv > {output_dir}/{sample_file}_allChr_iHS.csv")

def cleanup_intermediate_files(plink_file, sample_file):
    """Cleanup intermediate files"""
    # Remove VCF files, BIM files, and other intermediate files
    for i in range(1, 23):
        # Remove VCFs and GZ files
        os.remove(f"{plink_file}_chr{i}.vcf.gz")
        os.remove(f"{plink_file}_chr{i}_phased.vcf.gz")
        os.remove(f"{sample_file}_chr{i}_phased.bim")
        os.remove(f"chr{i}_gm.bim")
        os.remove(f"chr{i}_gm.map")
        # If needed, also remove other temporary files
        # os.remove(f"{plink_file}_chr{i}_recodedAA_recodedAA.vcf.gz")
        # os.remove(f"{plink_file}_chr{i}_phased_beagle.vcf.gz")
        # os.remove(f"{plink_file}_chr{i}_recodedAA.vcf.gz")
    
if __name__ == "__main__":
    import argparse

    # Argument parser
    parser = argparse.ArgumentParser(description="Selection Scan iHS Workflow")
    parser.add_argument("--plink_file", required=True, help="Base PLINK file prefix")
    parser.add_argument("--sample_file", required=True, help="Sample FID_IID file for population of interest (without .txt extension)")
    parser.add_argument("--ancestral_allele_file", required=True, help="Directory containing ancestral allele files")
    parser.add_argument("--beagle_jar", required=True, help="Path to the Beagle JAR file")
    parser.add_argument("--genetic_map_dir", required=True, help="Directory containing genetic map files for phasing")
    parser.add_argument("--genetic_map_dir_shapeit", required=True, help="Directory containing genetic map files for creating MAP files")
    parser.add_argument("--genome_version", type=int, choices=[19, 38], required=True, help="Genome version: 19 for hg19, 38 for hg38")
    args = parser.parse_args()

    # Create directory for output files
    output_dir = args.sample_file
    os.makedirs(output_dir, exist_ok=True)

    # Step-wise execution
    split_chromosomes(args.plink_file)
    polarize_ancestral_alleles(args.plink_file, args.ancestral_allele_file, args.genome_version)
    phase_genotypes(args.plink_file, args.beagle_jar, args.genetic_map_dir)
    split_population(args.sample_file, args.plink_file)
    create_map_files(args.sample_file, args.plink_file, args.genetic_map_dir_shapeit)
    run_selscan(args.sample_file, args.plink_file, output_dir)
    concatenate_csvs(args.sample_file, output_dir)

    # Cleanup intermediate files
    #cleanup_intermediate_files(args.plink_file, args.sample_file)
