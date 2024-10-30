#!/usr/bin/env python3

# This script by Obed A. Garcia recodes your VCF file by ancestral allele. Github/Twitter: obedaram
# The Ancestral allele files provided were queried from the 1000 genomes AA field in the VCF files. These sites are all biallelic SNPs and only contains sites that were determined that it can be successfully flipped, i.e. the AA field matches either the REF or the ALT.
# This script will automatically filter and only keep those sites intersecting the Ancestral Allele text file and print out a summary of what was done as well as a file of how each allele was treated.
# Any questions please email or message me at oag@umich.edu

import sys
import gzip
import os
import subprocess
from datetime import datetime

def normalize_chromosome(chrom):
    if (chrom.startswith('chr')):
        return chrom[3:]
    return chrom

def get_ancestral_allele_info(file):
    aa_dict = {}
    if file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    with opener(file, 'rt') as f:
        for line in f:
            if line.startswith('CHR'):
                continue
            cols = line.strip().split()
            chrom = normalize_chromosome(cols[0])
            pos = cols[1]
            ref = cols[3]
            alt = cols[4]
            aa = cols[5]
            key = (chrom, pos)
            aa_dict[key] = (ref, alt, aa)
    return aa_dict

def recode_vcf(vcf_file, aa_dict):
    if vcf_file.endswith('.gz'):
        vcf_opener = gzip.open
    else:
        vcf_opener = open

    output_file = vcf_file.replace('.vcf.gz', '_recodedAA.vcf').replace('.vcf', '_recodedAA.vcf')
    status_file = vcf_file.replace('.vcf.gz', '_status.txt').replace('.vcf', '_status.txt')
    summary_file = vcf_file.replace('.vcf.gz', '_summary.txt').replace('.vcf', '_summary.txt')

    total_variants = 0
    total_removed = 0
    total_flipped = 0
    total_same_ref = 0

    with vcf_opener(vcf_file, 'rt') as vcf, open(output_file, 'wt') as out_vcf, open(status_file, 'wt') as out_status:
        for line in vcf:
            if line.startswith('#'):
                out_vcf.write(line)
                continue
            total_variants += 1
            cols = line.strip().split('\t')
            chrom = normalize_chromosome(cols[0])
            pos = cols[1]
            ref = cols[3]
            alt = cols[4]
            key = (chrom, pos)
            if key not in aa_dict:
                out_status.write(f"{chrom}\t{pos}\t{ref}\t{alt}\tRemoved: No AA data\n")
                total_removed += 1
                continue
            aa_ref, aa_alt, aa = aa_dict[key]
            if ref == aa:
                out_vcf.write(line)
                out_status.write(f"{chrom}\t{pos}\t{ref}\t{alt}\tKept: Matches AA\n")
                total_same_ref += 1
            elif alt == aa:
                flipped_genotypes = flip_genotypes(cols[9:])
                out_vcf.write(f"{cols[0]}\t{cols[1]}\t{cols[2]}\t{alt}\t{ref}\t" + "\t".join(cols[5:9]) + "\t" + "\t".join(flipped_genotypes) + "\n")
                out_status.write(f"{chrom}\t{pos}\t{ref}\t{alt}\tFlipped\n")
                total_flipped += 1
            else:
                out_status.write(f"{chrom}\t{pos}\t{ref}\t{alt}\tRemoved: No AA match\n")
                total_removed += 1

    with open(summary_file, 'w') as summary:
        summary.write(f"Total number variants input: {total_variants}\n")
        summary.write(f"Total Variants Removed: {total_removed}\n")
        summary.write(f"Total Variants Flipped: {total_flipped}\n")
        summary.write(f"Total Variants same REF as AA: {total_same_ref}\n")

    # Compress the output VCF file using bgzip
    subprocess.run(["bgzip", output_file])

    return output_file + '.gz'


def flip_genotypes(genotypes):
    flipped = []
    for gt in genotypes:
        if gt == "./." or gt == ".|.":
            flipped.append(gt)  # Skip flipping for fully missing genotypes
        else:
            if '|' in gt:
                alleles = gt.split('|')
                separator = '|'
            else:
                alleles = gt.split('/')
                separator = '/'
            
            # Handle cases where one allele is missing
            if alleles[0] == '.':
                flipped_gt = f"./{1-int(alleles[1])}"
            elif alleles[1] == '.':
                flipped_gt = f"{1-int(alleles[0])}/."
            else:
                flipped_gt = f"{1-int(alleles[0])}{separator}{1-int(alleles[1])}"
            
            flipped.append(flipped_gt)
    return flipped

def main():
    if len(sys.argv) != 3:
        print("Usage: python recodeAA.py <ancestral_allele.txt> <vcf_file.gz>")
        sys.exit(1)

    aa_file = sys.argv[1]
    vcf_file = sys.argv[2]

    aa_dict = get_ancestral_allele_info(aa_file)
    output_file = recode_vcf(vcf_file, aa_dict)

if __name__ == "__main__":
    main()