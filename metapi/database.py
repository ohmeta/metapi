#!/usr/bin/env/python
import os
import pandas as pd


# parse NCBI assembly report, include genbank and refseq, and check the path exist on BGI cluster

def check_ncbi_assembly_report_on_bgi(asm_summary, database, domain, check=True, low_memory=False):
    asm_summary = pd.read_csv(asm_summary, sep='\t', skiprows=[0], low_memory=False)\
                    .astype({"ftp_path": str})\
                    .assign(database=database, domain=domain)
    asm_summary["fna_path_ncbi"] = asm_summary.apply(
        lambda x: os.path.join(x["ftp_path"], os.path.basename(x["ftp_path"]) + "_genomic.fna.gz"), axis=1)
    
    asm_summary["fna_path_bgi"] = asm_summary.apply(
        lambda x: x["fna_path_ncbi"].replace("ftp://ftp.ncbi.nlm.nih.gov", "/hwfssz1/pub/database/ftp.ncbi.nih.gov"), axis=1)
    
    if check:
        asm_summary["fna_path_bgi_exists"] = asm_summary.apply(
            lambda x: 1 if os.path.exists(x["fna_path_bgi"]) else 0, axis=1)
    return asm_summary


'''
mkdir -p ASSEMBLY_REPORTS

#wget -c -P ASSEMBLY_REPORTS https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ANI_report_bacteria.txt
wget -c -P ASSEMBLY_REPORTS https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ANI_report_prokaryotes.txt
wget -c -P ASSEMBLY_REPORTS https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_ANI_report_prokaryotes.txt
wget -c -P ASSEMBLY_REPORTS https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_assembly_summary.txt

wget -c -P ASSEMBLY_REPORTS https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_change_notice.txt
wget -c -P ASSEMBLY_REPORTS https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/prokaryote_type_strain_report.txt
wget -c -P ASSEMBLY_REPORTS https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/species_genome_size.txt.gz

wget -c -P ASSEMBLY_REPORTS https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
wget -c -P ASSEMBLY_REPORTS https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank_historical.txt
wget -c -P ASSEMBLY_REPORTS https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
wget -c -P ASSEMBLY_REPORTS https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq_historical.txt

mkdir -p refseq/{archaea,bacteria,fungi,viral}

wget -c -P refseq https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
wget -c -P refseq https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq_historical.txt

wget -c -P refseq/archaea https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt
wget -c -P refseq/archaea https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary_historical.txt

wget -c -P refseq/bacteria https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
wget -c -P refseq/bacteria https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary_historical.txt

wget -c -P refseq/fungi https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt
wget -c -P refseq/fungi https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary_historical.txt

wget -c -P refseq/viral https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
wget -c -P refseq/viral https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary_historical.txt

mkdir -p genbank/{archaea,bacteria,fungi,viral}

wget -c -P genbank https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
wget -c -P genbank https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank_historical.txt

wget -c -P genbank/archaea https://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt
wget -c -P genbank/archaea https://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary_historical.txt

wget -c -P genbank/bacteria https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
wget -c -P genbank/bacteria https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary_historical.txt

wget -c -P genbank/fungi https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt
wget -c -P genbank/fungi https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary_historical.txt

wget -c -P genbank/viral https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt
wget -c -P genbank/viral https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary_historical.txt
'''


asm_summary_genbank_archaea = check_ncbi_assembly_report_on_bgi("genbank/archaea/assembly_summary.txt", "genbank", "archaea")


asm_summary_genbank_bacteria = check_ncbi_assembly_report_on_bgi("genbank/bacteria/assembly_summary.txt", "genbank", "bacteria")


asm_summary_genbank_fungi = check_ncbi_assembly_report_on_bgi("genbank/fungi/assembly_summary.txt", "genbank", "fungi")


asm_summary_genbank_viral = check_ncbi_assembly_report_on_bgi("genbank/viral/assembly_summary.txt", "genbank", "viral")


asm_summary_refseq_archaea = check_ncbi_assembly_report_on_bgi("refseq/archaea/assembly_summary.txt", "refseq", "archaea")


asm_summary_refseq_bacteria = check_ncbi_assembly_report_on_bgi("refseq/bacteria/assembly_summary.txt", "refseq", "bacteria")


asm_summary_refseq_fungi = check_ncbi_assembly_report_on_bgi("refseq/fungi/assembly_summary.txt", "refseq", "fungi")


asm_summary_refseq_viral = check_ncbi_assembly_report_on_bgi("refseq/viral/assembly_summary.txt", "refseq", "viral")