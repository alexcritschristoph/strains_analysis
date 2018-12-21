import csv
import sys
import glob
import argparse
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict

def characterize_snp(gene_fasta, freq_file):
    ''' tests if SNPs are synonymous or non-synynomous with respect to a reference
    gene file
    '''
    print("analyzing fasta")
    gene_index = {}
    gene_starts = {}
    seqs = {}
    for record in SeqIO.parse(gene_fasta, 'fasta'):
        gene = str(record.id)
        gene_start = int(record.description.split("#")[1].strip())
        gene_end = int(record.description.split("#")[2].strip())
        gene_starts[gene] = gene_start
        seqs[gene] = record.seq
        for i in range(gene_start, gene_end+1):
            gene_index[gene_scaf + "_" + str(i)] = feature['feature']

    print("assigning snps to genes")
    freq = pd.read_table(freq_file)

    for snp in freq.itterrows():
        #get gene for this snp
        if snp['SNV'].split(":")[0] in gene_index:
            gene = gene_index[snp]
            gene_absolute_start = gene_starts['start']

            #calculate position of SNP within gene (0 based index)
            snp_start = int(snp[position].split("_")[-1]) - gene_absolute_start

            original_sequence = seqs[gene]
            new_sequence = original_sequence
            new_sequence[snp_start] = snp['SNV'].split(":")[1].strip()

            if new_sequence[snp_start] != original_sequence[snp_start]:
                old_aa_sequence = original_sequence.translate()
                new_aa_sequence = new_sequence.translate()
                mut = 'S:' + str(snp_start)
                for aa in range(0, len(old_aa_sequence)+1):
                    if new_aa_sequence[aa] != old_aa_sequence[aa]:
                        mut = 'N:' + str(old_aa_sequence[aa]) + str(snp_start) + str(new_aa_sequence[aa])
                        break
            else:
                mut = 'C'
        else:
            mut = 'I'
        freq['mutation'] = mut

    freq.head()


def create_gene_index(gene_file):
    ## Read in gene file
    Gdb = pd.read_table(gene_file, \
        names=['feature', 'scaffold', 'feature_number', 'scaff_length', 'scaff_GC',\
              'scaff_cov', 'idk', 'taxonomy', 'start', 'end', 'cu', 'anno1', 'anno2', 'anno3'])

    # Create index of scaf+positions to genes
    gene_index = {}
    gene_lengths = {}
    gene_starts = {}
    print("Building gene index")
    for index, feature in Gdb.iterrows():
        gene_scaf = feature['scaffold']
        gene_start = feature['start']
        gene_end = feature['end']
        gene_lengths[feature['feature']] = int(gene_end) - int(gene_start)
        gene_starts[feature['feature']] = int(gene_start)
        for i in range(gene_start, gene_end+1):
            gene_index[gene_scaf + "_" + str(i)] = feature['feature']

    return gene_index, gene_lengths, gene_starts

def get_gene_index(x, gene_index):
    try:
        return(gene_index[x])
    except:
        return None

def calculate_gene_clonality(gene_index, clonality_table, sample, min_cov = 10, min_gene_size = 120):
    '''
    gene_index = a dictionary of scafs and positions to gene IDs calculated by create_gene_index()
    clonality_table = filepath to a clonality .clonal table produced by main program
    min_cov = minimum coverage required for a position to calculate clonality
    min_gene_size = minimum number of positions > min_cov in a gene to calculate clonality for that gene.
    '''

    # Read in clonalities
    cl = pd.read_table(clonality_table)
    print("Calculating gene clonalities for " + clonality_table)

    # Create gene column
    cl['gene'] = cl['position'].apply(lambda x: get_gene_index(x, gene_index) )

    # Filter positions to at least min_cov, only positions within genes
    cl_filtered = cl[cl['coverage'] > min_cov]
    cl_filtered = cl_filtered[cl_filtered['gene'] != None]

    # calculate means of clonality and coverage per gene
    gene_table = cl_filtered.groupby('gene').agg({'clonality':'mean', 'coverage':'mean','position':'count'}).reset_index()
    
    # remove genes with fewer than min_gene_size positions
    gene_table = gene_table[gene_table['position'] >= min_gene_size]

    # add this sample name to the sample column
    gene_table['sample'] = sample
   
    return gene_table


def collate_gene_clonalities(dir, prefix, suffix, gene_file):
    '''
    Runs calculate_gene_clonality() on several .clonal files (has to be the SAME genome across different samples) 
    E.g., for files of the name 'genome:sample_0.98.clonal', prefix = 'genome', suffix = '0.98'. sample names will be read from 'sample'.
    DEFAULT NAMING SCHEME: genome:_sample_suffix.clonal
    ''' 
    gene_index, gene_lengths, gene_starts = create_gene_index(gene_file)
    gene_avg_clonality = defaultdict(dict)
    gene_avg_coverage = defaultdict(dict)

    i = 0
    for fn in glob.glob(dir.rstrip("/") + "/" + prefix + "*" + suffix + ".clonal"):
        sample = fn.split(":_")[1].replace("_" + suffix + ".clonal", "")  ### NOTE: CHANGE THIS LINE DEPENDING ON YOUR FILE NAMING SCHEME
        if i == 0:
            gene_table = calculate_gene_clonality(gene_index, fn, sample)
            i += 1
        else:
            gene_table = pd.concat([gene_table, calculate_gene_clonality(gene_index, fn, sample)])

    # write out gene table
    gene_table.to_csv(prefix + "_" + suffix + ".genes",sep='\t', quoting=csv.QUOTE_NONE)

if __name__ == '__main__':



    parser = argparse.ArgumentParser(description= """
        A script that runs two different gene-based analyses on output of the strain data calculation script:\n

        gene_clonality: This function calculates the average clonality of each gene in several samples mapped to the same reference genome. \n
        usage: python gene_analysis.py gene_clonality -d dir_of_files/ -p genome -s suffix -g gene_file.ql\n
        Where gene_file.ql is a 14 column tsv where each line is a gene, where the 1st column is the gene name,\n the 2nd is the scaffold name, the 9th is the gene's start position and the 10th is the gene's end position.

        \n\n

        gene_snps: This function determines whether SNPs are synonymous or non-synonymous with respect to a reference .fna gene file from prodigal\n
        usage: python gene_analysis.py gene_snps -f prodigal_fna.fna -n snps_file.freq
        """, formatter_class=argparse.RawTextHelpFormatter)

    # Required positional arguments
    parser.add_argument("analysis", help="which gene-level analysis to run: options: gene_clonality (for avg gene clonality) or gene_snps (to label snps as dn/ds)")

    # Argumemtns for gene_clonality
    parser.add_argument("-d", "--dir", action="store", default=None, \
        help='Directory of .clonal files')
    parser.add_argument("-p", "--prefix", action="store", default=None, \
        help='Prefix of genome .clonal file names - probably the name of your genome in these files. We assume a file naming structure similar to genome:_sample_suffix.clonal')
    parser.add_argument("-s", "--suffix", action="store", default=None, \
        help='Suffix of genome .clonal file names - usually the PID level at which the main script was run - such as 0.98, 0.96.')
    parser.add_argument("-g", "--gene_file", action="store", default=None, \
        help='Absolute minimum number of reads to confirm a SNV (both this AND -f must be true)')


    # Arguments for gene_snps

    parser.add_argument("-f", "--fasta", action="store", default=None, \
        help='Path to prodigal .fna fasta file of predicted genes (dna sequences of genes, not contigs or AA)')

    parser.add_argument("-n", "--freq", action="store", default=None, \
        help='Path to the .freq file from the main strain data calculation script.')

    args = parser.parse_args()

    if args.analysis == 'gene_clonality':
        if args.dir and args.prefix and args.suffix and args.gene_file:
            # usage: python gene_clonality.py files_directory prefix suffix gene_file
            collate_gene_clonalities(args.dir, args.prefix, args.suffix, args.gene_file )
        else:
            sys.exit("Error: You did not specific all four required arguments to calculate gene clonality")

    if args.analysis == 'gene_snps':
        if args.fasta and args.freq:
            characterize_snp(args.fasta, args.freq)
        else:
            sys.exit("Error: need to specific a fasta file and a frequency file")
