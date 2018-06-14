#!/usr/bin/env python3

#    DON'T CALL IT THAT
import pysam
from Bio import SeqIO
from collections import defaultdict
import math
import numpy as np
from tqdm import tqdm
import sys

def entropy2(counts):
    probs = []
    total = sum(counts)
    for c in counts:
        probs.append(float(c) / total)

    ent = 0
    for i in probs:
        if i != 0:
            ent -= i * math.log(i, math.e)
    return ent


def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts



def call_snv_site(counts, min_cov = 5, min_snp = 3):
    '''
    Determines whether a site has a variant based on its nucleotide count frequencies.
    '''
    P2C = {'A':0, 'C':1, 'T':2, 'G':3}
    C2P = {0:'A', 1:'C', 2:'T', 3:'G'}
    total =  sum(counts)
    if total >= min_cov:
        i = 0
        for c in counts:
            if c >= min_snp:
                i += 1
        if i > 1:
            #return consensus nucleotide and 2nd 
            return C2P[counts.index(max(counts))]
    else:
        return False


def _get_base_counts(pileupcolumn, minimum_mapq = 0, pair_mapqs = {}):
    '''
    From a pileupcolumn object, return a list with the counts of [A, C, T, G]
    '''
    P2C = {'A':0, 'C':1, 'T':2, 'G':3}
    C2P = {0:'A', 1:'C', 2:'T', 3:'G'}

    counts = [0,0,0,0]
    empty = [0,0,0,0]

    for pileupread in pileupcolumn.pileups:
        # print(pileupread.)
        if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_qualities[pileupread.query_position] >= 30:
            read_name = pileupread.alignment.query_name
            if pair_mapqs[read_name] >= minimum_mapq:
                try:
                    counts[P2C[pileupread.alignment.query_sequence[pileupread.query_position]]] += 1
                except KeyError: # This would be like an N or something not A/C/T/G
                    pass
    if counts != empty:
        return counts
    else:
        return False

def _get_scaffold_positions(gene_list = None, fasta_file = None):
    ''' Returns a list of windows to record SNVs in'''

    if not fasta_file and not gene_list:
        print("ERROR: REQUIRED TO SUPPLY FASTA or GENE FILE")
        sys.exit(1)
    
    if gene_list:
        f = open(gene_list)
        positions = []
        for line in f.readlines():
            positions.append([line.split(",")[1],int(line.split(",")[2]),int(line.split(",")[3])])
        f.close()
    else:
        positions = []
        for rec in SeqIO.parse(fasta_file, "fasta"):
            positions.append([str(rec.id),1,len(rec.seq)])

    return positions


def run_strain_profiler(bam, positions, min_coverage = 5, min_snp = 3):


    minimum_mapq = 2
    P2C = {'A':0, 'C':1, 'T':2, 'G':3}
    C2P = {0:'A', 1:'C', 2:'T', 3:'G'}
    global mate_pair_mappings

    # Where will we be looking at SNPs?

    #Set up variables
    raw_counts_data = defaultdict(dict) # Set up SNP table
    alpha_snvs = 0
    total_positions = 0
    read_to_snvs = defaultdict(list)
    snvs_to_reads = defaultdict(list)
    snvs_frequencies = defaultdict(int)
    total_read_length = 0
    ## ##################
    ## START READING BAM

    samfile = pysam.AlignmentFile(bam)
    sample = bam.split("/")[-1].split(".bam")[0]

    print("READING BAM: " + bam.split("/")[-1])

    #Get mapping quality for paired reads
    #assumes that paired reads have the same "query name"
    #Is this a good assumption?
    pair_mapqs = defaultdict(int)
    for gene in tqdm(positions, desc='Getting read pairs: '):
        for read in samfile.fetch(gene[0], gene[1], gene[2]):
            if pair_mapqs[read.query_name] < read.mapping_quality:
                pair_mapqs[read.query_name] = read.mapping_quality

    ## Start looping through each gene region
    for gene in tqdm(positions, desc='Finding SNVs ...'):
        scaff = gene[0]
        for pileupcolumn in samfile.pileup(scaff, gene[1], gene[2], stepper = 'nofilter'):
            #is this position an SNV?
            position = scaff + "_" + str(pileupcolumn.pos)
            counts = _get_base_counts(pileupcolumn, minimum_mapq = minimum_mapq, pair_mapqs = pair_mapqs)

            consensus = False
            if counts:
                total_positions += 1
                consensus = call_snv_site(counts, min_cov = min_coverage, min_snp = min_snp)
                total_read_length += sum(counts)

            if consensus:
                #there's an SNV at this site
                # add to SNV frequencies
                for pileupread in pileupcolumn.pileups:
                    read_name = pileupread.alignment.query_name
                    if not pileupread.is_del and not pileupread.is_refskip:
                        if pileupread.alignment.query_qualities[pileupread.query_position] >= 30 and pair_mapqs[read_name] >= minimum_mapq:
                            try:
                                val = pileupread.alignment.query_sequence[pileupread.query_position]
                                #if value is not the consensus value
                                if val != consensus and counts[P2C[val]] > min_snp:
                                    #this is a variant read!
                                    read_to_snvs[read_name].append(position + ":" + val)
                                    snvs_to_reads[position+":"+val].append(read_name)
                                    alpha_snvs += 1

                            except KeyError: # This would be like an N or something not A/C/T/G
                                pass


                freqs = {}
                nucl_count = 0 
                for nucl in counts:
                    if nucl > min_snp:
                        freq = float(nucl) / float(sum(counts))
                        snp = position + ":" + C2P[nucl_count]
                        snvs_frequencies[snp] = freq
                    nucl_count += 1
                    
    print("Total number of reads: " + str(len(read_to_snvs.keys())))
    print("Total SNVs: " + str(alpha_snvs))
    print("Total number of positions looked at: " + str(total_positions))
    return {'alpha_snvs': alpha_snvs,
            'total_read_length': total_read_length,
            'snvs_frequencies':snvs_frequencies,
            'read_to_snvs': read_to_snvs,
            'snvs_to_reads':snvs_to_reads
        }
    # TO DO
    # Link reads


def write_tables(fasta, bam, results):

    # Generate tables
    print("*** PRINTING TABLES ***")
    sample = bam.split("/")[-1].split(".bam")[0]
    genome = fasta.split("/")[-1].split(".")[0]
    #Number of alpha svs per read
    f = open('./' + genome + ".reads", "w+")
    f.write("Genome\tSample\tRead\tCount\n")
    for read in results['read_to_snvs']:
        f.write(genome + "\t" + sample + "\t" + read + "\t" + str(len(results['read_to_snvs'][read])) +"\n")
    f.close()

    f = open('./' + genome + ".freq", "w+")
    f.write("Genome\tSample\tSNV\tFrequency\n")
    for snp in results['snvs_frequencies']:
        f.write(genome + "\t" + sample + "\t" + snp + "\t" + str(results['snvs_frequencies'][snp]) +"\n")
    f.close()

    f = open('./' + genome + ".freq", "w+")
    f.write("Genome\tSNVs\tTotal Read Length\n")
    f.write(genome + "\t" + results['alpha_snvs'] + "\t" + results['total_read_length'])
    f.close()

def main():
    if len(sys.argv) > 1:
        bam = sys.argv[1]
    else:
        print("ERROR: FIRST ARGUMENT NEEDS TO BE THE BAM OR SAM")
        sys.exit(1)

    if len(sys.argv) > 2:
        fasta = sys.argv[2]
    else:
        print("ERROR: SECOND ARGUMENT NEEDS TO BE THE FASTA")
        sys.exit(1)

    myargs = getopts(sys.argv)
    
    if '-g' in myargs:
        genes = myargs['-g']
    else:
        genes = None
    if '-c' in myargs:  
        min_coverage = int(myargs['-c'])
    else:
        min_coverage = 5
    if '-s' in myargs:  
        min_snp = int(myargs['-s'])
    else:
        min_snp = 3
    positions = _get_scaffold_positions(genes, fasta)
    results = run_strain_profiler(bam, positions, min_coverage = min_coverage, min_snp = min_snp)
    
    write_tables(fasta, bam, results)

    # write_tables(genome, results)

if __name__ == '__main__':
    main()
