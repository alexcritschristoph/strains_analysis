#!/usr/bin/env python3

#    DON'T CALL IT THAT
import pysam
from Bio import SeqIO
from collections import defaultdict
import math
import numpy as np
from tqdm import tqdm
import sys
import pandas as pd
import cPickle
import csv

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

class SNVdata:

    # The class "constructor" - It's actually an initializer 
    def __init__(self):
        self.fasta = ""
        self.bam = ""
        self.results = False
        self.positions = []
        self.snv_table = None
        self.reads_table = None
        self.alpha_snvs = None
        self.read_to_snvs = None
        self.snvs_frequencies = None
        self.snvs_to_reads = None
        self.total_read_length = None
        self.total_positions = None
        self.total_snv_sites = None 
        self.non_consensus_snvs = None

    def save(self):
        if self.results:
            # Generate tables
            print("Printing to tables.")
            sample = self.bam.split("/")[-1].split(".bam")[0]
            genome = self.fasta.split("/")[-1].split(".")[0]
            print(genome)
            self.snv_table.to_csv(genome + ".reads",sep='\t', quoting=csv.QUOTE_NONE)
            self.reads_table.to_csv(genome + ".freq",sep='\t', quoting=csv.QUOTE_NONE)

            f = open('./' + genome + ".data", "w+")
            f.write("Genome\tSNVs\tTotal Read Length\n")
            f.write(genome + "\t" + str(self.alpha_snvs) + "\t" + str(self.total_read_length))
            f.close()

            f = open(genome + ".data", 'wb')
            cPickle.dump(self.__dict__, f, 2)
            f.close()
        else:
            print("No data to save.")

    def load(self, name):
        self.__dict__.clear()
        f = open(name + ".data", 'rb')
        tmp_dict = cPickle.load(f)
        f.close()          

        self.__dict__.update(tmp_dict) 

    def get_scaffold_positions(self, gene_list = None, fasta_file = None):
        ''' Returns a list of windows to record SNVs in'''
        if not fasta_file and not gene_list:
            print("ERROR: REQUIRED TO SUPPLY FASTA or GENE FILE")
            sys.exit(1)
        
        if gene_list:
            f = open(gene_list)
            for line in f.readlines():
                self.positions.append([line.split(",")[1],int(line.split(",")[2]),int(line.split(",")[3])])
            f.close()
        else:
            for rec in SeqIO.parse(fasta_file, "fasta"):
                self.positions.append([str(rec.id),1,len(rec.seq)])
        self.fasta = fasta_file



    def run_strain_profiler(self, bam, min_coverage = 5, min_snp = 3):


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
        total_snv_sites = 0
        non_consensus_snvs = []
        ## ##################
        ## START READING BAM

        samfile = pysam.AlignmentFile(bam)
        sample = bam.split("/")[-1].split(".bam")[0]

        print("READING BAM: " + bam.split("/")[-1])

        #Get mapping quality for paired reads
        #assumes that paired reads have the same "query name"
        #Is this a good assumption?
        pair_mapqs = defaultdict(int)
        for gene in tqdm(self.positions, desc='Getting read pairs: '):
            for read in samfile.fetch(gene[0], gene[1], gene[2]):
                if pair_mapqs[read.query_name] < read.mapping_quality:
                    pair_mapqs[read.query_name] = read.mapping_quality

        ## Start looping through each gene region
        for gene in tqdm(self.positions[0:10], desc='Finding SNVs ...'):
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
                    total_snv_sites += 1

                    # add to SNV frequencies
                    for pileupread in pileupcolumn.pileups:
                        read_name = pileupread.alignment.query_name
                        if not pileupread.is_del and not pileupread.is_refskip:
                            if pileupread.alignment.query_qualities[pileupread.query_position] >= 30 and pair_mapqs[read_name] >= minimum_mapq:
                                try:
                                    val = pileupread.alignment.query_sequence[pileupread.query_position]
                                    #if value is not the consensus value
                                    if counts[P2C[val]] >= min_snp:
                                        #this is a variant read!
                                        read_to_snvs[read_name].append(position + ":" + val)
                                        snvs_to_reads[position+":"+val].append(read_name)
                                        if val != consensus:
                                            alpha_snvs += 1
                                            if position + ":" + val not in non_consensus_snvs:
                                                non_consensus_snvs.append(position + ":" + val)

                                except KeyError: # This would be like an N or something not A/C/T/G
                                    pass


                    freqs = {}
                    nucl_count = 0 
                    for nucl in counts:
                        if nucl >= min_snp:
                            freq = float(nucl) / float(sum(counts))
                            snp = position + ":" + C2P[nucl_count]
                            snvs_frequencies[snp] = freq
                        nucl_count += 1
                        
        #Calculate SNP per read, Frequencies intersection 
        print("Calculating frequency - SNVs per read intersection...")
        snv_table = defaultdict(list)
        for snv in snvs_frequencies:
            reads = snvs_to_reads[snv]
            if len(reads) != 0:
                snvs_per_read_mean = 0.0
                for read in reads:
                    snvs_per_read_mean += len(read_to_snvs[read])
                snvs_per_read_mean = snvs_per_read_mean / len(reads)
                snv_table['SNV'].append(snv)
                snv_table['freq'].append(snvs_frequencies[snv])
                snv_table['SNVs-per-read'].append(snvs_per_read_mean)
            else:
                print("ERROR AT SNV: " + snv)
        snv_table = pd.DataFrame(snv_table)
        #Convert SNPs per read dictionary to Pandas DF
        reads_table = defaultdict(list)
        for read in read_to_snvs:
            reads_table['read'].append(read)
            reads_table['snvs'].append(str(len(read_to_snvs[read])))
        reads_table = pd.DataFrame(reads_table)


        # Calculate SNP linkage network
        # print("Calculating SNV linkage network...")
        # for snv in snvs_frequencies:
        #     reads = snvs_to_reads[snv]
        #     for read in reads:


        # Final statistics
        print("Total SNVs-sites: " + str(total_snv_sites))
        print("Total SNV-bases: " + str(alpha_snvs))
        print("Total sites: " + str(total_positions))
        print("Total number of bases: " + str(total_read_length))

        self.total_read_length = total_read_length
        self.alpha_snvs = alpha_snvs
        self.total_snv_sites = total_snv_sites

        self.snv_table = snv_table
        self.reads_table = reads_table
        self.read_to_snvs = read_to_snvs
        self.snvs_to_reads = snvs_to_reads
        self.total_positions = total_positions
        self.non_consensus_snvs = non_consensus_snvs

        self.results = True

        test = 0
        for snv in snvs_to_reads:
            # print(len(snvs_to_reads[snv]))
            if snv in non_consensus_snvs:
                test += len(snvs_to_reads[snv])
        print(test)
        print(len(non_consensus_snvs))


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


    strains = SNVdata()

    strains.get_scaffold_positions(genes, fasta)
    strains.run_strain_profiler(bam, min_coverage = min_coverage, min_snp = min_snp)
    strains.save()

    # write_tables(genome, results)

if __name__ == '__main__':
    main()