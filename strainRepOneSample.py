#!/usr/bin/env python

'''
Possible names for this program:
* strainRep
* blockStrain
* destrain
* instrain
'''

# Get the version
from _version import __version__

# Import
import csv
import sys
import math
import pysam
import pickle
import argparse
import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict
from sklearn.decomposition import PCA

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
        self.snv_net = None
        self.graph = None
        self.graph_model = None
        self.output = None
        self.testing = False

    def save(self, size=0):
        if self.results:
            # Generate tables
            if size == 0:
                self.read_to_snvs = None
                self.snvs_frequencies = None
                self.snvs_to_reads = None

            print("Printing to tables.")
            sample = self.bam.split("/")[-1].split(".bam")[0]

            genome = self.output
            print(genome)
            self.snv_table.to_csv(genome + ".freq",sep='\t', quoting=csv.QUOTE_NONE)
            self.reads_table.to_csv(genome + ".reads",sep='\t', quoting=csv.QUOTE_NONE)

            f = open(genome + ".data", "w+")
            f.write("Genome\tSNVs\tTotal Read Length\n")
            f.write(genome + "\t" + str(self.alpha_snvs) + "\t" + str(self.total_read_length))
            f.close()

            f = open(genome + ".data", 'wb')
            pickle.dump(self.__dict__, f, 2)
            f.close()
        else:
            print("No data to save.")

    def load(self, name):
        self.__dict__.clear()
        f = open(name + ".data", 'rb')
        tmp_dict = pickle.load(f)
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

    def calc_linkage_network(self):
        ''' Calculates the SNV linkage network - saves it as a dictionary of edges in self.snv_net. 
        Writes it out to a file, genome.net for reading in through other programs like node2vec
        '''
        snv_net = defaultdict(int)
        print("Calculating SNV linkage network...")

        for snv in self.snvs_to_reads:
            reads = self.snvs_to_reads[snv]
            for read in reads:
                for snv2 in self.read_to_snvs[read]:
                    if snv2 != snv:
                        snv_pair = frozenset([snv, snv2])
                        snv_net[snv_pair] += 1
        self.snv_net = snv_net

        print("There were " + str(len(snv_net.keys())) + " edges in the network")
        print("Writing graph to file...")
        f = open(self.output + '.net', 'w+')
        for edge in snv_net:
            nodes = list(edge)
            f.write(nodes[0] + '\t' + nodes[1] + "\t" + str(snv_net[edge]) + "\n")
        f.close()

        # Write to file


    def calc_network_structure(self, genome):
        '''runs node2vec on a graph'''
        pass


    def plot_graph(self):
        ''' Plots the output of node2vec models as a 2d histogram '''
        if self.model:
            print("Running PCA...")
            my_pca = PCA(n_components=d)
            pca_output = my_pca.fit_transform(model)

        else:
            print("ERROR: NO GRAPH GENERATED.")

    def plot(self, viz_type = None):
        pass

    def calc_graph(self, fasta = None):
        ''' Takes the snv_net object and creates a networkx network object from it.
        Not called automatically, for use in jupyter notebooks to play around with net.
        '''
        if self.snv_net:
            #create networkx graph
            G=nx.Graph()
            #add nodes
            nodes = set()
            for edge in self.snv_net:
                nodes_set = list(edge)
                node1 = nodes_set[0]
                node2 = nodes_set[1]
                if node1 not in nodes:
                    G.add_node(node1)
                    nodes.add(node1)
                if node2 not in nodes:
                    G.add_node(node2) 
                    nodes.add(node2)
                G.add_edge(node1,node2, weight=self.snv_net[edge])
            self.graph = G
        else:
            print("ERROR: NO LINKAGE NET GENERATED")
            return False


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
        non_consensus_snvs = set()
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

        if self.testing:
            self.positions = self.positions[0:10]

        ## Start looping through each gene region            
        for gene in tqdm(self.positions, desc='Finding SNVs ...'):
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
                                                non_consensus_snvs.add(position + ":" + val)

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

        # Final statistics
        print("Total SNVs-sites: " + str(total_snv_sites))
        print("Total SNV-bases: " + str(alpha_snvs))
        print("Total sites: " + str(total_positions))
        print("Total number of bases: " + str(total_read_length))
        print("Non-consensus SNVs: " + str(len(non_consensus_snvs)))

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

def main(args):
    '''
    Main entry point
    '''

    strains = SNVdata()

    if args.testing:
        strains.testing = True

    if not args.output:
        strains.output = args.fasta.split("/")[-1].split(".")[0]
    else:
        strains.output = args.output

    strains.get_scaffold_positions(args.genes, args.fasta)
    strains.run_strain_profiler(args.bam, min_coverage = int(args.min_coverage), min_snp = int(args.min_snp))
    strains.calc_linkage_network()
    # strains.calc_graph(args.output)
    strains.save(args.output)


    # write_tables(genome, results)

if __name__ == '__main__':
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional arguments
    parser.add_argument("bam", help="Sorted .bam file")
    parser.add_argument("fasta", help="Fasta file the bam is mapped to")

    # Optional arguments
    parser.add_argument("-o", "--output", action="store", default=None, \
        help='Output prefix')
    parser.add_argument("-g", "--genes", action="store", default=None, \
        help='Optional genes file')
    parser.add_argument("-c", "--min_coverage", action="store", default=5, \
        help='Minimum SNV coverage')
    parser.add_argument("-s", "--min_snp", action="store", default=3, \
        help='Minimum number of reads to confirm a SNV')

    parser.add_argument('--testing', action='store_true', default=False, \
        help ="Testing command runs only on first 10 contigs to run faster.")

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    # Parse
    args = parser.parse_args()

    main(args)
