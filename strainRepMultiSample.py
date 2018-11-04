#!/usr/bin/env python

'''
Possible names for this program:
* strainRep - not a legit acronym
* blockStrain - not a legit acronym
* destrain - not a legit acronym
* instrain - not a legit acronym
* MIPS - metagenome interference of population structure - a legt acronym
'''

# Get the version
from _version import __version__

# Import
import csv
import sys
import glob
import math
import pysam
import pickle
import argparse
import itertools
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from Bio import SeqIO
from scipy import stats
from collections import defaultdict
from sklearn.decomposition import PCA

def major_allele(snv, counts):
    d = {'A': counts[0], 'C': counts[1], 'T': counts[2], 'G': counts[3] }
    nucl = sorted(d, key=d.get, reverse=True)[0]
    return [snv + ":" + nucl, d[nucl]]

def minor_allele(snv, counts):
    d = {'A': counts[0], 'C': counts[1], 'T': counts[2], 'G': counts[3] }
    nucl = sorted(d, key=d.get, reverse=True)[1]
    return [snv + ":" + nucl, d[nucl]]


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


def calculate_clonality(counts, min_cov = 5):
    '''
    Calculates the probability that two reads have the same allele at this position
    '''
    prob = (float(counts[0]) / sum(counts)) * (float(counts[0]) / sum(counts)) + (float(counts[1]) / sum(counts)) * (float(counts[1]) / sum(counts)) + (float(counts[2]) / sum(counts)) * (float(counts[2]) / sum(counts)) + (float(counts[3]) / sum(counts)) * (float(counts[3]) / sum(counts))
    return prob


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


def _get_base_counts(pileupcolumn, filtered_reads):
    '''
    From a pileupcolumn object, return a list with the counts of [A, C, T, G]
    '''
    P2C = {'A':0, 'C':1, 'T':2, 'G':3}
    C2P = {0:'A', 1:'C', 2:'T', 3:'G'}

    counts = [0,0,0,0]
    empty = [0,0,0,0]

    for pileupread in pileupcolumn.pileups:
        # print(pileupread.)
        if not pileupread.is_del and not pileupread.is_refskip:
            read_name = pileupread.alignment.query_name
            if read_name in filtered_reads:
                try:
                    counts[P2C[pileupread.alignment.query_sequence[pileupread.query_position]]] += 1
                except KeyError: # This would be like an N or something not A/C/T/G
                    pass
    if counts != empty:
        return counts
    else:
        return False

########################
#### START SNVData class

class SNVdata:

    # The class "constructor" - It's actually an initializer
    def __init__(self):

        # Parameters
        self.fasta = ""
        self.bam = ""
        self.results = False
        self.output = None
        self.testing = False        
        self.positions = []


        # Data structures
        self.snv_table = None         #
        self.read_to_snvs = None      #
        self.windows_to_snvs = None   # 
        self.snv_counts = None        # Counts of AGCT for each SNV position
        self.clonality_table = None         # Dict of clonality histograms (lists) by window
        self.snv_graph = None         # Weighted networkx graph that tracks SNVs( 100:A, 100:T are diff nodes) that are linked
        self.position_graph = None    # Unweighted networkx graph that tracks positions that are linked
        self.r2linkage_table = None         # Dict of r2 histograms (list of lists) by window

        # General statistics
        self.coverages = None
        self.total_positions = None
        self.total_snv_sites = None
        self.alpha_snvs = None


        

    def save(self, size=0):
        if self.results:
            # Generate tables
            if size == 0:
                self.read_to_snvs = None

            print("Printing to tables.")
            sample = self.bam.split("/")[-1].split(".bam")[0]

            genome = self.output
            print(genome)
            self.snv_table.to_csv(genome + ".freq",sep='\t', quoting=csv.QUOTE_NONE)
            self.clonality_table.to_csv(genome + ".clonal",sep='\t', quoting=csv.QUOTE_NONE)

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
#            print("ERROR: REQUIRED TO SUPPLY FASTA or GENE FILE")
            sys.exit(1)

        if gene_list:
            f = open(gene_list)
            for line in f.readlines():
                self.positions.append([line.split(",")[1],int(line.split(",")[2]),int(line.split(",")[3])])
            f.close()
        else:
            for rec in SeqIO.parse(fasta_file, "fasta"):
                start = 0
                while True:
                    chunk = start + 15000
                    if chunk <= len(rec.seq):
                        self.positions.append([str(rec.id),start,start+10000])
                    else:
                        self.positions.append([str(rec.id),start,len(rec.seq)])
                        break
                    start += 10000
                    start += 1

        self.fasta = fasta_file

    def get_snps_per_gene(self, gene_file = None):
        pass

    def calc_linkage_network(self):
        ''' Calculates the SNV linkage network - saves it as a dictionary of edges in self.snv_net. 
        Writes it out to a file, genome.net for reading in through other programs like node2vec
        '''
        G=nx.Graph()
        G_pos = nx.Graph()

        print("Calculating SNV linkage network...")
        for read in self.read_to_snvs:
            snvs = self.read_to_snvs[read]
            for pair in itertools.combinations(snvs, 2):
                if G.has_edge(pair[0], pair[1]):
                    G[pair[0]][pair[1]]['weight'] += 1
                else:
                    G.add_edge(pair[0], pair[1], weight = 1)
                    if not G_pos.has_edge(pair[0].split(":")[0], pair[1].split(":")[0]):
                        G_pos.add_edge(pair[0].split(":")[0], pair[1].split(":")[0])

        print(str("Of " + str(self.total_snv_sites) + " SNP-sites, there were " + str(len(G_pos)) + " SNPs that could be linked to at least one other SNP."))
        print("The average SNP was linked to " + str(int(np.mean([x[1] for x in list(G_pos.degree())]))) + " other SNPs.")
        self.snv_graph = G
        self.position_graph = G_pos


    def calc_ld(self, snp_a, snp_b, min_snp):
        '''
        A function that calculates the LD between two SNPs
        '''

        #distance between snps
        distance = abs(int(snp_a.split("_")[-1]) - int(snp_b.split("_")[-1]))

        #calculate allele frequencies
        allele_A = major_allele(snp_a,self.snv_counts[snp_a]) #get major allele
        allele_a = minor_allele(snp_a,self.snv_counts[snp_a]) # get minor allele
        freq_A = float(allele_A[1]) / (allele_A[1] + allele_a[1])
        freq_a = float(allele_a[1]) / (allele_A[1] + allele_a[1])

        # get snp B frequencies
        allele_B = major_allele(snp_b,self.snv_counts[snp_b])
        allele_b = minor_allele(snp_b,self.snv_counts[snp_b])
        freq_B = float(allele_B[1]) / (allele_B[1] + allele_b[1])
        freq_b = float(allele_b[1]) / (allele_B[1] + allele_b[1])

        # Get frequencies of linkages
        countAB, countAb, countaB, countab  = 0,0,0,0
        if self.snv_graph.has_edge(allele_A[0], allele_B[0]):
            countAB = self.snv_graph[allele_A[0]][allele_B[0]]['weight']
        if self.snv_graph.has_edge(allele_A[0], allele_b[0]):
            countAb = self.snv_graph[allele_A[0]][allele_b[0]]['weight']
        if self.snv_graph.has_edge(allele_a[0], allele_B[0]):
            countaB = self.snv_graph[allele_a[0]][allele_B[0]]['weight']
        if self.snv_graph.has_edge(allele_a[0], allele_b[0]):
            countab = self.snv_graph[allele_a[0]][allele_b[0]]['weight']

        total = countAB + countAb + countaB + countab
        
        #Requires at least 10 linkages
        if total > min_snp and total > 10:

            linkage_points_x = []
            linkage_points_y = []
            for point in range(0,countAB):
                linkage_points_x.append(1)
                linkage_points_y.append(1)
            for point in range(0,countAb):
                linkage_points_x.append(1)
                linkage_points_y.append(0)
            for point in range(0,countaB):
                linkage_points_x.append(0)
                linkage_points_y.append(1)
            for point in range(0,countab):
                linkage_points_x.append(0)
                linkage_points_y.append(0)

            freq_AB = float(countAB) / total
            freq_Ab = float(countAb) / total
            freq_aB = float(countaB) / total
            freq_ab = float(countab) / total

            # Calculate linkage_D and linkage_d (major and minor alleles)
            linkD = freq_AB - freq_A * freq_B
            linkd = freq_ab - freq_a * freq_b
            
            # Calculate r^2
            r2 = stats.pearsonr(linkage_points_x, linkage_points_y)[0]
            if str(r2) != 'nan':
                # print("**** DEBUG ****")
                # print("freq_A " + str(freq_A))
                # print("freq_a " + str(freq_a))
                # print("freq_B " + str(freq_B))
                # print("freq_b " + str(freq_b))

                # print("countAB " + str(countAB))
                # print("countAb " + str(countAb))
                # print("countaB " + str(countaB))
                # print("countab " + str(countab))
                # print("total " + str(total))

                # print("linkD " +str(linkD))
                # print("linkd " + str(linkd))
                # print("r2 " + str(r2))
                # print(linkage_points_x)
                # print(linkage_points_y)
                return([distance, r2, linkD, linkd])
        else:
            return False

    def calc_ld_all_sites(self, min_snp):
        '''
        Calculates Linkage Disequilibrium for all SNVs in a window.
        '''
        r2_total = {}
        sigma_total = {}
        print("Calculating Linkage Disequilibrium")

        for window in self.positions:
            r2_spectrum = []
            sigma_spectrum = []

            window_name = window[0] + ":" + str(window[1]) + ":" + str(window[2])
            window_snvs = self.windows_to_snvs[window_name]

            for edge in self.position_graph.edges(window_snvs):
                snp_a = edge[0]
                snp_b = edge[1]
                ld_result = self.calc_ld(snp_a, snp_b, min_snp)
                if ld_result:
                    r2_spectrum.append(ld_result)
            
            r2_total[window_name] = r2_spectrum
        
        # Create r2 linkage table
        r2linkage_table = defaultdict(list)

        for window in r2_total:
            for datum in r2_total[window]:
                r2linkage_table['Window'].append(window)
                r2linkage_table['Distance'].append(datum[0])
                r2linkage_table['r2'].append(datum[1])
                r2linkage_table['linkD'].append(datum[2])
                r2linkage_table['linkd'].append(datum[3])

        self.r2linkage_table = pd.DataFrame(r2linkage_table)




    def calc_network_structure(self, genome):
        '''runs node2vec on a graph'''
        pass

    def plot(self, viz_type = None):
        pass


    def run_strain_profiler(self, bams, min_coverage = 5, min_snp = 3, filter_cutoff = 0):
        ''' 
        Main class for finding SNVs and generating data profile for a genome.
        '''

        minimum_mapq = 2
        P2C = {'A':0, 'C':1, 'T':2, 'G':3}
        C2P = {0:'A', 1:'C', 2:'T', 3:'G'}
        global mate_pair_mappings


        ### *************************************
        ### STEP 0: SET UP VARIABLES
        ### *************************************
        raw_counts_data = defaultdict(dict) # Set up SNP table
        alpha_snvs = 0
        total_positions = 0
        read_to_snvs = defaultdict(list)
        snvs_frequencies = defaultdict(dict)
        clonality_by_window = defaultdict(dict)
        windows_to_snvs = defaultdict(list)
        snv_counts = {}

        coverages = defaultdict(dict)
        total_snv_sites = 0

        # FOR TESTING: calculate insert sizes

        if self.testing:
            self.positions = self.positions[0:10]


        counts_sum_by_bam = {}
        counts_data_by_bam = defaultdict(dict)
        filtered_reads_by_bam = {}
        ### *************************************
        ### Start reading BAMs to get read counts
        ### *************************************
        for bam in bams:
            samfile = pysam.AlignmentFile(bam)
            sample = bam.split("/")[-1].split(".bam")[0]

            print("READING BAM: " + bam.split("/")[-1])

            print("Using reads with >" + str(filter_cutoff) + "% PID to consensus reference.")


            ### *****************************************************
            ### STEP 1: GET READ INFORMATION AND FILTER READS
            ### *****************************************************

            #Get mapping quality for paired reads
            #assumes that paired reads have the same "query name"
            pair_mapqs = defaultdict(int)
            insert_sizes = defaultdict(int)
            insert_sizes_r1 = defaultdict(lambda: -1)
            found_pairs = set()
            observed_read1s = set()
            observed_read2s = set()

            # mismatches per read / read pair work
            subset_reads = set()
            read_pair_mismatches = {} 
            read_pair_pid = {}

            total_read_count = 0
            for gene in tqdm(self.positions, desc='Getting read pairs: '):
                
                for read in samfile.fetch(gene[0], gene[1], gene[2]):
                    total_read_count += 1

                    #second read in pair
                    if (read.is_read2 and read.query_name in observed_read1s) or (read.is_read1 and read.query_name in observed_read2s):
                        if read.query_name not in found_pairs:
                            if read.get_reference_positions() != [] and insert_sizes_r1[read.query_name] != -1:
                                found_pairs.add(read.query_name)
                                read_pair_pid[read.query_name] = 1-(float(read_pair_mismatches[read.query_name][0]) + float(read.get_tag('NM'))) / ( float(read_pair_mismatches[read.query_name][1]) + read.infer_query_length())
                                insert_sizes[read.query_name] = read.get_reference_positions()[-1] - insert_sizes_r1[read.query_name]
                    
                    #this is the first read in a pair
                    else:


                        if read.get_reference_positions() != []:
                                insert_sizes_r1[read.query_name] = read.get_reference_positions()[0]
                                if read.is_read1:
                                    observed_read1s.add(read.query_name)
                                else:
                                    observed_read2s.add(read.query_name)

                                read_pair_mismatches[read.query_name] = [read.get_tag('NM'), read.infer_query_length()]

                                if pair_mapqs[read.query_name] < read.mapping_quality:
                                    pair_mapqs[read.query_name] = read.mapping_quality


            min_insert = 50 # paired reads must be 50 bp apart
            max_insert = np.median(list(insert_sizes.values())) * 2 # they can't be more than 2 * apart as the average


            too_short = 0
            too_long = 0
            good_length = 0
            for read_pair in found_pairs:
                if insert_sizes[read_pair] > min_insert:
                    if insert_sizes[read_pair] < max_insert:
                        if pair_mapqs[read_pair] > minimum_mapq:
                            good_length += 1

                            # Which set does this read go into?
                            if read_pair_pid[read_pair] > filter_cutoff:
                                subset_reads.add(read_pair)
                    else:
                        too_long += 2
                else:
                    too_short += 2

            print("median insert size: " + str(max_insert / 2))
            print("total reads found: " + str(total_read_count))
            print("reads with pair found: " + str(len(found_pairs) * 2))
            print("paired reads < 50 bp apart: " + str(too_short))
            print("paired reads > " + str(max_insert) + " apart: " + str(too_long))
            print("reads which pass pair insert size filter: " + str(good_length*2))
            print("reads which pass read pair PID >" + str(filter_cutoff) + "%: " + str(len(subset_reads)*2))

            filtered_reads_by_bam[sample] = subset_reads

            ### ***********************************
            ### STEP 2: SCAN BAMs for counts columns
            ### ***********************************

            ## Start looping through each region  
            for gene in tqdm(self.positions, desc='Getting all counts ...'):
                scaff = gene[0]
                window = gene[0] + ":" + str(gene[1]) + ":" + str(gene[2])
                for pileupcolumn in samfile.pileup(scaff, gene[1], gene[2], truncate = True, stepper = 'samtools', compute_baq= True, ignore_orphans = True, ignore_overlaps = True,  min_base_quality = 30):
                    ## Step 1: Are there any reads at this position?
                   
                    position = scaff + "_" + str(pileupcolumn.pos)
                    counts = _get_base_counts(pileupcolumn, filtered_reads = subset_reads)
                    if counts:
                        if position in counts_sum_by_bam:
                            # print(counts_sum_by_bam[position])
                            counts_sum_by_bam[position] = [sum(x) for x in zip(counts, counts_sum_by_bam[position])]
                        else:
                            counts_sum_by_bam[position] = counts
                    counts_data_by_bam[position][sample] = counts

        ### ****************************
        ### STEP 3: FINISHED first BAM loop, loop again and CALL SNVs
        ### ****************************        
        for bam in bams:
            samfile = pysam.AlignmentFile(bam)
            sample = bam.split("/")[-1].split(".bam")[0]

            print("READING BAM: " + bam.split("/")[-1])

            for gene in tqdm(self.positions, desc='Calling SNVs ...'):
                scaff = gene[0]
                window = gene[0] + ":" + str(gene[1]) + ":" + str(gene[2]) 
                for pileupcolumn in samfile.pileup(scaff, gene[1], gene[2], truncate = True, stepper = 'samtools', compute_baq= True, ignore_orphans = True, ignore_overlaps = True,  min_base_quality = 30):
                    ## Step 1: Are there any reads at this position?
                   
                    position = scaff + "_" + str(pileupcolumn.pos)
                    counts = counts_data_by_bam[position][sample]

                    consensus = False
                    # Yes there were reads at this position
                    if counts:
                        total_positions += 1
                        if sum(counts) > min_coverage:
                            pos_clonality = calculate_clonality(counts)
                            if sample in clonality_by_window[window]:
                                clonality_by_window[window][sample].append([position, pos_clonality])
                            else:
                                clonality_by_window[window][sample] = [[position, pos_clonality]]
                            consensus = call_snv_site(counts_sum_by_bam[position], min_cov = min_coverage, min_snp = min_snp)
                            coverages[position][sample] = sum(counts)


                    ## Strep 2: Is there an SNV at this position?
                    if consensus:
                        #there's an SNV at this site
                        total_snv_sites += 1
                        windows_to_snvs[window].append(position)
                        snv_counts[position] = counts_sum_by_bam[position]
                        # Get reads to snvs

                        for pileupread in pileupcolumn.pileups:
                            read_name = pileupread.alignment.query_name
                            if not pileupread.is_del and not pileupread.is_refskip:
                                if read_name in subset_reads:
                                    try:
                                        val = pileupread.alignment.query_sequence[pileupread.query_position]
                                        #if value is not the consensus value
                                        if counts[P2C[val]] >= min_snp:
                                            #this is a variant read!
                                            read_to_snvs[read_name].append(position + ":" + val)
                                            if val != consensus:
                                                alpha_snvs += 1
                                    except KeyError: # This would be like an N or something not A/C/T/G
                                        pass

                        # Add to frequencies
                        nucl_count = 0
                        for nucl in counts:
                            if nucl >= min_snp:
                                freq = float(nucl) / float(sum(counts))
                                snp = position + ":" + C2P[nucl_count]
                                snvs_frequencies[snp][sample] = [freq, window]

                            nucl_count += 1 

        #### SNP FINDING COMPLETE
        # Create SNV frequency table
        snv_table = defaultdict(list)
        for snv in snvs_frequencies:
            for sample in snvs_frequencies[snv]:
                snv_table['SNV'].append(snv)
                snv_table['sample'].append(sample)
                for bam in bams:
                    sample = bam.split("/")[-1].split(".bam")[0]
                    if sample in snvs_frequencies[snv]:
                        snv_table[sample].append(snvs_frequencies[snv][sample][0])
                    else:
                        snv_table[sample].append(0)
                snv_table['Window'].append(snvs_frequencies[snp][sample][1])

        snv_table = pd.DataFrame(snv_table)

        # Create clonality table
        clonality_table = defaultdict(list)
        for window in clonality_by_window:
            for position_pair in clonality_by_window[window]:
                clonality_table['scaffold'].append(window.split(":")[0])
                clonality_table['window_name'].append(window)
                clonality_table['position'].append(position_pair[0])
                for bam in bams:
                    sample = bam.split("/")[-1].split(".bam")[0]
                    if sample in clonality_by_window[window]:
                        clonality_table[sample].append(position_pair[1])
                        clonality_table[sample + '_coverage'].append(coverages[position_pair[0]])
                    else:
                        clonality_table[sample].append(np.nan)
                        clonality_table[sample + '_coverage'].append(0)


        clonality_table_final = pd.DataFrame(clonality_table)

        #calc average coverage
        avg_coverage = 0.0
        for counts in counts_sum_by_bam:
            avg_coverage += sum(counts_sum_by_bam[counts]) 

        avg_coverage = avg_coverage / len(counts_sum_by_bam)

        # Final statistics
        print("Total SNVs-sites (across all BAMs): " + str(total_snv_sites))
        print("Total SNV-bases: (across all BAMs): " + str(alpha_snvs))
        # print("Mean clonality: " + str(float(sum(clonality_table['clonality'])) / float(len(clonality_table['clonality'])) ))
        print("Total sites (across all BAMs): " + str(total_positions))
        print("Average coverage: " + str(avg_coverage))
        # print("Total number of bases: " + str(sum(coverages.values())))

        self.coverages = coverages
        self.alpha_snvs = alpha_snvs
        self.total_snv_sites = total_snv_sites
        self.clonality_table = clonality_table_final
        self.snv_table = snv_table
        self.read_to_snvs = read_to_snvs
        self.total_positions = total_positions
        self.windows_to_snvs = windows_to_snvs
        self.snv_counts = snv_counts

        self.results = True


####### END STRAINPROFILER CLASS
################################
def strain_pipeline(args, filter_cutoff):
    strains = SNVdata()

    if args.testing:
        strains.testing = True

    if not args.output:
        strains.output = args.fasta.split("/")[-1].split(".")[0] + "_" + str(filter_cutoff)
    else:
        strains.output = args.output + "_" + str(filter_cutoff)

    if '*' in args.bams:
        bams = glob.glob(args.bams) 
    else:
        bams = args.bams.split(",")
    
    print("Running at resolution: (>" + str(filter_cutoff) + "%)")
    strains.get_scaffold_positions(args.genes, args.fasta)
    strains.run_strain_profiler(bams, min_coverage = int(args.min_coverage), min_snp = int(args.min_snp), filter_cutoff = filter_cutoff)
    strains.calc_linkage_network()
    strains.calc_ld_all_sites(int(args.min_snp))
    strains.save()

def main(args):
    '''
    Main entry point
    '''

    strain_pipeline(args, 0.98)
    strain_pipeline(args, 0.96)
    strain_pipeline(args, 0.94)
    strain_pipeline(args, 0.92)
    strain_pipeline(args, 0.90)
    strain_pipeline(args, 0)



if __name__ == '__main__':
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional arguments
    parser.add_argument("bams", help="Sorted .bams file,, comma separated list.")
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


