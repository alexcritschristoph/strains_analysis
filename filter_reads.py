import sys
import pysam
import argparse
import numpy as np
from tqdm import tqdm
from Bio import SeqIO


def get_fasta(fasta_file = None):
    positions = []
    for rec in SeqIO.parse(fasta_file, "fasta"):
        start = 0
        
        positions.append([str(rec.id),start, len(rec.seq)])

    return positions

def filter_reads(bam, positions, filter_cutoff = 0.97, max_insert_relative = 3, min_insert = 50, min_mapq = 2, write_data = None, write_bam = False):

    # read sets
    observed_read1s = set()
    observed_read2s = set()
    mapped_pairs = set()
    final_reads = set()
    # counters
    total_read_count = 0
    total_read_pairs = 0
    total_mapped_pairs = 0

    # storing data
    read_data = {}
    pair_mapqs = {}
    pair_lengths = {}
    pair_mismatch = {}
    pair_inserts = {}
    pair_lengths = {}


    samfile = pysam.AlignmentFile(bam)
    print("READING BAM: " + bam.split("/")[-1])

    print("Using reads with >" + str(filter_cutoff) + "% PID to consensus reference.")

    ## STEP 1: collect paired reads and their information
    for gene in tqdm(positions, desc='Getting read pairs: '):
        for read in samfile.fetch(gene[0], gene[1], gene[2]):
            total_read_count += 1

            ## If we've seen this read's pair before
            if (read.is_read2 and read.query_name in observed_read1s) or (read.is_read1 and read.query_name in observed_read2s):

                #But if we haven't already seen this complete pair, then we can complete the pair and store the information
                #Also check that the pair is on the same scaffold
                if read.query_name not in mapped_pairs and gene[0] == read_data[read.query_name]['scaf']:
                    total_read_pairs += 1

                    if read.get_reference_positions() != []:
                        total_mapped_pairs += 1
                        mapped_pairs.add(read.query_name) #add to found 

                        #set mismatch percentage
                        pair_mismatch[read.query_name] = 1- ( ( float(read_data[read.query_name]['nm']) + float(read.get_tag('NM')) ) / ( float(read_data[read.query_name]['len']) + read.infer_query_length()) )
                        #set insert size
                        if read.get_reference_positions()[-1] > read_data[read.query_name]['start']:
                            pair_inserts[read.query_name] = read.get_reference_positions()[-1] - read_data[read.query_name]['start']
                        else:
                            pair_inserts[read.query_name] = read_data[read.query_name]['stop'] - read.get_reference_positions()[0]
                        #set mapq
                        pair_mapqs[read.query_name] = read.mapping_quality
                        if read_data[read.query_name]['mapq'] > read.mapping_quality:
                            pair_mapqs[read.query_name] = read_data[read.query_name]['mapq'] 

            #this is the first time we see a read from this pair and don't double count 
            elif (read.is_read1 and read.query_name not in observed_read1s) or (read.is_read2 and read.query_name not in observed_read2s):
                if read.get_reference_positions() != []: # don't use unmapped reads
                        if read.is_read1:
                            observed_read1s.add(read.query_name)
                        else:
                            observed_read2s.add(read.query_name)
                        #record the data for this read
                        read_data[read.query_name] = {"nm": read.get_tag('NM'), "len": read.infer_query_length(), "mapq": read.mapping_quality, "start": read.get_reference_positions()[0], 'stop': read.get_reference_positions()[-1], 'scaf': gene[0]}


    ## STEP 2: INSERT SIZE CUTOFF, MAPQ CUTOFF, AND MISMATCH CUTOFF
    max_insert = np.median(list(pair_inserts.values())) * max_insert_relative #insert size should be less than max_insert_relative * median value
    too_short = 0.0
    too_long = 0.0
    good_length = 0.0
    mapq_good = 0.0
    filter_cutoff_good = 0.0
    

    for read_pair in mapped_pairs:
        if pair_inserts[read_pair] > min_insert:
            if pair_inserts[read_pair] < max_insert:
                good_length += 2
                if pair_mapqs[read_pair] > min_mapq:
                    mapq_good += 2

                    # Which set does this read go into?
                    if pair_mismatch[read_pair] > filter_cutoff:
                        filter_cutoff_good += 2
                        final_reads.add(read_pair)
            else:
                too_long += 2
        else:
            too_short += 2
    
    print("**READ STATSTICS**")
    print("total reads found: " + str(total_read_count))
    print("total paired reads: " + str(total_read_pairs*2) + " (" + str(int(100*total_read_pairs*2.0 / total_read_count)) + "%)")
    print("total same scaffold mapped paired reads: " + str(total_mapped_pairs*2) + " (" + str(int(100*total_mapped_pairs*2.0 / total_read_count)) + "%)")
    print("")
    print("median insert size: " + str(max_insert / max_insert_relative))
    print("paired reads < 50 bp apart: " + str(too_short))
    print("paired reads > " + str(max_insert) + " apart: " + str(too_long))
    print("reads which also pass both pair insert size filters: " + str(good_length) + " (" + str(int(100*float(good_length) / total_read_count)) + "%)")
    print("(final) reads which also pass read pair PID >" + str(filter_cutoff) + "%: " + str(filter_cutoff_good) + " (" + str(int(100*float(filter_cutoff_good) / total_read_count)) + "%)")

    ## STEP 3: WRITE DATA IF NEEDED
    if write_data:
        f = open(write_data, 'w+')
        for read_pair in mapped_pairs:
            f.write(read_pair + "\t" + scafs[read_pair] + "\t" + str(pair_inserts[read_pair]) + "\t" + str(pair_mapqs[read_pair]) + "\t" + str(pair_mismatch[read_pair]) + "\n")
        f.close()
    ## STEP 4: WRITE NEW BAM IF NEEDED (TODO)

    return final_reads

        
            


if __name__ == '__main__':

    """ This is executed when run from the command line """

    parser = argparse.ArgumentParser(description="""Reports read statistics for a BAM mapping file. Do you have a SAM? Not to fear, just run these three commands:\n 
samtools view -S -b sample.sam > sample.bam\n
samtools sort sample.bam -o sample.sorted.bam\n
samtools index sample.sorted.bam\n in that order!""", formatter_class=argparse.RawTextHelpFormatter)


    # Required positional arguments
    parser.add_argument("bam", help="Sorted .bam file")
    parser.add_argument("fasta", help="Fasta file the bam is mapped to")


    # Optional arguments
    parser.add_argument("-m", "--mismatch_threshold", action="store", default=0.97, \
        help='Minimum percent identity of read pairs to consensus to use the reads - default is to run at 0.97.')

    parser.add_argument("-q", "--min_mapq", action="store", default=2, \
        help='Minimum mapq score of EITHER read in a pair to use that pair. Default: 2.')

    parser.add_argument("-l", "--max_insert_length", action="store", default=3, \
        help='Maximum insert size between two reads - default is to use 3x median insert size.')

    parser.add_argument("-u", "--min_insert_length", action="store", default=50, \
        help='Minimum insert size between two reads - default is to use 50 bp.')

    parser.add_argument("-w", "--write", action="store", default=None, \
        help='File name to write read statistics to.')



    # Parse
    args = parser.parse_args()
    positions = get_fasta(args.fasta)
    filter_reads(args.bam, positions, args.mismatch_threshold, args.max_insert_length, args.min_insert_length, args.min_mapq, write_data = args.write)