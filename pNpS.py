import csv
import sys
import glob
import argparse
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio.codonalign import CodonSeq
from collections import defaultdict
from Bio.codonalign.codonalphabet import default_codon_alphabet, default_codon_table 


def calc_pn_ps(gene_fasta, freq_file):
    """ calculates pN/pS for all genes, writes to table
    """
    all_genes = []
    for record in SeqIO.parse(gene_fasta, 'fasta'):
      all_genes.append(record.id)

    site_data = calc_S_nS_sites(gene_fasta, freq_file)

    freq = pd.read_table(freq_file)
    
    #calculate data per gene
    snps_nS = defaultdict(int)
    snps_S = defaultdict(int)
    for gene in freq.gene.unique():
      for index, snp in freq.loc[freq['gene'] == gene].iterrows():
        if 'N:' in snp['mutation']:
          snps_nS[gene] += 1
        elif 'S:' in snp['mutation']:
          snps_S[gene] += 1


    print("Gene\tpN/pS\tN\tS\tsites_N\tsites_S")
    for gene in all_genes:
      if gene in site_data:
        if snps_S[gene] == 0 or site_data[gene][1] == 0 or site_data[gene][0] == 0:
          pnps = 'inf'
        else:
          pnps = round(float(snps_nS[gene]) / float(site_data[gene][1]) / (float(snps_S[gene]) / float(site_data[gene][0])), 3)
        print(gene + "\t" + str(pnps) + "\t" + str(round(float(snps_nS[gene]),3)) + "\t" + str(round(float(snps_S[gene]),3)) + "\t" + str(round(site_data[gene][1],3)) + "\t" + str(round(site_data[gene][0],3)) )




def calc_S_nS_sites(gene_fasta, freq):
    """calculates the number of S and nS sites per gene.
    """
    
    genes = {}
    for record in SeqIO.parse(gene_fasta, 'fasta'):
        my_seq = record.seq.transcribe()
        if (len(my_seq) % 3 != 0):
            print(feature.qualifiers['product'][0] + " is not divisible by 3")
        elif 'N' not in str(my_seq):
            cod_seq = CodonSeq(str(my_seq))
            listt = get_codon_list(cod_seq)
            # get rid of stop codon
            listt.pop()            
            tS, tN = count(listt)
            genes[record.id] = [tS, tN]
        elif 'N' in str(my_seq):
          print("Note: gene " + str(record.id) +" was skipped because it has N's in it.")
    return genes

def get_codon_list(codonseq): 
      """List of codons according to full_rf_table for counting (PRIVATE).""" 
      full_rf_table = codonseq.get_full_rf_table() 
      codon_lst = [] 
      for i, k in enumerate(full_rf_table): 
          if isinstance(k, int): 
              start = k 
              try: 
                  end = int(full_rf_table[i+1]) 
              except IndexError: 
                  end = start+3 
              this_codon = str(codonseq[start:end]) 
              if len(this_codon) == 3: 
                  codon_lst.append(this_codon) 
              else: 
                  codon_lst.append(str(this_codon.ungap())) 
          elif str(codonseq[int(k):int(k)+3]) == "---": 
              codon_lst.append("---") 
          else: 
              # this may be problematic, as normally no codon shoud 
              # fall into this condition 
              codon_lst.append(codonseq[int(k):int(k)+3]) 
      return codon_lst 

def count(codon_lst, k=1, codon_table=default_codon_table):
      S_site = 0.0  # synonymous sites 
      N_site = 0.0  # non-synonymous sites 
      purine = ('A', 'G') 
      pyrimidine = ('T', 'C') 
      base_tuple = ('A', 'T', 'C', 'G') 
      for codon in codon_lst: 
          neighbor_codon = {'transition': [], 'transversion': []} 
          # classify neighbor codons 
          codon = codon.replace('U', 'T') 
          if codon == '---': 
              continue 
          for n, i in enumerate(codon): 
              for j in base_tuple: 
                  if i == j: 
                      pass 
                  elif i in purine and j in purine: 
                      codon_chars = [c for c in codon] 
                      codon_chars[n] = j 
                      this_codon = ''.join(codon_chars) 
                      neighbor_codon['transition'].append(this_codon) 
                  elif i in pyrimidine and j in pyrimidine: 
                      codon_chars = [c for c in codon] 
                      codon_chars[n] = j 
                      this_codon = ''.join(codon_chars) 
                      neighbor_codon['transition'].append(this_codon) 
                  else: 
                      codon_chars = [c for c in codon] 
                      codon_chars[n] = j 
                      this_codon = ''.join(codon_chars) 
                      neighbor_codon['transversion'].append(this_codon) 
          # count synonymous and non-synonymous sites 
          #codon = codon.replace('T', 'U')
          if (codon == 'TAG'):
           print("STOP DETECTED")
           continue
          aa = codon_table.forward_table[codon]
          this_codon_N_site = this_codon_S_site = 0 
          for neighbor in neighbor_codon['transition']: 
              if neighbor in codon_table.stop_codons: 
                  this_codon_N_site += 1 
              elif codon_table.forward_table[neighbor] == aa: 
                  this_codon_S_site += 1 
              else: 
                  this_codon_N_site += 1 
          for neighbor in neighbor_codon['transversion']: 
              if neighbor in codon_table.stop_codons: 
                  this_codon_N_site += k 
              elif codon_table.forward_table[neighbor] == aa: 
                  this_codon_S_site += k 
              else: 
                  this_codon_N_site += k 
          norm_const = (this_codon_N_site + this_codon_S_site)/3 
          S_site += float(this_codon_S_site) / float(norm_const) 
          N_site += float(this_codon_N_site) / float(norm_const) 
      return (S_site, N_site) 



if __name__ == '__main__':



    parser = argparse.ArgumentParser(description= """
        Calculates gene-level pN/pS (a population dN/dS without complete genotypes) from an allele .freq file and a prodigal .fna file.\n
        The calculation here is: (# of non-synonymous SNPs / # of non-synonymous sites) / (# of synonymous SNPs / # of synonymous sites).\n
        Requires your .freq file to have already identified SNPs as S or N using the gene_analysis.py script (will end in _aa.).\n  
        Output: a new file 
        """, formatter_class=argparse.RawTextHelpFormatter)

    # Arguments for gene_snps
    parser.add_argument("-f", "--fasta", action="store", default=None, \
        help='Path to prodigal .fna fasta file of predicted genes (dna sequences of genes, not contigs or AA)')

    parser.add_argument("-n", "--freq", action="store", default=None, \
        help='Path to the .freq file from the main strain data calculation script.')

    args = parser.parse_args()

    if args.fasta and args.freq:
        calc_pn_ps(args.fasta, args.freq)
    else:
        sys.exit("Error: need to specific a fasta file and a frequency file")
