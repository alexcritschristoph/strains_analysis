import sys
import glob

sys.path.append('/data6/Angelo/alexcc/AngeloStrainsPaper/strainrep/strains_analysis/')
import strainRep2
null_model = strainRep2.generate_snp_model('./combined_null1000000.txt')

# first parameter is directory of .data files
files = sys.argv[1]
# second parameter is the clonality level to test (0.98, etc.)
level = sys.argv[2]

# third parameter is the genome
genome = sys.argv[3]

from collections import defaultdict

print(files)
print("Getting data...")

f = open('../genomes_above_cutoffs.txt')
above_cutoff = []

for line in f.readlines():
    above_cutoff.append(":".join(line.split()))
f.close()



counts = defaultdict(lambda: defaultdict(list))

samples = []
## Calculate the combined counts
for fn in glob.glob(files + genome + "*" + level + ".data"):

    if genome + ":" + sample in above_cutoff:
        samples.append()
        sample = fn.split(":")[1].split("_" + level)[0]
        print(sample)

        s = strainRep2.SNVdata()
        s.load(fn.replace(".data",""))

        #add counts to all 
        for position in s.all_counts:
            counts[genome + ":" + sample][position] = list(s.all_counts[position])

            if position in counts['all']:
                counts['all'][position][0] += s.all_counts[position][0]
                counts['all'][position][1] += s.all_counts[position][1]
                counts['all'][position][2] += s.all_counts[position][2]
                counts['all'][position][3] += s.all_counts[position][3]
            else:
                counts['all'][position] = list(s.all_counts[position])


print("Calling SNPs")

f = open('')
snp_sites = []
print("Contig,Position", end="")
for sample in samples:
    print("," + sample + "-A", end="")
    print("," + sample + "-C", end="")
    print("," + sample + "-G", end="")
    print("," + sample + "-T", end="")

print("")
for pos in counts['all']:
    consensus = strainRep2.call_snv_site(counts['all'][pos], null_snp_model, min_cov = 20, min_freq = 0.05)
    if consensus:
    	print("")
        #there is a SNP here
        print(pos.split(":")[0], end='')
        print("," + pos.split(":")[1], end='')

        for sample in samples:

            print("," + str(counts[genome + sample][0]), end="") #A
            print("," + str(counts[genome + sample][1]), end="") #C
            print("," + str(counts[genome + sample][3]), end="") #G
            print("," + str(counts[genome + sample][2]), end="") #T (note reversal bc counts is ACTG)