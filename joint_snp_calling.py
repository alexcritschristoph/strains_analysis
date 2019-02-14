import sys
import glob

sys.path.append('/data6/Angelo/alexcc/AngeloStrainsPaper/strainrep/strains_analysis/')
import strainRep2
null_snp_model = strainRep2.generate_snp_model('./combined_null1000000.txt')

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
print(files + genome + "*" + level + ".data")
for fn in glob.glob(files + genome + "*" + level + ".data"):
    sample = fn.split(":")[1].split("_" + level)[0]
    if genome + ":" + sample in above_cutoff:
        samples.append(sample)
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

snp_sites = []

f = open(genome + ".freqs", "w+")

f.write("Contig,Position")
for sample in samples:
    f.write("," + sample + "-A")
    f.write("," + sample + "-C")
    f.write("," + sample + "-G")
    f.write("," + sample + "-T")

f.write("\n")

for pos in counts['all']:
    consensus = strainRep2.call_snv_site(counts['all'][pos], null_snp_model, min_cov = 20, min_freq = 0.05)
    if consensus:
        #there is a SNP here
        f.write("_".join(pos.split("_")[:-1]))
        f.write("," + pos.split("_")[-1])

        for sample in samples:
            if len(counts[genome + ":" + sample][pos]) == 0:
                counts[genome + ":" + sample][pos] = [0,0,0,0]

            f.write("," + str(counts[genome + ":" + sample][pos][0])) #A
            f.write("," + str(counts[genome + ":" + sample][pos][1])) #C
            f.write("," + str(counts[genome + ":" + sample][pos][3])) #G
            f.write("," + str(counts[genome + ":" + sample][pos][2])) #T (note reversal bc counts is ACTG)
        f.write("\n")
f.close()
