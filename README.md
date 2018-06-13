# strains_analysis
Code for analyzing population genomics in genome-resolved metagenomes

### Usage

```
python strainRepOneSample.py bam_file.bam genome.fasta [gene_list.txt] -c min_coverage_int -s min_snp_int
```

Where `gene_list.txt` is optional. If you do make this file, it should look like the below. The `-c` defines the minimum coverage required at a site to call SNVs on it (default 5). The `-s` site defines the minimum number of reads that have a given nucleotide at a position for that variant to be called as an SNV (default: 3). 

`Gene_list.txt`:
```
COG0001,14_0903_02_20cm_scaffold_13826,7320,7544
COG0002,14_0903_02_20cm_scaffold_5262,27,1703
COG0003,14_0903_02_20cm_scaffold_5308,3879,5147
COG0004,14_0903_02_20cm_scaffold_5308,5156,5650
COG0005,14_0903_02_20cm_scaffold_5308,6075,6737
COG0006,14_0903_02_20cm_scaffold_5308,6827,7180
```

### TODO
0. Upload code to GitHub
1. Add max coverage?
2. SNP-frequency x SNVs per read matrix
3. Does this actually work?
