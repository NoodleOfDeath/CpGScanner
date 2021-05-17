# CpG-Scanner
Simple scripts for identifying CpG islands in a genomic sequence string.

# OOTB CLI Usage (Python)

```
$ python3 scan.py
Took 0.00 second(s) to generate 1024 character sequence
Searching 1024 character sequence for CpG islands using 2 threads, threshold of 0.60, min length of 8, and slicing into chunks of size 4 or less
Sliced 1024 character sequence into 256 chunks
Found 19 CpG islands matching the criteria
Took 12.89 second(s) to find 19 CpG islands in 1024 character sequence
[(GGGGTGGCGTGC, 32, 12), (CTGGTGCG, 72, 8), (CACCACCGAGCC, 152, 12), (CTCGCACC, 180, 8), (CGACGGCG, 236, 8), (CGCTGTGC, 260, 8), (GACGAGCC, 276, 8), (AGGCCGGAGCCG, 320, 12), (ACGGGACC, 368, 8), (GCCTGACG, 384, 8), (ACCCCGCG, 396, 8), (CCGGGCCA, 476, 8), (GGTCCCAG, 548, 8), (CCGGTGCG, 680, 8), (GGGGTCGCTGCGAGGC, 696, 16), (CCGGGCAGCACG, 716, 12), (CGCAGCTCGCGAGGCCGGAG, 764, 20), (TGGGCGCAGTGG, 788, 12), (CAGCCGGCGGCCCACC, 956, 16)]
```
