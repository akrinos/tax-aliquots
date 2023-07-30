# Code for `tax-aliquots` (Assigning Lineage to Queries Over Two Steps)

This code has been made available for testing and reproducibility purposes in conjunction with the manuscript "Missing microbial eukaryotes and misleading meta-omic conclusions". In the future, tax-aliquots will be incorporated into the `EUKulele` taxonomic annotation tool (https://github.com/AlexanderLabWHOI/EUKulele).

### Arguments

* `-i` : input file from DeepClust (https://github.com/bbuchfink/diamond)
* `-c` : FASTA file that contains all sequences present in the `DeepClust` file, both reference and query
* `--output_file` : prefix for output file names
* `--output_dir` : directory in which to save output
* `-t` : taxonomy table, in the style of `EUKulele`
* `--prot_map` : protein name mapping, in the style of `EUKulele`
* `--kmer_len` : amino acid _k_-mer length to use
* `-s` : sample delimiter in sequence file
* `--level_of_interest` : taxonomic level of interest (e.g. phylum, family; should match EUKulele
* `--list_of_interest` : underscore-separated list of taxonomic labels to pull from the DeepClust file

### Invocation
`tax-aliquots` can be invoked using the `process_clusters.py` file in `src/tax-aliquots-scripts`. `tax-aliquots` uses `faSomeRecords` (https://github.com/santiagosnchez/faSomeRecords) to read in records from fasta files

### Example command

```
python process_clusters.py \
    -i deepclust_contigname.mad.50.out \
    -c combined_seqs.fasta \
    --output_file="test" \
    --output_dir="test_dir" \
    -t tax-table.txt \
    --prot_map=prot-map.json \
    --kmer_len=3 -s test --level_of_interest family \
    --list_of_interest Hemiaulaceae_Rhizosoleniaceae_Thalassiosiraceae_Skeletonemataceae
```

### Citations

Buchfink, Benjamin, et al. "Sensitive clustering of protein sequences at tree-of-life scale using DIAMOND DeepClust." bioRxiv (2023): 2023-01.
Sanchez-Ramirez, Santiago, et al. "faSomeRecords." https://github.com/santiagosnchez/faSomeRecords