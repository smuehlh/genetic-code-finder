# genetic-code-finder

genetic-code-finder is a Ruby library to identify genetic codes in yeasts.

It creates input for the proteomics software package [MaxQuant](https://www.maxquant.org/) and parses MaxQuant's output for codon translations.

Original publications based these scripts are
 * https://doi.org/10.1016/j.cub.2018.04.085
 * https://doi.org/10.1101/gr.200931.115

Copyright (c) 2021, by Göttingen University

## Installation

No installation required other than Ruby 3.0 or later.


## Usage

A complete work flow will consist of the following steps

 1. Create FASTA file to be used as MaxQuant database, containing a codon (e.g. CTG) translated into each amino acid.
 2. [Start MaxQuant]
 3. Combine and enrich MaxQuant output tables evidence.txt and msms.txt for further analyses.
 4. Collect statistics about found codon translations.

This library contains scripts for each of the above mentioned steps, named in their order of execution.

### Script usage

#### 1. Prepare database input for MaxQuant search

```bash
ruby 01_create_maxquant_db.rb --input sample_data/Clavispora_cDNA_excerpt.fasta --output Clavispora_maxquant_db.fas --map Clavispora_maxquant_db_map.csv [--codon CTG] [--cleavage K,R] [--ile]
```

 - Split genomic FASTA (parameter `--input`) at cleavage sites (parameter `--cleavage`, defaults to K, R)
 - Translate *in-silico* cleaved sequences, translating designated codon (parameter `--codon`, defaults to CTG) into all amino acids (resulting in 19  peptides all differing in their CTG translation)
 - Per default, translation into isoleucine will be omitted (parameter `--ile` to create peptides for both isoleucine and leucine)
- Standardise original sequence names (headers in input FASTA) for use in MaxQuant, storing the mapping between original and standardised names (parameter `--map`)
- Save resulting peptides as FASTA (parameter `--output`) to use them as input in MaxQuant search

#### 2. Combine and enrich MaxQuant output files `evidence.txt` and `msms.txt` for downstream analysis

```bash
ruby 02_combine_maxquant_tables.rb --evidence sample_data/evidence.txt --msms sample_data/msms.txt --map Clavispora_maxquant_db_map.csv --cdna sample_data/Clavispora_cDNA_excerpt.fasta --output Clavispora_enriched_evidence.txt [--codon CTG]
```

 - Map peptide spectrum matches in evidence file (parameter `--evidence`) proteins (parameter `--cdna`, optional parameter `--codon`) and original protein names (parameter `--map`)
 - Map matched peptide fragments from msms file (parameter `--msms`) onto corresponding peptide spectrum matches
 - Remove decoy hits
 - Save resulting data (parameter `--output`) to use it as input for downstream analysis

#### 3. Analyse MaxQuant output and collect statistics about found codon (e.g. CTG) translations

```bash
ruby 03_make_statistics.rb --input Clavispora_enriched_evidence.txt --cdna sample_data/Clavispora_cDNA_excerpt.fasta --output Clavispora_statistics.txt --genes Clavispora_gene_statistics.csv
```

 - Collect generic counts such as number of peptide spectrum matches and recovered proteins and CTG positions (parameters `--input`, `--cdna`)
 - Collect CTG translation counts, based on supported CTG positions
 - Save generic counts as flat text file (parameter `--output`)
 - Save generic counts per gene (parameter `--genes`)

#### 4. [Optional] plot statistics

Rscripts are meant as starting point for custom plots.
```bash
Rscript 04_plot_translation.R Clavispora_statistics.txt [CTG]
```
- Plot CTG translation counts (parameter `input.txt`)
    - Fraction of PSMs containing CTG with a given translation

```bash
Rscript 04_plot_genes.R Clavispora_gene_statistics.csv [CTG]
```

- Plot counts per gene (parameter `input.csv`)
    - Sequence coverage
    - Number of PSMs/ non-redundant peptides
    - Number of PSMs containing CTG
    - CTG coverage

## Test
### Unit tests

```bash
ruby test/sequences.rb  # test lib/sequence.rb
ruby test/evidence.rb  # test lib.parse_evidence.rb
ruby test/msms.rb  # test lib/parse_msms.rb
ruby test/enriched_evidence.rb  # test lib/parse_enriched_evidence.rb
ruby test/statistics.rb  # test lib/statistics.rb
ruby test/file_helper.rb  # test lib/file_helper.rb
```

### Integration tests

```bash
ruby test/01_create.rb  # test script 01_create_maxquant_db.rb
ruby test/02_combine_tables.rb  # test script 02_combine_maxquant_tables.rb
ruby test/03_statistics.rb  # test script 03_make_statistics.rb
```

## Test data
For playing around with this library, some sample data is provided. It is an excerpt of [ProteomeXchange](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD009494-1&test=no), Clavispora processed data.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Author
Stefanie Mühlhausen

## License
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.de.html)
