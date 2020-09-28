# genetic-code-finder

Genetic-code-finder is a Ruby library to identify genetic codes in yeasts.

It creates input for the proteomics software package [MaxQuant](https://www.maxquant.org/) and parses MaxQuant's output for codon translations.

Original publications based these scripts are
 * https://doi.org/10.1016/j.cub.2018.04.085
 * https://doi.org/10.1101/gr.200931.115

Copyright (c) 2020, by Göttingen University

## Installation

No installation required other than Ruby 2.7 or later.


## Usage

A complete work flow will consist of the following steps

 1. Create FASTA file to be used as MaxQuant database, containing a codon (e.g. CTG) translated into each amino acid.
 2. [Run MaxQuant]
 3. Combine and enrich MaxQuant output tables evidence.txt and msms.txt for further analyses.

This library contains scripts for each of the above mentioned steps, named in their order of execution.

```bash
ruby 01_create_maxquant_db.rb -h
ruby 02_combine_maxquant_tables.rb -h
```

## Test
### Unit tests

```bash
ruby test/sequences.rb # test Sequence
ruby test/evidence.rb # test ParseEvidence
ruby test/msms.rb # test ParseMsms
```
### Integration tests

```bash
ruby test/01_create.rb # test 01_create_maxquant_db.rb
ruby test/02_combine_tables.rb # test 02_combine_maxquant_tables.rb
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
