# !/usr/bin/env ruby
require "test/unit/assertions"
include Test::Unit::Assertions
require "tempfile"

=begin
    Test lib/parse_evidence.rb:

    TODO

        - read FASTA, fail if file contains duplicate headers
        - write FASTA
        - split sequence into codons
        - translate sequence, omit stop codon
        - translate sequence with translation for single codon codon set
        - identify stop codon
        - # amino acids: 20
        - # codons: 64

=end

# require .rb files in library (including all subfolders)
Dir[File.join(__dir__, "..", "lib", "**", "*.rb")].each do |file|
    require File.absolute_path(file)
end

header = "Sequence	Length	Modifications	Modified sequence	Oxidation (M) Probabilities	Oxidation (M) Score Diffs	Acetyl (Protein N-term)	Oxidation (M)	Missed cleavages	Proteins	Leading Proteins	Leading Razor Protein	Type	Raw file	Fraction	Experiment	MS/MS m/z	Charge	m/z	Mass	Resolution	Uncalibrated - Calibrated m/z [ppm]	Uncalibrated - Calibrated m/z [Da]	Mass Error [ppm]	Mass Error [Da]	Uncalibrated Mass Error [ppm]	Uncalibrated Mass Error [Da]	Max intensity m/z 0	Retention time	Retention length	Calibrated retention time	Calibrated retention time start	Calibrated retention time finish	Retention time calibration	Match time difference	Match m/z difference	Match q-value	Match score	Number of data points	Number of scans	Number of isotopic peaks	PIF	Fraction of total spectrum	Base peak fraction	PEP	MS/MS Count	MS/MS Scan Number	Score	Delta score	Combinatorics	Intensity	Reverse	Potential contaminant	id	Protein group IDs	Peptide ID	Mod. peptide ID	MS/MS IDs	Best MS/MS	AIF MS/MS IDs	Oxidation (M) site IDs"


mq_data = ParseEvidence.new("CTG")

# TODO
# - parse header:
# only once.
# - return values correctly:
#      - protein region
#      - codon transle
#      - is a decoy match??
#      - get evid
#      - get peptide,
#         - get peptide with unnown AA
#     - get prot
#     - get score
#     - get mass err
#     - is mass error not specified
#     - get codon pos

# ParseEvidence.is_header
# - treat (only) first line as header line
# - identify Sequence, Proteins, id, Score, Mass error [ppm], Reverse and Potential Contaminant rows
assert_true mq_data.is_header
mq_data.parse_header(header)
assert_false mq_data.is_header

mq_data.parse_line(line)
mq_data.get_peptide

require "byebug"
debugger
assert_equal 64, Sequence.codons.size

# ParseEvidence.parse_header
# - ???
assert_equal 64, Sequence.codons.size

# ParseEvidence.parse_line
#  - ???
