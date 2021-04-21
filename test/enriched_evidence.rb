# !/usr/bin/env ruby
require "test/unit/assertions"
include Test::Unit::Assertions
require "tempfile"

=begin
    Test lib/parse_enriched_evidence.rb:

        - parse header and subsequent lines
        - parse original protein name and simplified name used in MaxQuant tables
        - parse peptide positions in protein
        - parse peptide codons
        - parse CTG position in peptide
        - parse supported positions
=end

# require .rb files in library (including all subfolders)
Dir[File.join(__dir__, "..", "lib", "**", "*.rb")].each do |file|
    require File.absolute_path(file)
end

header = "cDNA\tStart pos in protein\tEnd pos in protein\tCTG codon pos\tOriginal protein name\tb/y-ion supported pos\tCorresponding scan numbers\tSequence\tLength\tModifications\tModified sequence\tOxidation (M) Probabilities\tOxidation (M) Score Diffs\tAcetyl (Protein N-term)\tOxidation (M)\tMissed cleavages\tProteins\tLeading Proteins\tLeading Razor Protein\tType\tRaw file\tFraction\tExperiment\tMS/MS m/z\tCharge\tm/z\tMass\tResolution\tUncalibrated - Calibrated m/z [ppm]\tUncalibrated - Calibrated m/z [Da]\tMass Error [ppm]\tMass Error [Da]\tUncalibrated Mass Error [ppm]\tUncalibrated Mass Error [Da]\tMax intensity m/z 0\tRetention time\tRetention length\tCalibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\tRetention time calibration\tMatch time difference\tMatch m/z difference\tMatch q-value\tMatch score\tNumber of data points\tNumber of scans\tNumber of isotopic peaks\tPIF\tFraction of total spectrum\tBase peak fraction\tPEP\tMS/MS Count\tMS/MS Scan Number\tScore\tDelta score\tCombinatorics\tIntensity\tReverse\tPotential contaminant\tid\tProtein group IDs\tPeptide ID\tMod. peptide ID\tMS/MS IDs\tBest MS/MS\tAIF MS/MS IDs\tOxidation (M) site IDs"
line = "GCCGCCGCAGCACTGGGGGCTGCGCTTGCGCCGCAGCGC\t110\t122\t4\tEEQ37961 cdna supercontig:ASM383v1:CH408077:1767990:1769492:-1 gene:CLUG_02084 gene_biotype:protein_coding transcript_biotype:protein_coding description:hypothetical protein\t0;1;2;3;4;5;6;7;8;9;10;11;12\t8270\tAAAASGAALAPQR\t13\tUnmodified\t_AAAASGAALAPQR_\t0\t0\t0\tg1_102_127-CTG-S\tg1_102_127-CTG-S\tg1_102_127-CTG-S\tMULTI-MSMS\tH_Schmitt_211216_190117_L1_R1_14\t14\tL1\t5.778.184\t2\t577.817.461\t115.362.037\t42112.35\t-0.42427\t-0.00024515\t0.91792\t0.00053039\t0.49364\t0.00028524\t577.818.033.203.664\t16.334\t0.21083\t16.472\t16.339\t16.549\t0.13784\t62\t35\t2\t0\t0\t0\t1,28E-04\t1\t8270\t116.24\t99.661\t1\t16324000\t63\t1435\t11\t11\t51\t51"


mq_data = ParseEnrichedEvidence.new()

# ParseEnrichedEvidence.is_header & parse_line
# - treat (only) first line as header line
# - identify Sequence, Proteins, id rows
# no use testing for
#   - Score and Mass error [ppm], as those are converted to float
#   - Reverse and Potential Contaminant, as those are tested elsewhere
assert_true mq_data.is_header
assert_equal "CTG", mq_data.parse_header(header)
assert_false mq_data.is_header

mq_data.parse_line(header)
assert_equal "Original protein name", mq_data.get_protein
assert_equal "Proteins", mq_data.get_protein_name_used_in_db

# ParseEnrichedEvidence.get_*
# - codons are ["GCC", "GCC", "GCA", "GCA", "CTG", "GGG", "GCT", "GCG", "CTT", "GCG", "CCG", "CAG", "CGC"]
# - CTG pos is 4
# - peptide start is 110
# - peptide stop is 122
# - protein pos is 110
# - supported pos are 0 - 12
# - scan numbers is 8270
# - orignal protein name is "EEQ37961 cdna supercontig:ASM383v1:CH408077:1767990:1769492:-1 gene:CLUG_02084 gene_biotype:protein_coding transcript_biotype:protein_coding description:hypothetical protein"
# - protein name as listed in MaxQuant tables is "g1"
mq_data.parse_line(line)
assert_equal ["GCC", "GCC", "GCA", "GCA", "CTG", "GGG", "GCT", "GCG", "CTT", "GCG", "CCG", "CAG", "CGC"], mq_data.get_codons
assert_equal [4], mq_data.get_codon_pos
assert_equal 110, mq_data.get_peptide_start
assert_equal 122, mq_data.get_peptide_stop
assert_equal 110, mq_data.convert_peptide_to_protein_pos(0)
assert_equal "EEQ37961 cdna supercontig:ASM383v1:CH408077:1767990:1769492:-1 gene:CLUG_02084 gene_biotype:protein_coding transcript_biotype:protein_coding description:hypothetical protein", mq_data.get_protein
assert_equal "g1", mq_data.get_protein_name_used_in_db
assert_equal ["8270"], mq_data.get_scannumbers
assert_equal (0..12).to_a, mq_data.get_supported_pos
