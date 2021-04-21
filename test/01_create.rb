# !/usr/bin/env ruby
require "test/unit/assertions"
include Test::Unit::Assertions
require "tempfile"

=begin
    Integration test 01_create_maxquant_db:

        - Input: FASTA file containing 1 CTG
        - expected output:
            - subsequence before part with CTG
            - subsequence with CTG translated into 19 amino acids (omitting Ile)
            - subsequence with CTG allows for 2 missed cleavages at N-term and C-term
            - subsequence after part with CTG
            - all generated subsequences are unique

=end

# require .rb files in library (including all subfolders)
Dir[File.join(__dir__, "..", "lib", "**", "*.rb")].each do |file|
    require File.absolute_path(file)
end

# input contains codon in question
input = Tempfile.new("gcf")
# full translation: MSDTLYSQREHQNNQKFEQLASTLHQFRTTVDHDIHNNVQQENSLLDSLNDSFNSLMVLVKQTSGELRTVMNRNASLTRIVGMILLGFFIIWMLYKLI
input.write(">h1\nATGTCAGATACATTATATTCTCAAAGAGAACATCAGAACAATCAGAAATTTGAACAATTAGCTTCAACTTTACACCAATTCAGAACCACCGTAGATCATGATATTCACAATAATGTTCAACAAGAAAACTCCCTACTAGATTCATTGAATGATAGTTTCAATTCTTTAATGGTGCTGGTTAAACAAACTTCTGGAGAATTAAGAACAGTTATGAATAGAAATGCTAGTCTTACTCGAATAGTAGGAATGATATTGCTTGGATTCTTTATAATTTGGATGTTGTATAAATTAATATAA")
input.close
output = Tempfile.new("gcf")
output_csv = Tempfile.new("gcf")

original_stdout = $stdout.clone
$stdout.reopen(File.new('/dev/null', 'w'))
system("ruby", File.join(__dir__, "..", "01_create_maxquant_db.rb"), "-i", input.path, "-o", output.path, "-m", output_csv.path)
$stdout.reopen(original_stdout)

# sequence for MaxQuant DB
# - contains one CTG
# - substr before CTG = MSDTLYSQR
# - unique sequences
# - allowing for two missed cleavages before CTG -> substring starts with EHQNN
# - allowign for two missed cleavages after CTG -> substring ends with TVMNR
# -> EHQNNQKFEQLASTLHQFRTTVDHDIHNNVQQENSLLDSLNDSFNSLMV_L_VKQTSGELRTVMNR
# - substr after CTG = NASLTRIVGMILLGFFIIWMLYKLI
# - expect total of 21 sequences (part before CTG, part after CTG, 19x part with CTG)
# - expect no translation of CTG into Ile
headers, seqs = Sequence.read_fasta(output.path)
assert_equal headers.size, 21
assert_equal seqs.size, 21
assert_equal seqs.uniq.size, 21
assert_include seqs, "MSDTLYSQR"
assert_include seqs, "NASLTRIVGMILLGFFIIWMLYKLI"
# exemplary tests for the expected 19 CTG translations
assert_include seqs, "EHQNNQKFEQLASTLHQFRTTVDHDIHNNVQQENSLLDSLNDSFNSLMVLVKQTSGELRTVMNR" # Leu
assert_include seqs, "EHQNNQKFEQLASTLHQFRTTVDHDIHNNVQQENSLLDSLNDSFNSLMVRVKQTSGELRTVMNR" # Arg
assert_not_include seqs, "EHQNNQKFEQLASTLHQFRTTVDHDIHNNVQQENSLLDSLNDSFNSLMVIVKQTSGELRTVMNR" # Ile

input.unlink
output.unlink
output_csv.unlink
