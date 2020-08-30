# !/usr/bin/env ruby
require "test/unit/assertions"
include Test::Unit::Assertions
require "tempfile"

=begin
    Test lib/sequences:

        - read FASTA, fail if file contains duplicate headers
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

# Sequence.read_fasta
# - returns headers and sequences
# - aborts on duplicate headers
file = Tempfile.new("gcf")
file.write(">h1\nSEQ")
file.close
assert_equal [["h1"], ["SEQ"]], Sequence.read_fasta(file.path)
file.unlink

file = Tempfile.new("gcf")
file.write(">h1\nSEQ\n>h1\nSEQ")
file.close
# temporarily redirect stderr to suppress expected stderr- message
original_stderr = $stderr.clone
$stderr.reopen(File.new('/dev/null', 'w'))
assert_raise SystemExit do
    Sequence.read_fasta(file.path)
end
file.unlink
$stderr.reopen(original_stderr)

# Sequence.split_cdna_into_codons
# - split each 3. char
assert_equal ["ATG"], Sequence.split_cdna_into_codons("ATG")
assert_equal ["ATG", "A"], Sequence.split_cdna_into_codons("ATGA")

# Sequence.is_stopcodon
# - identify stop codons in standard genetic code
assert_true Sequence.is_stopcodon("TAA")
assert_false Sequence.is_stopcodon("ATG")

# Sequence.translate
# - translate accordig to standard genetic code
# - omit trailing stop codon
# - 'translate' invalid codons as "X"
assert_equal "ML", Sequence.translate("ATGCTG")
assert_equal "M", Sequence.translate("ATGTAA")
assert_equal "MX", Sequence.translate("ATGNAA")

# Sequence.translate_codons_with_one_codon_set
# - translate all but <codon> to standard genetic code
assert_equal "AL", Sequence.translate_codons_with_one_codon_set(["ATG", "CTG"], "ATG", "A")
assert_not_equal "AL", Sequence.translate("ATGCTG")

# Sequence.codons
# - 64 values
assert_equal 64, Sequence.codons.size

# Sequence.amino_acids
# - 20 values
assert_equal 20, Sequence.amino_acids.size
