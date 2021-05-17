# !/usr/bin/env ruby
require "test/unit/assertions"
include Test::Unit::Assertions
require "tempfile"

=begin
    Integration test 02_combine_maxquant_tables:
        - Input: mockup data
        - expected output:
            evidence file enriched with
            - peptide information
            - msms data

=end

# require .rb files in library
Dir[File.join(__dir__, "..", "lib", "*.rb")].each do |file|
    require file
end

# dictionary mapping original gene name onto standardised one
orig_header = "EEQ37961 cdna supercontig:ASM383v1:CH408077:1767990:1769492:-1 gene:CLUG_02084 gene_biotype:protein_coding transcript_biotype:protein_coding description:hypothetical protein"
mockup_data_map = Tempfile.new("gcf")
mockup_data_map.write("#{orig_header};g1")
mockup_data_map.close
mockup_data_evidence = File.join(__dir__, "..", "sample_data", "evidence.txt")
mockup_data_msms = File.join(__dir__, "..", "sample_data", "msms.txt")
mockup_data_cdna = File.join(__dir__, "..", "sample_data", "Clavispora_cDNA_excerpt.fasta")
output = Tempfile.new("gcf")

original_stdout = $stdout.clone
$stdout.reopen(File.new('/dev/null', 'w'))

system("ruby", File.join(__dir__, "..", "02_combine_maxquant_tables.rb"), "-e", mockup_data_evidence, "--msms", mockup_data_msms, "--cdna", mockup_data_cdna, "--map", mockup_data_map.path, "-o", output.path)
$stdout.reopen(original_stdout)

ind_cdna, ind_ctg_pos, ind_orig_header, ind_msms = nil, nil, nil, nil
is_first_line = true

IO.foreach(output.path) do |line|
    parts = line.split("\t")
    if is_first_line
        # header line
        ind_cdna = parts.index("cDNA")
        ind_ctg_pos = parts.index("CTG codon pos")
        ind_orig_header = parts.index("Original protein name")
        ind_msms = parts.index("Corresponding scan numbers")
        is_first_line = false
    else
        # test first line after header
        # peptide AAAALGAALAPQR
        # contains one CTG at index 4
        # corresponding Scan number is 8270
        assert_true ind_cdna.is_a? Integer
        assert_true ind_orig_header.is_a? Integer
        assert_true ind_ctg_pos.is_a? Integer
        assert_true ind_msms.is_a? Integer
        assert_equal "AAAALGAALAPQR", Sequence.translate(parts[ind_cdna])
        assert_equal orig_header, parts[ind_orig_header]
        assert_equal "4", parts[ind_ctg_pos]
        assert_equal "8270", parts[ind_msms]
        break
    end
end
assert_false is_first_line

mockup_data_map.unlink
output.unlink
