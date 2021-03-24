# !/usr/bin/env ruby
require "test/unit/assertions"
include Test::Unit::Assertions
require "tempfile"

=begin
    Integration test 03_make_statistics:

        - Input: mockup data
        - expected output:
            -

        TODO
=end

# require .rb files in library (including all subfolders)
Dir[File.join(__dir__, "..", "lib", "**", "*.rb")].each do |file|
    require File.absolute_path(file)
end

# combined evidence and msms tables
header = "cDNA\tStart pos in protein\tEnd pos in protein\tCodon pos\tOriginal protein name\tb/y-ion supported pos\tCorresponding scan numbers\tSequence\tLength\tModifications\tModified sequence\tOxidation (M) Probabilities\tOxidation (M) Score Diffs\tAcetyl (Protein N-term)\tOxidation (M)\tMissed cleavages\tProteins\tLeading Proteins\tLeading Razor Protein\tType\tRaw file\tFraction\tExperiment\tMS/MS m/z\tCharge\tm/z\tMass\tResolution\tUncalibrated - Calibrated m/z [ppm]\tUncalibrated - Calibrated m/z [Da]\tMass Error [ppm]\tMass Error [Da]\tUncalibrated Mass Error [ppm]\tUncalibrated Mass Error [Da]\tMax intensity m/z 0\tRetention time\tRetention length\tCalibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\tRetention time calibration\tMatch time difference\tMatch m/z difference\tMatch q-value\tMatch score\tNumber of data points\tNumber of scans\tNumber of isotopic peaks\tPIF\tFraction of total spectrum\tBase peak fraction\tPEP\tMS/MS Count\tMS/MS Scan Number\tScore\tDelta score\tCombinatorics\tIntensity\tReverse\tPotential contaminant\tid\tProtein group IDs\tPeptide ID\tMod. peptide ID\tMS/MS IDs\tBest MS/MS\tAIF MS/MS IDs\tOxidation (M) site IDs"
line = "GCCGCCGCAGCACTGGGGGCTGCGCTTGCGCCGCAGCGC\t110\t122\t4\tEEQ37961 cdna supercontig:ASM383v1:CH408077:1767990:1769492:-1 gene:CLUG_02084 gene_biotype:protein_coding transcript_biotype:protein_coding description:hypothetical protein\t0;1;2;3;4;5;6;7;8;9;10;11;12\t8270\tAAAASGAALAPQR\t13\tUnmodified\t_AAAASGAALAPQR_\t0\t0\t0\tg1_102_127-CTG-S\tg1_102_127-CTG-S\tg1_102_127-CTG-S\tMULTI-MSMS\tH_Schmitt_211216_190117_L1_R1_14\t14\tL1\t5.778.184\t2\t577.817.461\t115.362.037\t42112.35\t-0.42427\t-0.00024515\t0.91792\t0.00053039\t0.49364\t0.00028524\t577.818.033.203.664\t16.334\t0.21083\t16.472\t16.339\t16.549\t0.13784\t62\t35\t2\t0\t0\t0\t1,28E-04\t1\t8270\t116.24\t99.661\t1\t16324000\t63\t1435\t11\t11\t51\t51"
mockup_data_enriched_evidence = Tempfile.new("gcf")
mockup_data_enriched_evidence.write("#{header}\n#{line}")
mockup_data_enriched_evidence.close

mockup_data_cdna = File.join(__dir__, "..", "sample_data", "Clavispora_cDNA_excerpt.fasta")
output = Tempfile.new("gcf")
output_psms = Tempfile.new("gcf")

original_stdout = $stdout.clone
$stdout.reopen(File.new('/dev/null', 'w'))

system("ruby", File.join(__dir__, "..", "03_make_statistics.rb"), "-i", mockup_data_enriched_evidence.path, "--cdna", mockup_data_cdna, "-o", output.path, "-p", output_psms.path, "--codon", "CTG")
$stdout.reopen(original_stdout)


require "byebug"

is_first_line = true
is_translation_part = false
IO.foreach(output.path) do |line|
    desc, val = line.split(":")
    # debugger
    # puts "??"
    if is_first_line
        assert_equal "Total number of PSMs", desc
        # assert_equal "1", val
        is_first_line = false
    end

    # # translation table sub-part of file
    # if is_translation_part
    #     # either serine line or all 0 lines
    #     if desc.starts_with?("S")
    #         assert_equal "S 1   1   1   0", desc
    #     else
    #         assert_match /[A-Z\/](\t0){4}/, desc
    #     end
    # end
    # if desc.starts_with?("Translation")
    #     is_translation_part
    # end
end
# assert_true is_translation_part # should have seen translation part in file

ind_psm, ind_masserr, ind_transl, ind_flag = nil, nil, nil, nil
is_first_line_after_header = false
IO.foreach(output.path) do |line|
    parts = line.split(",")
    is_first_line_after_header = false

    if is_first_line_after_header
        # assert equals on first line.
        assert_equal "AAAASGAALAPQR", parts[ind_psm]
        assert_equal "0.49364", parts[ind_masserr]
        assert_equal "S", parts[ind_transl]
        debugger
        puts "??"
        assert_equal "true", parts[ind_flag]
    else
        # header line
        ind_psm = parts.index("PSM")
        ind_masserr = parts.index("Mass error [ppm]")
        ind_transl = parts.index("CTG translation(s)")
        ind_flag = parts.index("Has b/y supported CTG position?")
        is_first_line_after_header = true
    end
end

mockup_data_enriched_evidence.unlink
output.unlink
output_psms.unlink
