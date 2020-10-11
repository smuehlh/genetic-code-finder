# !/usr/bin/env ruby
require "test/unit/assertions"
include Test::Unit::Assertions
require "tempfile"

=begin
    Test lib/parse_evidence.rb:

        - parse header and subsequent lines
        - identify decoy match
        - parse evidenceid
        - parse protein, protein region and applied special translation
        - parse score
        - parse mass error

=end

# require .rb files in library (including all subfolders)
Dir[File.join(__dir__, "..", "lib", "**", "*.rb")].each do |file|
    require File.absolute_path(file)
end

header = "Sequence\tLength\tModifications\tModified sequence\tOxidation (M) Probabilities\tOxidation (M) Score Diffs\tAcetyl (Protein N-term)\tOxidation (M)\tMissed cleavages\tProteins\tLeading Proteins\tLeading Razor Protein\tType\tRaw file\tFraction\tExperiment\tMS/MS m/z\tCharge\tm/z\tMass\tResolution\tUncalibrated - Calibrated m/z [ppm]\tUncalibrated - Calibrated m/z [Da]\tMass Error [ppm]\tMass Error [Da]\tUncalibrated Mass Error [ppm]\tUncalibrated Mass Error [Da]\tMax intensity m/z 0\tRetention time\tRetention length\tCalibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\tRetention time calibration\tMatch time difference\tMatch m/z difference\tMatch q-value\tMatch score\tNumber of data points\tNumber of scans\tNumber of isotopic peaks\tPIF\tFraction of total spectrum\tBase peak fraction\tPEP\tMS/MS Count\tMS/MS Scan Number\tScore\tDelta score\tCombinatorics\tIntensity\tReverse\tPotential contaminant\tid\tProtein group IDs\tPeptide ID\tMod. peptide ID\tMS/MS IDs\tBest MS/MS\tAIF MS/MS IDs\tOxidation (M) site IDs"
decoy1 = "AACLLPK\t7\tUnmodified\t_AACLLPK_\t\t\t0\t0\t0\tCON__P02768-1\tCON__P02768-1\tCON__P02768-1\tMULTI-MSMS\tH_Schmitt_211216_190117_L1_R1_13\t13\tL1\t3.867.227\t2\t386.722.924\t771.431.295\t47187.57\t-0.054016\t-2,09E-01\t0.47997\t0.00018562\t0.42595\t0.00016473\t386.722.233.403.272\t22.474\t12.711\t22.575\t21.736\t23.007\t0.1014\t\t\t\t\t229\t142\t3\t0\t0\t0\t0.020385\t1\t12770\t84.359\t72.638\t1\t164120000\t\t+\t187\t14\t45\t45\t159\t159\t\t"
decoy2 = "AEANLAANK\t9\tUnmodified\t_AEANLAANK_\t\t\t0\t0\t0\t\tREV__g5027_0_1011\tREV__g5027_0_1011\tMULTI-MSMS\tH_Schmitt_211216_190117_L1_R1_16\t16\tL1\t4.512.402\t2\t451.240.524\t900.466.495\t43165.19\t-0.076729\t-3,46E-01\t0.4042\t0.00018239\t0.32748\t0.00014777\t451.240.307.402.434\t19.142\t0.37056\t19.346\t19.218\t19.589\t0.20444\t\t\t\t\t88\t46\t3\t0\t0\t0\t0.050821\t1\t10200\t48.794\t21.491\t1\t66807000\t+\t\t2190\t9522\t479\t496\t1751\t1751\t\t"
line = "AAAASGAALAPQR\t13\tUnmodified\t_AAAASGAALAPQR_\t\t\t0\t0\t0\tg1_102_127-CTG-S\tg1_102_127-CTG-S\tg1_102_127-CTG-S\tMULTI-MSMS\tH_Schmitt_211216_190117_L1_R1_14\t14\tL1\t5.778.184\t2\t577.817.461\t115.362.037\t42112.35\t-0.42427\t-0.00024515\t0.91792\t0.00053039\t0.49364\t0.00028524\t577.818.033.203.664\t16.334\t0.21083\t16.472\t16.339\t16.549\t0.13784\t\t\t\t\t62\t35\t2\t0\t0\t0\t1,28E-04\t1\t8270\t116.24\t99.661\t1\t16324000\t\t\t63\t1435\t11\t11\t51\t51\t\t"


mq_data = ParseEvidence.new("CTG")

# ParseEvidence.is_header & parse_line
# - treat (only) first line as header line
# - identify Sequence, Proteins, id rows
# no use testing for
#   - Score and Mass error [ppm], as those are converted to float
#   - Reverse and Potential Contaminant, as those are tested elsewhere
assert_true mq_data.is_header
mq_data.parse_header(header)
assert_false mq_data.is_header

mq_data.parse_line(header)
assert_equal "Sequence", mq_data.get_peptide
assert_equal "Proteins", mq_data.get_protein
assert_equal "id", mq_data.get_evidenceid

# ParseEvidence.get_*
# - peptide is AAAASGAALAPQR
# - evidence id is 63
# g1_102_127-CTG-S
# - protein is g1
# - protein region is 102 - 127
# - codon and special translation is CTG / S
# - score is 116.24
# - mass error is specified and is 0.91792
mq_data.parse_line(line)
assert_equal "AAAASGAALAPQR", mq_data.get_peptide
assert_equal "63", mq_data.get_evidenceid
assert_equal "g1", mq_data.get_protein
assert_equal [102, 127], mq_data.get_protein_region
assert_equal ["CTG", "S"], mq_data.get_codon_transl_applied_to_protein
assert_equal 116.24, mq_data.get_score
assert_equal 0.91792, mq_data.get_masserr
assert_false mq_data.is_masserr_unspecified?
assert_equal mq_data.get_peptide, mq_data.get_peptide_with_unknown_aminoacids_denoted_as_X

# - denote unknown amino acids as X rather than U
mq_data.parse_line(line.sub("AAAASGAALAPQR", "AAAASGUALAPQR"))
assert_equal "AAAASGUALAPQR", mq_data.get_peptide
assert_equal "AAAASGXALAPQR", mq_data.get_peptide_with_unknown_aminoacids_denoted_as_X
assert_not_equal mq_data.get_peptide, mq_data.get_peptide_with_unknown_aminoacids_denoted_as_X

# ParseEvidence.is_decoy_match
# - identity (only) matches in Reverse and Contaminant DBs as decoy
mq_data.parse_line(decoy1)
assert_true mq_data.is_decoy_match
mq_data.parse_line(decoy2)
assert_true mq_data.is_decoy_match
mq_data.parse_line(line)
assert_false mq_data.is_decoy_match
