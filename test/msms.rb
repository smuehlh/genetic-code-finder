# !/usr/bin/env ruby
require "test/unit/assertions"
include Test::Unit::Assertions
require "tempfile"

=begin
    Test lib/parse_msms.rb:

        - parse header and subsequent lines
        - parse evidenceid
        - parse peptide
        - parse matched ions
        - parse scan number
        - identify b/y-ion supported positions

=end

# require .rb files in library (including all subfolders)
Dir[File.join(__dir__, "..", "lib", "**", "*.rb")].each do |file|
    require File.absolute_path(file)
end

header = "Raw file\tScan number\tScan index\tSequence\tLength\tMissed cleavages\tModifications\tModified sequence\tOxidation (M) Probabilities\tOxidation (M) Score Diffs\tAcetyl (Protein N-term)\tOxidation (M)\tProteins\tCharge\tFragmentation\tMass analyzer\tType\tScan event number\tIsotope index\tm/z\tMass\tMass Error [ppm]\tSimple Mass Error [ppm]\tRetention time\tPEP\tScore\tDelta score\tScore diff\tLocalization prob\tCombinatorics\tPIF\tFraction of total spectrum\tBase peak fraction\tPrecursor Full ScanNumber\tPrecursor Intensity\tPrecursor Apex Fraction\tPrecursor Apex Offset\tPrecursor Apex Offset Time\tMatches\tIntensities\tMass Deviations [Da]\tMass Deviations [ppm]\tMasses\tNumber of Matches\tIntensity coverage\tPeak coverage\tNeutral loss level\tETD identification type\tReverse\tAll scores\tAll sequences\tAll modified sequences\tid\tProtein group IDs\tPeptide ID\tMod. peptide ID\tEvidence ID\tOxidation (M) site IDs"

line = "H_Schmitt_211216_190117_L1_R1_14\t8270\t2904\tAAAASGAALAPQR\t13\t0\tUnmodified\t_AAAASGAALAPQR_\t\t\t0\t0\tg1_102_127-CTG-S\t2\tHCD\tFTMS\tMULTI-MSMS\t3\t0\t577.81746\t1153.6204\t0.91792\t1.6273642\t16.285\t1.2802E-08\t116.24\t99.661\tNaN\tNaN\t1\t0\t0\t0\t-1\t0\t0\t0\t0\ty1;y3;y4;y5;y6;y7;y8;y9;y10;y11;y1-NH3;y3-NH3;y4-NH3;y9-H2O;y9-NH3;y10-H2O;y10(2+);y12(2+);a2;b2;b3;b4;b6;b6-H2O;b7;b7-H2O;b8;b11\t92159;124870.9;137909.8;120648.3;140553.2;52794.8;179005.9;241070.8;83605.9;37491.4;16749.2;92401.2;20138;3179.3;10447.9;4430.2;137909.8;4627.5;86763.5;210341.4;215262.9;40058;24701.5;30206.8;17482.8;17531.7;9207.6;3179.3\t-0.0004631022;0.0002068901;0.0002963727;-0.0002496136;0.0002206965;-0.001487416;-0.000554896;-0.0007182999;0.0006209789;-0.001324622;-0.0002788053;-0.0005707426;-0.001878586;-0.004762337;-0.005732103;-0.001462025;-0.005320321;0.002578742;-0.0004543421;-0.0001915144;-0.0002481368;0.001437485;0.0007886927;-0.0006650063;0.002619525;0.0007264661;-0.001460782;-0.01599572\t-2.644494;0.5169278;0.6288847;-0.4271633;0.3367416;-2.047578;-0.7082745;-0.825177;0.6595524;-1.308198;-1.763556;-1.489395;-4.135643;-5.586493;-6.716324;-1.583123;-11.2894;4.755228;-3.947813;-1.338497;-1.158874;5.041078;1.837552;-1.617236;5.236497;1.506456;-2.557012;-18.7639\t175.119415283203;400.230086654317;471.267110959512;584.351720926237;655.388364403903;726.427186304559;783.447717507831;870.479909321653;941.515683830579;1012.55474321948;158.092681884766;383.204315185547;454.242736816406;852.473388671875;853.458374023438;923.507202148438;471.267110959512;542.296325683594;115.087043762207;143.081695556641;214.118865966797;285.154294132757;429.208435058594;411.199324071272;500.243718014252;482.235046386719;571.284912109375;852.473388671875\t28\t0.4874856\t0.1313131\tNone\tUnknown\t\t116.2378;16.5771;15.6294\tAAAASGAALAPQR;RGAEEQLVPR;AANPRARWGR\t_AAAASGAALAPQR_;_RGAEEQLVPR_;_AANPRARWGR_\t51\t1435\t11\t11\t63\t"


mq_data = ParseMsms.new()

# ParseMsms.is_header & parse_line
# - treat (only) first line as header line
# - identify Sequence, evidence id, Matches and Scan number rows
assert_true mq_data.is_header
mq_data.parse_header(header)
assert_false mq_data.is_header

mq_data.parse_line(header)
assert_equal "Sequence", mq_data.get_peptide
assert_equal "Evidence ID", mq_data.get_evidenceid
assert_equal ["Matches"], mq_data.get_matched_ions # splits str at ';'
assert_equal "Scan number", mq_data.get_scannumber

# ParseMsms.get_*
# - peptide is AAAASGAALAPQR
# - evidence id is 63
# - matched ions are 'y1;y3;y4;y5;y6;y7;y8;y9;y10;y11;y1-NH3;y3-NH3;y4-NH3;y9-H2O;y9-NH3;y10-H2O;y10(2+);y12(2+);a2;b2;b3;b4;b6;b6-H2O;b7;b7-H2O;b8;b11' split at ';'
# - scannumber is 8270
# - all positions are supported
mq_data.parse_line(line)
assert_equal "AAAASGAALAPQR", mq_data.get_peptide
assert_equal "63", mq_data.get_evidenceid
assert_equal ["y1", "y3", "y4", "y5", "y6", "y7", "y8", "y9", "y10", "y11", "y1-NH3", "y3-NH3", "y4-NH3", "y9-H2O", "y9-NH3", "y10-H2O", "y10(2+)", "y12(2+)", "a2", "b2", "b3", "b4", "b6", "b6-H2O", "b7", "b7-H2O", "b8", "b11"], mq_data.get_matched_ions
assert_equal "8270", mq_data.get_scannumber
assert_equal (0...mq_data.get_peptide.size).to_a, mq_data.get_positions_with_b_y_type_support
