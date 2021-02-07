# !/usr/bin/env ruby
require "byebug"
require "optparse"

=begin
    Make statistics about codon translation.

    Collect statistics:
        - number of PSMs
        - number of non-redundant peptides
        - mean/ median mass error per PSM
        - number and percentage of identified proteins
        both overall and based on those PSMs containing b/y supported <codon> only

        - number and percentage of recovered <codon> positions
        - mean/ median number of PSMs per <codon> position

        Counts per amino acid:
        - number of PSMs with <codon> translated into this amino acid
        - number of non-redundant peptides with <codon> translated into this amino acid
        - number of <codon> positions translated into this amino acid
        - number of <codon> positions translated into this and another amino acid (= ambiguous translation)

    <Codon> counts are based on those PSMs where <codon> position is supported by b/y-ion fragments. General counts such as total number of PSMs are not, as b/y support gets meaningless when inspecting the PSM as a whole (each PSM has at least one supported position by definition).


    Args:
        input (str): path to input file (enriched evidence file; output from script 02_combine_maxquant_tables)
        codon (str): codon translated by script 01_create_maxquant_dbs into each amino acid
        cdna (str): path input FASTA (cDNA sequences; used as input for 01_create_maxquant_dbs)
        output (str): path to output TXT (statistics described above)
        psm (str): path to output CSV (PSM subset for generating plots)

    Returns
        statistics about the dataset in plain text
        a subset of those PSMs not containing codon at all or containing supported codon

=end

# require .rb files in library (including all subfolders)
Dir[File.join(File.dirname(__FILE__), "lib", "**", "*.rb")].each do |file|
    require File.absolute_path(file)
end

class OptParser
    def self.parse(args)

        options = Hash.new

        # mandatory parameters
        options[:input] = nil
        options[:output] = nil
        options[:codon] = nil
        options[:cdna] = nil

        opt_parser = OptionParser.new do |opts|
            opts.banner = "Make statistics about codon translation."
            opts.separator ""
            opts.separator "Copyright (c) 2020, by Göttingen University"
            opts.separator "Author: Stefanie Mühlhausen"
            opts.separator "This program comes with ABSOLUTELY NO WARRANTY"

            opts.separator ""
            opts.separator "Usage: ruby #{File.basename($PROGRAM_NAME)} -i combined-evidence-msms -o output --codon codon --cdna cdna"

            opts.on("-i", "--input FILE",
                "Path to combinded MaxQuant tables.") do |path|
                FileHelper.file_exist_or_die(path)
                options[:input] = path
            end
            opts.on("--codon CODON",
            "Codon translated into all amino acids." ) do |codon|
                options[:codon] = codon
                unless Sequence.codons.include?(codon)
                    abort "#{codon} not a valid codon."
                end
            end
            opts.on("--cdna FILE",
                "Path to cDNA file (input to 01_create_maxquant_dbs), ",
                "in FASTA format.") do |path|
                FileHelper.file_exist_or_die(path)
                options[:cdna] = path
            end
            opts.on("-o", "--output FILE",
                "Path to statistics output file, in TXT format.") do |path|
                options[:output] = path
            end
            opts.on("-p", "--psm FILE",
                "Path to auxiliary output file containing non-decoy PSM and ",
                "selected information, in CSV format.") do |path|
                options[:psm] = path
            end

            opts.separator ""
            opts.on_tail("-h", "--help", "Show this message") do
                puts opts
                exit
            end

        end

        if args.empty? then
            # display help and exit if program is called without any argument
            abort opt_parser.help
        end

        opt_parser.parse(args)

        # ensure mandatory arguments are present
        abort "Missing mandatory argument: --input" unless options[:input]
        abort "Missing mandatory argument: --codon" unless options[:codon]
        abort "Missing mandatory argument: --cdna" unless options[:cdna]
        abort "Missing mandatory argument: --output" unless options[:output]

        return options
    end
end

options = OptParser.parse(ARGV)

# read in cDNA as reference
all_prots, all_seqs = Sequence.read_fasta(options[:cdna])
all_protein_subset_with_codon = []
all_codon_pos = [] # values: protein-pos
all_prots.each_with_index do |prot, ind|
    seq = all_seqs[ind]
    codons = Sequence.split_cdna_into_codons(seq)

    codons.each_with_index do |codon, ind|
        if codon == options[:codon]
            # collect codon position and protein
            key = prot + "-" + ind.to_s
            all_codon_pos.push(key)

            all_protein_subset_with_codon.push(prot)
        end
    end

end
all_protein_subset_with_codon = all_protein_subset_with_codon.uniq

# collect MaxQuant data
# counts for all PSMs/ proteins, irrepective if they contain codon or not
psms = [] #  reduce to unique values for non-redundant peptide count
mass_errors = []
proteins = []
# counts for PSMs/proteins with codon only
# NOTE - can't use codon_pos to derive PSM counts, as PSM might contain more than one codon pos
psm_subset_with_codon = []
mass_error_subset_with_codon = []
protein_subset_with_codon = []
codon_pos = {} # keys: codon positions, values: number of PSMs per position

fh_psms = File.open(options[:psm], "w")
fh_psms.puts "PSM,Mass error [ppm],#{options[:codon]} translation(s),Has b/y supported #{options[:codon]} position?"
mq_data = ParseEnrichedEvidence.new(options[:codon])
IO.foreach(options[:input]) do |line|
    if mq_data.is_header
        mq_data.parse_header(line)
        next
    end
    mq_data.parse_line(line)

    # update generic counts
    psms.push(mq_data.get_peptide)
    proteins.push(mq_data.get_protein_name_used_in_db)
    unless mq_data.is_masserr_unspecified?
        mass_errors.push(mq_data.get_masserr)
    end

    found_translations = []
    has_supported_codon_pos = false

    # update codon-specific counts, require codon pos to be b/y supported
    if (mq_data.get_supported_pos && mq_data.get_codon_pos).any?
        psm_subset_with_codon.push(mq_data.get_peptide)
        protein_subset_with_codon.push(mq_data.get_protein_name_used_in_db)
        unless mq_data.is_masserr_unspecified?
            mass_error_subset_with_codon.push(mq_data.get_masserr)
        end
        (mq_data.get_supported_pos && mq_data.get_codon_pos).each do |pos|
            pos_in_prot = mq_data.convert_peptide_to_protein_pos(pos)
            key = mq_data.get_protein_name_used_in_db + "-" + pos_in_prot.to_s
            transl = mq_data.get_peptide[pos]
            unless codon_pos[key]
                codon_pos[key] = {}
            end
            unless codon_pos[key][transl]
                codon_pos[key][transl] = []
            end
            codon_pos[key][transl].push(mq_data.get_peptide)

            found_translations.push(transl)
        end

        has_supported_codon_pos = true
    end

    # output PSM
    fh_psms.print "#{mq_data.get_peptide},#{mq_data.get_masserr},"
    fh_psms.print "#{found_translations.join("/")},"
    fh_psms.print "#{has_supported_codon_pos}\n"
end

fh_stats = File.open(options[:output], "w")

fh_stats.print "Total number of PSMs: "
fh_stats.print psms.size.to_s + "\n"

fh_stats.print "Total number of non-redundant peptides: "
fh_stats.print psms.uniq.size.to_s + "\n"

fh_stats.print "Mass error (mean / median): "
fh_stats.print Statistics.mean(mass_errors).round(4).to_s + " / "
fh_stats.print Statistics.median(mass_errors).round(4).to_s + "\n"

fh_stats.print "Number of identified proteins: "
fh_stats.print proteins.uniq.size.to_s + "\n"

fh_stats.print "Percentage of identified proteins: "
fh_stats.print Statistics.percentage(proteins.uniq.size, all_prots.size).round(2).to_s + "\n"
fh_stats.puts ""

fh_stats.puts "The following counts are based on those PSMs where #{options[:codon]} "
fh_stats.puts "position is supported by b/y-type ions"
fh_stats.puts ""

fh_stats.print "Total number of PSMs containing #{options[:codon]}: "
fh_stats.print psm_subset_with_codon.size.to_s + "\n"

fh_stats.print "Total number of non-redundant peptides containing #{options[:codon]}: "
fh_stats.print psm_subset_with_codon.uniq.size.to_s + "\n"

fh_stats.print "Mass error (mean / median) of PSMs containing #{options[:codon]}: "
fh_stats.print Statistics.mean(mass_error_subset_with_codon).round(4).to_s + " / "
fh_stats.print Statistics.median(mass_error_subset_with_codon).round(4).to_s + "\n"

fh_stats.print "Number of identified proteins containing #{options[:codon]}: "
fh_stats.print protein_subset_with_codon.uniq.size.to_s + "\n"

fh_stats.print "Percentage of identified proteins containing #{options[:codon]}: "
fh_stats.print Statistics.percentage(protein_subset_with_codon.uniq.size, all_protein_subset_with_codon.size).round(2).to_s + "\n"
fh_stats.puts ""

fh_stats.print "Number of #{options[:codon]} positions covered by PSMs: "
fh_stats.print codon_pos.keys.uniq.size.to_s + "\n"

fh_stats.print "Total number of #{options[:codon]} positions: "
fh_stats.print all_codon_pos.size.to_s + "\n"

fh_stats.print "Percentage of recovered #{options[:codon]} positions: "
fh_stats.print Statistics.percentage(codon_pos.keys.uniq.size, all_codon_pos.size).round(2).to_s + "\n"
fh_stats.puts ""

fh_stats.print "Number of PSMs per recovered #{options[:codon]} position (mean / median): "
psms_per_pos = codon_pos.collect{|_,v| v.values.flatten.size}
fh_stats.print Statistics.mean(psms_per_pos).round(4).to_s + " / "
fh_stats.print Statistics.median(psms_per_pos).round(4).to_s + "\n"
fh_stats.puts ""

fh_stats.puts "The following counts are broken down to individual translations"
fh_stats.puts ["Translation", "#PSMs", "# non-redundant peptides", "#{options["codon"]} positions", "#{options["codon"]} positions found with additional translations"].join("\t")
Sequence.amino_acids.each do |aa|
    if aa == "I"
        # skip as isoleucine and leucine can't be discriminated
        next
    end
    n_psms = 0
    n_peptides = 0
    n_pos = 0
    n_ambig_transl_pos = 0
    codon_pos.each do |key, val|
        n_psms += val[aa].to_a.size # number of PSMs
        n_peptides += val[aa].to_a.uniq.size # number of non-redundant peptides
        if val[aa].to_a.size > 0
            n_pos += 1 # number of codon positions
            if val.keys.size > 1
                n_ambig_transl_pos += 1 # number of codon positions with multiple translations (including aa)
            end
        end
    end
    if aa == "L"
        # adjust for output
        aa = "L/I"
    end
    fh_stats.puts [aa, n_psms, n_peptides, n_pos, n_ambig_transl_pos].join("\t")
end
fh_stats.close
fh_psms.close
