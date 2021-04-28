# !/usr/bin/env ruby
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

        Counts per gene:
        - sequence coverage
        - number of PSMs
        - number of non-redundant peptides
        - number of PSMs with <codon>
        - <codon> coverage

    <Codon> counts are based on those PSMs where <codon> position is supported by b/y-ion fragments. General counts such as total number of PSMs are not, as b/y support gets meaningless when inspecting the PSM as a whole (each PSM has at least one supported position by definition).

    Retrieve designated <codon> from enriched evidence file, it should be noted there in header as "<codon> codon pos".

    To plot found translations, call Rscript 03_plot_translation.R on statistical output.

    Args:
        input (str): path to input file (enriched evidence file; output of script 02_combine_maxquant_tables)
        cdna (str): path input FASTA (cDNA sequences; used as input for 01_create_maxquant_dbs)
        output (str): path to output TXT (statistics described above)
        genes (str): path to output CSV (gene statistics described above)

    Returns
        statistics about the dataset in plain text
        statistics about recovered genes in CSV format

=end

# require .rb files in library (including all subfolders)
Dir[File.join(File.dirname(__FILE__), "lib", "**", "*.rb")].each do |file|
    require File.absolute_path(file)
end
# also require 01_create_maxquant_db to use method simplify_header()
path_to_01_script = Dir[File.join(File.dirname(__FILE__), "01*.rb")][0]
begin
    require File.absolute_path(path_to_01_script)
rescue OptionParser::InvalidOption, OptionParser::AmbiguousOption
end

class OptParser
    def self.parse(args)

        options = Hash.new

        # mandatory parameters
        options[:input] = nil
        options[:output] = nil
        options[:genes] = nil
        options[:cdna] = nil

        opt_parser = OptionParser.new do |opts|
            opts.banner = "Make statistics about codon translation."
            opts.separator ""
            opts.separator "Copyright (c) 2020-2021, by Göttingen University"
            opts.separator "Author: Stefanie Mühlhausen"
            opts.separator "This program comes with ABSOLUTELY NO WARRANTY"

            opts.separator ""
            opts.separator "Usage: ruby #{File.basename($PROGRAM_NAME)} -i combined-evidence-msms -c cdna -o stats-output -g gene-stats-output"

            opts.on("-i", "--input FILE",
                "Path to input file enriched evidence",
                "(output of 02_combine_maxquant_tables).") do |path|
                FileHelper.file_exist_or_die(path)
                options[:input] = path
            end
            opts.on("-c", "--cdna FILE",
                "Path to cDNA file (input to 01_create_maxquant_dbs), ",
                "in FASTA format.") do |path|
                FileHelper.file_exist_or_die(path)
                options[:cdna] = path
            end
            opts.on("-o", "--output FILE",
                "Path to statistics output file, in TXT format.") do |path|
                options[:output] = path
            end
            opts.on("-g", "--genes FILE",
                "Path to auxiliary output file containing recovered genes, ",
                "in CSV format.") do |path|
                options[:genes] = path
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
        abort "Missing mandatory argument: --cdna" unless options[:cdna]
        abort "Missing mandatory argument: --output" unless options[:output]
        abort "Missing mandatory argument: --genes" unless options[:genes]

        return options
    end
end

def get_codon_db_was_prepared_for(path)
    codon = ParseEnrichedEvidence.codon_used_for_db_preparation(path)
    unless Sequence.codons.include?(codon)
        abort "Codon uses as basis for DB (#{codon}) not a valid codon."
    end
    codon
end

options = OptParser.parse(ARGV)
designated_codon = get_codon_db_was_prepared_for(options[:input])

# read in cDNA as reference
all_prots, all_seqs = Sequence.read_fasta(options[:cdna])
proteins_with_coverage_data = {} # use to collect PSMs, too
total_num_prots_with_codon = 0 # all proteins with designated codon
total_num_codon_pos = 0 # total number of codon positions
all_prots.each_with_index do |prot, ind|
    patched_name = simplify_header(ind+1)
    seq = all_seqs[ind]
    codons = Sequence.split_cdna_into_codons(seq)

    proteins_with_coverage_data[patched_name] = {
        len: codons.size,  # total length
        covered_pos: [], # found positions, any support
        codon_pos: 0, # total codon pos
        b_y_covered_codon_pos: [], # found codon pos, only b/y supported
        psms: [], # found PSMs, any support
        psms_codon: [], # found PSMs covering codon pos, only b/y supported
    }
    codons.each_with_index do |codon, ind|
        if codon == designated_codon
            # collect codon position
            total_num_codon_pos += 1
            proteins_with_coverage_data[patched_name][:codon_pos] += 1
        end
    end
    if proteins_with_coverage_data[patched_name][:codon_pos] > 0
        total_num_prots_with_codon += 1
    end
end

# collect MaxQuant data
# counts for all PSMs/ proteins, irrepective if they contain codon or not
mass_errors = []

# counts for PSMs/proteins with codon only, codon pos must be b/y supported
# NOTE - can't use codon_transl to derive PSM counts, as PSM might contain more than one codon pos
mass_error_subset_with_codon = []
codon_transl = {} # keys: codon positions, values: PSMs per translation

mq_data = ParseEnrichedEvidence.new()
IO.foreach(options[:input]) do |line|
    if mq_data.is_header
        mq_data.parse_header(line)
        next
    end
    mq_data.parse_line(line)

    protein = mq_data.get_protein_name_used_in_db
    psm = mq_data.get_peptide
    all_pos = (mq_data.get_peptide_start..mq_data.get_peptide_stop).to_a

    # update generic counts, irrespective of b/y support
    proteins_with_coverage_data[protein][:psms].push(psm)
    proteins_with_coverage_data[protein][:covered_pos].push(*all_pos)
    unless mq_data.is_masserr_unspecified?
        mass_errors.push(mq_data.get_masserr)
    end

    # update codon-specific counts, require codon pos to be b/y supported
    supported_codon_pos = mq_data.get_codon_pos & mq_data.get_supported_pos
    supported_codon_pos.each do |peptide_pos|
        pos = mq_data.convert_peptide_to_protein_pos(peptide_pos)
        key = protein + "-" + pos.to_s
        transl = psm[peptide_pos]

        unless codon_transl[key]
            codon_transl[key] = {}
        end
        unless codon_transl[key][transl]
            codon_transl[key][transl] = []
        end
        codon_transl[key][transl].push(psm)

        proteins_with_coverage_data[protein][:b_y_covered_codon_pos].push(pos)
    end
    if supported_codon_pos.any?
        # update the following counts only once, even if PSM contains multiple supported codon pos

        proteins_with_coverage_data[protein][:psms_codon].push(psm)
        unless mq_data.is_masserr_unspecified?
            mass_error_subset_with_codon.push(mq_data.get_masserr)
        end
    end
end

# prepare for output
proteins = proteins_with_coverage_data.filter_map{|k,v| k if v[:psms].any?}
psms = proteins_with_coverage_data.filter_map{|k,v| v[:psms]}.flatten
protein_subset_with_codon = proteins_with_coverage_data.filter_map{|k,v| k if v[:psms_codon].any?}
psm_subset_with_codon = proteins_with_coverage_data.filter_map{|k,v| v[:psms_codon]}.flatten
num_covered_pos = codon_transl.keys.uniq.size

# write statistics output
fh_stats = File.open(options[:output], "w")

fh_stats.print "Total number of PSMs: "
fh_stats.print psms.size.to_s + "\n"

fh_stats.print "Total number of non-redundant peptides: "
fh_stats.print psms.uniq.size.to_s + "\n"

fh_stats.print "Mass error (mean / median): "
fh_stats.print Statistics.mean(mass_errors).round(4).to_s + " / "
fh_stats.print Statistics.median(mass_errors).round(4).to_s + "\n"

fh_stats.print "Number of identified proteins: "
fh_stats.print proteins.size.to_s + "\n"

fh_stats.print "Percentage of identified proteins: "
fh_stats.print Statistics.percentage(proteins.size, all_prots.size).round(2).to_s + "\n"
fh_stats.puts ""

fh_stats.puts "The following counts are based on those PSMs where #{designated_codon} "
fh_stats.puts "position is supported by b/y-type ions"
fh_stats.puts ""

fh_stats.print "Total number of PSMs containing #{designated_codon}: "
fh_stats.print psm_subset_with_codon.size.to_s + "\n"

fh_stats.print "Total number of non-redundant peptides containing #{designated_codon}: "
fh_stats.print psm_subset_with_codon.uniq.size.to_s + "\n"

fh_stats.print "Mass error (mean / median) of PSMs containing #{designated_codon}: "
fh_stats.print Statistics.mean(mass_error_subset_with_codon).round(4).to_s + " / "
fh_stats.print Statistics.median(mass_error_subset_with_codon).round(4).to_s + "\n"

fh_stats.print "Number of identified proteins containing #{designated_codon}: "
fh_stats.print protein_subset_with_codon.size.to_s + "\n"

fh_stats.print "Percentage of identified proteins containing #{designated_codon}: "
fh_stats.print Statistics.percentage(protein_subset_with_codon.size, total_num_prots_with_codon).round(2).to_s + "\n"
fh_stats.puts ""

fh_stats.print "Number of #{designated_codon} positions covered by PSMs: "
fh_stats.print num_covered_pos.to_s + "\n"

fh_stats.print "Total number of #{designated_codon} positions: "
fh_stats.print total_num_codon_pos.to_s + "\n"

fh_stats.print "Percentage of recovered #{designated_codon} positions: "
fh_stats.print Statistics.percentage(num_covered_pos, total_num_codon_pos).round(2).to_s + "\n"
fh_stats.puts ""

fh_stats.print "Number of PSMs per recovered #{designated_codon} position (mean / median): "
psms_per_pos = codon_transl.collect{|_,v| v.values.flatten.size}
fh_stats.print Statistics.mean(psms_per_pos).round(4).to_s + " / "
fh_stats.print Statistics.median(psms_per_pos).round(4).to_s + "\n"
fh_stats.puts ""

fh_stats.puts "The following counts are broken down to individual translations"
fh_stats.puts ["Translation", "#PSMs", "# non-redundant peptides", "#{designated_codon} positions", "#{designated_codon} positions found with additional translations"].join("\t")
Sequence.amino_acids.each do |aa|
    if aa == "I"
        # skip as isoleucine and leucine can't be discriminated
        next
    end
    n_psms = 0
    n_peptides = 0
    n_pos = 0
    n_ambig_transl_pos = 0
    codon_transl.each do |key, val|
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

# write genes output
fh_genes = File.open(options[:genes], "w")
fh_genes.puts ["Gene", "Sequence coverage [%]", "#PSMs", "#non-redundant peptides", "#PSMs containing #{designated_codon}", "#{designated_codon} coverage [%]"].join(",")

proteins_with_coverage_data.each do |prot, data|
    next if data[:psms].empty?

    seq_coverage = Statistics.percentage(data[:covered_pos].uniq.size, data[:len])
    codon_coverage = Statistics.percentage(data[:b_y_covered_codon_pos].uniq.size, data[:codon_pos])
    n_psms = data[:psms].size
    n_non_red_peptides = data[:psms].uniq.size
    n_psms_codon = data[:psms_codon].size

    fh_genes.puts [prot, seq_coverage, n_psms, n_non_red_peptides, n_psms_codon, codon_coverage].join(",")
end

fh_stats.close
fh_genes.close
