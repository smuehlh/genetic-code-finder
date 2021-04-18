# !/usr/bin/env ruby
require "optparse"

=begin
    Combine MaxQuant evidence and msms tables.

    Enrich evidence table
        - skip lines matching the decoy database (which MaxQuant automatically generates and searches against)
        - add information about b/y-ion support
        - add original gene name (was simplified for DB creation)
        - add cDNA of peptide sequnce
        - add start and end position of peptide in gene
        - add position of codon <codon> in peptide

    Args:
        evidence (str): path to input file (MaxQuant evidence.txt)
        msms (str): path to input file (MaxQuant msms.txt)
        map (str): path to input file (mapping between original and shortened FASTA headers)
        output (str): path to output TSV (enriched evidence file)
        codon (str): codon translated by script 01_create_maxquant_dbs into each amino acid
        cdna (str): path input FASTA (cDNA sequences; used as input for 01_create_maxquant_dbs)

    Returns
        an enriched evidence table

=end

# require .rb files in library (including all subfolders)
Dir[File.join(File.dirname(__FILE__), "lib", "**", "*.rb")].each do |file|
    require File.absolute_path(file)
end

class OptParser
    def self.parse(args)

        options = Hash.new

        # mandatory parameters
        options[:evidence] = nil
        options[:msms] = nil
        options[:map] = nil
        options[:output] = nil
        options[:codon] = nil
        options[:cdna] = nil

        opt_parser = OptionParser.new do |opts|
            opts.banner = "Combine MaxQuant tables evidence and msms."
            opts.separator ""
            opts.separator "Copyright (c) 2020-2021, by Göttingen University"
            opts.separator "Author: Stefanie Mühlhausen"
            opts.separator "This program comes with ABSOLUTELY NO WARRANTY"

            opts.separator ""
            opts.separator "Usage: ruby #{File.basename($PROGRAM_NAME)} -e evidence -m msms -o output --map dict --codon codon --cdna cdna"

            opts.on("-e", "--evidence FILE",
                "Path to input file evidence.txt.") do |path|
                FileHelper.file_exist_or_die(path)
                options[:evidence] = path
            end
            opts.on("-m", "--msms FILE",
                "Path to input file msms.txt.") do |path|
                FileHelper.file_exist_or_die(path)
                options[:msms] = path
            end
            opts.on("--map FILE",
                "Path to input file mapping original ",
                "to shortened FASTA headers.") do |path|
                options[:map] = path
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
                "Path to output file, in TSV format.") do |path|
                options[:output] = path
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
        abort "Missing mandatory argument: --evidence" unless options[:evidence]
        abort "Missing mandatory argument: --msms" unless options[:msms]
        abort "Missing mandatory argument: --map" unless options[:map]
        abort "Missing mandatory argument: --codon" unless options[:codon]
        abort "Missing mandatory argument: --cdna" unless options[:cdna]
        abort "Missing mandatory argument: --output" unless options[:output]

        return options
    end
end

def read_msms(path)
    evids_with_msms_data = {}
    mq_data = ParseMsms.new()
    IO.foreach(path) do |line|
        if mq_data.is_header
            mq_data.parse_header(line)
            next
        end
        mq_data.parse_line(line)
        evid = mq_data.get_evidenceid
        scannr = mq_data.get_scannumber
        supported_pos = mq_data.get_positions_with_b_y_type_support

        unless evids_with_msms_data[evid]
            evids_with_msms_data[evid] = {scannrs: [], supported: []}
        end

        evids_with_msms_data[evid][:scannrs].push scannr
        evids_with_msms_data[evid][:supported] |= supported_pos
    end
    evids_with_msms_data
end

def read_gene_names_map(path)
    map = {}
    IO.foreach(path) do |line|
        line = line.chomp
        orig, _, patch = line.rpartition(";")
        map[patch] = orig
    end
    map
end

def get_corresponding_msms_data(mq_data, evids_with_msms_data)
    evid = mq_data.get_evidenceid
    supported = evids_with_msms_data[evid][:supported] || []
    scannrs = evids_with_msms_data[evid][:scannrs] || []

    [supported, scannrs]
end

def get_corresponding_original_genename(mq_data, prots_with_orig_names)
    prots_with_orig_names[mq_data.get_protein]
end

def get_protein_cdna(original_name, headers, seqs)
    ind = headers.index(original_name)
    seqs[ind]
end

def map_peptide_onto_cdna(mq_data, cdna)
    peptide = mq_data.get_peptide_with_unknown_aminoacids_denoted_as_X
    region_start, region_stop = mq_data.get_protein_region
    codon, special_transl = mq_data.get_codon_transl_applied_to_protein

    # map cdna to codons and reduce to rough protein region
    codons = Sequence.split_cdna_into_codons(cdna)
    region = codons[region_start..region_stop]
    transl =
        if codon
            Sequence.translate_codons_with_one_codon_set(region, codon, special_transl)
        else
            Sequence.translate(region.join(""))
        end

    # map rough region to actual peptide position in protein
    peptide_start_in_region = transl.index(peptide)
    peptide_start_in_protein = region_start + peptide_start_in_region
    peptide_stop_in_protein = peptide_start_in_protein + peptide.size - 1
    peptide_codons = codons[peptide_start_in_protein..peptide_stop_in_protein]

    [peptide_codons, peptide_start_in_protein, peptide_stop_in_protein]
end

def get_codon_pos(codons, special_codon)
    codons.each_with_index.collect do |this_codon, ind|
        ind if this_codon == special_codon
    end.compact
end

def get_corresponding_msms_data(mq_data, evids_with_msms_data)
    msms = evids_with_msms_data[mq_data.get_evidenceid] || {}
    supported = msms[:supported] || []
    scannrs = msms[:scannrs] || []

    [supported, scannrs]
end


options = OptParser.parse(ARGV)
msms_data = read_msms(options[:msms])
proteins_with_original_names = read_gene_names_map(options[:map])
headers, seqs = Sequence.read_fasta(options[:cdna])

fh = File.open(options[:output], "w")
mq_data = ParseEvidence.new(options[:codon])
IO.foreach(options[:evidence]) do |line|
    if mq_data.is_header
        mq_data.parse_header(line)

        fh.puts "cDNA\tStart pos in protein\tEnd pos in protein\tCodon pos\tOriginal protein name\tb/y-ion supported pos\tCorresponding scan numbers\t#{line}"
        next
    end
    mq_data.parse_line(line)
    next if mq_data.is_decoy_match

    # map to MSMS data
    supported_pos, scannrs = get_corresponding_msms_data(mq_data, msms_data)
    # map to original gene names for enhanced output
    original_name = get_corresponding_original_genename(mq_data, proteins_with_original_names)

    # map peptide onto cDNA to get corresponding codon sequence
    cdna = get_protein_cdna(original_name, headers, seqs)
    codons, startpos, stoppos = map_peptide_onto_cdna(mq_data, cdna)
    codon_pos = get_codon_pos(codons, options[:codon])

    additional_data = "#{codons.join("")}\t#{startpos}\t#{stoppos}\t#{codon_pos.join(";")}\t#{original_name}\t#{supported_pos.join(";")}\t#{scannrs.join(";")}"
    fh.puts "#{additional_data}\t#{line}"
end
fh.close
