# !/usr/bin/env ruby
require "byebug"
require "optparse"

=begin
    Combine MaxQuant evidence and msms tables.


    Args:
        evidence (str): path to input file (MaxQuant evidence.txt)
        msms (str): path to input file (MaxQuant msms.txt)
        map (str): path to input file (mapping between original and shortened FASTA headers)
        output (str): path to output TSV (enriched evidence file)

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

        opt_parser = OptionParser.new do |opts|
            opts.banner = "Combine MaxQuant tables evidence and msms."
            opts.separator ""
            opts.separator "Copyright (c) 2020, by Göttingen University"
            opts.separator "Author: Stefanie Mühlhausen"
            opts.separator "This program comes with ABSOLUTELY NO WARRANTY"

            opts.separator ""
            opts.separator "Usage: ruby #{File.basename($PROGRAM_NAME)} -e evidence -m msms -o output --map dict -c codon"

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
            opts.on("-c", "--codon CODON",
            "Codon translated into all amino acids." ) do |codon|
                options[:codon] = codon
                unless Sequence.codons.include?(codon)
                    abort "#{codon} not a valid codon."
                end
            end
            opts.on("-o", "--output FILE",
                "Path to output file, in FASTA format.") do |path|
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

def get_corresponding_msms_data(mq_data, evids_with_msms_data)
    evid = mq_data.get_evidenceid
    supported = evids_with_msms_data[evid][:supported] || []
    scannrs = evids_with_msms_data[evid][:scannrs] || []

    [supported, scannrs]
end

options = OptParser.parse(ARGV)
msms_data = read_msms(options[:msms])

fh = File.open(options[:output], "w")
mq_data = ParseEvidence.new(options[:codon])
IO.foreach(options[:evidence]) do |line|
    if mq_data.is_header
        mq_data.parse_header(line)
        fh.puts "b/y-ion supported pos\tCorresponding scan numbers\t#{line}"
        next
    end
    mq_data.parse_line(line)
    next if mq_data.is_decoy_match

    supported_pos, scannrs = get_corresponding_msms_data(mq_data, msms_data)
    additional_data = "#{supported_pos.join(";")}\t#{scannrs.join(";")}"
    fh.puts "#{additional_data}\t#{line}"
end
fh.close

