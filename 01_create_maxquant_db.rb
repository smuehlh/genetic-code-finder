# !/usr/bin/env ruby
require "byebug"
require "optparse"

=begin
    Create a FASTA file to be used as MaxQuant database.

    Skip all sequences that cannot be splitted into codons (whose length is
    not a multiple of 3)


    Args:
        input (str): path input FASTA (cDNA sequences)
        output (str): path to output FASTA (protein sequences)
        map (str): path to output CSV (dictionary of simplified gene names)
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
        options[:map] = nil
        options[:codon] = nil

        # optional parameters
        options[:omit_ile] = true
        options[:cleavage] = ["K", "R"]

        opt_parser = OptionParser.new do |opts|
            opts.banner = "Generate FASTA file to be used as MaxQuant database."
            opts.separator ""
            opts.separator "Copyright (c) 2020, by Göttingen University"
            opts.separator "Author: Stefanie Mühlhausen"
            opts.separator "This program comes with ABSOLUTELY NO WARRANTY"

            opts.separator ""
            opts.separator "Usage: ruby #{File.basename($PROGRAM_NAME)} -i input -o output -m dict"

            opts.on("-i", "--input FILE",
                "Path to input file, in FASTA format.") do |path|
                FileHelper.file_exist_or_die(path)
                options[:input] = path
            end
            opts.on("-o", "--output FILE",
                "Path to output file, in FASTA format.") do |path|
                options[:output] = path
            end
            opts.on("-m", "--map FILE",
                "Path to auxially output file mapping original ",
                "to simplified FASTA headers, in CSV format.") do |path|
                options[:map] = path
            end
            opts.on("-c", "--codon CODON",
            "Codon to translate into all amino acids." ) do |codon|
                options[:codon] = codon
                unless Sequence.genetic_code.key?(codon)
                    abort "#{codon} not a valid codon."
                end
            end
            opts.on("--ile",
                "Optional, specify to translate into both Ile and Leu",
                "NOTE - MaxQuant will pick up only Leu translation") do |opt|
                options[:omit_ile] = false
            end
            opts.on("--cleavage K,R", Array,
                "Optional, specify to cleave after amino acids other than 'K' and 'R'") do |aminoacids|
                invalid = aminoacids - Sequence.genetic_code.values
                if invalid.any?
                    abort "#{invalid} invalid amino acids"
                end
                options[:cleavage] = aminoacids
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
        abort "Missing mandatory argument: --output" unless options[:output]
        abort "Missing mandatory argument: --map" unless options[:map]
        abort "Missing mandatory argument: --codon" unless options[:codon]

        if options[:omit_ile]
            puts "INFO - do not translate #{options[:codon]} into Ile (recommended)"
        else
            puts "INFO - translate #{options[:codon]} into both Ile and Leu (not recommended)"
        end
        puts "INFO - prepare DB for cleavage after #{options[:cleavage].join(",")}"
        return options
    end
end

def simplify_header(num)
    "g#{num}"
end

options = OptParser.parse(ARGV)

headers, seqs = Sequence.read_fasta(options[:input])
fh_map = File.open(options[:map], "w")
seqs.each_with_index do |seq, ind|
    orig_header = headers[ind]
    patched_header = simplify_header(ind+1)
    fh_map.puts "#{orig_header};#{patched_header}"

    codons = Sequence.split_cdna_into_codons(seq)
    if codons.last.size != 3
        puts "Skip #{orig_header}, cannot split cDNA into codons"
        next
    end

end
fh_map.close
