# !/usr/bin/env ruby
require "optparse"

=begin
    Create a FASTA file to be used as MaxQuant database.

    Translate a given codon into all amino acids (with the exception of Ile, for which no separate translation will be created by default as Ile and Leu are indistinguishable using mass spectrometry).

    Size of the resulting database is reduced by separating sequences at cleavage sizes into parts with and without the specified codon and adding the parts without only once to the DB. Thus, only in sequence parts that could be detected as cleaved peptides the specified codon will be translated into all amino acids (save Ile). Cutting sequences a cleavage sites will allow for two missed cleaves. Resulting subsequences must be of minimum length 7 to reflect MaxQuant's min length setting.

    Skip all sequences that cannot be split into codons (whose length is
    not a multiple of 3)

    Args:
        input (str): path input FASTA (cDNA sequences)
        output (str): path to output FASTA (protein sequences)
        map (str): path to output CSV (dictionary of simplified gene names)
        codon (str): optional argument; if present codon other than CTG will be translated into each amino acid
        ile: optional argument; if present codon will be translated into both
            Ile and Leu
        cleavage (array): optional argument; amino acids to cleave at


    Returns:
        a FASTA file to be used as MaxQuant database
        a CSV file mapping shortened onto original FASTA headers
=end

# require .rb files in library
Dir[File.join(__dir__, "lib", "**", "*.rb")].each do |file|
    require file
end

class OptParser
    def self.parse(args)

        options = Hash.new

        # mandatory parameters
        options[:input] = nil
        options[:output] = nil
        options[:map] = nil

        # optional parameters
        options[:codon] = "CTG"
        options[:omit_ile] = true
        options[:cleavage] = ["K", "R"]

        opt_parser = OptionParser.new do |opts|
            opts.banner = "Generate FASTA file to be used as MaxQuant database."
            opts.separator ""
            opts.separator "Copyright (c) 2020-2021, by Göttingen University"
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
                "Path to auxiliary output file mapping original ",
                "to shortened FASTA headers, in CSV format.") do |path|
                options[:map] = path
            end
            opts.on("--codon CODON",
            "Optional, specify to translate codon other than CTG into all amino acids." ) do |codon|
                options[:codon] = codon
                unless Sequence.codons.include?(codon)
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
                invalid = aminoacids - Sequence.amino_acids
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

        if options[:omit_ile]
            puts "INFO - do not translate #{options[:codon]} into Ile (recommended)"
        else
            puts "INFO - translate #{options[:codon]} into both Ile and Leu (not recommended)"
        end
        puts "INFO - prepare DB for cleavage after #{options[:cleavage].join(",")}"
        puts "INFO - prepare DB for codon #{options[:codon]}"
        return options
    end
end

def simplify_header(num)
    "g#{num}"
end

def append_header(header, start, stop, *codon_transl)
    if codon_transl.any?
        "#{header}_#{start}_#{stop}-#{codon_transl[0]}-#{codon_transl[1]}"
    else
        "#{header}_#{start}_#{stop}"
    end
end

def divide_codons_into_parts_wo_and_parts_with_codon_of_interest(codons, cleavage_sites, codon_of_interest)
    is_collecting_part_with_codon = false
    cleavage_indices = []
    parts_wo_codon = [] # array of [start, stop]
    parts_with_codon = [] # array of [start, stop]
    start = 0

    codons.each_with_index do |codon, codon_ind|
        if cleavage_sites.include?(Sequence.translate(codon))
            cleavage_indices.push codon_ind
            possible_stop = cleavage_indices.last
            if is_collecting_part_with_codon && cleavage_indices.size >= 3 && codons[start..possible_stop].size >= 6
                parts_with_codon.push [start, possible_stop]

                start = possible_stop + 1
                is_collecting_part_with_codon = false
                cleavage_indices = []
            end
        elsif codon == codon_of_interest
            possible_stop = cleavage_indices[-3]
            if cleavage_indices.size >= 3 && codons[start..possible_stop].size >= 6
                parts_wo_codon.push [start, possible_stop]

                start = possible_stop + 1
            end
            is_collecting_part_with_codon = true
            cleavage_indices = []
        end
    end

    stop = codons.size - 1
    # prepare to add remaining sequence to previous part if it would be too short otherwise
    if codons[start..stop].size < 6
        last_start,_ =
            if parts_wo_codon.empty?
                is_collecting_part_with_codon = true
                parts_with_codon.pop
            elsif parts_with_codon.empty?
                parts_wo_codon.pop
            else
                last_start_wo,_ = parts_wo_codon.last
                last_start_w,_ = parts_with_codon.last
                if last_start_w > last_start_wo
                    is_collecting_part_with_codon = true
                    parts_with_codon.pop
                else
                    parts_wo_codon.pop
                end
            end
        start = last_start
    end

    if is_collecting_part_with_codon
        parts_with_codon.push [start, stop]
    else
        parts_wo_codon.push [start, stop]
    end

    [parts_wo_codon, parts_with_codon]
end

options = OptParser.parse(ARGV)

headers, seqs = Sequence.read_fasta(options[:input])
fh_map = File.open(options[:map], "w")
fh_seq = File.open(options[:output], "w")

seqs.each_with_index do |seq, ind|
    orig_header = headers[ind]
    patched_header = simplify_header(ind+1)
    fh_map.puts "#{orig_header};#{patched_header}"

    codons = Sequence.split_cdna_into_codons(seq)
    if codons.last.size != 3
        puts "Skip #{orig_header}, cannot split cDNA into codons"
        next
    end

    if codons.include?(options[:codon])
        # 'dynamic part': translate codon into each amino acid
        parts_wo_codon, parts_with_codon = divide_codons_into_parts_wo_and_parts_with_codon_of_interest(
                codons, options[:cleavage], options[:codon]
            )

        # 'static sub-part'
        parts_wo_codon.each do |start, stop|
            codons_part = codons[start..stop]
            header = append_header(patched_header, start, stop)
            Sequence.write_fasta(fh_seq, header, Sequence.translate(codons_part.join("")))
        end

        # 'dynamic sub-part': translate into all amino acids
        parts_with_codon.each do |start, stop|
            codons_part = codons[start..stop]
            Sequence.amino_acids.each do |aa|
                next if aa == "I" && options[:omit_ile]
                transl = Sequence.translate_codons_with_one_codon_set(codons_part, options[:codon], aa)
                header = append_header(patched_header, start, stop, options[:codon], aa)
                Sequence.write_fasta(fh_seq, header, transl)
            end
        end
    else
        # 'static part': translate sequence only once
        header = append_header(patched_header, 0, codons.size-1)
        Sequence.write_fasta(fh_seq, header, Sequence.translate(codons.join("")))
    end

end
fh_map.close
fh_seq.close
