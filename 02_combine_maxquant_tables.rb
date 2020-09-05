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

        opt_parser = OptionParser.new do |opts|
            opts.banner = "Combine MaxQuant tables evidence and msms."
            opts.separator ""
            opts.separator "Copyright (c) 2020, by Göttingen University"
            opts.separator "Author: Stefanie Mühlhausen"
            opts.separator "This program comes with ABSOLUTELY NO WARRANTY"

            opts.separator ""
            opts.separator "Usage: ruby #{File.basename($PROGRAM_NAME)} -i input -o output -m dict"

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
        abort "Missing mandatory argument: --output" unless options[:output]

        return options
    end
end

options = OptParser.parse(ARGV)

fh = File.open(options[:output], "w")
IO.foreach(options[:evidence]) do |line|
    fh.puts line

end
fh.close

