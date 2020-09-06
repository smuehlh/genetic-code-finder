class ParseEvidence

    def initialize(codon = "CTG")
        @parsed_line = []
        @is_headerline = true
        @codon = codon
    end

    def is_header
        @is_headerline
    end

    def parse_header(line)
        split_line(line)
        @is_headerline = false

        @ind_peptide = @parsed_line.index("Sequence")
        @ind_protein = @parsed_line.index("Proteins")
        @ind_id = @parsed_line.index("id")
        @ind_score = @parsed_line.index("Score")
        @ind_masserr = @parsed_line.index("Mass error [ppm]") || @parsed_line.index("Mass Error [ppm]")
        @ind_reverse = @parsed_line.index("Reverse")
        @ind_contaminant = @parsed_line.index("Potential contaminant")
    end

    def parse_line(line)
        split_line(line)
    end

    def get_approximate_protein_region
        # the protein-region the peptide falls into. no exact peptide start-/stop-positions
        _, protein_data = split_matched_proteins
        [protein_data[:region_start], protein_data[:region_stop]]
    end

    def get_codon_transl_applied_to_protein
        _, protein_data = split_matched_proteins
        [protein_data[:codon], protein_data[:transl]]
    end

    def is_decoy_match
        @parsed_line[@ind_reverse] == "+" || @parsed_line[@ind_contaminant] == "+"
    end

    def get_peptide_with_unknown_aminoacids_denoted_as_X
        # maxquant replaces "X" in geneprediction (= db) by "U"
        get_peptide.gsub("U","X")
    end

    def get_evidenceid
        @parsed_line[@ind_id]
    end

    def get_peptide
        @parsed_line[@ind_peptide]
    end

    def get_protein
        split_matched_proteins.first
    end

    def get_score
        @parsed_line[@ind_score].to_f
    end

    def is_masserr_unspecified?
        @parsed_line[@ind_masserr] == "NaN"
    end

    def get_masserr
        @parsed_line[@ind_masserr].to_f
    end

    def get_codon_pos(codon)
        get_codons.each_with_index.collect do |this, ind|
            ind if this == codon
        end.compact
    end

    private

    def split_line(line)
        line = line.chomp
        @parsed_line = line.split("\t")
    end

    def split_matched_proteins
        # NOTE this returns unique proteins if multiple variants of the same protein were matched
        raw_names = @parsed_line[@ind_protein].split(";")
        parsed =
            raw_names.map do |str|
                gene_region, codon, transl = str.split("-")
                gene, start, stop = gene_region.split("_")
                data = {
                    region_start: start.to_i,
                    region_stop: stop.to_i
                }
                data[:codon] = codon if codon
                data[:transl] = transl if transl

                [gene, data]
            end.to_h
        parsed.first
    end
end