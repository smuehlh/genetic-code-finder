class ParseEnrichedEvidence < ParseEvidence

    alias_method :parent_get_protein, :get_protein

    def initialize(codon = "CTG")
        super(codon)
    end

    def parse_header(line)

        super(line)

        @ind_cdna = @parsed_line.index("cDNA")
        @ind_codon = @parsed_line.index("Codon pos")
        @ind_startpos = @parsed_line.index("Start pos in protein")
        @ind_stoppos = @parsed_line.index("End pos in protein")
        @ind_orig_name = @parsed_line.index("Original protein name")
        @ind_supportedpos = @parsed_line.index("b/y-ion supported pos")
        @ind_scan_nrs = @parsed_line.index("Corresponding scan numbers")
    end

    def get_codon_pos
        # position(s) in peptide
        @parsed_line[@ind_codon].split(";").map(&:to_i)
    end

    def get_codons
        Sequence.split_cdna_into_codons(@parsed_line[@ind_cdna])
    end

    def get_protein
        @parsed_line[@ind_orig_name]
    end

    def get_protein_name_used_in_db
        parent_get_protein
    end

    def get_peptide_start
        @parsed_line[@ind_startpos].to_i
    end

    def get_peptide_stop
        @parsed_line[@ind_stoppos].to_i
    end

    def convert_peptide_to_protein_pos(peptide_pos)
        get_peptide_start + peptide_pos
    end

    def get_supported_pos
        @parsed_line[@ind_supportedpos].split(";").map(&:to_i)
    end

    def get_scannumbers
        @parsed_line[@ind_scan_nrs].split(";")
    end
end