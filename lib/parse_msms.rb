class ParseMsms

    def initialize
        @parsed_line = []
        @is_headerline = true
    end

    def is_header
        @is_headerline
    end

    def parse_header(line)
        split_line(line)
        @is_headerline = false

        @ind_peptide = @parsed_line.index("Sequence")
        @ind_protein = @parsed_line.index("Proteins")
        @ind_evidenceid = @parsed_line.index("Evidence ID")
        @ind_matched_ions = @parsed_line.index("Matches")
        @ind_scannumber = @parsed_line.index("Scan number")
    end

    def parse_line(line)
        split_line(line)
    end

    def get_evidenceid
        @parsed_line[@ind_evidenceid]
    end

    def get_peptide
        @parsed_line[@ind_peptide]
    end

    def get_matched_ions
        @parsed_line[@ind_matched_ions].split(";")
    end

    def get_positions_with_b_y_type_support
        peptide_inds = (0..get_peptide.size-1).to_a
        peptide_inds.collect do |ind|
            ind if is_pos_sufficiently_supported?(ind)
        end.compact
    end

    def get_scannumber
        @parsed_line[@ind_scannumber]
    end

    private

    def split_line(line)
        line = line.chomp
        @parsed_line = line.split("\t")
    end

    def is_pos_sufficiently_supported?(pos)
        # NOTE: this does not include that measured and theoretical mass equal
        b, b_minus1, y, y_minus1 = get_support_for(pos)

        # imaginary fragments for first/last peptide pos
        if pos == 0
            b_minus1 = nil # b0 does not exist
            y = true # Mass of y(n) = MW(peptide)
        end
        if pos == get_peptide.size - 1
            y_minus1 = nil # y0 does not exist
            b = true # Mass of b(n) = MW(peptide)
        end

        (y_minus1 && y) || (b_minus1 && b) || (b && y) || (y_minus1 && b_minus1)
    end

    def get_support_for(pos)
        peptide_pos_human_couting = pos + 1
        reverse_peptide_pos_human_counting = get_peptide.size - pos

        b = is_ion_matched?("b#{peptide_pos_human_couting}")
        b_minus1 = is_ion_matched?("b#{peptide_pos_human_couting-1}")
        y = is_ion_matched?("y#{reverse_peptide_pos_human_counting}")
        y_minus1 = is_ion_matched?("y#{reverse_peptide_pos_human_counting-1}")

        [b, b_minus1, y, y_minus1]
    end

    def is_ion_matched?(ion)
        get_matched_ions.index(ion) ||
        get_matched_ions.index("#{ion}-NH3") ||
        get_matched_ions.index("#{ion}-H2O") ||
        get_matched_ions.index("#{ion}(2+)")
    end
end