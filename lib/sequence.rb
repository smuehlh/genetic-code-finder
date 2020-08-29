module Sequence
    extend self

    def genetic_code
        {
            "TTT" => "F", "TTC" => "F",
            "TTA" => "L", "TTG" => "L", "CTT" => "L", "CTC" => "L",
                "CTA" => "L", "CTG" => "L",
            "ATT" => "I", "ATC" => "I", "ATA" => "I",
            "ATG" => "M",
            "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
            "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
                "AGT" => "S", "AGC" => "S",
            "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
            "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
            "TAT" => "Y", "TAC" => "Y",
            "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
            "CAT" => "H", "CAC" => "H",
            "CAA" => "Q", "CAG" => "Q",
            "AAT" => "N", "AAC" => "N",
            "AAA" => "K", "AAG" => "K",
            "GAT" => "D", "GAC" => "D",
            "GAA" => "E", "GAG" => "E",
            "TGT" => "C", "TGC" => "C",
            "TGG" => "W",
            "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
                "AGA" => "R", "AGG" => "R",
            "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G",
            "TAA" => "*", "TAG" => "*", "TGA" => "*"
        }
    end

    def codons
        genetic_code.keys
    end

    def amino_acids
        genetic_code.values.uniq - ["*"]
    end

    def is_stopcodon(codon)
        genetic_code.fetch(codon, "") == "*"
    end

    def translate(cdna)
        codons = split_cdna_into_codons(cdna)
        delete_trailing_stopcodon_if_present(codons)
        codons.map do |codon|
            genetic_code.fetch(codon).upcase
        end.join("")
    end

    def split_cdna_into_codons(cdna)
        cdna.scan(/.{1,3}/)
    end

    def delete_trailing_stopcodon_if_present(codons)
        if is_stopcodon(codons.last)
            codons.pop
        end
    end

    def read_fasta(path)
        headers_with_seqs = {}
        current_header = ""
        IO.foreach(path) do |line|
            line.chomp!
            if line.start_with?(">")
                # fasta header
                current_header = line[1..-1]
                if headers_with_seqs[current_header]
                    abort "duplicate FASTA header #{current_header}"
                end
                headers_with_seqs[current_header] = ""
            else
                # fasta sequence
                headers_with_seqs[current_header] += line
            end
        end
        [headers_with_seqs.keys, headers_with_seqs.values]
    end

    def write_fasta(fh, header, seq)
        header = ">#{header}" unless header.start_with?(">")
        fh.puts header
        fh.puts seq.scan(/.{1,80}/).join("\n")
    end

end