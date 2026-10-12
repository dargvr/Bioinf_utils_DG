def convert_fasta(input_fasta: str) -> list[tuple[str, str]]:
    with open(input_fasta) as fasta_file:
        reads = []
        info_line = ""
        sequence_total = ""
        for line in fasta_file:
            if line.startswith(">"):
                if sequence_total != "":
                    reads.append((info_line, sequence_total))
                    sequence_total = ""
                info_line = line.strip()
            if not line.startswith(">"):
                sequence = line.strip()
                sequence_total += sequence
        reads.append((info_line, sequence_total))
        return reads


def make_output_file(output_fasta: str, convert_reads: list[tuple[str, str]]) -> None:
    with open(output_fasta, mode="w") as fasta_file:
        for info_line, sequence_total in convert_reads:
            fasta_file.write(info_line + "\n")
            fasta_file.write(sequence_total + "\n")
