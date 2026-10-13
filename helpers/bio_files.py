def convert_fasta(input_fasta: str) -> list[tuple[str, str]]:
    """
    Read input FASTA file and join multiple sequence lines into one line

    Arguments:
    input_fasta: str: name of file or path to file

    Returns:
    list[tuple[str, str]]: list of information about sequence and sequence in one line
    """
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
    """
    Make output file with reads including read's information and sequence

    Arguments:
    output_fasta: str: name of output file
    convert_reads: list[tuple[str, str]]: list of information about sequence and sequence in one line
    """
    with open(output_fasta, mode="w") as fasta_file:
        for info_line, sequence_total in convert_reads:
            fasta_file.write(info_line + "\n")
            fasta_file.write(sequence_total + "\n")


def separate_blast_table(input_file: str) -> list[str]:
    """
    Read a BLAST output file and select the first protein
    in the "Description" column for each QUERY

    Arguments:
    input_file: str: name of file or path to file

    Returns:
    list[str]: list of the first proteins from the "Description"
    """
    first_proteins = []
    with open(input_file) as blast_file:
        for line in blast_file:
            if line.startswith("Description"):
                first_protein = blast_file.readline().strip()
                first_protein = first_protein.split("  ")[0]
                if "..." in first_protein:
                    first_protein = first_protein.split("...")[0]
                first_proteins.append(first_protein)
        return first_proteins


def make_sorted_file(output_file: str, first_proteins: list[str]) -> None:
    """
    Make file with protein names, one per line, in alphabetical order

    Arguments:
    output_file: str: name of output file
    first_proteins: list[str]: list of proteins
    """
    first_proteins.sort()
    with open(output_file, mode="w") as text_file:
        for protein in first_proteins:
            text_file.write(protein + "\n")
