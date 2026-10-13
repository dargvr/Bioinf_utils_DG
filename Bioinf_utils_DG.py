import os
from helpers.fastq import read_fastq, check_gc, check_length, check_quality, write_fastq


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple[float, float] | float = (0, 100),
    length_bounds: tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Filter FASTQ input file by length, quality and GC content to output file in directory "filtered"

    Arguments:
    input_fastq: str (path to file)
    output_fastq: str (name of output file)
    gc_bounds: tuple[float, float] | float = (0, 100)) (lower and upper bounds
    or upper bound (then lower bound = 0) or by default (0, 100))
    length_bounds: tuple[int, int] | int = (0, 2**32)) (lower and upper bounds
    or upper bound (then lower bound = 0) or by default (0, 2**32))
    quality_threshold: float (threshold quality level, threshold by default = 0)
    """
    reads = read_fastq(input_fastq)
    filtered_reads = []
    for read_id, seq, plus_info, qual_read in reads:
        if (
            check_gc(seq, gc_bounds)
            and check_length(seq, length_bounds)
            and check_quality(qual_read, quality_threshold)
        ) == True:
            filtered_reads.append((read_id, seq, plus_info, qual_read))
    write_fastq(output_fastq, filtered_reads)


from helpers.bio_files import convert_fasta, make_output_file


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str | None = None
) -> None:
    if output_fasta is None:
        path_file = os.path.split(input_fasta)
        file_name = path_file[1]
        output_fasta = "oneline_" + file_name
    reads = convert_fasta(input_fasta)
    make_output_file(output_fasta, reads)


from helpers.bio_files import separate_blast_table, make_sorted_file


def parse_blast_output(input_file: str, output_file: str | None = None) -> None:
    if output_file is None:
        path_file = os.path.split(input_file)
        file_name = path_file[1]
        output_file = "sorted_proteins_" + file_name
    first_proteins_blast = separate_blast_table(input_file)
    make_sorted_file(output_file, first_proteins_blast)
