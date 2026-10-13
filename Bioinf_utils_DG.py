import os

from helpers.bio_files import (
    convert_fasta,
    make_output_file,
    make_sorted_file,
    separate_blast_table,
)
from helpers.fastq import check_gc, check_length, check_quality, read_fastq, write_fastq


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
    input_fastq: str: path to file
    output_fastq: str: name of output file
    gc_bounds: tuple[float, float] | float = (0, 100): lower and upper bounds
    or upper bound (then lower bound = 0) or by default (0, 100)
    length_bounds: tuple[int, int] | int = (0, 2**32): lower and upper bounds
    or upper bound (then lower bound = 0) or by default (0, 2**32)
    quality_threshold: float = 0: threshold quality level, threshold by default = 0
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


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str | None = None
) -> None:
    """
    Read input FASTA file and create FASTA output file
    with each sequence on one line
    If name of output file is not given,
    create name of output file ("oneline_" + input file name)

    Arguments:
    input_fasta: str: name of file or path to input file
    output_fasta: str | None: name of output file (optional)
    """
    if output_fasta is None:
        path_file = os.path.split(input_fasta)
        file_name = path_file[1]
        output_fasta = "oneline_" + file_name
    reads = convert_fasta(input_fasta)
    make_output_file(output_fasta, reads)


def parse_blast_output(input_file: str, output_file: str | None = None) -> None:
    """
    Read a BLAST output file and select the first protein
    in the "Description" column for each QUERY
    and make file with protein names, one per line, in alphabetical order
    If name of output file is not given,
    create name of output file ("sorted_proteins_" + input file name)

    Arguments:
    input_file: str: name of file or path to input file
    output_file: str | None: name of output file (optional)
    """
    if output_file is None:
        path_file = os.path.split(input_file)
        file_name = path_file[1]
        output_file = "sorted_proteins_" + file_name
    first_proteins_blast = separate_blast_table(input_file)
    make_sorted_file(output_file, first_proteins_blast)
