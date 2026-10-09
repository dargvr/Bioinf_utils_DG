from helpers.fastq import read_fastq, check_gc, check_length, check_quality, write_fastq


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple[float, float] | float = (0, 100),
    length_bounds: tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
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
