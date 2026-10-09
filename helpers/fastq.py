import os


def read_fastq(input_fastq: str) -> list[tuple[str, str, str, str]]:
    """
    Read FASTQ and returns reads

    Arguments:
    input_fastq: str (path to file)

    Returns:
    reads: list[tuple[str, str, str, str]]
    4 specifications of 1 read (read_id, seq, plus_info, qual_read)
    """
    with open(input_fastq) as fastq_file:
        reads = []
        for line in fastq_file:
            if line.startswith("@"):
                read_id = line.strip()
                seq = fastq_file.readline().strip()
                plus_info = fastq_file.readline().strip()
                qual_read = fastq_file.readline().strip()
                reads.append((read_id, seq, plus_info, qual_read))
        return reads


def gc_content(seq: str) -> float:
    """
    Calculate GC-content (percent %) of a nucleotide sequence

    Arguments:
    seq: str (nucleotide sequence)

    Returns:
    gc_procent: float (percentage of GC in nucleotide sequence)
    """
    seq = seq.upper()  # на всякий случай
    gc_count = seq.count("G") + seq.count("C")
    len_seq = len(seq)
    gc_procent = gc_count / len_seq * 100
    return gc_procent


def quality_score(qual_read: str) -> float:
    """
    Calculate mean quality of nucleotide sequence

    Arguments:
    qual_read: str (quality code (phred33) of each nucleotide)

    Returns:
    mean_qual_score: float (mean quality score of nucleotide sequence)
    """
    sum_score = 0
    for i in qual_read:
        score = ord(i) - 33
        sum_score += score
    len_read = len(qual_read)
    mean_qual_score = sum_score / len_read
    return mean_qual_score


def check_gc(seq: str, gc_bounds: tuple[float, float] | float = (0, 100)) -> bool:
    """
    Check GC content within upper and lower bounds (inclusive)

    Arguments:
    seq: str (nucleotide sequence)
    gc_bounds: tuple[float, float] | float = (0, 100)) (lower and upper bounds
    or upper bound (then lower bound = 0) or by default (0, 100))

    Returns:
    True - if GC-content within bounds
    False - if GC-content out of bounds
    """
    if type(gc_bounds) == tuple:
        low = gc_bounds[0]
        high = gc_bounds[1]
    else:
        low = 0
        high = gc_bounds
    return low <= gc_content(seq) <= high


def check_length(seq: str, length_bounds: tuple[int, int] | int = (0, 2**32)) -> bool:
    """
    Check length of sequence within upper and lower bounds (inclusive)

    Arguments:
    seq: str (nucleotide sequence)
    length_bounds: tuple[int, int] | int = (0, 2**32)) (lower and upper bounds
    or upper bound (then lower bound = 0) or by default (0, 2**32))

    Returns:
    True - if length of sequence within bounds
    False - if length of sequence out of bounds
    """
    if type(length_bounds) == tuple:
        low = length_bounds[0]
        high = length_bounds[1]
    else:
        low = 0
        high = length_bounds
    return low <= len(seq) <= high


def check_quality(qual_read: str, quality_threshold: float = 0) -> bool:
    """
    Check mean quality of sequence againts a threshold

    Arguments:
    qual_read: str (quality code (phred33) of each nucleotide)
    quality_threshold: float (threshold quality level)

    Returns:
    True - if mean quality of sequence equal or higher than threshold
    False - if mean quality of sequence less than threshold
    """
    return quality_score(qual_read) >= quality_threshold


def write_fastq(
    output_fastq: str, filtered_reads: list[tuple[str, str, str, str]]
) -> None:
    """
    Write fastq file in directory "filtered"
    and make directory "filtered" if it doesn't exist

    Arguments:
    output_fastq: str (name of output file)
    filtered_reads: list[tuple[str, str, str, str]]
    (4 specifications of 1 read (read_id, seq, plus_info, qual_read))
    """
    data_dir = "filtered"
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)
    with open(os.path.join(data_dir, output_fastq), mode="w") as fastq_file:
        for read_id, seq, plus_info, qual_read in filtered_reads:
            fastq_file.write(read_id + "\n")
            fastq_file.write(seq + "\n")
            fastq_file.write(plus_info + "\n")
            fastq_file.write(qual_read + "\n")
