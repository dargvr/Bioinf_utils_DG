import os


def read_fastq(input_fastq: str) -> list[tuple[str, str, str, str]]:
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
    seq = seq.upper()  # на всякий случай
    gc_count = seq.count("G") + seq.count("C")
    len_seq = len(seq)
    gc_procent = gc_count / len_seq * 100
    return gc_procent


def quality_score(qual_read: str) -> float:
    sum_score = 0
    for i in qual_read:
        score = ord(i) - 33
        sum_score += score
    len_read = len(qual_read)
    mean_qual_score = sum_score / len_read
    return mean_qual_score


def check_gc(seq: str, gc_bounds: tuple[float, float] | float = (0, 100)) -> bool:
    if type(gc_bounds) == tuple:
        low = gc_bounds[0]
        high = gc_bounds[1]
    else:
        low = 0
        high = gc_bounds
    return low <= gc_content(seq) <= high


def check_length(seq: str, length_bounds: tuple[int, int] | int = (0, 2**32)) -> bool:
    if type(length_bounds) == tuple:
        low = length_bounds[0]
        high = length_bounds[1]
    else:
        low = 0
        high = length_bounds
    return low <= len(seq) <= high


def check_quality(qual_read: str, quality_threshold: float = 0) -> bool:
    return quality_score(qual_read) >= quality_threshold


def write_fastq(
    output_fastq: str, filtered_reads: list[tuple[str, str, str, str]]
) -> None:
    data_dir = "filtered"
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)
    with open(os.path.join(data_dir, output_fastq), mode="w") as fastq_file:
        for read_id, seq, plus_info, qual_read in filtered_reads:
            fastq_file.write(read_id + "\n")
            fastq_file.write(seq + "\n")
            fastq_file.write(plus_info + "\n")
            fastq_file.write(qual_read + "\n")
