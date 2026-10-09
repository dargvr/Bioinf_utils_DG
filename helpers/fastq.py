def read_fastq(input_fastq: str) -> list[tuple[str, str, str]]:
    with open(input_fastq) as fastq_file:
        reads = []
        for line in fastq_file:
            if line.startswith("@"):
                read_id = line.strip()
                seq = fastq_file.readline().strip()
                plus_info = fastq_file.readline().strip()
                qual_read = fastq_file.readline().strip()
                reads.append((read_id, seq, qual_read))
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
