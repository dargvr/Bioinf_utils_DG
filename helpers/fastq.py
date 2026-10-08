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
