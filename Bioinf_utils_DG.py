def filter_fastq(
    input_fastq: str, 
    output_fastq:str, 
    gc_bounds: tuple[float, float] | float = (0,100), 
    length_bounds: tuple[int, int] | int = (0, 2**32), 
    quality_threshold: float = 0) -> None:

