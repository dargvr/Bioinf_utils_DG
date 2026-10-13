# Bioinf_utils_DG
Утилиты на python для работы с файлами FASTA, FASTQ и результатами BLAST

 ![Пример работы Bioinf_utils_DG](image_for_readme/Bioinf_utils_DG.jpg)

## Структура проекта
```python
-/
|- README.md
|- Bioinf_utils_DG.py # утилиты (импорты, filter_fastq и функции из части B: convert_multiline_fasta_to_oneline и parse_blast_output)
|- helpers/
    |- fastq.py       # вспомогательные функции для filter_fastq
    |- bio_files.py   # вспомогательные функции для части B
```

## Использование утилит 
```python
from Bioinf_utils_DG import (
    filter_fastq,
    convert_multiline_fasta_to_oneline,
    parse_blast_output,
)
```
## Поддерживаемые функции
### filter_fastq
Фильтрует риды из FASTQ-файла по GC-составу, длине и качеству (phred33)
Отфильтрованные риды сохраняются в папку filtered 
```python
filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple[float, float] | float = (0, 100),
    length_bounds: tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0,
)
```

### convert_multiline_fasta_to_oneline
Читает FASTA-файл, в котором последовательность нуклеотидов/аминокислот может быть записана на нескольких строках и сохраняет новый FASTA-файл, в котором последовательность записана одной строкой. 
Если имя выходного файла не указано в функции, функция создаст его автоматически. 
```python
convert_multiline_fasta_to_oneline(
    input_fasta: str, 
    output_fasta: str | None = None
)
```

### parse_blast_output
Читает файл с результатами BLAST и выбирает наилучшее совпадение - первую строку в столбце Description и сохраняет названия выбранных белков в выходной файл в алфавитном порядке по одному на строку.
Если имя выходного файла не указано в функции, функция создаст его автоматически. 
```python
parse_blast_output(
    input_file: str, 
    output_file: str | None = None
)
```