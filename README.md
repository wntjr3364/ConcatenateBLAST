# ConcatenateBLAST

A Python tool for concatenating fragmented BLAST local alignments into continuous query-subject alignments.

## Features

- Merges fragmented BLAST alignments based on coverage and seed size
- Calculates normalized coverage between query and subject sequences

## Installation
```bash
git clone https://github.com/yourusername/ConcatenateBLAST.git
cd ConcatenateBLAST
```

## Usage
```bash
python ConcatenateBLAST.py -i <blast_result_file> -o <output_file> [options]
```

## Arguments

- `-i, --input`: BLAST output file in tabular format (outfmt 6)
- `-o, --output`: Output file name
- `-c, --coverage`: Coverage threshold (default: 0.7)
- `-s, --seed`: Seed size for alignment merging (default: 23)
- `-f, --format`: Output format ('result' or 'link', default: 'result')

### Example
```bash
python ConcatenateBLAST.py -i example/blastn.result.txt -o example/concat.result.txt -c 0.8 -s 25
```

## Output Formats

### Result Format (default)

| Index | Field Name | Description |
|-------|------------|-------------|
| 0 | qseq | Query sequence ID |
| 1 | sseq | Subject sequence ID |
| 2 | qlen | Query sequence length |
| 3 | slen | Subject sequence length |
| 4 | nident | Number of identical matches |
| 5 | qstart | Query start position |
| 6 | qend | Query end position |
| 7 | sstart | Subject start position |
| 8 | send | Subject end position |
| 9 | qcov | Query coverage |
| 10 | scov | Subject coverage |
| 11 | normCov | Normalized coverage |
| 12 | qpos | Query alignment positions |
| 13 | spos | Subject alignment positions |
| 14 | match | Number of matches per alignment section |

### Link Format

| Index | Field Name | Description |
|-------|------------|-------------|
| 0 | qseq | Query sequence ID |
| 1 | qstart | Query start position |
| 2 | qend | Query end position |
| 3 | qcov | Query coverage |
| 4 | sseq | Subject sequence ID |
| 5 | sstart | Subject start position |
| 6 | send | Subject end position |
| 7 | scov | Subject coverage |
| 8 | match | Number of matches for this alignment |
| 9 | pmatch | Percentage identity for this alignment |
