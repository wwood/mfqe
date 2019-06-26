```
Usage: zcat my.fastq.gz |mfqe --fastq-read-name-lists <LIST1> .. --output-fastq-files <OUTPUT1> ..

Purpose is to extract one or more sets of reads from a FASTQ (or FASTA) file by specifying their read names.

Read name files are uncompressed text files with read names(without comments).
Output is gzip-compressed, input may or may not be.

Other FASTQ options:

--input-fastq <PATH>: Use this file as input FASTQ [default: Use STDIN]

An analogous set of options is implemented for FASTA:

--fasta-read-name-lists <LIST1> ..
--output-fasta-files <OUTPUT1> ..
--input-fasta <PATH>
```
