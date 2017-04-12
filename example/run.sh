ATROPOS_HOME=~/projects/atropos
atropos trim \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
  -pe1 test.1.fq.gz -pe2 test.2.fq.gz -o trimmed.1.fq.gz -p trimmed.2.fq.gz \
  --report-file summary --report-formats txt json --stats both:tiles
  #-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG \
  #-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
  #-pe1 $ATROPOS_HOME/tests/data/big.1.fq -pe2 $ATROPOS_HOME/tests/data/big.2.fq \
# Run FastQC for comparison
fastqc -f fastq test.1.fq.gz test.2.fq.gz
fastqc -f fastq trimmed.1.fq.gz trimmed.2.fq.gz