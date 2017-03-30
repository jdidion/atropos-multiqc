ATROPOS_HOME=~/projects/atropos
atropos trim --threads 2 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -pe1 $ATROPOS_HOME/tests/data/big.1.fq -pe2 $ATROPOS_HOME/tests/data/big.2.fq -o test1.fq -p test2.fq --report-file summary --report-formats txt json --stats both
