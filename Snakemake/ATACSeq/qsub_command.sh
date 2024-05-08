snakemake -s ATACSeq.snake --sdm conda --executor cluster-generic --cluster-generic-submit-cmd "qsub -P bf528 -pe omp {threads}" --jobs 5
