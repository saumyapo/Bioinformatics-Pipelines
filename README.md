# Bioinformatics-Pipelines
Bioinformatics pipelines developed either using Nextflow or Snakemake.

## Snakemake Pipelines
- The `Snakemake` directory has different Snakemake pipelines developed for different data types. 
- Each directory contains the Snakemake file for that respective workflow and also a Discussion.Rmd which contains information about what references and methods were used to recreate the results via the pipeline and what the conclusions/learning outcomes were.
- The directory also includes other scripts written outside of Snakemake (using either Python or R) which were used for further analysis of the data.
- Where relevant there is an "images" directory present which holds the images generated during downstream analysis of the data.
- The conda env used to run each workflow is present and additionally there is an `envs` repository with the yml used for each tool during analysis.

### Links to Snakemake files
1. [Snakemake RNASeq:](Snakemake/RNASeq)
   *  [Methods and Results](Snakemake/RNASeq/RNASeq_Discussion.Rmd)
2. [Snakemake ChIPSeq:](Snakemake/ChIPSeq)
   *  [Methods and Results](Snakemake/ChIPSeq/ChIPSeq_Discussion.Rmd)
   *  [Images obtained from downstream analysis](Snakemake/ChIPSeq/images)
3. [Snakemake ATACSeq:](Snakemake/ATACSeq)
   *  [Methods and Results](Snakemake/ATACSeq/ATACSeq_Discussion.md)
   *  [Images obtained from downstream analysis](Snakemake/ATACSeq/images)


## Nextflow Pipelines
The Nextflow pipelines were developed using the same workflow to additionally practice Groovy and are stored in `Nextflow`. Since the same methods and data were used, the results are identical. The links for the pipelines can be found below:
1. [NF RNASeq](Nextflow/RNASeq.nf)
2. [NF ChIPSeq](Nextflow/ChIPSeq.nf)

   

