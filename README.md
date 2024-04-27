# Bioinformatics-Pipelines
Bioinformatics pipelines developed either using Nextflow or Snakemake.

## Snakemake Pipelines
The `Snakemake` directory has different Snakemake pipelines developed for different data types. 

Each directory contains the Snakemake file for that respective workflow and also a Discussion.Rmd which contains information about what references and methods were used to recreate the results via the pipeline and what the conclusions/learning outcomes were. The directory also includes other scripts written outside of Snakemake (using either Python or R) which were used for further analysis of the data.
Where relevant there is an "images" directory present which holds the images generated during downstream analysis of the data.

### Links to relevant files can be found below:
1. [RNASeq:](RNASeq)
- [Methods and Results](RNASeq/RNASeq_Discussion.Rmd)
2. [ChIPSeq:](/ChIPSeq)
- [Methods and Results](ChIPSeq/ChIPSeq_Discussion.Rmd)
- [Images obtained from downstream analysis](ChIPSeq/images)


## Nextflow Pipelines
<strong> Note: Nextflow pipelines of the same workflow were developed to additionally practice using Nextflow. Since the same methods and data were used, the results are identical. The links for the pipelines can be found below:</strong>
1. [RNASeq](https://github.com/saumyapo/Nextflow/blob/fe14420bd0fb7950e446df1020833cf1dcf0398a/RNASeq.nf)
2. [ChIPSeq](https://github.com/saumyapo/Nextflow/blob/fe14420bd0fb7950e446df1020833cf1dcf0398a/ChIPSeq.nf)
   

