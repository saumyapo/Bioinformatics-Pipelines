names = ['ADrep1_subsample_R1', 'ADrep1_subsample_R2', 'ADrep2_subsample_R1', 'ADrep2_subsample_R2',
         'P0rep1_subsample_R1', 'P0rep1_subsample_R2', 'P0rep2_subsample_R1', 'P0rep2_subsample_R2',
         'P4rep1_subsample_R1', 'P4rep1_subsample_R2', 'P4rep2_subsample_R1', 'P4rep2_subsample_R2',
         'P7rep1_subsample_R1', 'P7rep1_subsample_R2', 'P7rep2_subsample_R1', 'P7rep2_subsample_R2']

output_dir = "results"

process fastqc {
    input:
    val names from names

    output:
    file("results/fastqc/${names}_fastqc.html") into fastqc_html,
    file("results/fastqc/${names}_fastqc.zip") into fastqc_zip

    script:
    """
    fastqc samples/${names}.fastq.gz -o results/fastqc
    """
}

process multiqc {
    input:
    file(fastqc_zip) from fastqc_zip.collect()

    output:
    file("results/multiqc_report.html")

    script:
    """
    multiqc results/fastqc/ -o results/
    """
}

process get_m39 {
    output:
    file("reference/GRCm39.primary_assembly.genome.fa.gz")

    script:
    """
    wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/GRCm39.primary_assembly.genome.fa.gz" -P reference/
    """
}

process get_m39_gtf {
    output:
    file("reference/gencode.vM34.primary_assembly.annotation.gtf")

    script:
    """
    wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.primary_assembly.annotation.gtf.gz" -P reference/
    gzip -d reference/*gtf*
    """
}

process star {
    input:
    val names from names

    output:
    file("results/star/${names}.Aligned.out.bam")

    script:
    """
    STAR --genomeDir reference/m39_subset_star/ --readFilesIn samples/${names}_R1.fastq.gz samples/${names}_R2.fastq.gz --readFilesCommand zcat --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix results/star/${names}.
    """
}

process samtools_flagstat {
    input:
    val names from names

    output:
    file("results/samtools/${names}_flagstat.txt")

    script:
    """
    samtools flagstat results/star/${names}.Aligned.out.bam > results/samtools/${names}_flagstat.txt
    """
}

process verse {
    input:
    val names from names

    output:
    file("results/verse/${names}.exon.txt")

    script:
    """
    verse -S -a reference/gencode.vM34.primary_assembly.annotation.gtf -o results/verse/${names} results/star/${names}.Aligned.out.bam
    """
}

process concat_verse {
    input:
    val names from names

    output:
    file("results/verse/verse_concat.csv")

    script:
    """
    python concat_df.py -i $(find results/verse/ -name "${names}.exon.txt") -o results/verse/verse_concat.csv
    """
}

process filter_cts {
    input:
    file("results/verse/verse_concat.csv")

    output:
    file("results/verse/verse_concat_filtered.csv")

    script:
    """
    python filter_cts_mat.py -i results/verse/verse_concat.csv -o results/verse/verse_concat_filtered.csv
    """
}

process txn_mapping {
    output:
    file("results/id2gene.txt")

    script:
    """
    python parse_gtf.py -i reference/gencode.vM34.primary_assembly.annotation.gtf -o results/id2gene.txt
    """
}

workflow {
    input:
    set val names from names.collect { it.substring(0, -9) } , file("samples/${names}.fastq.gz") from names.collect { it.substring(0, -9) }

    // Define the workflow steps
    fastqc()
    multiqc()
    get_m39()
    get_m39_gtf()
    star()
    samtools_flagstat()
    verse()
    concat_verse()
    filter_cts()
    txn_mapping()
}
