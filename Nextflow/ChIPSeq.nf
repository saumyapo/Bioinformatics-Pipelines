
import groovy.transform.Field

sample_csv = file('sample_sheet.csv').text.readLines().drop(1).collect { it.split(",") as [String] }.collectEntries { [(it[0]): it[1]] }

CONDITIONS = sample_csv.values().collect { it[1] }.toSet()
REPS = sample_csv.values().collect { it[2] }.toSet()


@Field def wildcards = [
    condition: CONDITIONS.join('|'),
    rep: REPS.join('|')
]

output_dir = "results"

process wget_files {
    input:
    val condition, rep from sample_csv

    output:
    file("samples/${condition}_${rep}.fastq.gz")

    script:
    """
    wget -O samples/${condition}_${rep}.fastq.gz ${sample_csv[condition, rep]}
    """
}

process bowtie2_build_gencode {
    output:
    file("reference/bowtie_index_GRCh38/human_reference.4.bt2")

    script:
    """
    bowtie2-build --threads 16 reference/GRCh38.primary_assembly.genome.fa.gz reference/index/human_reference
    """
}

process fastqc {
    input:
    file(fastq) from glob("samples/*_{wildcards.rep}.fastq.gz")

    output:
    file("results/fastqc/{condition}_{rep}_fastqc.html") into fastqc_out

    script:
    """
    fastqc ${fastq} -o results/fastqc --threads 4
    """
}

process trimmomatic {
    input:
    file(fastq) from glob("samples/*_{wildcards.rep}.fastq.gz")

    output:
    file("results/trimmomatic/{condition}_{rep}_trim.fastq.gz")

    script:
    """
    trimmomatic SE ${fastq} results/trimmomatic/{condition}_{rep}_trim.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
    """
}

process bowtie2_align {
    input:
    file(fastq) from glob("results/trimmomatic/*_{wildcards.rep}_trim.fastq.gz")

    output:
    file("results/bam/{condition}_{rep}.bam")

    script:
    """
    bowtie2 -x reference/index/human_reference -U ${fastq} | samtools view -b > results/bam/{condition}_{rep}.bam
    """
}

process samtools_sort {
    input:
    file(bam) from glob("results/bam/*_{wildcards.rep}.bam")

    output:
    file("results/bam/{condition}_{rep}_sorted.bam")

    script:
    """
    samtools sort -@ 16 -o results/bam/{condition}_{rep}_sorted.bam ${bam}
    """
}

process samtools_idx {
    input:
    file(bam) from glob("results/bam/*_{wildcards.rep}_sorted.bam")

    output:
    file("results/bam/{condition}_{rep}_sorted.bam.bai")

    script:
    """
    samtools index ${bam}
    """
}

process samtools_flagstats {
    input:
    file(bam) from glob("results/bam/*_{wildcards.rep}_sorted.bam")

    output:
    file("results/flagstats/{condition}_{rep}_flagstats.txt")

    script:
    """
    samtools flagstat ${bam} > results/flagstats/{condition}_{rep}_flagstats.txt
    """
}

process multiqc {
    input:
    file(flagstats) from glob("results/flagstats/*_{wildcards.rep}_flagstats.txt"),
    file(fastqc) from glob("results/fastqc/*_{wildcards.rep}_fastqc.html")

    output:
    file("results/multiqc_report.html")

    script:
    """
    multiqc -o results/ ${flagstats} ${fastqc}
    """
}

process bamCoverage {
    input:
    file(bam) from glob("results/bam/*_{wildcards.rep}_sorted.bam")

    output:
    file("results/bigwig/{condition}_{rep}.bw")

    script:
    """
    bamCoverage -b ${bam} -o results/bigwig/{condition}_{rep}.bw
    """
}

process multiBwSummary {
    input:
    file(bigwigs) from glob("results/bigwig/*_{wildcards.rep}.bw")

    output:
    file("results/multiBwSummary_matrix.txt")

    script:
    """
    multiBigwigSummary bins -b ${bigwigs} -out results/multiBwSummary_matrix.txt
    """
}

process plotCorrelation {
    input:
    file(matrix) from file("results/multiBwSummary_matrix.txt")

    output:
    file(plot) into correlation_plot,
    file(matrix) into correlation_matrix

    script:
    """
    plotCorrelation --corData ${matrix} --corMethod pearson --skipZeros --whatToPlot heatmap --plotFile ${plot} --outFileCorMatrix ${matrix} --plotNumbers
    """
}

process make_tag_dir {
    input:
    file(bam) from glob("results/bam/*_{wildcards.rep}_sorted.bam")

    output:
    directory("results/tagdir/{condition}_{rep}_tagDir")

    script:
    """
    makeTagDirectory results/tagdir/{condition}_{rep}_tagDir ${bam}
    """
}

process findPeaks {
    input:
    file(control) from glob("results/tagdir/INP_{wildcards.rep}_tagDir"),
    file(exp) from glob("results/tagdir/RUNX1_{wildcards.rep}_tagDir")

    output:
    file("results/peaks/{rep}.peaks")

    script:
    """
    findPeaks ${exp} -style factor -o results/peaks/{rep}.peaks -i ${control}
    """
}

process convertPeakFiles {
    input:
    file(peaks) from glob("results/peaks/*_{wildcards.rep}.peaks")

    output:
    file("results/bed/{rep}.bed")

    script:
    """
    pos2bed.pl ${peaks} > results/bed/{rep}.bed
    """
}

process intersect_peaks {
    input:
    file(peaks1) from file("results/bed/rep1.bed"),
    file(peaks2) from file("results/bed/rep2.bed")

    output:
    file("results/all_intersected.bed")

    script:
    """
    bedtools intersect -a ${peaks1} -b ${peaks2} -f 0.5 -r > results/all_intersected.bed
    """
}

process filter_blacklist {
    input:
    file(peaks) from file("results/all_intersected.bed"),
    file(blacklist) from file("reference/hg38-blacklist.v2.bed")

    output:
    file("results/filtered.bed")

    script:
    """
    bedtools intersect -a ${peaks} -b ${blacklist} -v > results/filtered.bed
    """
}

process unzip_gtf {
    input:
    file(gtf) from file("reference/gencode.v45.primary_assembly.annotation.gtf.gz")

    output:
    file("reference/gencode.v45.primary_assembly.annotation.gtf")

    script:
    """
    gunzip ${gtf}
    """
}

process annotate_peaks {
    input:
    file(peaks) from glob("results/filtered.bed"),
    file(gtf) from file("reference/gencode.v45.primary_assembly.annotation.gtf"),
    file(genome) from file("reference/GRCh38.primary_assembly.genome.fa")

    output:
    file("results/annotated.txt")

    script:
    """
    annotatePeaks.pl ${peaks} hg38 -gtf ${gtf} > results/annotated.txt
    """
}

process unzip_genome {
    input:
    file(genome) from file("reference/GRCh38.primary_assembly.genome.fa.gz")

    output:
    file("reference/GRCh38.primary_assembly.genome.fa")

    script:
    """
    gunzip ${genome}
    """
}

process motifs {
    input:
    file(peaks) from glob("results/filtered.bed"),
    file(genome) from file("reference/GRCh38.primary_assembly.genome.fa")

    output:
    directory("results/motifs"),
    file("results/motifs/knownResults.html")

    script:
    """
    findMotifsGenome.pl ${peaks} ${genome} results/motifs -size 200 -mask
    """
}

process computeMatrix {
    input:
    file(bigwig) from file("results/bigwig/RUNX1_{wildcards.rep}.bw"),
    file(bed) from file("reference/hg38_genes.bed")

    output:
    file("results/RUNX1_{rep}_matrix.gz")

    script:
    """
    computeMatrix scale-regions -R ${bed} -S ${bigwig} -b 2000 -a 2000 --outFileName results/RUNX1_{rep}_matrix.gz
    """
}

process plotMatrix {
    input:
    file(matrix) from file("results/RUNX1_{wildcards.rep}_matrix.gz")

    output:
    file("results/RUNX1_{rep}_coverage_plot.png")

    script:
    """
    plotProfile -m ${matrix} -o results/RUNX1_{rep}_coverage_plot.png
    """
}


workflow {
    sample_csv.collect { condition, rep ->
        [condition: condition, rep: rep]
    }.into { samples ->
        input:
        set condition, file(rep) from samples.collectEntries { [it.condition, file("samples/${it.condition}_${it.rep}.fastq.gz")] }
    }

    wget_files()
    bowtie2_build_gencode()
    fastqc()
    trimmomatic()
    bowtie2_align()
    samtools_sort()
    samtools_idx()
    samtools_flagstats()
    multiqc()
    bamCoverage()
    multiBwSummary()
    plotCorrelation()
    make_tag_dir()
    findPeaks()
    convertPeakFiles()
    intersect_peaks()
    filter_blacklist()
    unzip_gtf()
    annotate_peaks()
    unzip_genome()
    motifs()
    computeMatrix()
    plotMatrix()
}
