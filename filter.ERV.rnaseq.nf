// Filter RNA-seq for ERV detection
// 
// Input: BAM files, coordinate sorted
// Output: BAM files, coordinated sorted

// Default
params.mapq_threshold = 1

Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .set { ch_input_bam  }  

 /*
 * PREPROCESSING - Prepare genome intervals for filtering
 */
if (params.erv_bed) { ch_erv_bed = file(params.erv_bed, checkIfExists: true) } else { exit 1, "ERV BED file not specified!" }
if (params.genome_file) { ch_genome_file = file(params.genome_file, checkIfExists: true) } else { exit 1, "Genome file not specified!" }


/*
 ERVmap filter
*/
process ERVFilter {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".flagstat")) "samtools_stats/ERV_filter/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/ERV_filter/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/ERV_filter/$filename"
                          else filename
                }

    input:
    set val(name), file(bam) from ch_input_bam

    output:
    set val(name), file("*.sorted.bam") into ch_filter_bam, ch_picard_bam
    file "*.{flagstat,idxstats,stats}" into ch_filter_bam_stats_mqc

    script:
    prefix = "${name}.ERVflt"
    filter_params = "-F 1804"
    """
    samtools view \\
        $filter_params \\
        -@ $task.cpus \\
        -b -o ${prefix}.unsorted.bam ${bam[0]}

    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam ${prefix}.unsorted.bam    

    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.idxstats
    samtools stats ${prefix}.sorted.bam > ${prefix}.stats
    """
}


/*
 Count with bedtools coverage
*/
process Count {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    set val(name), file(bam) from ch_filter_bam
    file erv_bed from ch_erv_bed
    file genome_file from ch_genome_file

    output:
    set val(name), file("*.tab") into ch_count

    script:
    prefix = "${name}.ERVflt"
    filter_params = "-F 1804"
    """
    bedtools coverage -b ${bam[0]} -a ${erv_bed} -counts -sorted -g ${genome_file} > ${prefix}.tab
    bedtools coverage -b ${bam[0]} -a ${erv_bed} -counts -sorted -g ${genome_file} -s > ${prefix}.stranded.tab
    """
}