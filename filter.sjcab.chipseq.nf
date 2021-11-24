// Filter 
// Using the same parameters as SJ CAB ChIP-seq (2021/11/19)
// https://wiki.stjude.org/display/CAB/ChIPseq+QC+and+peak+calling
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
if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) } else { exit 1, "Fasta file not specified!" }


/*
 SJ CAB ChIP-seq filter
*/
process SJCABChIPFilter {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".flagstat")) "samtools_stats/SJCAB_filter/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/SJCAB_filter/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/SJCAB_filter/$filename"
                          else filename
                }

    input:
    set val(name), file(bam) from ch_input_bam

    output:
    set val(name), file("*.sorted.bam") into ch_filter_bam, ch_picard_bam
    file "*.{flagstat,idxstats,stats}" into ch_filter_bam_stats_mqc

    script:
    prefix = "${name}.flt"
    filter_params = "-F 1024 -b"
    multimap_params = params.keep_multi_map ? "" : "-q ${params.mapq_threshold}"
    """
    samtools view \\
        $filter_params \\
        $multimap_params \\
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
 * MultiQC for Filter
 */
process MultiQCFilter {
    label 'multiqc'
    publishDir "${params.outdir}/multiqc/SJCAB_filter", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file ('samtools_stats/SJCAB_filter/*') from ch_filter_bam_stats_mqc.collect()

    output:
    file "*multiqc_report.html" into ch_multiqc_report_filter1
    file "*_data"
    file "multiqc_plots"

    script:
    """
    multiqc . -f -p -m samtools
    """
}


/*
 * Picard Collect Multiple Metrics
 */
process PicardCollectMultipleMetrics {
    tag "$name"
    label "picard"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                      if (filename.endsWith("_metrics")) "picard_metrics/SJCAB_filter/$filename"
                      else if (filename.endsWith(".pdf")) "picard_metrics/SJCAB_filter/pdf/$filename"
                      else null
                }

    when:
    !params.skip_picard_metrics

    input:
    set val(name), file(bam), file(bai) from ch_picard_bam
    file fasta from ch_fasta

    output:
    file "*_metrics" into ch_collectmetrics_mqc
    file "*.pdf" into ch_collectmetrics_pdf

    script:
    prefix = "${name}.flt"
    avail_mem = params.picard_mem
    """
    picard -Xmx${avail_mem}g CollectMultipleMetrics \\
        INPUT=${bam[0]} \\
        OUTPUT=${prefix}.CollectMultipleMetrics \\
        REFERENCE_SEQUENCE=$fasta \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=tmp
    """
}


