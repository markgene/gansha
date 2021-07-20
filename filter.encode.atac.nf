// Filter
// The filtering is based on ENCODE ATAC-seq pipeline v1 specifications (2019/09/18).
//
// Link of ENCODE ATACSeq Pipeline v1 specifications (2019/09/18)
// https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
// 
// Input: sorted BAMs
// Output: filtered BAMs.

Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .set { ch_input_bam  }  

/*
 ATAC-seq BAM filter based on ENCODE ATAC-seq pipeline for paired-end
*/
process AtacBamEncodeFilterRound1 {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".flagstat")) "samtools_stats/filter_round1/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/filter_round1/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/filter_round1/$filename"
                          else null

                }

    input:
    set val(name), file(bam) from ch_input_bam

    output:
    set val(name), file("*.unsorted.bam") into ch_filter1_bam
    file "*.{flagstat,idxstats,stats}" into ch_filter1_bam_stats_mqc

    script:
    prefix = "${name}.flt1"
    filter_params = "-F 1804 -f 2"
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
 * MultiQC for Filter round 1
 */
process MultiQCFilterRound1 {
    label 'multiqc'
    publishDir "${params.outdir}/multiqc/filter_round1", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file ('samtools_stats/filter_round1/*') from ch_filter1_bam_stats_mqc.collect()

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
 Fixmate
*/
process SamtoolsFixmate {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".flagstat")) "samtools_stats/fixmate/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/fixmate/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/fixmate/$filename"
                          else null

                }

    input:
    set val(name), file(bam) from ch_filter1_bam

    output:
    set val(name), file("*.fixmate.bam") into ch_fixmate_bam
    file "*.{flagstat,idxstats,stats}" into ch_fixmate_bam_stats_mqc

    script:
    prefix = "${name}.fixmate"
    """
    samtools sort -n -@ $task.cpus -o ${name}.sortByName.bam ${bam[0]}
    samtools fixmate -@ $task.cpus -r ${name}.sortByName.bam ${prefix}.bam
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam ${prefix}.bam

    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.idxstats
    samtools stats ${prefix}.sorted.bam > ${prefix}.stats
    """
}

/*
 * MultiQC for Fixmat
 */
process MultiQCFixmate {
    label 'multiqc'
    publishDir "${params.outdir}/multiqc/fixmate", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file ('samtools_stats/fixmate/*') from ch_fixmate_bam_stats_mqc.collect()

    output:
    file "*multiqc_report.html" into ch_multiqc_report_fixmate
    file "*_data"
    file "multiqc_plots"

    script:
    """
    multiqc . -f -p -m samtools
    """
}


/*
 Remove orphan
*/
process RemoveOrphan {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".flagstat")) "samtools_stats/orphan_rm/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/orphan_rm/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/orphan_rm/$filename"
                          else null

                }

    input:
    set val(name), file(bam) from ch_fixmate_bam

    output:
    set val(name), file("*.orphan_rm.bam") into ch_orphan_rm_bam
    file "*.{flagstat,idxstats,stats}" into ch_orphan_rm_bam_stats_mqc

    script:
    prefix = "${name}.orphan_rm"
    filter_params = "-F 1804 -f 2"
    multimap_params = params.keep_multi_map ? "" : "-q ${params.mapq_threshold}"
    """
    samtools view $filter_params $multimap_params -@ $task.cpus -u ${bam[0]} | samtools sort -@ $task.cpus -o ${prefix}.bam -

    samtools index ${prefix}.bam
    samtools flagstat ${prefix}.bam > ${prefix}.flagstat
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
    samtools stats ${prefix}.bam > ${prefix}.stats
    """
}

/*
 * MultiQC for Orphan removal
 */
process MultiQCOrphanRemoval {
    label 'multiqc'
    publishDir "${params.outdir}/multiqc/orphan_rm", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file ('samtools_stats/orphan_rm/*') from ch_orphan_rm_bam_stats_mqc.collect()

    output:
    file "*multiqc_report.html" into ch_multiqc_report_orphan_rm
    file "*_data"
    file "multiqc_plots"

    script:
    """
    multiqc . -f -p -m samtools
    """
}

/*
 Picard mark duplicates
*/
process PicardMarkDuplicates {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".flagstat")) "samtools_stats/mark_dups/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/mark_dups/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/mark_dups/$filename"
                          else if (filename.endsWith('.metrics.txt')) "picard_metrics/$filename"
                          else null
                }

    input:
    set val(name), file(bam) from ch_orphan_rm_bam

    output:
    set val(name), file("*.mark_dups.bam") into ch_mark_dups_bam
    file "*.{flagstat,idxstats,stats}" into ch_mark_dups_bam_stats_mqc
    file "*.metrics.txt" into ch_mark_dups_metrics_mqc

    script:
    // ENCODE ATAC-seq pipeline parameters:
    // INPUT=${FILT_BAM_FILE} 
    // OUTPUT=${TMP_FILT_BAM_FILE} 
    // METRICS_FILE=${DUP_FILE_QC} 
    // VALIDATION_STRINGENCY=LENIENT 
    // ASSUME_SORTED=true 
    // REMOVE_DUPLICATES=false
    prefix = "${name}.mark_dups"
    """
    picard -Xmx5g MarkDuplicates \\
        INPUT=${bam[0]} \\
        OUTPUT=${prefix}.bam \\
        ASSUME_SORTED=true \\
        REMOVE_DUPLICATES=false \\
        METRICS_FILE=${prefix}.metrics.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=tmp

    samtools index ${prefix}.bam
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
    samtools flagstat ${prefix}.bam > ${prefix}.flagstat
    samtools stats ${prefix}.bam > ${prefix}.stats
    """
}

/*
 * MultiQC for Picard MarkDuplicates
 */
process MultiQCPicardMarkDuplicates {
    label 'multiqc'
    publishDir "${params.outdir}/multiqc/mark_dups", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file ('samtools_stats/mark_dups/*') from ch_mark_dups_bam_stats_mqc.collect()
    file ('picard_metrics/*') from ch_mark_dups_metrics_mqc.collect()

    output:
    file "*multiqc_report.html" into ch_multiqc_report_mark_dups
    file "*_data"
    file "multiqc_plots"

    script:
    """
    multiqc . -f -p
    """
}