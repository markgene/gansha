// Filter
// The filtering is based on ENCODE ATAC-seq pipeline v1 specifications (2019/09/18).
//
// Link of ENCODE ATACSeq Pipeline v1 specifications (2019/09/18)
// https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
// 
// Input: BAM files, coordinate sorted
// Output: BAM files, coordinated sorted

Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .set { ch_input_bam  }  

 /*
 * PREPROCESSING - Prepare genome intervals for filtering
 */
if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) } else { exit 1, "Fasta file not specified!" }   

/*
 ATAC-seq BAM filter based on ENCODE ATAC-seq pipeline for paired-end

# =============================
# Remove  unmapped, mate unmapped
# not primary alignment, reads failing platform
# Only keep properly paired reads
# Obtain name sorted BAM file
# ==================
FILT_BAM_PREFIX="${OFPREFIX}.filt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
TMP_FILT_BAM_PREFIX="tmp.${FILT_BAM_PREFIX}.nmsrt"
TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam"

## for multimapping > 0:
samtools view -F 524 -f 2 -u ${RAW_BAM_FILE} | samtools -n /dev/stdin -o ${TMP_FILT_BAM_FILE}

TMP_FILT_FIXMATE_BAM_FILE="${TMP_FILT_BAM_PREFIX}.fixmate.bam"

samtools view -h ${TMP_FILT_BAM_FILE} | assign_multimappers.py -k $multimapping --paired-end | samtools fixmate -r /dev/stdin ${TMP_FILT_FIXMATE_BAM_FILE}

## for multimapping==0:
MAPQ_THRESH=30
samtools view -F 1804 -f 2 -q ${MAPQ_THRESH} -u ${RAW_BAM_FILE} | samtools sort -n /dev/stdin -o ${TMP_FILT_BAM_FILE} 
samtools fixmate -r ${TMP_FILT_BAM_FILE} ${TMP_FILT_FIXMATE_BAM_FILE}
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
    filter_params = "-F 524 -f 2"
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
                          else if (filename.endsWith(".fixmate.bam") && params.filter_save_int_bam) "bam/fixmate/$filename"
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

# Remove orphan reads (pair was removed)
# and read pairs mapping to different chromosomes
# Obtain position sorted BAM


samtools view -F 1804 -f 2 -u ${TMP_FILT_FIXMATE_BAM_FILE} | samtools sort /dev/stdin -o ${FILT_BAM_FILE}

rm  ${TMP_FILT_FIXMATE_BAM_FILE}

rm ${TMP_FILT_BAM_FILE}

*/
process RemoveOrphan {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".flagstat")) "samtools_stats/orphan_rm/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/orphan_rm/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/orphan_rm/$filename"
                          else if (filename.endsWith(".orphan_rm.bam") && params.filter_save_int_bam) "bam/orphan_rm/$filename"
                          else null
                }

    input:
    set val(name), file(bam) from ch_fixmate_bam

    output:
    set val(name), file("*.orphan_rm.bam") into ch_orphan_rm_bam // The BAM file is already sorted
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

# =============
# Mark duplicates
# =============
module add picard-tools/1.126

TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
MARKDUP="/srv/gs1/software/picard-tools/1.126/MarkDuplicates.jar"
DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc"

java -Xmx4G -jar ${MARKDUP} INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}

*/
process PicardMarkDuplicates {
    tag "$name"
    label "picard"

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
    set val(name), file("*.mark_dups.bam") into ch_mark_dups_bam // The BAM file is already sorted
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
    avail_mem = params.picard_mem
    """
    picard -Xmx${avail_mem}g MarkDuplicates \\
        -INPUT ${bam[0]} \\
        -OUTPUT ${prefix}.bam \\
        -ASSUME_SORTED true \\
        -REMOVE_DUPLICATES false \\
        -METRICS_FILE ${prefix}.metrics.txt \\
        -VALIDATION_STRINGENCY LENIENT \\
        -MAX_RECORDS_IN_RAM 250000 \\
        -TMP_DIR tmp

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

/*
 ATAC-seq BAM filter based on ENCODE ATAC-seq pipeline for paired-end
*/
process AtacBamEncodeFilterRound2 {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".flagstat")) "samtools_stats/filter_round2/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/filter_round2/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/filter_round2/$filename"
                          else filename
                }

    input:
    set val(name), file(bam) from ch_mark_dups_bam // The BAM file is already sorted

    output:
    set val(name), file("*.flt.bam"), file("*.flt.bam.bai") into ch_filter2_bam, ch_in_lib_complexity // The BAM file is already sorted
    file "*.{flagstat,idxstats,stats}" into ch_filter2_bam_stats_mqc

    script:
    prefix = "${name}.flt"
    filter_params = "-F 1804 -f 2"
    multimap_params = params.keep_multi_map ? "" : "-q ${params.mapq_threshold}"
    """
    samtools view \\
        $filter_params \\
        $multimap_params \\
        -@ $task.cpus \\
        -b -o ${prefix}.bam ${bam[0]}

    samtools index ${prefix}.bam
    samtools flagstat ${prefix}.bam > ${prefix}.flagstat
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
    samtools stats ${prefix}.bam > ${prefix}.stats
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
                      if (filename.endsWith("_metrics")) "picard_metrics/filter_round2/$filename"
                      else if (filename.endsWith(".pdf")) "picard_metrics/filter_round2/pdf/$filename"
                      else null
                }

    when:
    !params.skip_picard_metrics

    input:
    set val(name), file(bam), file(bai) from ch_filter2_bam
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


/*
 * MultiQC for Filter round 2

# ============================
# Remove duplicates
# Index final position sorted BAM
# Create final name sorted BAM
# ============================
FINAL_BAM_PREFIX="${OFPREFIX}.nodup"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_FILE}.bai"
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}


# Index Final BAM file
samtools index ${FINAL_BAM_FILE}

samtools sort -n --threads 10 ${FINAL_BAM_FILE} -O SAM  | SAMstats --sorted_sam_file -  --outf ${FINAL_BAM_FILE_MAPSTATS}

 */
process MultiQCFilterRound2 {
    label 'multiqc'
    publishDir "${params.outdir}/multiqc/filter_round2", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file ('samtools_stats/filter_round2/*') from ch_filter2_bam_stats_mqc.collect()
    file ('picard_metrics/filter_round2/*_metrics') from ch_collectmetrics_mqc.collect()

    output:
    file "*multiqc_report.html" into ch_multiqc_report_filter2
    file "*_data"
    file "multiqc_plots"

    script:
    """
    multiqc . -f -p
    """
}


/*
 * Library complexity metrics based on ENCODE Guidelines

# =============================
# Compute library complexity
# =============================
# Sort by name
# convert to bedPE and obtain fragment coordinates
# sort by position and strand
# Obtain unique count statistics

module add bedtools/2.26

PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair


samtools sort -n ${FILT_BAM_FILE} -o ${OFPREFIX}.srt.tmp.bam
bedtools bamtobed -bedpe -i ${OFPREFIX}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}
rm ${OFPREFIX}.srt.tmp.bam


rm ${FILT_BAM_FILE}


 */
process LibraryComplexity {
    tag "$name"

    publishDir "${params.outdir}/library_complexity", mode: 'copy'

    input:
    set val(name), file(bam), file(bai) from ch_in_lib_complexity

    output:
    file "*.nrf_pbc.txt" into ch_lib_complexity

    script:
    prefix = "${name}.flt"
    mito_filter = params.spikein_mito_name ? "${params.spikein_mito_name}\\|${params.mito_name}" : params.mito_name
    // Format:
    // TotalReadPairs [tab]: mt
    // DistinctReadPairs [tab]: m0
    // OneReadPair [tab]: m1
    // TwoReadPairs [tab]: m2
    // NRF=Distinct/Total [tab]: m0/mt 
    // PBC1=OnePair/Distinct [tab]: m1/m0 
    // PBC2=OnePair/TwoPair: m1/m2
    //
    // Explain:
    // 1. The bamtobed command and awk output chrom1, start1, chrom2, end2, strand1, strand2.
    // 2. Mito reads are skipped.
    // 3. Unique read pairs defined by the field above.
    // 4. Count.
    //
    // nameSrt suffix is consistent with Samtools fixmate manual
    """
    samtools sort -@ $task.cpus -n -o ${prefix}.nameSrt.bam ${bam[0]}
    bedtools bamtobed -bedpe -i ${prefix}.nameSrt.bam \\
      | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$4,\$6,\$9,\$10}' \\
      | grep -v "${mito_filter}" | sort | uniq -c \\
      | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{m2 == 0 ? pbc2 = -1 : pbc2 = m1/m2; printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,pbc2}' > ${prefix}.nrf_pbc.txt
    """
}
