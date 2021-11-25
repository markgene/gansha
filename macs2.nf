// Default MACS peak calling.
// Narrow peak: NOT filter by insert size.
// Broad peak

// Header files for MultiQC
ch_peak_count_header = file("$baseDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_frip_score_header = file("$baseDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_blacklist = file(params.blacklist, checkIfExists: true)


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     PARSE DESIGN FILE                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Channel.fromPath(params.design, checkIfExists: true)
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.group, row.ip_name, row.control_name, row.peak_type, file(row.ip_bam, checkIfExists: true), file(row.control_bam, checkIfExists: true) ]  }
    .set { ch_design_macs }

// For debugging:
// ch_design_macs.view { row -> "${row[0]} - ${row[1]} - ${row[2]} - ${row[5]}"}    


/*
 * Call peaks with MACS2 and calculate FRiP score
 */
process Macs2 {
    tag "${ip_name} vs ${control_name}, ${peak_type}"
    label 'macs2'
    publishDir "${params.outdir}/macs/${ip_name}/${control_name}/${peak_type}", mode: 'copy',
        saveAs: { filename ->
                      if (filename.endsWith(".tsv")) "qc/$filename"
                      else if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                      else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                      else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                      else filename
                }

    when:
    params.macs_gsize

    input:
    set val(group_name), val(ip_name), val(control_name), val(peak_type), file(ip_bam), file(control_bam) from ch_design_macs
    file blacklist from ch_blacklist
    file peak_count_header from ch_peak_count_header
    file frip_score_header from ch_frip_score_header

    output:
    set val(ip_name), file("*{bed,xls,Peak,bdg}") into ch_macs_output
    file "*_mqc.tsv" into ch_macs_mqc
    file "*.{flagstat,idxstats,stats}" into ch_macs_bam_stats_mqc
    // set val(antibody), val(replicatesExist), val(multipleGroups), val(ip), val(control), file("*.$PEAK_TYPE") into ch_macs_homer, ch_macs_qc, ch_macs_consensus
    
    script:
    broad = peak_type == "narrow" ? '' : "--broad --broad-cutoff ${params.macs2_broad_cutoff}"
    format = params.single_end ? "BAM" : "BAMPE"
    pileup = params.macs2_save_pileup ? "-B --SPMR" : ""
    macs2_peak_file = "${ip_name}_peaks.${peak_type}Peak"
    macs2_blacklist_filtered = "${ip_name}.blFlt.${peak_type}Peak"
    """
    samtools index ${ip_bam}
    samtools idxstats ${ip_bam} > ${ip_name}.idxstats
    samtools flagstat ${ip_bam} > ${ip_name}.flagstat
    samtools stats ${ip_bam} > ${ip_name}.stats

    samtools index ${control_bam}
    samtools idxstats ${control_bam} > ${control_name}.idxstats
    samtools flagstat ${control_bam} > ${control_name}.flagstat
    samtools stats ${control_bam} > ${control_name}.stats

    macs2 callpeak \\
        -t ${ip_bam} \\
        -c ${control_bam} \\
        $broad \\
        -f $format \\
        -g $params.macs_gsize \\
        -n ${ip_name} \\
        $pileup \\
        --keep-dup all

    bedtools intersect -v -a $macs2_peak_file -b $blacklist \\
            | awk 'BEGIN{OFS="\\t"} {if (\$5>1000) \$5=1000; print \$0}' \\
            | grep -P 'chr[\\dXY]+[ \\t]' > $macs2_blacklist_filtered

    cat $macs2_blacklist_filtered | wc -l | awk -v OFS='\\t' '{ print "${ip_name}", "No", "${peak_type}", \$1 }' | cat $peak_count_header - > ${ip_name}_blFlt_peaks.count_mqc.tsv

    READS_IN_PEAKS=\$(intersectBed -a ${ip_bam} -b $macs2_blacklist_filtered -bed -c -f ${params.frip_min_overlap} | awk -F '\\t' '{sum += \$NF} END {print sum}')
    samtools view -c ${ip_bam} | awk -v a="\$READS_IN_PEAKS" -v OFS='\\t' '{print "${ip_name}", "No", "${peak_type}", a/\$1, a, \$1}' | cat $frip_score_header - > ${ip_name}_blFlt_peaks.FRiP_mqc.tsv
    """
}

