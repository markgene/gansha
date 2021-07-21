/*
* Run MACS2 for ATAC-seq peak calling. 
* The parameters (especially extsize) is setup according to author's 
* recommendation on GitHub, SJ CAB and ENCODE ATAC-seq pipeline.
*
* Input: BAM files, coordinate sorted
* Output: MACS2 output
*/

// Check input
Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .set { ch_bam  }  
ch_peak_count_header = file("$baseDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_frip_score_header = file("$baseDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_blacklist = file(params.blacklist, checkIfExists: true)

/*
* ATAC-seq peak calling with MACS2
*/
process Macs2Atac {
    tag "${name}"
    label 'macs2'

    publishDir "${params.outdir}/macs2_atac/${name}", mode: 'copy',
        saveAs: { filename ->
                      if (filename.endsWith(".tsv")) "qc/$filename"
                      else filename
                }

    when:
    params.macs2_gsize

    input:
    set val(name), val(bam) from ch_bam
    file blacklist from ch_blacklist
    file peak_count_header from ch_peak_count_header
    file frip_score_header from ch_frip_score_header

    output:
    set val(name), file("*{bed,xls,Peak,bdg}") into ch_macs_output
    file "*_mqc.tsv" into ch_macs_mqc
    
    script:
    peak_type = params.macs2_peak_type
    broad = peak_type == "narrow" ? '' : "--broad --broad-cutoff ${params.macs2_broad_cutoff}"
    format = params.single_end ? "BAM" : "BAMPE"
    pileup = params.macs2_save_pileup ? "-B --SPMR" : ""
    extsize = params.macs2_extsize > 0 ? "--extsize ${params.macs2_extsize}" : ""
    macs2_peak_file = "${name}_peaks.${peak_type}Peak"
    macs2_blacklist_filtered = "${name}.blFlt.${peak_type}Peak"
    // Explain:
    // 1. Run MACS2
    // 2. Filter blacklist (command line borrowed by ENCODE ATAC-seq)
    // 3. Calculate FRiP
    """
    macs2 callpeak \\
        -t ${bam[0]} \\
        -q ${params.macs2_q} \\
        $broad \\
        -f $format \\
        -g $params.macs2_gsize \\
        -n ${name} \\
        $extsize \\
        $pileup \\
        --keep-dup all \\
        --nomodel

    bedtools intersect -v -a $macs2_peak_file -b $blacklist \\
            | awk 'BEGIN{OFS="\\t"} {if (\$5>1000) \$5=1000; print \$0}' \\
            | grep -P 'chr[\\dXY]+[ \\t]' > $macs2_blacklist_filtered

    cat $macs2_blacklist_filtered | wc -l | awk -v OFS='\\t' '{ print "${name}", "No", "${peak_type}", \$1 }' | cat $peak_count_header - > ${name}_blFlt_peaks.count_mqc.tsv

    READS_IN_PEAKS=\$(intersectBed -a ${bam[0]} -b $macs2_blacklist_filtered -bed -c -f ${params.frip_min_overlap} | awk -F '\\t' '{sum += \$NF} END {print sum}')
    samtools view -c ${bam[0]} | awk -v a="\$READS_IN_PEAKS" -v OFS='\\t' '{print "${name}", "No", "${peak_type}", a/\$1, a, \$1}' | cat $frip_score_header - > ${name}_blFlt_peaks.FRiP_mqc.tsv
    """
}