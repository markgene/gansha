/*
Input: BAM files after filtering
Output: 
    1. BAM files split into target and spike-in
    2. BAM files split into target and spike-in nucleosome free regions (NFRs)
    3. BigWig files normalizing target by spike-in using all reads
    4. BigWig files normalizing target by spike-in using NFR reads
    5. Library complexity of target and spike-in
    6. Samtools summary by MultiQC

Steps

1. Split by placed chromosomes (excluding mito chromosome)
2. Library complexity of targets and spike-ins
3. Normalize target by spike-in library size.
    1. All reads.
    2. Nucleosome free regions.
*/

// params.split_save_align_intermeds = true

Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .into { ch_bam1; ch_bam2  }

if (params.target_chr) { ch_target_chr   = file(params.target_chr,  checkIfExists: true) } else { exit 1, "Chr list file of target genome not specified!" }
if (params.spikein_chr) { ch_spikein_chr = file(params.spikein_chr, checkIfExists: true) } else { exit 1, "Chr list file of spike-in genome not specified!" }
if (params.chrom_sizes) { ch_genome_sizes_bigwig = file(params.chrom_sizes, checkIfExists: true) } else { exit 1, "Chromosome sizes file not specified!" }

// ch_bam1.view()

process SplitTargetBAM {
    tag "Target BAM of $name"
    label "SplitTargetSpikeinBam"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                      if (params.single_end || params.split_save_align_intermeds) {
                          if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                          else if (filename.endsWith(".target.bam")) "bam/$filename"
                          else if (filename.endsWith(".target.bam.bai")) "bam/$filename"
                          else null
                      }
                }

    input:
    set val(name), file(bam) from ch_bam1
    file target_chr from ch_target_chr

    output:
    set val(name), file("*.target.bam") into ch_target_bam
    set val(name), file("*.target.bam.bai") into ch_target_bai
    set val(name), file("*.flagstat") into ch_target_bam_flagstat
    set val(name), file("*.target.bam"), file("*.target.bam.bai") into ch_target_lib_complexity // The BAM file is already sorted
    file "*.{idxstats,stats}" into ch_target_bam_stats_mqc

    script:
    """
    samtools index ${bam}
    samtools view ${bam} \$(cat ${target_chr} | tr '\\n' ' ') -@ $task.cpus -b -o ${name}.0.bam
    samtools sort -@ $task.cpus -o ${name}.target.bam ${name}.0.bam
    samtools index ${name}.target.bam
    samtools flagstat ${name}.target.bam > ${name}.target.bam.flagstat
    samtools idxstats ${name}.target.bam > ${name}.target.bam.idxstats
    samtools stats ${name}.target.bam > ${name}.target.bam.stats
    """
}

process SplitSpikeinBAM {
    tag "Spike-in BAM of $name"
    label "SplitTargetSpikeinBam"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                      if (params.single_end || params.split_save_align_intermeds) {
                          if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                          else if (filename.endsWith(".spikein.bam")) "bam/$filename"
                          else if (filename.endsWith(".spikein.bam.bai")) "bam/$filename"
                          else null
                      }
                }

    input:
    set val(name), file(bam) from ch_bam2
    file spikein_chr from ch_spikein_chr

    output:
    set val(name), file("*.spikein.bam") into ch_spikein_bam
    set val(name), file("*.spikein.bam.bai") into ch_spikein_bai
    set val(name), file("*.flagstat") into ch_spikein_bam_flagstat
    file "*.{idxstats,stats}" into ch_spikein_bam_stats_mqc
    set val(name), file("*.spikein.bam"), file("*.spikein.bam.bai") into ch_spikein_lib_complexity // The BAM file is already sorted

    script:
    """
    samtools index ${bam}
    samtools view ${bam} \$(cat ${spikein_chr} | tr '\\n' ' ') -@ $task.cpus -b -o ${name}.0.bam
    samtools sort -@ $task.cpus -o ${name}.spikein.bam ${name}.0.bam
    samtools index ${name}.spikein.bam
    samtools flagstat ${name}.spikein.bam > ${name}.spikein.bam.flagstat
    samtools idxstats ${name}.spikein.bam > ${name}.spikein.bam.idxstats
    samtools stats ${name}.spikein.bam > ${name}.spikein.bam.stats
    """
}

ch_target_bam
    .combine(ch_spikein_bam, by: 0)
    .into { ch_bigwig_bams; ch_bigwig_bams_nuc_free }

// ch_bigwig_bams.view()

process BigWig {
    tag "$name"
    label 'bigwig'
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                      if (filename.endsWith("scale_factor.txt")) "scale/$filename"
                      else if (filename.endsWith(".bigWig")) "bigwig/$filename"
                      else null
                }

    input:
    set val(name), file(target_bam), file(spikein_bam) from ch_bigwig_bams
    file sizes from ch_genome_sizes_bigwig

    output:
    set val(name), file("*.bigWig") into ch_bigwig_plotprofile
    file "*scale_factor.txt" into ch_bigwig_scale

    script:
    prefix = "${name}"
    pe_fragment = params.single_end ? "" : "-pc"
    extend = (params.single_end && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
    """
    samtools flagstat ${spikein_bam} > ${prefix}.spikein.flagstat
    SCALE_FACTOR=\$(grep 'mapped (' ${prefix}.spikein.flagstat | head -1 | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.per1mSpikein.scale_factor.txt
    genomeCoverageBed -ibam ${target_bam} -bg -scale \$SCALE_FACTOR $pe_fragment $extend | LC_ALL=C sort -k1,1 -k2,2n >  ${prefix}.per1mSpikein.bedGraph

    bedGraphToBigWig ${prefix}.per1mSpikein.bedGraph $sizes ${prefix}.per1mSpikein.bigWig
    """
}


/*
* BAM and BigWig file of reads generated from putative nucleosome regions
*/
process NucFreeBigWig {
    tag "${name}"

    publishDir path: "${params.outdir}/nuc_free", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats//$filename"
                          else if (filename.endsWith(".bigWig")) "bigwig/$filename"
                          else if (filename.endsWith(".flagstat.scale_factor.txt")) "bigwig/scale/$filename"
                          else if (filename.endsWith(".bam")) "bam/$filename"
                          else if (filename.endsWith(".bai")) "bam/$filename"
                          else null
                }

    input:
    set val(name), file(target_bam), file(spikein_bam) from ch_bigwig_bams_nuc_free
    file sizes from ch_genome_sizes_bigwig

    output:
    set val(name), file("${prefix}.bam"), file("${prefix}.bam.bai") into ch_nuc_free_macs2
    set val(name), file("*.bigWig") into ch_nuc_free_bw
    set val(name), file("*.flagstat.scale_factor.txt") into ch_nuc_free_scale_factor
    file "*.{flagstat,idxstats,stats}" into ch_nuc_free_mqc
    
    script:
    prefix = "${name}.nuc_free${params.nuc_free_max_len}"
    target_prefix = "${prefix}.target"
    spikein_prefix = "${prefix}.spikein"
    scale_mb = params.scale_to / 1000000
    prefix_scale_mb = "${prefix}.scale${scale_mb}mb"
    pe_fragment = params.single_end ? "" : "-pc"
    extend = (params.single_end && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
    """
    samtools view -h ${target_bam} -@ $task.cpus \\
        | awk 'BEGIN { FS="\\t"; SIZE=${params.nuc_free_max_len}; S2=SIZE*SIZE }  /^@/ { print \$0; next } { if (\$9*\$9 < S2) print \$0}' \\
        | samtools view -@ $task.cpus -Sb - > ${target_prefix}.bam
    samtools view -h ${spikein_bam} -@ $task.cpus \\
        | awk 'BEGIN { FS="\\t"; SIZE=${params.nuc_free_max_len}; S2=SIZE*SIZE }  /^@/ { print \$0; next } { if (\$9*\$9 < S2) print \$0}' \\
        | samtools view -@ $task.cpus -Sb - > ${spikein_prefix}.bam

    samtools index ${target_prefix}.bam
    samtools index ${spikein_prefix}.bam
    samtools idxstats ${spikein_prefix}.bam > ${spikein_prefix}.idxstats
    samtools flagstat ${spikein_prefix}.bam > ${spikein_prefix}.flagstat
    samtools stats ${spikein_prefix}.bam > ${spikein_prefix}.stats

    SCALE_FACTOR=\$(grep 'mapped (' ${spikein_prefix}.flagstat | head -1 | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.per1mSpikein.scale_factor.txt
    genomeCoverageBed -ibam ${target_prefix}.bam -bg -scale \$SCALE_FACTOR $pe_fragment $extend | LC_ALL=C sort -k1,1 -k2,2n >  ${prefix}.per1mSpikein.bedGraph
    bedGraphToBigWig ${prefix}.per1mSpikein.bedGraph $sizes ${prefix}.per1mSpikein.bigWig
    """
}


/*
 * MultiQC for nucleosome free regions
 */
process MultiQCFilterNucFree {
    label 'multiqc'
    publishDir "${params.outdir}/nuc_free/multiqc", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file ('nuc_free/samtools_stats/*') from ch_nuc_free_mqc.collect()

    output:
    file "*multiqc_report.html" into ch_nuc_free_multiqc_report
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
process LibraryComplexityTarget {
    tag "$name"

    publishDir "${params.outdir}/library_complexity", mode: 'copy'

    input:
    set val(name), file(bam), file(bai) from ch_target_lib_complexity

    output:
    file "*.nrf_pbc.txt" into ch_target_lib_complexity

    script:
    prefix = "${name}.target"
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
    samtools sort -n -o ${prefix}.nameSrt.bam ${bam[0]}
    bedtools bamtobed -bedpe -i ${prefix}.nameSrt.bam \\
      | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$4,\$6,\$9,\$10}' \\
      | grep -v "${params.mito_name}" | sort | uniq -c \\
      | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{m2 == 0 ? pbc2 = -1 : pbc2 = m1/m2; printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,pbc2}' > ${prefix}.nrf_pbc.txt
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
process LibraryComplexitySpikein {
    tag "$name"

    publishDir "${params.outdir}/library_complexity", mode: 'copy'

    input:
    set val(name), file(bam), file(bai) from ch_spikein_lib_complexity

    output:
    file "*.nrf_pbc.txt" into ch_spikein_lib_complexity

    script:
    prefix = "${name}.spikein"
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
    samtools sort -n -o ${prefix}.nameSrt.bam ${bam[0]}
    bedtools bamtobed -bedpe -i ${prefix}.nameSrt.bam \\
      | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$4,\$6,\$9,\$10}' \\
      | grep -v "${params.spikein_mito_name}" | sort | uniq -c \\
      | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{m2 == 0 ? pbc2 = -1 : pbc2 = m1/m2; printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,pbc2}' > ${prefix}.nrf_pbc.txt
    """
}