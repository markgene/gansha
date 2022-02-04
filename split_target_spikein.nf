// Input: BAMs.
// Output: Pairs of split BAMs.

// params.save_align_intermeds = true

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
                      if (params.single_end || params.save_align_intermeds) {
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
                      if (params.single_end || params.save_align_intermeds) {
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
    .set { ch_bigwig_bams }

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
    SCALE_FACTOR=\$(grep 'mapped (' ${prefix}.spikein.flagstat | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.per1mSpikein.scale_factor.txt
    genomeCoverageBed -ibam ${target_bam} -bg -scale \$SCALE_FACTOR $pe_fragment $extend | LC_ALL=C sort -k1,1 -k2,2n >  ${prefix}.per1mSpikein.bedGraph

    bedGraphToBigWig ${prefix}.per1mSpikein.bedGraph $sizes ${prefix}.per1mSpikein.bigWig
    """
}