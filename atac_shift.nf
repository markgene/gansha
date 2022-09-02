/*
* Shift ATAC-seq alignments
*/


// Check input
Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .into { ch_input_bam }  

/*
* BigWig file of all reads
*/
process AlignmentSieve {
    label 'atac_shift'
    tag "${name}"

    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                          else filename
                }

    input:
    set val(name), file(bam) from ch_input_bam

    output:
    set val(name), file("*.shift.bam"), file("*.shift.bam.bai") into ch_shift_bam // The BAM file is already sorted
    file "*.{flagstat,idxstats,stats}" into ch_filter2_bam_stats_mqc

    script:
    prefix = "${name}.shift"
    """
    samtools index ${bam[0]}
    alignmentSieve -b ${bam[0]} -o ${name}.tmp.bam -p $task.cpus --ATACshift
    samtools sort -@ $task.cpus -O bam -o ${prefix}.bam ${name}.tmp.bam
    samtools index ${prefix}.bam
    samtools flagstat ${prefix}.bam > ${prefix}.flagstat
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
    samtools stats ${prefix}.bam > ${prefix}.stats
    """
}

