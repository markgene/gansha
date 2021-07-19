// Bowtie2 
// Using the same parameters as ENCODE ATACSeq Pipeline v1 specifications (2019/09/18)
// https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit

// Check parameters
if (params.fasta) {
    lastPath = params.fasta.lastIndexOf(File.separator)
    bowtie2_base = params.fasta.substring(lastPath+1)
    ch_fasta = file(params.fasta, checkIfExists: true)
} else {
    exit 1, "Fasta file not specified!"
}

if (params.bowtie2_index) {
    lastPath = params.bowtie2_index.lastIndexOf(File.separator)
    bowtie2_dir =  params.bowtie2_index.substring(0,lastPath+1)
    bowtie2_base = params.bowtie2_index.substring(lastPath+1)
    Channel
        .fromPath(bowtie2_dir, checkIfExists: true)
        .set { ch_bowtie2_index }
}

Channel
    .fromFilePairs( params.reads, size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .set { ch_reads }

 /*
 * Align read Bowtie2
 */
process Bowtie2Align {
    tag "$name"

    input:
    set val(name), file(reads) from ch_reads
    file index from ch_bowtie2_index.collect()

    output:
    set val(name), file("*.bam") into ch_bowtie2_bam

    script:
    prefix = "${name}.bowtie2.unsorted"
    // sm = "SM:${name.split('_')[0..-2].join('_')}"
    sm = "SM:${name}"
    pl = "PL:ILLUMINA"
    lb = "LB:${name}"
    rg = "\'@RG\\tID:${name}\\tSM:${sm}\\tPL:ILLUMINA\\tLB:${name}\\tPU:1\'"
    if (params.seq_center) {
        rg = "\'@RG\\tID:${name}\\tSM:${sm}\\tPL:ILLUMINA\\tLB:${name}\\tPU:1\\tCN:${params.seq_center}\'"
    }
    // Using the same parameters as ENCODE ATACSeq Pipeline v1 specifications (2019/09/18)
    // https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
    // No multimapping allowed
    """
    bowtie2 --local -q \\
    -p $task.cpus \\
    --end-to-end --very-sensitive -X ${params.bowtie2_encode_maxins} \\
    --rg-id $name \\
    --rg $sm --rg $pl --rg $lb --rg PU:1 \\
    -x ${index}/${bowtie2_base} \\
    -1 ${reads[0]} -2 ${reads[1]} \\
    | samtools view -@ $task.cpus -b -h -O BAM -o ${prefix}.bam -
    """
}


/*
 * Convert .bam to coordinate sorted .bam
 */
process SortBAM {
    tag "$name"

    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                        if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                        else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                        else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                        else filename
                }
    
    input:
    set val(name), file(bam) from ch_bowtie2_bam

    output:
    set val(name), file("*.sorted.{bam,bam.bai}") into ch_sort_bam_merge
    file "*.{flagstat,idxstats,stats}" into ch_sort_bam_flagstat_mqc

    script:
    prefix = "${name}.bowtie2.sorted"
    """
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $name $bam
    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
    """
}


/*
 * MultiQC for Trim Galore
 */
process MultiQCBowtie2 {
    label 'multiqc'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file ('samtools_stats/*') from ch_sort_bam_flagstat_mqc.collect()

    output:
    file "*multiqc_report.html" into ch_multiqc_report_bowtie2
    file "*_data"
    file "multiqc_plots"

    script:
    """
    multiqc . -f -p -m samtools
    """
}