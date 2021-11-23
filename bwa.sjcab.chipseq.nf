// BWA 
// Using the same parameters as SJ CAB ChIP-seq (2021/11/19)
// https://wiki.stjude.org/display/CAB/ChIPseq+QC+and+peak+calling
//
// Input: FASTQ files
// Output: BAM files, coordinated sorted

// Check parameters
if (!params.read_len) {
    exit 1, "Read length not specified!"
}

if (params.fasta) {
    ch_fasta = file(params.fasta, checkIfExists: true)
} else {
    exit 1, "Fasta file not specified!"
}

if (params.bwa_index) {
    lastPath = params.bwa_index.lastIndexOf(File.separator)
    bwa_dir =  params.bwa_index.substring(0,lastPath+1)
    bwa_base = params.bwa_index.substring(lastPath+1)
    Channel
        .fromPath(bwa_dir, checkIfExists: true)
        .set { ch_bwa_index }
}

Channel
    .fromFilePairs( params.reads, size: 1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .set { ch_reads }

/*
* Align read BWA
*/
process BwaAlign {
    tag "$name"

    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                        if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                        else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                        else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                        else filename
                }

    input:
    set val(name), file(reads) from ch_reads
    file index from ch_bwa_index.collect()

    output:
    set val(name), file("*.marked.bam") into ch_bwa_bam

    script:
    prefix = "${name}.bwa.unsorted"
    script:
    if (params.read_len >= 75) {
        """
        bwa mem -t $task.cpus ${index}/${bwa_base} ${reads[0]} | samtools view -b - > ${prefix}.unmark.bam
        cat ${prefix}.unmark.bam | bamsormadup > ${prefix}.marked.bam
        """
    } else {
        """
        bwa aln -t $task.cpus ${index}/${bwa_base} ${reads[0]} > ${name}.fq.gz.sai
        bwa samse ${index}/${bwa_base} ${name}.fq.gz.sai ${reads[0]} | samtools view -b - > ${prefix}.unmark.bam
        cat ${prefix}.unmark.bam | bamsormadup > ${prefix}.marked.bam
        """
    }
}


/*
 * Convert .bam to coordinate sorted .bam
 */
process SamtoolsSortDefault {
    tag "$name"

    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                        if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                        else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                        else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                        else filename
                }
    
    input:
    set val(name), file(bam) from ch_bwa_bam

    output:
    set val(name), file("*.sorted.{bam,bam.bai}") into ch_sort_bam_merge
    file "*.{flagstat,idxstats,stats}" into ch_sort_bam_flagstat_mqc

    script:
    prefix = "${name}.bwa.sorted.marked"
    """
    samtools sort -@ $task.cpus -o ${prefix}.bam -T $name $bam
    samtools index ${prefix}.bam
    samtools flagstat ${prefix}.bam > ${prefix}.flagstat
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
    samtools stats ${prefix}.bam > ${prefix}.stats
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
    file "*multiqc_report.html" into ch_multiqc_report_bwa
    file "*_data"
    file "multiqc_plots"

    script:
    """
    multiqc . -f -p -m samtools
    """
}