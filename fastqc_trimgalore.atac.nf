// Validate inputs
if (!params.outdir) {
    exit 1, "Output directory not specified!"
}

Channel
    .fromFilePairs( params.reads, size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .into { ch_raw_reads_fastqc; ch_raw_reads_trimgalore }

/*
 * FastQC
 */
process FastQC {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: { filename ->
                      filename.endsWith(".zip") ? "zips/$filename" : "$filename"
                }

    input:
    set val(name), file(reads) from ch_raw_reads_fastqc

    output:
    file "*.{zip,html}" into ch_fastqc_reports_mqc

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    if (params.single_end) {
        """
        [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
        fastqc -q -t $task.cpus ${name}.fastq.gz
        """
    } else {
        """
        [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
        [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
        fastqc -q -t $task.cpus ${name}_1.fastq.gz
        fastqc -q -t $task.cpus ${name}_2.fastq.gz
        """
    }
}

/*
 * STEP 2 - Trim Galore!
 */
process TrimGalore {
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: { filename ->
                        if (filename.endsWith(".html")) "fastqc/$filename"
                        else if (filename.endsWith(".zip")) "fastqc/zips/$filename"
                        else if (filename.endsWith("trimming_report.txt")) "logs/$filename"
                        else params.save_trimmed ? filename : null
                }

    input:
    set val(name), file(reads) from ch_raw_reads_trimgalore

    output:
    set val(name), file("*.fq.gz") into ch_trimmed_reads
    file "*.txt" into ch_trimgalore_results_mqc
    file "*.{zip,html}" into ch_trimgalore_fastqc_reports_mqc

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
    if (params.single_end) {
        """
        [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
        trim_galore --fastqc --gzip --length ${params.trim_min_len} $c_r1 $tpc_r1 $nextseq ${name}.fastq.gz -j $task.cpus
        """
    } else {
        """
        [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
        [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
        trim_galore --paired --fastqc --gzip --length ${params.trim_min_len} $c_r1 $c_r2 $tpc_r1 $tpc_r2 $nextseq ${name}_1.fastq.gz ${name}_2.fastq.gz -j $task.cpus
        mv ${name}_1_val_1.fq.gz ${name}_R1_trimmed.fq.gz
        mv ${name}_2_val_2.fq.gz ${name}_R2_trimmed.fq.gz
        """
    }
}

/*
 * MultiQC for FastQC
 */
process MultiQCFastQC {
    label 'multiqc'
    publishDir "${params.outdir}/multiqc/fastqc", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file ('fastqc/*') from ch_fastqc_reports_mqc.collect().ifEmpty([])

    output:
    file "*multiqc_report.html" into ch_multiqc_report_fastqc
    file "*_data"
    file "multiqc_plots"

    script:
    """
    multiqc . -f -p -m fastqc
    """
}

/*
 * MultiQC for Trim Galore
 */
process MultiQCTrimGalore {
    label 'multiqc'
    publishDir "${params.outdir}/multiqc/trimgalore", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file ('trimgalore/*') from ch_trimgalore_results_mqc.collect().ifEmpty([])
    file ('trimgalore/fastqc/*') from ch_trimgalore_fastqc_reports_mqc.collect().ifEmpty([])

    output:
    file "*multiqc_report.html" into ch_multiqc_report_trimgalore
    file "*_data"
    file "multiqc_plots"

    script:
    """
    multiqc . -f -p -m cutadapt -m fastqc
    """
}
