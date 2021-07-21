/*
* The workflow has multiple functions:
* 
* 1. Partition alignments into reads generated from putative nucleosome free 
* regions of DNA, and reads likely derived from nucleosome associated DNA 
* (in the future). 
* 2. Create BigWig track files.
* 3. Prepare the nucleosome free regions BAM file for peak calling. 
*
*
* Input: BAM files, coordinate sorted
* Output: 
*   1) BigWig files of all reads, scaled, 
*   2) BAM files of nucleosome free regions, coordinate sorted
*   3) BigWig files of nucleosome free regions, scaled
*/


// Examine input
ch_chrom_size = file(params.chrom_sizes, checkIfExists: true)
Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .into { ch_input_bam; ch_input_nuc_free  }  

/*
* BigWig file of all reads
*/
process BigWig {
    tag "${name}"

    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".bigWig")) "bigwig/$filename"
                          else if (filename.endsWith(".flagstat.scale_factor.txt")) "bigwig/scale/$filename"
                          else null
                }

    input:
    set val(name), file(bam) from ch_input_bam
    file sizes from ch_chrom_size

    output:
    set val(name), file("*.bigWig") into ch_bw
    set val(name), file("*.flagstat.scale_factor.txt") into ch_scale_factor
    
    script:
    prefix = "${name}"
    scale_mb = params.scale_to / 1000000
    prefix_scale_mb = "${prefix}.scale${scale_mb}mb"
    pe_fragment = params.single_end ? "" : "-pc"
    """
    samtools index ${bam[0]}
    samtools flagstat ${bam[0]} > ${prefix}.tmp.flagstat

    SCALE_FACTOR=\$(grep 'mapped (' ${prefix}.tmp.flagstat | head -1 | awk '{print ${params.scale_to}/\$1}')
    echo \$SCALE_FACTOR > ${prefix_scale_mb}.flagstat.scale_factor.txt
    genomeCoverageBed -ibam ${bam[0]} -bg -scale \$SCALE_FACTOR $pe_fragment | LC_ALL=C sort -k1,1 -k2,2n >  ${prefix_scale_mb}.bedGraph
    bedGraphToBigWig ${prefix_scale_mb}.bedGraph $sizes ${prefix_scale_mb}.bigWig
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
                          else filename
                }

    input:
    set val(name), file(bam) from ch_input_nuc_free
    file sizes from ch_chrom_size

    output:
    set val(name), file("${prefix}.bam") into ch_nuc_free_macs2
    set val(name), file("*.bigWig") into ch_nuc_free_bw
    set val(name), file("*.flagstat.scale_factor.txt") into ch_nuc_free_scale_factor
    file "*.{flagstat,idxstats,stats}" into ch_nuc_free_mqc
    
    script:
    prefix = "${name}.nuc_free${params.nuc_free_max_len}"
    scale_mb = params.scale_to / 1000000
    prefix_scale_mb = "${prefix}.scale${scale_mb}mb"
    pe_fragment = params.single_end ? "" : "-pc"
    """
    samtools view -h ${bam[0]} -@ $task.cpus \\
        | awk 'BEGIN { FS="\\t"; SIZE=${params.nuc_free_max_len}; S2=SIZE*SIZE }  /^@/ { print \$0; next } { if (\$9*\$9 < S2) print \$0}' \\
        | samtools view -@ $task.cpus -Sb - > ${prefix}.bam

    samtools index ${prefix}.bam
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
    samtools flagstat ${prefix}.bam > ${prefix}.flagstat
    samtools stats ${prefix}.bam > ${prefix}.stats

    SCALE_FACTOR=\$(grep 'mapped (' ${prefix}.flagstat | head -1 | awk '{print ${params.scale_to}/\$1}')
    echo \$SCALE_FACTOR > ${prefix_scale_mb}.flagstat.scale_factor.txt
    genomeCoverageBed -ibam ${prefix}.bam -bg -scale \$SCALE_FACTOR $pe_fragment | LC_ALL=C sort -k1,1 -k2,2n >  ${prefix_scale_mb}.bedGraph
    bedGraphToBigWig ${prefix_scale_mb}.bedGraph $sizes ${prefix_scale_mb}.bigWig
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