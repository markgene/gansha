/*
* The workflow has a simple function: create BigWig file from sorted BAM file
*
* Input: BAM files, coordinate sorted
* Output: BigWig files
*/


// Check input
ch_chrom_size = file(params.chrom_sizes, checkIfExists: true)
Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .set { ch_input_bam  }  

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