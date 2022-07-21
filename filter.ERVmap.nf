// IMPORTANT: I JUST STARTED WORKING ON IT. ALSO SEE PERL SCRIPT parse_bam.pl
// Filter RNA-seq for ERV detection
// Borrow from ERVmap (2022/07/19)
// https://github.com/mtokuyama/ERVmap/blob/master/scripts/erv_genome.pl
// /bin/bash -c '$bwa mem -t 8 -p $genome ${btrimout} | $samtools view -Sh -F4 - | $filter | $samtools view -bSh - > $bwabam'
// 
// Input: BAM files, coordinate sorted
// Output: 
//   1. BAM files, coordinated sorted
//   2. Tab files of read count per ERV region taking account strand info or not.
//
//   The table has seven columns. The last column is the read count. 
//   By default, overlap of 1 bp is counted.

// Default
params.mapq_threshold = 1

Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .set { ch_input_bam  }  

 /*
 * PREPROCESSING - Prepare genome intervals for filtering
 */
if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) } else { exit 1, "Fasta file not specified!" }


/*
 ERVmap filter
*/
process ERVmapFilter {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
                          if (filename.endsWith(".flagstat")) "samtools_stats/ERVmap_filter/$filename"
                          else if (filename.endsWith(".idxstats")) "samtools_stats/ERVmap_filter/$filename"
                          else if (filename.endsWith(".stats")) "samtools_stats/ERVmap_filter/$filename"
                          else filename
                }

    input:
    set val(name), file(bam) from ch_input_bam

    output:
    set val(name), file("*.sorted.bam") into ch_filter_bam, ch_picard_bam
    file "*.{flagstat,idxstats,stats}" into ch_filter_bam_stats_mqc

    script:
    prefix = "${name}.ERVmapFlt"
    filter_params = "-F 4"
    """
    samtools view \\
        $filter_params \\
        $multimap_params \\
        -@ $task.cpus \\
        -b -o ${prefix}.unsorted.bam ${bam[0]}

    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam ${prefix}.unsorted.bam    

    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.idxstats
    samtools stats ${prefix}.sorted.bam > ${prefix}.stats
    """
}