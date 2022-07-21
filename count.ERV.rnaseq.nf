// Filter RNA-seq for ERV detection
// 
// Input: BAM files, coordinate sorted
// Output: 
//   1. BAM files, coordinated sorted
//   2. Tab files of read count per ERV region taking account strand info or not.
//
//   The table has seven columns. The last column is the read count. 
//   By default, overlap of 1 bp is counted.

Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .set { ch_filter_bam  }  

 /*
 * PREPROCESSING - Prepare genome intervals for filtering
 */
if (params.erv_bed) { ch_erv_bed = file(params.erv_bed, checkIfExists: true) } else { exit 1, "ERV BED file not specified!" }
if (params.genome_file) { ch_genome_file = file(params.genome_file, checkIfExists: true) } else { exit 1, "Genome file not specified!" }

/*
 Count with bedtools coverage
*/
process Count {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    set val(name), file(bam) from ch_filter_bam
    file erv_bed from ch_erv_bed
    file genome_file from ch_genome_file

    output:
    set val(name), file("*.tab") into ch_count

    script:
    prefix = "${name}.ERVflt.B${params.overlap_read}"
    filter_params = "-F 1804"
    """
    bedtools coverage -b ${bam[0]} -a ${erv_bed} -counts -sorted -g ${genome_file} -F ${params.overlap_read} > ${prefix}.tab
    bedtools coverage -b ${bam[0]} -a ${erv_bed} -counts -sorted -g ${genome_file} -F ${params.overlap_read} -s > ${prefix}.stranded.tab
    """
}