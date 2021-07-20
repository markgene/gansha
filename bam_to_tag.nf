// Convert PE BAM to tagAlign (BED 3+3 format)
// Inspired by ENCODE ATAC-seq Pipeline v1 section 2.a
//
// Link of ENCODE ATACSeq Pipeline v1 specifications (2019/09/18)
// https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit

Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .set { ch_flt_bam  }


/*
* Samtools sort by read name
*/
process SamtoolsSortByReadName {
    tag "$name"

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename -> params.save_align_intermeds ? filename : null }

    input:
    set val(name), file(bam) from ch_flt_bam

    output:
    set val(name), file("*.nameSrt.bam") into ch_name_sort_bam

    script:
    prefix = "${name}.nameSrt"
    """
    samtools sort -n -@ $task.cpus -o ${prefix}.bam -T $name $bam
    """
}


/*
* Create bedpe file (zipped)
*/
process BamToBedpe {
    tag "$name"

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename -> params.save_bedpe_intermeds ? filename : null }

    input:
    set val(name), file(bam) from ch_name_sort_bam

    output:
    set val(name), file("*.bedpe.gz") into ch_bedpe_gz

    script:
    prefix = "${name}"
    """
    bedtools bamtobed -bedpe -mate1 -i ${bam[0]} | gzip -nc > ${prefix}.bedpe.gz
    """
}