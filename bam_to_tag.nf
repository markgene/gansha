// Convert PE BAM to tagAlign (BED 3+3 format)
// Inspired by ENCODE ATAC-seq Pipeline v1 section 2.a
//
// Link of ENCODE ATACSeq Pipeline v1 specifications (2019/09/18)
// https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
//
// Input: filtered BAM files
// Output: to be added

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
* Create BEDPE file (zipped)
*/
process BamToBedpe {
    tag "$name"

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename -> params.save_bedpe_intermeds ? filename : null }

    input:
    set val(name), file(bam) from ch_name_sort_bam

    output:
    set val(name), file("*.bedpe.gz") into ch_bedpe_gz, ch_subsample_bedpe_gz

    script:
    prefix = "${name}"
    """
    bedtools bamtobed -bedpe -mate1 -i ${bam[0]} | gzip -nc > ${prefix}.bedpe.gz
    """
}

/*
* Create tagAlign file (zipped) from BEDPE file (zipped)
*/
process BedpeToTagAlign {
    tag "$name"

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename -> params.save_ta_intermeds ? filename : null }

    input:
    set val(name), file(bedpe) from ch_bedpe_gz

    output:
    set val(name), file("${prefix}.ta.gz") into ch_ta_gz
    set val(name), file("${subsample_ta_file}") into ch_subsampled_ta_gz

    script:
    prefix = "${name}"
    subsample_mb = params.xcor_subsample_size / 1000000
    subsample_ta_file = "${prefix}.sample${subsample_mb}mb.mate1.ta.gz"
    // The so-called tagAlign format is BED6:
    // 1. Chrom1
    // 2. Start1
    // 3. End1
    // 4. Name: all set to N
    // 5. Score: all set to 1000
    // 6. Strand1
    // 
    // Explain: the read pair on each line of BEDPE file is split into two rows in the output.
    """
    zcat ${bedpe} | awk 'BEGIN{OFS="\\t"}{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' | gzip -nc > ${prefix}.ta.gz

    zcat ${bedpe} \\
      | grep -v "${params.mito_name}" \\
      | shuf -n ${params.xcor_subsample_size} --random-source=<(openssl enc -aes-256-ctr -pass pass:\$(zcat -f ${prefix}.ta.gz | wc -c) -nosalt </dev/zero 2>/dev/null) \\
      | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3,"N","1000",\$9}' | gzip -nc > ${subsample_ta_file}
    """
}
