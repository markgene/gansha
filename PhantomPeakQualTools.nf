/*
* Cross-correlation analysis
*
* Input: BAM files, coordinate sorted
* Output: 
*/

Channel
    .fromFilePairs( params.bams, size: -1 )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .set { ch_input_bam  }

ch_spp_correlation_header = file("$baseDir/assets/multiqc/spp_correlation_header.txt", checkIfExists: true)
ch_spp_nsc_header = file("$baseDir/assets/multiqc/spp_nsc_header.txt", checkIfExists: true)
ch_spp_rsc_header = file("$baseDir/assets/multiqc/spp_rsc_header.txt", checkIfExists: true)
params.skip_spp = false

process PhantomPeakQualTools {
    tag "$name"
    publishDir path: "${params.outdir}", mode: 'copy'

    when:
    !params.skip_spp

    input:
    set val(name), file(bam) from ch_input_bam
    file spp_correlation_header from ch_spp_correlation_header
    file spp_nsc_header from ch_spp_nsc_header
    file spp_rsc_header from ch_spp_rsc_header

    output:
    file '*.pdf' into ch_spp_plot
    file '*.spp.out' into ch_spp_out,
                          ch_spp_out_mqc
    file '*_mqc.tsv' into ch_spp_csv_mqc

    script:
    """
    RUN_SPP=`which run_spp.R`
    Rscript -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="${bam[0]}" -savp="${name}.spp.pdf" -savd="${name}.spp.Rdata" -out="${name}.spp.out" -p=$task.cpus
    cp $spp_correlation_header ${name}_spp_correlation_mqc.tsv
    Rscript -e "load('${name}.spp.Rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${name}_spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"

    awk -v OFS='\\t' '{print "${name}", \$9}' ${name}.spp.out | cat $spp_nsc_header - > ${name}_spp_nsc_mqc.tsv
    awk -v OFS='\\t' '{print "${name}", \$10}' ${name}.spp.out | cat $spp_rsc_header - > ${name}_spp_rsc_mqc.tsv
    """
}

