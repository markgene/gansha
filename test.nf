mito_filter = params.spikein_mito_name ? "${params.spikein_mito_name}\\|${params.mito_name}" : params.mito_name
println "${mito_filter}"

if (params.filter_save_int_bam) {
    println "save intermediate bam"
} else {
    println "Do not save intermediate bam"
}
