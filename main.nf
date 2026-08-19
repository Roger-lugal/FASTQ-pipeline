#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samplesheet = null
params.input = null
params.selected = null
params.output = 'results'
params.variant_start = 0
params.variant_end = null
params.min_length = 1
params.publish_mode = 'copy'

if (!params.samplesheet && (!params.input || !params.selected)) {
    exit 1, "Please provide either --samplesheet OR both --input and --selected."
}

Channel
    .fromPath(params.samplesheet ?: 'NOFILE', checkIfExists: params.samplesheet ? true : false)
    .ifEmpty {
        if (params.samplesheet) {
            exit 1, "Samplesheet not found: ${params.samplesheet}"
        }
    }
    .map { file ->
        def rows = file.readLines()
        rows = rows.findAll { it && !it.trim().startsWith('#') }
        if (!rows) {
            exit 1, "Samplesheet is empty: ${file}"
        }
        def header = rows[0].split(',')
        def idx = [:]
        header.eachWithIndex { name, i -> idx[name.trim()] = i }
        def data = rows.tail()
        return data.collect { line ->
            def values = line.split(',')
            def sample_id = values[idx['sample_id']].trim()
            def input_fastq = values[idx['input']].trim()
            def selected_fastq = values[idx['selected']].trim()
            tuple(sample_id, file(input_fastq), file(selected_fastq))
        }
    }
    .flatMap { it -> it }
    .set { samples }

if (!params.samplesheet) {
    Channel
        .of(tuple('sample', file(params.input), file(params.selected)))
        .set { samples }
}

process run_dms {
    tag { sample_id }
    publishDir "${params.output}", mode: params.publish_mode
    cpus 2
    memory '4 GB'
    input:
        tuple val(sample_id), path(input_fastq), path(selected_fastq)
    output:
        path "${sample_id}.dms.tsv", emit: dms_table
    script:
        def end_arg = params.variant_end ? "--variant-end ${params.variant_end}" : ""
        def out_name = "${sample_id}.dms.tsv"
        """
        python3 ${projectDir}/dms_analysis.py \
            --input ${input_fastq} \
            --selected ${selected_fastq} \
            --variant-start ${params.variant_start} \
            ${end_arg} \
            --min-length ${params.min_length} \
            --output ${out_name}
        """
}

workflow {
    run_dms(samples)
}
