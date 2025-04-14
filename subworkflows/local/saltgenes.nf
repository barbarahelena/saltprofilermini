/*
 * Salt gene profiler
 */

include { SALTGENES_FILTER                  } from '../../modules/local/saltgenes/filter'
include { SALTGENES_CATPERASSEMBLY          } from '../../modules/local/saltgenes/catperassembly'
include { SALTGENES_MERGE                   } from '../../modules/local/saltgenes/merge'

workflow SALTGENES {

    take:
    genes
    bakta_output

    main:

    ch_versions = Channel.empty()
    ch_saltgenes = bakta_output.combine(genes)

    // Get fastas
    SALTGENES_FILTER ( ch_saltgenes )
    ch_versions = ch_versions.mix(SALTGENES_FILTER.out.versions.first())

    // Group and concatenate the fasta/gff per gene per sample
    ch_seqs = SALTGENES_FILTER.out.seqs.groupTuple(by: 0)
    SALTGENES_CATPERASSEMBLY( ch_seqs )
    ch_versions = ch_versions.mix(SALTGENES_CATPERASSEMBLY.out.versions.first())

    ch_alltab = SALTGENES_CATPERASSEMBLY.out.merged.map{ _meta, _gene, csv_path -> csv_path }.collect()
    SALTGENES_MERGE( ch_alltab )
    ch_versions = ch_versions.mix(SALTGENES_MERGE.out.versions.first())

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}
