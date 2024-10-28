/*
 * Salt gene profiler
 */


include { SALTGENES_FILTER          } from '../../modules/local/saltgenes/filter'
include { SALTGENES_CATPERGENE      } from '../../modules/local/saltgenes/catpergene'
include { SALTGENES_CATPERSAMPLE    } from '../../modules/local/saltgenes/catpersample'
include { SALTGENES_MAP             } from '../../modules/local/saltgenes/map'
include { SALTGENES_BOWTIE2BUILD    } from '../../modules/local/saltgenes/bowtie2build'
include { SALTGENES_BOWTIE2ALIGN    } from '../../modules/local/saltgenes/bowtie2align'
include { SALTGENES_COUNT           } from '../../modules/local/saltgenes/count'

workflow SALTGENES {

    take:
    genes
    prokka_output
    reads

    main:

    ch_versions = Channel.empty()

    ch_saltgenes = prokka_output.combine(genes)

    SALTGENES_FILTER ( 
        ch_saltgenes
    )
    ch_versions = ch_versions.mix(SALTGENES_FILTER.out.versions.first())
   
    // Concatenate the seqs and gffs per gene per sample
    ch_fapergene = SALTGENES_FILTER.out.seqs
                            .map { metadata, gene, fasta -> 
                                        def lastDashIndex = metadata.id.lastIndexOf('-')
                                        def dotIndex = metadata.id.indexOf('.', lastDashIndex)
                                        def subject_id = metadata.id.substring(lastDashIndex + 1, dotIndex)
                                        tuple( subject_id, metadata.id, gene, fasta )
                                    }
                            .groupTuple(by: [0,2])
    ch_gffpergene = SALTGENES_FILTER.out.gff
                            .map { metadata, gene, gff -> 
                                        def lastDashIndex = metadata.id.lastIndexOf('-')
                                        def dotIndex = metadata.id.indexOf('.', lastDashIndex)
                                        def subject_id = metadata.id.substring(lastDashIndex + 1, dotIndex)
                                        tuple( subject_id, metadata.id, gene, gff )
                                    }
                            .groupTuple(by: [0,2])
    SALTGENES_CATPERGENE( ch_fapergene, ch_gffpergene )
    ch_versions = ch_versions.mix(SALTGENES_CAT.out.versions.first())

    // Concatenate the seqs and gffs per sample
    ch_fapersample = ch_fapergene.groupTuple(by: 0)
    ch_gffpersample = ch_gffpergene.groupTuple(by: 0)
    SALTGENES_CATPERSAMPLE( ch_fapersample, ch_gffpersample )

    // Mapping with Bowtie2: now per sample because seemed more efficient
    SALTGENES_BOWTIE2BUILD( SALTGENES_CATPERSAMPLE.out.mergedfa )
    ch_versions = ch_versions.mix(SALTGENES_BOWTIE2BUILD.out.versions.first())

    // Combine the index channel, reads channel and gff channel (annotation)
    ch_reads = reads.map { metadata, reads -> tuple( metadata.id, reads ) }
    ch_indexgffreads = SALTGENES_BOWTIE2BUILD.out.index
                            .combine(ch_reads, by: 0)
                            .combine(SALTGENES_CATPERSAMPLE.out.mergedgff, by: 0)

    SALTGENES_BOWTIE2ALIGN( ch_indexgffreads )
    ch_versions = ch_versions.mix(SALTGENES_BOWTIE2ALIGN.out.versions.first())

    SALTGENES_COUNT(SALTGENES_BOWTIE2ALIGN.out.bam, ch_indexgffreads )
    ch_versions = ch_versions.mix(SALTGENES_COUNT.out.versions.first())

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}
