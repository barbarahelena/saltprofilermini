/*
 * Salt gene profiler
 */

include { SALTGENES_FILTER                  } from '../../modules/local/saltgenes/filter'
include { SALTGENES_CATPERGENE              } from '../../modules/local/saltgenes/catpergene'
include { SALTGENES_CATPERSAMPLE            } from '../../modules/local/saltgenes/catpersample'
include { SALTGENES_MUSCLE                  } from '../../modules/local/saltgenes/muscle'
include { SALTGENES_FASTTREE                } from '../../modules/local/saltgenes/fasttree'
include { SALTGENES_BWAMEM2BUILD            } from '../../modules/local/saltgenes/bwamem2build'
include { SALTGENES_BWAMEM2ALIGN            } from '../../modules/local/saltgenes/bwamem2align'
include { SALTGENES_ADDTAX                  } from '../../modules/local/saltgenes/addtax'
include { SALTGENES_PLOTGENES               } from '../../modules/local/saltgenes/plotgenes'
include { COMBINE_TSV as COMBINE_SALTGENES  } from '../modules/local/combine_tsv'
include { MERGE_ALLTAB                      } from '../modules/local/merge_alltab'

workflow SALTGENES {

    take:
    genes
    prokka_output
    reads
    tax

    main:

    ch_versions = Channel.empty()

    ch_saltgenes = prokka_output.combine(genes)

    // Get fastas
    SALTGENES_FILTER ( ch_saltgenes )
    ch_versions = ch_versions.mix(SALTGENES_FILTER.out.versions.first())

    // Group and concatenate the fasta/gff per gene per sample
    ch_seqs = SALTGENES_FILTER.out.seqs
                    .map { metadata, gene, fasta, gff -> 
                                        def lastDashIndex = metadata.id.lastIndexOf('-')
                                        def dotIndex = metadata.id.indexOf('.', lastDashIndex)
                                        def subject_id = metadata.id.substring(lastDashIndex + 1, dotIndex)
                                        tuple( subject_id, metadata.id, gene, fasta, gff )
                                    }
                            .groupTuple(by: [0,2])
    SALTGENES_CATPERGENE( ch_seqs )
    ch_versions = ch_versions.mix(SALTGENES_CATPERGENE.out.versions.first())

    // MSA and tree of sequences per gene
    // ch_fastapergene = SALTGENES_CATPERGENE.out.mergedseqs 
    //                         .map { metadata, gene, fasta, gff -> 
    //                                 def fastalist = fasta.collect()
    //                                 def gfflist = gff.collect()
    //                                 tuple( metadata, gene, fasta, gff )
    //                                 }
    
    // SALTGENES_MUSCLE( ch_fastapergene )
    // SALTGENES_FASTTREE( SALTGENES_MUSCLE.out.msa )

    // Concatenate the fasta/gffs per sample
    ch_genespersample = SALTGENES_CATPERGENE.out.mergedseqs.groupTuple(by: 0)
    SALTGENES_CATPERSAMPLE( ch_genespersample )

    // Mapping with BWA-mem2 per sample
    SALTGENES_BWAMEM2BUILD( SALTGENES_CATPERSAMPLE.out.mergedgenes )
    ch_versions = ch_versions.mix(SALTGENES_BWAMEM2BUILD.out.versions.first())

    // Combine the index channel, reads channel and gff channel (annotation)
    ch_reads = reads.map { metadata, reads -> tuple( metadata.id, reads ) }
    ch_indexgffreads = SALTGENES_BWAMEM2BUILD.out.index
                            .combine(ch_reads, by: 0)

    SALTGENES_BWAMEM2ALIGN( ch_indexgffreads )
    ch_versions = ch_versions.mix(SALTGENES_BWAMEM2ALIGN.out.versions.first())

    // Merge saltgene count tables
    COMBINE_SALTGENES ( SALTGENES.out.counts.map{ it[1] }.collect(), "saltgenes_summary" )

    // Merge saltgene table with depth-and-tax table
    MERGE_ALLTAB ( COMBINE_SALTGENES.out.combined, tax )

    // MSA and tree per gene?

    emit:
    counts = SALTGENES_BWAMEM2ALIGN.out.stats
    versions = ch_versions                     // channel: [ versions.yml ]
}
