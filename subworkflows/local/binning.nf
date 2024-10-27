/*
 * Binning with MetaBAT2
 */

include { METABAT2_METABAT2                                            } from '../../modules/nf-core/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS                         } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { GUNZIP as GUNZIP_BINS                                        } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_UNBINS                                      } from '../../modules/nf-core/gunzip/main'

include { CONVERT_DEPTHS                        } from '../../modules/local/convert_depths'
include { SPLIT_FASTA                           } from '../../modules/local/split_fasta'

workflow BINNING {
    take:
    assemblies           // channel: [ val(meta), path(assembly), path(bams), path(bais) ]
    reads                // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()

    // generate coverage depths for each contig
    ch_summarizedepth_input = assemblies
                                .map { meta, assembly, bams, bais ->
                                    [ meta, bams, bais ]
                                }

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_summarizedepth_input )

    ch_metabat_depths = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .map { meta, depths ->
            def meta_new = meta + [binner: 'MetaBAT2']
            [ meta_new, depths ]
        }

    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first())

    // combine depths back with assemblies
    ch_metabat2_input = assemblies
        .map { meta, assembly, bams, bais ->
            def meta_new = meta + [binner: 'MetaBAT2']
            [ meta_new, assembly, bams, bais ]
        }
        .join( ch_metabat_depths, by: 0 )
        .map { meta, assembly, bams, bais, depths ->
            [ meta, assembly, depths ]
        }

    // main bins for decompressing for MAG_DEPTHS
    ch_final_bins_for_gunzip = Channel.empty()

    // final gzipped bins
    ch_binning_results_gzipped_final = Channel.empty()

    // run binning
    if ( !params.skip_metabat2 ) {
        METABAT2_METABAT2 ( ch_metabat2_input )
        // before decompressing first have to separate and re-group due to limitation of GUNZIP module
        ch_final_bins_for_gunzip = ch_final_bins_for_gunzip.mix( METABAT2_METABAT2.out.fasta.transpose() )
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix( METABAT2_METABAT2.out.fasta )
        ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions.first())
    }

    // decide which unbinned fasta files to further filter, depending on which binners selected
    // NOTE: CONCOCT does not produce 'unbins' itself, therefore not included here.
    if ( !params.skip_metabat2) {
        ch_input_splitfasta = METABAT2_METABAT2.out.unbinned
    } else {
        ch_input_splitfasta = Channel.empty()
    }

    SPLIT_FASTA ( ch_input_splitfasta )
    // large unbinned contigs from SPLIT_FASTA for decompressing for MAG_DEPTHS,
    // first have to separate and re-group due to limitation of GUNZIP module
    ch_split_fasta_results_transposed = SPLIT_FASTA.out.unbinned.transpose()
    ch_versions = ch_versions.mix(SPLIT_FASTA.out.versions)

    GUNZIP_BINS ( ch_final_bins_for_gunzip )
    ch_binning_results_gunzipped = GUNZIP_BINS.out.gunzip
        .groupTuple(by: 0)

    GUNZIP_UNBINS ( ch_split_fasta_results_transposed )
    ch_splitfasta_results_gunzipped = GUNZIP_UNBINS.out.gunzip
        .groupTuple(by: 0)

    ch_versions = ch_versions.mix(GUNZIP_BINS.out.versions.first())
    ch_versions = ch_versions.mix(GUNZIP_UNBINS.out.versions.first())

    emit:
    bins                                         = ch_binning_results_gunzipped
    bins_gz                                      = ch_binning_results_gzipped_final
    unbinned                                     = ch_splitfasta_results_gunzipped
    unbinned_gz                                  = SPLIT_FASTA.out.unbinned
    metabat2depths                               = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    versions                                     = ch_versions
}
