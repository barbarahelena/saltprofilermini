/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_saltprofiler'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { SALTGENES                       } from '../subworkflows/local/saltgenes'

//
// MODULE: Installed directly from nf-core/modules
//
include { PRODIGAL                        } from '../modules/nf-core/prodigal/main'
include { PROKKA                          } from '../modules/nf-core/prokka/main'


////////////////////////////////////////////////////
/* --  Create channel for reference databases  -- */
////////////////////////////////////////////////////

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SALTPROFILERMINI {

    take:
    input_assemblies
    input_genes

    main:

    ch_versions = Channel.empty()

    /*
    ================================================================================
                                    Predict proteins
    ================================================================================
    */

    if (!params.skip_prodigal){
        PRODIGAL (
            input_assemblies,
            'gff'
        )
        ch_versions = ch_versions.mix(PRODIGAL.out.versions.first())
    }
            
    /*
        * Prokka: Genome annotation
    */

    if (!params.skip_prokka){
        PROKKA (
            input_assemblies,
            [],
            []
        )
        ch_versions = ch_versions.mix(PROKKA.out.versions.first())

        /*
        * Overview of salt tolerance genes
        */

        if ( !params.skip_saltgenes ) {
            ch_prokka_output = PROKKA.out.gff.combine(PROKKA.out.fna, by: 0)
            SALTGENES(
                input_genes,
                ch_prokka_output
            )
            ch_versions = ch_versions.mix(SALTGENES.out.versions.first())

        }
    }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/