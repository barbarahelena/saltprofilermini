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
include { SALTGENES                } from '../subworkflows/local/saltgenes'

//
// MODULE: Installed directly from nf-core/modules
//
include { PRODIGAL                        } from '../modules/nf-core/prodigal/main'
include { BAKTA_BAKTA                     } from '../modules/nf-core/bakta/bakta/main'
include { BAKTA_BAKTADBDOWNLOAD           } from '../modules/nf-core/bakta/baktadbdownload/main'

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
    bakta_db = params.bakta_database ? Channel.fromPath( params.bakta_database ).first() : []

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
        * Bakta: Genome annotation
    */

    if (!params.skip_bakta){
        if ( ! bakta_db ){
            BAKTA_BAKTADBDOWNLOAD()
            bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db
        }         
        BAKTA_BAKTA( 
            input_assemblies, 
            bakta_db,
            [],
            []
        )
        ch_annotation = BAKTA_BAKTA.out.fna.join(BAKTA_BAKTA.out.gff)
        ch_versions = ch_versions.mix( BAKTA_BAKTA.out.versions.first() )

        /*
        * Overview of salt tolerance genes
        */

        if ( !params.skip_saltgenes ) {

            SALTGENES(input_genes, ch_annotation)

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