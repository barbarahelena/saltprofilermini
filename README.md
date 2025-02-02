<h1>
  <!-- <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/mag_logo_mascot_dark.png">
    <img alt="nf-core/mag" src="docs/images/mag_logo_mascot_light.png">
  </picture>
  -->
</h1>

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**saltprofilermini** is a very short bioinformatics pipeline for the annotation of metagenomes focused on salt tolerance genes.

<!-- <p align="center">
    <img src="docs/images/mag_workflow.png" alt="nf-core/mag workflow overview" width="90%">
</p> -->

## Pipeline summary

This pipeline uses assemblies as input. The pipeline then:
- predicts protein-coding genes for the assemblies using [Prodigal](https://github.com/hyattpd/Prodigal), and bins with [Prokka](https://github.com/tseemann/prokka)
- filters genes for genes of interest

Furthermore, the pipeline creates various reports in the results directory specified, including a [MultiQC](https://multiqc.info/) report summarizing some of the findings and software versions.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

For more details and further functionality, please refer to the [usage documentation] and the [parameter documentation].

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/mag/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation]().

## Credits

The salt tolerance-specific modules were written by Barbara, and the prokka module is part of the nf-core repository.

## Citations

If you use this pipeline, I advice to cite Prokka, since Prokka is crucial for the downstream steps in this pipeline.
