[![Build Status](https://travis-ci.com/broadinstitute/viral-pipelines.svg?branch=master)](https://travis-ci.com/broadinstitute/viral-pipelines)
[![Documentation Status](https://readthedocs.org/projects/viral-pipelines/badge/?version=latest)](http://viral-pipelines.readthedocs.io/en/latest/?badge=latest)

viral-pipelines
===============

A set of scripts and tools for the analysis of viral NGS data.

Workflows are written in [WDL](https://github.com/openwdl/wdl) format. This is a portable workflow language that allows for easy execution on a wide variety of platforms:
 - on individual machines (using miniWDL or Cromwell to execute)
 - on commercial cloud platforms like GCP, AWS, or Azure (using Cromwell or CromwellOnAzure)
 - on institutional HPC systems (using Cromwell)
 - on commercial platform as a service vendors (like DNAnexus)
 - on academic cloud platforms (like Terra)

Currently, all workflows are regularly deployed to a GCS bucket: [gs://viral-ngs-wdl](https://console.cloud.google.com/storage/browser/viral-ngs-wdl?forceOnBucketsSortingFiltering=false&organizationId=548622027621&project=gcid-viral-seq). 

Workflows are also available in the [Terra featured workspace](https://app.terra.bio/#workspaces/pathogen-genomic-surveillance/COVID-19).

Workflows are continuously deployed to a [DNAnexus CI project](https://platform.dnanexus.com/projects/F8PQ6380xf5bK0Qk0YPjB17P).

Continuous deploy to [Dockstore](https://dockstore.org/) is pending.

Basic execution
---------------

The easiest way to get started is on a single, Docker-capable machine (your laptop, shared workstation, or virtual machine) using [miniWDL](https://github.com/chanzuckerberg/miniwdl). MiniWDL can be installed via `pip` or `conda` (via conda-forge). After confirming that it works (`miniwdl run_self_test`, you can use [miniwdl run](https://github.com/chanzuckerberg/miniwdl#miniwdl-run) to invoke WDL workflows from this repository. For example: `miniwdl run https://storage.googleapis.com/viral-ngs-wdl/quay.io/broadinstitute/viral-pipelines/2.0.20.3/assemble_refbased.wdl` will execute the reference-based assembly pipeline, when provided with the appropriate inputs.

Available workflows
-------------------

The workflows provided here are more fully documented at our [ReadTheDocs](https://viral-pipelines.readthedocs.io/) page.
