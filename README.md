[![Build Status](https://travis-ci.com/broadinstitute/viral-pipelines.svg?branch=master)](https://travis-ci.com/broadinstitute/viral-pipelines)
[![Documentation Status](https://readthedocs.org/projects/viral-pipelines/badge/?version=latest)](http://viral-pipelines.readthedocs.io/en/latest/?badge=latest)

# viral-pipelines

A set of scripts and tools for the analysis of viral NGS data.

Workflows are written in [WDL](https://github.com/openwdl/wdl) format. This is a portable workflow language that allows for easy execution on a wide variety of platforms:
 - on individual machines (using [miniWDL](https://github.com/chanzuckerberg/miniwdl) or [Cromwell](https://github.com/broadinstitute/cromwell) to execute)
 - on commercial cloud platforms like GCP, AWS, or Azure (using [Cromwell](https://github.com/broadinstitute/cromwell) or [CromwellOnAzure](https://github.com/microsoft/CromwellOnAzure))
 - on institutional HPC systems (using [Cromwell](https://github.com/broadinstitute/cromwell))
 - on commercial platform as a service vendors (like [DNAnexus](https://dnanexus.com/))
 - on academic cloud platforms (like [Terra](https://app.terra.bio/))


## Obtaining the latest WDL workflows

Workflows from this repository are continuously deployed to [Dockstore](https://dockstore.org/organizations/BroadInstitute/collections/pgs), a GA4GH Tool Registry Service. They can then be easily imported to any bioinformatic compute platform that utilizes the TRS API and understands WDL (this includes Terra, DNAnexus, DNAstack, etc).

Workflows are also available in the [Terra featured workspace](https://app.terra.bio/#workspaces/pathogen-genomic-surveillance/COVID-19).

Workflows are continuously deployed to a [DNAnexus CI project](https://platform.dnanexus.com/projects/F8PQ6380xf5bK0Qk0YPjB17P).


## Basic execution

<img src="https://raw.githubusercontent.com/openwdl/learn-wdl/master/images/miniwdl-dev.png" width=600>

The easiest way to get started is on a single, Python & Docker-capable machine (your laptop, shared workstation, or virtual machine) using [miniWDL](https://github.com/chanzuckerberg/miniwdl) as shown above. MiniWDL can be installed either via `pip` or `conda` (via conda-forge). After confirming that it works (`miniwdl run_self_test`, you can use [miniwdl run](https://miniwdl.readthedocs.io/en/latest/getting_started.html) to invoke WDL workflows from this repository.

For example, to list the inputs for the assemble_refbased workflow:

```
miniwdl run https://raw.githubusercontent.com/broadinstitute/viral-pipelines/v2.1.8.0/pipes/WDL/workflows/assemble_refbased.wdl
```

This will emit:
```
missing required inputs for assemble_refbased: reads_unmapped_bams, reference_fasta

required inputs:
  Array[File]+ reads_unmapped_bams
  File reference_fasta

optional inputs:
  <really long list>

outputs:
  <really long list>
```

To then execute this workflow on your local machine, invoke it with like this:
```
miniwdl run \
  https://raw.githubusercontent.com/broadinstitute/viral-ngs-staging/master/pipes/WDL/workflows/assemble_refbased.wdl \
  reads_unmapped_bams=PatientA_library1.bam \
  reads_unmapped_bams=PatientA_library2.bam \
  reference_fasta=/refs/NC_045512.2.fasta \
  trim_coords_bed=/refs/NC_045512.2-artic_primers-3.bed \
  sample_name=PatientA
```

In the above example, reads from two sequencing runs are aligned and merged together before consensus calling. The optional bed file provided turns on primer trimming at the given coordinates.


## Available workflows

The workflows provided here are more fully documented at our [ReadTheDocs](https://viral-pipelines.readthedocs.io/en/latest/workflows.html) page.
