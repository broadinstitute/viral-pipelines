[![Build Status](https://travis-ci.com/broadinstitute/viral-pipelines.svg?branch=master)](https://travis-ci.com/broadinstitute/viral-pipelines)
[![Documentation Status](https://readthedocs.org/projects/viral-pipelines/badge/?version=latest)](http://viral-pipelines.readthedocs.io/en/latest/?badge=latest)

viral-pipelines
===============

A set of scripts and tools for the analysis of viral NGS data.

Workflows are written in [WDL](https://github.com/openwdl/wdl) format, currently in a mix of WDL "draft-2" and "1.0". This is a portable workflow language that allows for easy execution on a wide variety of platforms:
 - on individual machines (using miniWDL or Cromwell to execute)
 - on commercial cloud platforms like GCP, AWS, or Azure (using Cromwell or CromwellOnAzure)
 - on institutional HPC systems (using Cromwell)
 - on commercial platform as a service vendors (like DNAnexus)
 - on academic cloud platforms (like Terra)

(TO DO: describe how these are automatically pushed to dockstore, etc)


Basic execution
---------------

The easiest way to get started is on a single, Docker-capable machine (your laptop, shared workstation, or virtual machine) using [miniWDL](https://github.com/chanzuckerberg/miniwdl). MiniWDL can be installed via `pip` or `conda` (via conda-forge). After confirming that it works, you can use [miniwdl run](https://github.com/chanzuckerberg/miniwdl#miniwdl-run) to invoke WDL workflows from this repository.


Available workflows
-------------------

(TO DO)
