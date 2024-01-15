Using the WDL pipelines
=======================

Rather than chaining together viral-ngs pipeline steps as series of tool commands called in isolation, it is possible to execute them as a complete automated pipeline, from processing raw sequencer output to creating files suitable for GenBank submission. This utilizes the Workflow Description Language, which is documented at:
https://github.com/openwdl/wdl

There are various methods for executing these workflows on your infrastructure which are more thoroughly documented in our `README <https://github.com/broadinstitute/viral-pipelines/blob/master/README.md#viral-pipelines>`_.
