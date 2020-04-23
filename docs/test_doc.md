# assemble_refbased
Reference-based microbial consensus calling. Aligns short reads to a singular reference genome, calls a new consensus sequence, and emits: new assembly, reads aligned to provided reference, reads aligned to new assembly, various figures of merit, plots, and QC metrics. The user may provide unaligned reads spread across multiple input files and this workflow will parallelize alignment per input file before merging results prior to consensus calling.

## Inputs


### Required inputs
<p name="assemble_refbased.reads_unmapped_bams">
        <b>assemble_refbased.reads_unmapped_bams</b><br />
        <i>Array[File]+ &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.reference_fasta">
        <b>assemble_refbased.reference_fasta</b><br />
        <i>File &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.sample_name">
        <b>assemble_refbased.sample_name</b><br />
        <i>String &mdash; Default: None</i><br />
        ???
</p>



### Other inputs
<details>
<summary> Show/Hide </summary>
<p name="assemble_refbased.align_to_ref.aligner">
        <b>assemble_refbased.align_to_ref.aligner</b><br />
        <i>String? &mdash; Default: "novoalign"</i><br />
        ???
</p>
<p name="assemble_refbased.align_to_ref.docker">
        <b>assemble_refbased.align_to_ref.docker</b><br />
        <i>String? &mdash; Default: "quay.io/broadinstitute/viral-core"</i><br />
        ???
</p>
<p name="assemble_refbased.align_to_ref.machine_mem_gb">
        <b>assemble_refbased.align_to_ref.machine_mem_gb</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.align_to_ref.sample_name">
        <b>assemble_refbased.align_to_ref.sample_name</b><br />
        <i>String &mdash; Default: basename(basename(basename(reads_unmapped_bam,".bam"),".taxfilt"),".clean")</i><br />
        ???
</p>
<p name="assemble_refbased.align_to_self.aligner">
        <b>assemble_refbased.align_to_self.aligner</b><br />
        <i>String? &mdash; Default: "novoalign"</i><br />
        ???
</p>
<p name="assemble_refbased.align_to_self.docker">
        <b>assemble_refbased.align_to_self.docker</b><br />
        <i>String? &mdash; Default: "quay.io/broadinstitute/viral-core"</i><br />
        ???
</p>
<p name="assemble_refbased.align_to_self.machine_mem_gb">
        <b>assemble_refbased.align_to_self.machine_mem_gb</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.align_to_self.sample_name">
        <b>assemble_refbased.align_to_self.sample_name</b><br />
        <i>String &mdash; Default: basename(basename(basename(reads_unmapped_bam,".bam"),".taxfilt"),".clean")</i><br />
        ???
</p>
<p name="assemble_refbased.call_consensus.docker">
        <b>assemble_refbased.call_consensus.docker</b><br />
        <i>String? &mdash; Default: "quay.io/broadinstitute/viral-assemble"</i><br />
        ???
</p>
<p name="assemble_refbased.call_consensus.machine_mem_gb">
        <b>assemble_refbased.call_consensus.machine_mem_gb</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.call_consensus.major_cutoff">
        <b>assemble_refbased.call_consensus.major_cutoff</b><br />
        <i>Float? &mdash; Default: 0.5</i><br />
        ???
</p>
<p name="assemble_refbased.call_consensus.mark_duplicates">
        <b>assemble_refbased.call_consensus.mark_duplicates</b><br />
        <i>Boolean? &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.call_consensus.min_coverage">
        <b>assemble_refbased.call_consensus.min_coverage</b><br />
        <i>Int? &mdash; Default: 2</i><br />
        ???
</p>
<p name="assemble_refbased.ivar_trim.bam_basename">
        <b>assemble_refbased.ivar_trim.bam_basename</b><br />
        <i>String &mdash; Default: basename(aligned_bam,".bam")</i><br />
        ???
</p>
<p name="assemble_refbased.ivar_trim.docker">
        <b>assemble_refbased.ivar_trim.docker</b><br />
        <i>String? &mdash; Default: "andersenlabapps/ivar:1.2.1"</i><br />
        ???
</p>
<p name="assemble_refbased.ivar_trim.machine_mem_gb">
        <b>assemble_refbased.ivar_trim.machine_mem_gb</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.ivar_trim.min_keep_length">
        <b>assemble_refbased.ivar_trim.min_keep_length</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.ivar_trim.min_quality">
        <b>assemble_refbased.ivar_trim.min_quality</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.ivar_trim.sliding_window">
        <b>assemble_refbased.ivar_trim.sliding_window</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.merge_align_to_ref.docker">
        <b>assemble_refbased.merge_align_to_ref.docker</b><br />
        <i>String? &mdash; Default: "quay.io/broadinstitute/viral-core"</i><br />
        ???
</p>
<p name="assemble_refbased.merge_align_to_ref.machine_mem_gb">
        <b>assemble_refbased.merge_align_to_ref.machine_mem_gb</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.merge_align_to_ref.reheader_table">
        <b>assemble_refbased.merge_align_to_ref.reheader_table</b><br />
        <i>File? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.merge_align_to_self.docker">
        <b>assemble_refbased.merge_align_to_self.docker</b><br />
        <i>String? &mdash; Default: "quay.io/broadinstitute/viral-core"</i><br />
        ???
</p>
<p name="assemble_refbased.merge_align_to_self.machine_mem_gb">
        <b>assemble_refbased.merge_align_to_self.machine_mem_gb</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.merge_align_to_self.reheader_table">
        <b>assemble_refbased.merge_align_to_self.reheader_table</b><br />
        <i>File? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.comment">
        <b>assemble_refbased.multiqc_align_to_ref.comment</b><br />
        <i>String? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.config">
        <b>assemble_refbased.multiqc_align_to_ref.config</b><br />
        <i>File? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.config_yaml">
        <b>assemble_refbased.multiqc_align_to_ref.config_yaml</b><br />
        <i>String? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.data_dir">
        <b>assemble_refbased.multiqc_align_to_ref.data_dir</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.dirs">
        <b>assemble_refbased.multiqc_align_to_ref.dirs</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.dirs_depth">
        <b>assemble_refbased.multiqc_align_to_ref.dirs_depth</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.docker">
        <b>assemble_refbased.multiqc_align_to_ref.docker</b><br />
        <i>String &mdash; Default: "ewels/multiqc:latest"</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.exclude_modules">
        <b>assemble_refbased.multiqc_align_to_ref.exclude_modules</b><br />
        <i>Array[String]+? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.export">
        <b>assemble_refbased.multiqc_align_to_ref.export</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.file_name">
        <b>assemble_refbased.multiqc_align_to_ref.file_name</b><br />
        <i>String? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.file_with_list_of_input_paths">
        <b>assemble_refbased.multiqc_align_to_ref.file_with_list_of_input_paths</b><br />
        <i>File? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.flat">
        <b>assemble_refbased.multiqc_align_to_ref.flat</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.force">
        <b>assemble_refbased.multiqc_align_to_ref.force</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.full_names">
        <b>assemble_refbased.multiqc_align_to_ref.full_names</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.ignore_analysis_files">
        <b>assemble_refbased.multiqc_align_to_ref.ignore_analysis_files</b><br />
        <i>String? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.ignore_sample_names">
        <b>assemble_refbased.multiqc_align_to_ref.ignore_sample_names</b><br />
        <i>String? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.input_directory">
        <b>assemble_refbased.multiqc_align_to_ref.input_directory</b><br />
        <i>String &mdash; Default: "multiqc-input"</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.interactive">
        <b>assemble_refbased.multiqc_align_to_ref.interactive</b><br />
        <i>Boolean &mdash; Default: true</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.lint">
        <b>assemble_refbased.multiqc_align_to_ref.lint</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.machine_mem_gb">
        <b>assemble_refbased.multiqc_align_to_ref.machine_mem_gb</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.megaQC_upload">
        <b>assemble_refbased.multiqc_align_to_ref.megaQC_upload</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.module_to_use">
        <b>assemble_refbased.multiqc_align_to_ref.module_to_use</b><br />
        <i>Array[String]+? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.no_data_dir">
        <b>assemble_refbased.multiqc_align_to_ref.no_data_dir</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.out_dir">
        <b>assemble_refbased.multiqc_align_to_ref.out_dir</b><br />
        <i>String &mdash; Default: "./multiqc-output"</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.output_data_format">
        <b>assemble_refbased.multiqc_align_to_ref.output_data_format</b><br />
        <i>String? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.pdf">
        <b>assemble_refbased.multiqc_align_to_ref.pdf</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.report_filename">
        <b>assemble_refbased.multiqc_align_to_ref.report_filename</b><br />
        <i>String &mdash; Default: if defined(file_name) then basename(select_first([file_name]),".html") else "multiqc"</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.sample_names">
        <b>assemble_refbased.multiqc_align_to_ref.sample_names</b><br />
        <i>File? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.tag">
        <b>assemble_refbased.multiqc_align_to_ref.tag</b><br />
        <i>String? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.template">
        <b>assemble_refbased.multiqc_align_to_ref.template</b><br />
        <i>String? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.title">
        <b>assemble_refbased.multiqc_align_to_ref.title</b><br />
        <i>String? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.multiqc_align_to_ref.zip_data_dir">
        <b>assemble_refbased.multiqc_align_to_ref.zip_data_dir</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.novocraft_license">
        <b>assemble_refbased.novocraft_license</b><br />
        <i>File? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.plot_ref_coverage.bin_large_plots">
        <b>assemble_refbased.plot_ref_coverage.bin_large_plots</b><br />
        <i>Boolean? &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.plot_ref_coverage.binning_summary_statistic">
        <b>assemble_refbased.plot_ref_coverage.binning_summary_statistic</b><br />
        <i>String? &mdash; Default: "max"</i><br />
        ???
</p>
<p name="assemble_refbased.plot_ref_coverage.docker">
        <b>assemble_refbased.plot_ref_coverage.docker</b><br />
        <i>String? &mdash; Default: "quay.io/broadinstitute/viral-core"</i><br />
        ???
</p>
<p name="assemble_refbased.plot_ref_coverage.machine_mem_gb">
        <b>assemble_refbased.plot_ref_coverage.machine_mem_gb</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.plot_ref_coverage.plot_only_non_duplicates">
        <b>assemble_refbased.plot_ref_coverage.plot_only_non_duplicates</b><br />
        <i>Boolean? &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.plot_ref_coverage.skip_mark_dupes">
        <b>assemble_refbased.plot_ref_coverage.skip_mark_dupes</b><br />
        <i>Boolean? &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.plot_self_coverage.bin_large_plots">
        <b>assemble_refbased.plot_self_coverage.bin_large_plots</b><br />
        <i>Boolean? &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.plot_self_coverage.binning_summary_statistic">
        <b>assemble_refbased.plot_self_coverage.binning_summary_statistic</b><br />
        <i>String? &mdash; Default: "max"</i><br />
        ???
</p>
<p name="assemble_refbased.plot_self_coverage.docker">
        <b>assemble_refbased.plot_self_coverage.docker</b><br />
        <i>String? &mdash; Default: "quay.io/broadinstitute/viral-core"</i><br />
        ???
</p>
<p name="assemble_refbased.plot_self_coverage.machine_mem_gb">
        <b>assemble_refbased.plot_self_coverage.machine_mem_gb</b><br />
        <i>Int? &mdash; Default: None</i><br />
        ???
</p>
<p name="assemble_refbased.plot_self_coverage.plot_only_non_duplicates">
        <b>assemble_refbased.plot_self_coverage.plot_only_non_duplicates</b><br />
        <i>Boolean? &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.plot_self_coverage.skip_mark_dupes">
        <b>assemble_refbased.plot_self_coverage.skip_mark_dupes</b><br />
        <i>Boolean? &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.skip_mark_dupes">
        <b>assemble_refbased.skip_mark_dupes</b><br />
        <i>Boolean? &mdash; Default: false</i><br />
        ???
</p>
<p name="assemble_refbased.trim_coords_bed">
        <b>assemble_refbased.trim_coords_bed</b><br />
        <i>File? &mdash; Default: None</i><br />
        ???
</p>
</details>






<hr />

> Generated using WDL AID (0.1.1)
