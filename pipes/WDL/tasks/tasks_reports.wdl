version 1.0

task plot_coverage {
  input {
    File     aligned_reads_bam
    String   sample_name

    Boolean? skip_mark_dupes=false
    Boolean? plot_only_non_duplicates=false
    Boolean? bin_large_plots=false
    String?  binning_summary_statistic="max" # max or min

    Int?     machine_mem_gb
    String?  docker="quay.io/broadinstitute/viral-core"
  }
  
  command {
    set -ex -o pipefail

    read_utils.py --version | tee VERSION

    samtools view -c ${aligned_reads_bam} | tee reads_aligned
    if [ "$(cat reads_aligned)" != "0" ]; then
      samtools index -@ "$(nproc)" "${aligned_reads_bam}"

      PLOT_DUPE_OPTION=""
      if [[ "${skip_mark_dupes}" != "true" ]]; then
        PLOT_DUPE_OPTION="${true='--plotOnlyNonDuplicates' false="" plot_only_non_duplicates}"
      fi
      
      BINNING_OPTION="${true='--binLargePlots' false="" bin_large_plots}"

      # plot coverage
      reports.py plot_coverage \
        "${aligned_reads_bam}" \
        "${sample_name}.coverage_plot.pdf" \
        --outSummary "${sample_name}.coverage_plot.txt" \
        --plotFormat pdf \
        --plotWidth 1100 \
        --plotHeight 850 \
        --plotDPI 100 \
        $PLOT_DUPE_OPTION \
        $BINNING_OPTION \
        --binningSummaryStatistic ${binning_summary_statistic} \
        --plotTitle "${sample_name} coverage plot" \
        --loglevel=DEBUG

    else
      touch ${sample_name}.coverage_plot.pdf ${sample_name}.coverage_plot.txt
    fi

    # collect figures of merit
    samtools view -H ${aligned_reads_bam} | perl -n -e'/^@SQ.*LN:(\d+)/ && print "$1\n"' |  python -c "import sys; print(sum(int(x) for x in sys.stdin))" | tee assembly_length
    # report only primary alignments 260=exclude unaligned reads and secondary mappings
    samtools view -h -F 260 ${aligned_reads_bam} | samtools flagstat - | tee ${sample_name}.flagstat.txt
    grep properly ${sample_name}.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
    samtools view ${aligned_reads_bam} | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
    python -c "print (float("`cat bases_aligned`")/"`cat assembly_length`") if "`cat assembly_length`">0 else 0" > mean_coverage
  }

  output {
    File   coverage_plot                 = "${sample_name}.coverage_plot.pdf"
    File   coverage_tsv                  = "${sample_name}.coverage_plot.txt"
    Int    assembly_length               = read_int("assembly_length")
    Int    reads_aligned                 = read_int("reads_aligned")
    Int    read_pairs_aligned            = read_int("read_pairs_aligned")
    Int    bases_aligned                 = read_int("bases_aligned")
    Float  mean_coverage                 = read_float("mean_coverage")
    String viralngs_version              = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu: 2
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x4"
    preemptible: 1
  }
}

task coverage_report {
  input {
    Array[File]+ mapped_bams
    Array[File]  mapped_bam_idx # optional.. speeds it up if you provide it, otherwise we auto-index
    String       out_report_name="coverage_report.txt"

    Int?         machine_mem_gb
    String?      docker="quay.io/broadinstitute/viral-core"
  }

  command {
    reports.py --version | tee VERSION
    reports.py coverage_only \
      ${sep=' ' mapped_bams} \
      ${out_report_name} \
      --loglevel DEBUG
  }

  output {
    File   coverage_report  = "${out_report_name}"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 2]) + " GB"
    cpu: 2
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd2_v2_x4"
  }
}


task fastqc {
  input {
    File     reads_bam

    Int?     machine_mem_gb
    String?  docker="quay.io/broadinstitute/viral-core"
  }

  String   reads_basename=basename(reads_bam, ".bam")

  command {
    set -ex -o pipefail
    reports.py --version | tee VERSION
    reports.py fastqc ${reads_bam} ${reads_basename}_fastqc.html --out_zip ${reads_basename}_fastqc.zip
  }

  output {
    File   fastqc_html      = "${reads_basename}_fastqc.html"
    File   fastqc_zip      = "${reads_basename}_fastqc.zip"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: select_first([machine_mem_gb, 2]) + " GB"
    cpu: 1
    docker: "${docker}"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}

task align_and_count {
  input {
    File    reads_bam
    File    ref_db
    Int?    minScoreToFilter = 60
    Int?    topNHits = 3

    Int?    machine_mem_gb
    String? docker="quay.io/broadinstitute/viral-core"
  }

  String  reads_basename=basename(reads_bam, ".bam")
  String  ref_basename=basename(ref_db, ".fasta")

  command {
    set -ex -o pipefail

    read_utils.py --version | tee VERSION

    ln -s ${reads_bam} ${reads_basename}.bam
    read_utils.py bwamem_idxstats \
      ${reads_basename}.bam \
      ${ref_db} \
      --outStats ${reads_basename}.count.${ref_basename}.txt.unsorted \
      --minScoreToFilter=${minScoreToFilter} \
      --loglevel=DEBUG

      sort -b -r -n -k3 ${reads_basename}.count.${ref_basename}.txt.unsorted > ${reads_basename}.count.${ref_basename}.txt
      head -n ${topNHits} ${reads_basename}.count.${ref_basename}.txt > ${reads_basename}.count.${ref_basename}.top_${topNHits}_hits.txt
  }

  output {
    File   report           = "${reads_basename}.count.${ref_basename}.txt"
    File   report_top_hits  = "${reads_basename}.count.${ref_basename}.top_${topNHits}_hits.txt"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu: 4
    docker: "${docker}"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}

task align_and_count_summary {
  input {
    Array[File]+  counts_txt

    Int?          machine_mem_gb
    String?       docker="quay.io/broadinstitute/viral-core"
  }

  command {
    set -ex -o pipefail

    mkdir spike_summaries
    cp ${sep=' ' counts_txt} spike_summaries/

    reports.py --version | tee VERSION
    reports.py aggregate_spike_count spike_summaries/ count_summary.tsv \
      --loglevel=DEBUG
  }

  output {
    File   count_summary    = "count_summary.tsv"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: select_first([machine_mem_gb, 3]) + " GB"
    cpu: 2
    docker: "${docker}"
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}

task aggregate_metagenomics_reports {
  input {
    Array[File]+ kraken_summary_reports 
    String       aggregate_taxon_heading_space_separated  = "Viruses"
    String       aggregate_taxlevel_focus                 = "species"
    Int?         aggregate_top_N_hits                     = 5

    Int?         machine_mem_gb
    String?      docker="quay.io/broadinstitute/viral-classify"
  }

  parameter_meta {
    aggregate_taxon_heading_space_separated: { description: "The taxonomic heading to analyze. More than one can be specified." }
    aggregate_taxlevel_focus:                { description: "species,genus,family,order,class,phylum,kingdom,superkingdom" }
    aggregate_top_N_hits:                    { description: "only include the top N hits from a given sample in the aggregate report" }
  }

  String       aggregate_taxon_heading = sub(aggregate_taxon_heading_space_separated, " ", "_") # replace spaces with underscores for use in filename

  command {
    set -ex -o pipefail

    metagenomics.py --version | tee VERSION
    metagenomics.py taxlevel_summary \
      ${sep=' ' kraken_summary_reports} \
      --csvOut aggregate_taxa_summary_${aggregate_taxon_heading}_by_${aggregate_taxlevel_focus}_top_${aggregate_top_N_hits}_by_sample.csv \
      --noHist \
      --taxHeading ${aggregate_taxon_heading_space_separated} \
      --taxlevelFocus ${aggregate_taxlevel_focus} \
      --zeroFill --includeRoot --topN ${aggregate_top_N_hits} \
      --loglevel=DEBUG
  }

  output {
    File   krakenuniq_aggregate_taxlevel_summary = "aggregate_taxa_summary_${aggregate_taxon_heading}_by_${aggregate_taxlevel_focus}_top_${aggregate_top_N_hits}_by_sample.csv"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 4]) + " GB"
    cpu: 1
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd2_v2_x2"
    preemptible: 0
  }
}

task MultiQC {
  input {
    Array[File]     input_files = []

    Boolean         force = false
    Boolean         full_names = false
    String?         title
    String?         comment
    String?         file_name
    String          out_dir = "./multiqc-output"
    String?         template
    String?         tag
    String?         ignore_analysis_files
    String?         ignore_sample_names
    File?           sample_names
    Array[String]+? exclude_modules
    Array[String]+? module_to_use
    Boolean         data_dir = false
    Boolean         no_data_dir = false
    String?         output_data_format
    Boolean         zip_data_dir = false
    Boolean         export = false
    Boolean         flat = false
    Boolean         interactive = true
    Boolean         lint = false
    Boolean         pdf = false
    Boolean         megaQC_upload = false # Upload generated report to MegaQC if MegaQC options are found
    File?           config  # directory
    String?         config_yaml

    Int?            machine_mem_gb
    String          docker = "quay.io/biocontainers/multiqc:1.8--py_2"
  }

  parameter_meta {
    output_data_format: { description: "[tsv|yaml|json] default:tsv" }
  }

  # get the basename in all wdl use the filename specified (sans ".html" extension, if specified)
  String report_filename = if (defined(file_name)) then basename(select_first([file_name]), ".html") else "multiqc"

  command {
      set -ex -o pipefail

      echo "${sep='\n' input_files}" > input-filenames.txt

      multiqc \
      --file-list input-filenames.txt \
      --dirs \
      --outdir "${out_dir}" \
      ${true="--force" false="" force} \
      ${true="--fullnames" false="" full_names} \
      ${"--title " + title} \
      ${"--comment " + comment} \
      ${"--filename " + file_name} \
      ${"--template " + template} \
      ${"--tag " + tag} \
      ${"--ignore " + ignore_analysis_files} \
      ${"--ignore-samples" + ignore_sample_names} \
      ${"--sample-names " + sample_names} \
      ${true="--exclude " false="" defined(exclude_modules)}${sep=" --exclude " exclude_modules} \
      ${true="--module " false="" defined(module_to_use)}${sep=" --module " module_to_use} \
      ${true="--data-dir" false="" data_dir} \
      ${true="--no-data-dir" false="" no_data_dir} \
      ${"--data-format " + output_data_format} \
      ${true="--zip-data-dir" false="" zip_data_dir} \
      ${true="--export" false="" export} \
      ${true="--flat" false="" flat} \
      ${true="--interactive" false="" interactive} \
      ${true="--lint" false="" lint} \
      ${true="--pdf" false="" pdf} \
      ${false="--no-megaqc-upload" true="" megaQC_upload} \
      ${"--config " + config} \
      ${"--cl-config " + config_yaml }

      if [ -z "${file_name}" ]; then
        mv "${out_dir}/${report_filename}_report.html" "${out_dir}/${report_filename}.html"
      fi

      "${out_dir}/${report_filename}_data"
      tar -c "${out_dir}/${report_filename}_data" | gzip -c > "${report_filename}_data.tar.gz"
  }

  output {
      File multiqc_report            = "${out_dir}/${report_filename}.html"
      File multiqc_data_dir_tarball  = "${report_filename}_data.tar.gz"
  }

  runtime {
    memory: select_first([machine_mem_gb, 2]) + " GB"
    cpu: 1
    docker: "${docker}"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}
