version 1.0

task subsample_by_cases {
    meta {
        description: "Run subsampler to get downsampled dataset and metadata proportional to epidemiological case counts."
    }
    input {
        File    metadata
        File    case_data
        File    keep_file
        File    remove_file
        File    filter_file

        String? start_date
        String? end_date
        String  id_column
        String  geo_column
        String  date_column
        String  baseline        =   "0.001"
        String  refgenome_size  =   "1"
        String  max_missing     =   "99"
        String  seed_num        =   "2007"
        String  unit            =   "week"

        String  docker = "quay.io/broadinstitute/subsampler"
    }
    command {
        # subsampler is a snakemake target and must be run from /opt/subsampler
        cd /opt/subsampler
        snakemake subsample --config metadata=~{metadata} \
                                         case_data=~{case_data} \
                                         keep_file=~{keep_file} \
                                         remove_file=~{remove_file} \
                                         filter_file=~{filter_file} \
                                         id_column=~{id_column} \
                                         geo_column=~{geo_column} \
                                         date_column=~{date_column} \
                                         baseline=~{baseline} \
                                         refgenome_size=~{refgenome_size} \
                                         max_missing=~{max_missing} \
                                         seed_num=~{seed_num} \
                                         unit=~{unit} \
                                         ~{"start_date=" + start_date} \
                                         ~{"end_date=" + end_date}
        # put outputs in original working directory for delocalization
        cd -
        ln -s /opt/subsampler/outputs/* .
    }
    runtime {
        docker: docker
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File    genome_matrix_days              =   "genome_matrix_days.tsv"
        File    matrix_genomes_unit             =   "matrix_genomes_unit.tsv"
        File    matrix_cases_unit               =   "matrix_cases_unit.tsv"
        File    weekly_sampling_proportions     =   "weekly_sampling_proportions.tsv"
        File    weekly_sampling_bias            =   "weekly_sampling_bias.tsv"
        File    matrix_genomes_unit_corrected   =   "matrix_genomes_unit_corrected.tsv"
        File    selected_sequences              =   "selected_sequences.txt"
        File    selected_metadata               =   "selected_metadata.tsv"
        File    sampling_stats                  =   "sampling_stats.txt"
    }
}

task multi_align_mafft_ref {
  input {
    File         reference_fasta
    Array[File]+ assemblies_fasta # fasta files, one per sample, multiple chrs per file okay
    Int?         mafft_maxIters
    Float?       mafft_ep
    Float?       mafft_gapOpeningPenalty

    Int?         machine_mem_gb
    String       docker = "quay.io/broadinstitute/viral-phylo:2.1.20.2"
  }

  String         fasta_basename = basename(reference_fasta, '.fasta')
  Int disk_size = 200

  command {
    interhost.py --version | tee VERSION
    interhost.py multichr_mafft \
      ${reference_fasta} ${sep=' ' assemblies_fasta} \
      . \
      ${'--ep=' + mafft_ep} \
      ${'--gapOpeningPenalty=' + mafft_gapOpeningPenalty} \
      ${'--maxiters=' + mafft_maxIters} \
      --outFilePrefix align_mafft-${fasta_basename} \
      --preservecase \
      --localpair \
      --sampleNameListFile align_mafft-${fasta_basename}-sample_names.txt \
      --loglevel DEBUG
  }

  output {
    #File         sampleNamesFile = "align_mafft-${fasta_basename}-sample_names.txt"
    Array[File]+ alignments_by_chr = glob("align_mafft-${fasta_basename}*.fasta")
    String       viralngs_version  = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 60]) + " GB"
    cpu: 8
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem3_ssd1_v2_x8"
    maxRetries: 2
  }
}

task multi_align_mafft {
  input {
    Array[File]+ assemblies_fasta # fasta files, one per sample, multiple chrs per file okay
    String       out_prefix = "aligned"
    Int?         mafft_maxIters
    Float?       mafft_ep
    Float?       mafft_gapOpeningPenalty

    Int?         machine_mem_gb
    String       docker = "quay.io/broadinstitute/viral-phylo:2.1.20.2"
  }

  Int disk_size = 200

  command {
    interhost.py --version | tee VERSION
    interhost.py multichr_mafft \
      ${sep=' ' assemblies_fasta} \
      . \
      ${'--ep=' + mafft_ep} \
      ${'--gapOpeningPenalty=' + mafft_gapOpeningPenalty} \
      ${'--maxiters=' + mafft_maxIters} \
      --outFilePrefix ${out_prefix} \
      --preservecase \
      --localpair \
      --sampleNameListFile ${out_prefix}-sample_names.txt \
      --loglevel DEBUG
  }

  output {
    File        sampleNamesFile   = "${out_prefix}-sample_names.txt"
    Array[File] alignments_by_chr = glob("${out_prefix}*.fasta")
    String      viralngs_version  = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 30]) + " GB"
    cpu: 8
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem2_ssd1_v2_x8"
    maxRetries: 2
  }
}

task beast {
  input {
    File    beauti_xml

    Boolean beagle_double_precision=true
    String? beagle_order
    
    String? accelerator_type
    Int?    accelerator_count
    String? gpu_type
    Int?    gpu_count
    String? vm_size

    String  docker = "quay.io/broadinstitute/beast-beagle-cuda:1.10.5pre"
  }

  meta {
    description: "Execute GPU-accelerated BEAST. For tips on performance, see https://beast.community/performance#gpu"
  }
  parameter_meta {
    beagle_double_precision: {
      description: "If beagle_double_precision=true, use double-precision calculation (perhaps set to false to gain execution speed if MCMC chain convergence is possible using single-precision calculation)."
    }
    beagle_order: {
      description: "The order of CPU(0) and GPU(1+) resources used to process partitioned data."
    }

    accelerator_type: {
      description: "[GCP] The model of GPU to use. For availability and pricing on GCP, see https://cloud.google.com/compute/gpus-pricing#gpus"
    }
    accelerator_count: {
      description: "[GCP] The number of GPUs of the specified type to use."
    }
    gpu_type: {
      description: "[Terra] The model of GPU to use. For availability and pricing on GCP, see https://support.terra.bio/hc/en-us/articles/4403006001947-Getting-started-with-GPUs-in-a-Cloud-Environment"
    }
    gpu_count: {
      description: "[Terra] The number of GPUs of the specified type to use."
    }
  }

  Int disk_size    = 300
  Int boot_disk    = 50
  Int disk_size_az = disk_size + boot_disk

  # platform-agnostic number of GPUs we're actually using
  Int gpu_count_used = select_first([accelerator_count, gpu_count, 1])

  command {
    set -e
    beast -beagle_info
    nvidia-smi

    # if beagle_order is not specified by the user, 
    # create an appropriate string based on the gpu count

    default_beagle_order="$(seq -s, ~{gpu_count_used})"
    beagle_order=~{default="$default_beagle_order" beagle_order}
    echo "beagle_order: $beagle_order"

    bash -c "sleep 60; nvidia-smi" &
    beast \
      -beagle_multipartition off \
      -beagle_GPU \
      -beagle_cuda \
      -beagle_SSE \
      ~{true="-beagle_double" false="-beagle_single" beagle_double_precision} \
      -beagle_scaling always \
      ~{'-beagle_order ' + beagle_order} \
      ~{beauti_xml}
  }

  output {
    File        beast_log    = glob("*.log")[0]
    Array[File] trees        = glob("*.trees")
    Array[File] ops          = glob("*.ops")
    Array[File] rates        = glob("*.rates")
    Array[File] root         = glob("*.root")
    File        beast_stdout = stdout()
  }

  runtime {
    docker: "${docker}"
    memory: "7 GB"
    cpu:    4
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size_az + " GB"
    vm_size: select_first([accelerator_type, "Standard_NC6"])  # TES Azure
    maxRetries: 1
    bootDiskSizeGb: boot_disk
    gpu:                 true                # dxWDL
    dx_timeout:          "40H"               # dxWDL
    dx_instance_type:    "mem1_ssd1_gpu2_x8" # dxWDL
    acceleratorType:     select_first([accelerator_type, "nvidia-tesla-p4"])  # GCP PAPIv2
    acceleratorCount:    select_first([accelerator_count, 1])  # GCP PAPIv2
    gpuType:             select_first([gpu_type, "nvidia-tesla-p4"])  # Terra
    gpuCount:            select_first([gpu_count, 1])  # Terra
    nvidiaDriverVersion: "410.79"
  }
}

task index_ref {
  input {
    File   referenceGenome
    File?  novocraft_license

    Int?   machine_mem_gb
    String docker = "quay.io/broadinstitute/viral-core:2.1.33"
  }

  Int disk_size = 100

  command {
    read_utils.py --version | tee VERSION
    read_utils.py novoindex \
    "${referenceGenome}" \
    ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license}
    
    read_utils.py index_fasta_samtools "${referenceGenome}"
    read_utils.py index_fasta_picard "${referenceGenome}"
  }

  output {
    File   referenceNix     = "*.nix"
    File   referenceFai     = "*.fasta.fai"
    File   referenceDict    = "*.dict"
    String viralngs_version = read_string("VERSION")
  }
  runtime {
    docker: "${docker}"
    cpu: 2
    memory: select_first([machine_mem_gb, 4]) + " GB"
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES

    maxRetries: 2
  }
}

task trimal_clean_msa {
  input {
    File   in_aligned_fasta

    Int?   machine_mem_gb
    String docker = "quay.io/biocontainers/trimal:1.4.1--h6bb024c_3"

    String input_basename = basename(basename(in_aligned_fasta, ".fasta"), ".fa")
  }

  Int disk_size = 100

  command {
    trimal -fasta -automated1 -in "${in_aligned_fasta}" -out "${input_basename}_trimal_cleaned.fasta"
  }

  output {
    File trimal_cleaned_fasta = "${input_basename}_trimal_cleaned.fasta"
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x8"
    maxRetries: 2
  }
}

task merge_vcfs_bcftools {
  input {
    Array[File] in_vcfs_gz

    Int?   machine_mem_gb
    String docker = "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"

    String output_prefix = "merged"
  }

  parameter_meta {
    in_vcfs_gz: {
      description: "VCF files to merged; should be (b)gzipped.",
      patterns: ["*.vcf.gz"] }
  }

  command {

    # copy files to CWD (required for tabix indexing)
    INFILES=~{write_lines(in_vcfs_gz)}
    ln -s $(cat $INFILES) .
    for FN in $(cat $INFILES); do basename $FN; done > vcf_filenames.txt

    # tabix index input vcfs (must be gzipped)
    for FN in $(cat vcf_filenames.txt); do
      tabix -p vcf $FN
    done

    # see: https://samtools.github.io/bcftools/bcftools.html#merge
    # --merge snps allows snps to be merged to multi-allelic (multi-ALT) records, all other records are listed separately
    bcftools merge 
    --missing-to-ref \
    --force-samples \
    --merge snps \
    --output ${output_prefix}.vcf.gz \
    --output-type z \
    --threads "$(nproc --all)" \
    --file-list vcf_filenames.txt

    # tabix index the vcf to create .tbi file
    tabix -p vcf ${output_prefix}.vcf.gz
  }

  output {
    File merged_vcf_gz     = "${output_prefix}.vcf.gz"
    File merged_vcf_gz_tbi = "${output_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 3]) + " GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task merge_vcfs_gatk {
  input {
    Array[File] in_vcfs_gz
    File        ref_fasta

    Int?        machine_mem_gb
    String      docker = "quay.io/broadinstitute/viral-phylo:2.1.20.2"

    String      output_prefix = "merged"
  }

  parameter_meta {
    in_vcfs_gz: {
      description: "VCF files to merged; should be (b)gzipped.",
      patterns: ["*.vcf.gz"] 
    }
    ref_fasta: {
      description: "fasta file of reference genome relative to which the input VCF sites were called",
      patterns: ["*.fasta",".fa"]
    }
  }

  command {

    # tabix index input vcfs (must be gzipped)
    parallel -I ,, \
      "tabix -p vcf ,," \
      ::: "${sep=' ' in_vcfs_gz}"

    # index reference to create .fai and .dict indices
    samtools faidx "${ref_fasta}"
    picard CreateSequenceDictionary R="${ref_fasta}" O=$(basename $(basename "${ref_fasta}" .fasta) .fa).dict

    # store input vcf file paths in file
    for invcf in $(echo "${sep=' ' in_vcfs_gz}"); do 
      echo "$invcf" > input_vcfs.list
    done

    # merge
    gatk3 -T CombineVariants -R "${ref_fasta}" -V input_vcfs.list -o "${output_prefix}.vcf" -genotypeMergeOptions UNIQUIFY
    
    # bgzip output
    bgzip "${output_prefix}.vcf"

    # tabix index the vcf to create .tbi file
    tabix -p vcf "${output_prefix}.vcf.gz"
  }

  output {
    File merged_vcf_gz     = "${output_prefix}.vcf.gz"
    File merged_vcf_gz_tbi = "${output_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 3]) + " GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}
