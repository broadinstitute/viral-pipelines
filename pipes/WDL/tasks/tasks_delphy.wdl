task delphy {
  input {
    File sequences_msa_fasta

    Int? num_steps
    Int? tree_every_n_steps
    Int? snapshot_every_n_steps
    Int? log_every_n_steps

    Int? num_threads # min($(nproc),(num_sequences / 100)) suggested

    Int? cpus
    Int? machine_mem_gb
    Int disk_dize = 300

    String docker = "quay.io/broadinstitute/delphy:0.999"
  }

  String out_basename = basename(basename(sequences_msa_fasta,".fa"),".fasta")

  meta {
    description: "Execute Delphy; see: https://github.com/broadinstitute/delphy"
  }
  parameter_meta {
    sequences_msa_fasta: {
      description: "multiple alignment of input sequences in fasta format"
    }
    num_steps: {
      description: "(5000000 * num_sequences) suggested"
    }
    tree_every_n_steps: {
      description: "(num_steps / 200) suggested"
    }
    snapshot_every_n_steps: {
      description: "(num_steps / 200) suggested"
    }
    log_every_n_steps: {
      description: "(num_steps / 10000) suggested"
    }
    num_threads: {
      description: "min($(nproc),(num_sequences / 100)) suggested"
    }
  }

  command <<<
    set -e
    delphy --version > DELPHY_VERSION

    num_sequences=$(grep -c ">" ~{sequences_msa_fasta})

    steps=~{if defined(num_steps) then "~{num_steps}" else "$((($num_sequences * 5000000)))"}
    tree_every=~{if defined(tree_every_n_steps) then "~{tree_every_n_steps}" else "$((($steps / 200)))"}
    snapshot_every=~{if defined(snapshot_every_n_steps) then "~{snapshot_every_n_steps}" else "$((($steps / 200)))"}
    log_every=~{if defined(log_every_n_steps) then "~{log_every_n_steps}" else "$((($steps / 10000)))"}

    threads=~{if defined(num_threads) then "~{num_threads}" else "$(nproc)"}

    delphy \
       --v0-in-fasta ~{sequences_msa_fasta} \
       --v0-threads $threads \
       #[--v0-site-rate-heterogeneity] \
       --v0-steps $steps \
       --v0-tree-every $tree_every \
       --v0-delphy-snapshot-every $snapshot_every \
       --v0-log-every $log_every \
       --v0-out-trees-file  "~{out_basename}.dphy.trees" \
       --v0-out-log-file    "~{out_basename}.dphy.log" \
       --v0-out-delphy-file "~{out_basename}.dphy" \
       --v0-out-beast-xml   "~{out_basename}.dphy.for_beast.xml"
  >>>

  output {
    File   delphy_trees            = "~{out_basename}.dphy.trees"
    File   delphy_log              = "~{out_basename}.dphy.log"
    File   delphy_for_web_ui       = "~{out_basename}.dphy"
    File   delphy_output_for_beast = "~{out_basename}.dphy.for_beast.xml"
    File   delphy_stdout           = stdout()
    
    String delphy_version          = read_string('DELPHY_VERSION')
  }

  runtime {
    docker:           docker
    memory:           select_first([machine_mem_gb, 15]) + " GB"
    cpu:              select_first([cpus, 8])
    disks:            "local-disk " + disk_size + " HDD"
    disk:             disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_gpu2_x8" # dxWDL

    maxRetries: 1
  }
}