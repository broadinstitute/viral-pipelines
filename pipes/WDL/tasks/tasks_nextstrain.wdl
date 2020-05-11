version 1.0

task concatenate {
    meta {
        description: "This is nothing more than unix cat."
    }
    input {
        Array[File] infiles
        String      output_name
    }
    command {
        cat ~{sep=" " infiles} > "${output_name}"
    }
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File combined = "${output_name}"
    }
}

task filter_segments {
    input {
        File  all_samples_fasta
        Int?  segment = 1
        File? pre_assembled_samples_fasta

        Int? machine_mem_gb
    }
    command <<<
    python3 <<CODE

        segment = "-"+'~{segment}'
        segment_fasta = ""

        with open('~{all_samples_fasta}', 'r') as fasta:
            records=fasta.read().split(">")

            for r in records:
                if len(r.split("\n")) > 1:
                    header = r.split("\n")[0]

                    if segment in header:
                        new_header = header.replace(segment, "")
                        contents = r.replace(header, new_header)
                        segment_fasta += ">"+contents

        if '~{pre_assembled_samples_fasta}':
            with open('~{pre_assembled_samples_fasta}', 'r') as pre_assembled:
                segment_fasta += pre_assembled.read()
        print(segment_fasta)

    CODE
    >>>
    runtime {
        docker : "python"
        memory : select_first([machine_mem_gb, 3]) + " GB"
        cpu :    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File assembly_of_segment = stdout()
    }
}

task filter_subsample_sequences {
    meta {
        description: "Filter and subsample a sequence set. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/filter.html"
    }
    input {
        File     sequences_fasta
        File     metadata_tsv

        Int?     sequences_per_group
        String?  group_by
        File?    include
        File?    exclude

        Boolean  non_nucleotide=true

        String?  min_date
        String?  max_date
        Int?     min_length
        File?    priority
        Int?     subsample_seed
        String?  exclude_where
        String?  include_where

        Int?     machine_mem_gb
        String   docker = "nextstrain/base"
    }
    String in_basename = basename(sequences_fasta, ".fasta")
    command {
        augur filter \
            --sequences ~{sequences_fasta} \
            --metadata ~{metadata_tsv} \
            ~{"--min-date " + min_date} \
            ~{"--max-date " + max_date} \
            ~{"--min-length " + min_length} \
            ~{true="--non-nucleotide " false=""  non_nucleotide} \
            ~{"--exclude " + exclude} \
            ~{"--include " + include} \
            ~{"--priority " + priority} \
            ~{"--sequences-per-group " + sequences_per_group} \
            ~{"--group-by " + group_by} \
            ~{"--subsample-seed " + subsample_seed} \
            ~{"--exclude-where " + exclude_where} \
            ~{"--include-where " + include_where} \
            --output "~{in_basename}.filtered.fasta"
    }
    runtime {
        docker: docker
        memory: "4 GB"
        cpu :   2
        disks:  "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
    output {
        File filtered_fasta = "~{in_basename}.filtered.fasta"
    }
}

task augur_mafft_align {
    meta {
        description: "Align multiple sequences from FASTA. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/align.html"
    }
    input {
        File     sequences
        File     ref_fasta
        String   basename

        File?    existing_alignment
        Boolean  fill_gaps = true
        Boolean  remove_reference = true

        Int?     machine_mem_gb
        String   docker = "nextstrain/base"
    }
    command {
        augur align --sequences ~{sequences} \
            --reference-sequence ~{ref_fasta} \
            --output ~{basename}_aligned.fasta \
            ~{true="--fill-gaps" false="" fill_gaps} \
            ~{"--existing-alignment " + existing_alignment} \
            ~{true="--remove-reference" false="" remove_reference} \
            --debug \
            --nthreads auto
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MAX_RAM
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 104]) + " GB"
        cpu :   16
        disks:  "local-disk 375 LOCAL"
        preemptible: 2
        dx_instance_type: "mem3_ssd2_v2_x16"
    }
    output {
        File aligned_sequences = "~{basename}_aligned.fasta"
        File align_troubleshoot = stdout()
        File max_ram_usage_in_bytes = "MAX_RAM"
    }
}

task draft_augur_tree {
    meta {
        description: "Build a tree using a variety of methods. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/tree.html"
    }
    input {
        File     aligned_fasta
        String   basename

        String   method = "iqtree"
        String   substitution_model = "GTR"
        File?    exclude_sites
        File?    vcf_reference
        String?  tree_builder_args

        Int?     machine_mem_gb
        String   docker = "nextstrain/base"
    }
    command {
        augur tree --alignment ~{aligned_fasta} \
            --output ~{basename}_raw_tree.nwk \
            --method ~{default="iqtree" method} \
            --substitution-model ~{default="GTR" substitution_model} \
            ~{"--exclude-sites " + exclude_sites} \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--tree-builder-args " + tree_builder_args} \
            --nthreads auto
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 30]) + " GB"
        cpu :   16
        disks:  "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x16"
        preemptible: 2
    }
    output {
        File aligned_tree = "~{basename}_raw_tree.nwk"
    }
}

task refine_augur_tree {
    meta {
        description: "Refine an initial tree using sequence metadata. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/refine.html"
    }
    input {
        File     raw_tree
        File     aligned_fasta
        File     metadata
        String   basename

        Int?     gen_per_year
        Float?   clock_rate
        Float?   clock_std_dev
        Boolean  keep_root = false
        String?  root
        Boolean? covariance
        Boolean  keep_polytomies = false
        Int?     precision
        Boolean  date_confidence = true
        String?  date_inference = "marginal"
        String?  branch_length_inference
        String?  coalescent
        Int?     clock_filter_iqd = 4
        String?  divergence_units
        File?    vcf_reference

        Int?     machine_mem_gb
        String   docker = "nextstrain/base"
    }
    command {
        augur refine \
            --tree ~{raw_tree} \
            --alignment ~{aligned_fasta} \
            --metadata ~{metadata} \
            --output-tree ~{basename}_refined_tree.nwk \
            --output-node-data ~{basename}_branch_lengths.json \
            --timetree \
            ~{"--clock-rate " + clock_rate} \
            ~{"--clock-std-dev " + clock_std_dev} \
            ~{"--coalescent " + coalescent} \
            ~{"--clock-filter-iqd " + clock_filter_iqd} \
            ~{"--gen-per-year " + gen_per_year} \
            ~{"--root " + root} \
            ~{"--precision " + precision} \
            ~{"--date-inference " + date_inference} \
            ~{"--branch-length-inference " + branch_length_inference} \
            ~{"--divergence-units " + divergence_units} \
            ~{true="--covariance" false="--no-covariance" covariance} \
            ~{true="--keep-root" false="" keep_root} \
            ~{true="--keep-polytomies" false="" keep_polytomies} \
            ~{true="--date-confidence" false="" date_confidence} \
            ~{"--vcf-reference " + vcf_reference}
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 30]) + " GB"
        cpu :   16
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x16"
        preemptible: 2
    }
    output {
        File tree_refined  = "~{basename}_refined_tree.nwk"
        File branch_lengths = "~{basename}_branch_lengths.json"
    }
}

task ancestral_traits {
    meta {
        description: "Infer ancestral traits based on a tree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/traits.html"
    }
    input {
        File           tree
        File           metadata
        Array[String]+ columns
        String         basename

        Boolean        confidence = true
        File?          weights
        Float?         sampling_bias_correction

        Int?     machine_mem_gb
        String   docker = "nextstrain/base"
    }
    command {
        augur traits \
            --tree ~{tree} \
            --metadata ~{metadata} \
            --columns ~{sep=" " columns} \
            --output-node-data "~{basename}_nodes.json" \
            ~{"--weights " + weights} \
            ~{true="--confidence" false="" confidence}
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu :   2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 2
    }
    output {
        File node_data_json = "~{basename}_nodes.json"
    }
}

task ancestral_tree {
    meta {
        description: "Infer ancestral sequences based on a tree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/ancestral.html"
    }
    input {
        File     refined_tree
        File     aligned_fasta
        String   basename

        String   inference = "joint"
        Boolean  keep_ambiguous = false
        Boolean  infer_ambiguous = false
        Boolean  keep_overhangs = false
        File?    vcf_reference
        File?    output_vcf

        Int?     machine_mem_gb
        String   docker = "nextstrain/base"
    }
    command {
        augur ancestral \
            --tree ~{refined_tree} \
            --alignment ~{aligned_fasta} \
            --output-node-data ~{basename}_nt_muts.json \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--output-vcf " + output_vcf} \
            --output-sequences ~{basename}_ancestral_sequences.fasta \
            ~{true="--keep-overhangs" false="" keep_overhangs} \
            --inference ~{default="joint" inference} \
            ~{true="--keep-ambiguous" false="" keep_ambiguous} \
            ~{true="--infer-ambiguous" false="" infer_ambiguous}
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 7]) + " GB"
        cpu :   4
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x4"
        preemptible: 2
    }
    output {
        File nt_muts_json = "~{basename}_nt_muts.json"
        File sequences    = "~{basename}_ancestral_sequences.fasta"
    }
}

task translate_augur_tree {
    meta {
        description: "export augur files to json suitable for auspice visualization. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/translate.html"
    }
    input {
        String basename
        File   refined_tree
        File   nt_muts
        File   genbank_gb

        File?  genes
        File?  vcf_reference_output
        File?  vcf_reference

        Int?   machine_mem_gb
        String docker = "nextstrain/base"
    }
    command {
        augur translate --tree ~{refined_tree} \
            --ancestral-sequences ~{nt_muts} \
            --reference-sequence ~{genbank_gb} \
            ~{"--vcf-reference-output " + vcf_reference_output} \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--genes " + genes} \
            --output-node-data ~{basename}_aa_muts.json
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu :   2
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 2
    }
    output {
        File aa_muts_json = "~{basename}_aa_muts.json"
    }
}

task export_auspice_json {
    meta {
        description: "export augur files to json suitable for auspice visualization. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/export.html"
    }
    input {
        File   auspice_config
        File   metadata
        File   refined_tree
        File?  branch_lengths
        File?  traits
        File?  nt_muts
        File?  aa_muts
        File?  lat_longs_tsv
        File?  colors_tsv
        String basename

        Int?   machine_mem_gb
        String docker = "nextstrain/base"
    }
    command {

        NODE_DATA_FLAG=""
        if [[ -n "~{branch_lengths}" || -n "~{traits}" || -n "~{nt_muts}" || -n "~{aa_muts}" ]]; then
          NODE_DATA_FLAG="--node-data "
        fi

        augur export v2 --tree ~{refined_tree} \
            --metadata ~{metadata} \
            $NODE_DATA_FLAG ~{sep=' ' select_all([branch_lengths,traits,nt_muts,aa_muts])}\
            --auspice-config ~{auspice_config} \
            ~{"--lat-longs " + lat_longs_tsv} \
            ~{"--colors " + colors_tsv} \
            --output ~{basename}_auspice.json
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu :   2
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 2
    }
    output {
        File virus_json = "~{basename}_auspice.json"
    }
}
