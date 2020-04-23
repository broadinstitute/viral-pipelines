version 1.0

task concatenate {
    # this is nothing more than unix cat
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

task augur_mafft_align {
    input {
        File     sequences
        File     ref_fasta
        String   basename

        Boolean? existing_alignment = false
        Boolean? debug = false
        Boolean? fill_gaps = true
        Boolean? remove_reference = true

        Int?     machine_mem_gb
        String   docker = "nextstrain/base"
    }
    command {
        augur align --sequences ~{sequences} \
            --reference-sequence ~{ref_fasta} \
            --output ~{basename}_aligned.fasta \
            ~{true="--fill-gaps" false="" fill_gaps} \
            ~{true="--existing-alignment " false="" existing_alignment} \
            ~{true="--remove-reference" false="" remove_reference} \
            ~{true="--debug" false="" debug} \
            --nthreads auto
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes | tee MAX_RAM
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
        Int  max_ram_usage_in_bytes = read_int("MAX_RAM")
    }
}

task draft_augur_tree {
    input {
        File     aligned_fasta
        String   basename

        String?  method              # default iqtree
        String?  substitution_model  # default GTR
        File?    exclude_sites
        File?    vcf_reference

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
    input {
        File     raw_tree
        File     aligned_fasta
        File     metadata
        String   basename

        Int?     gen_per_year
        Boolean? clock_rate = false
        Boolean? clock_std_dev = false
        Boolean? keep_root = false
        String?  root
        Boolean? covariance = false
        Boolean? no_covariance = false
        Boolean? keep_polytomies = false
        Int?     precision
        Boolean? date_confidence = false
        String?  date_inference 
        String?  branch_length_inference
        Boolean? clock_filter_iqd = false
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
            --timetree ~{true="--clock-rate" false="" clock_rate} \
            ~{true="--clock-std-dev" false="" clock_std_dev} \
            --gen-per-year ~{default=50 gen_per_year} \
            ~{true="--covariance" false="" covariance} \
            ~{true="--no-covariance" false="" no_covariance} \
            --root ~{default="best" root} \
            ~{true="--keep-root" false="" keep_root} \
            --precision ~{default=1 precision} \
            ~{true="--keep-polytomies" false="" keep_polytomies} \
            --date-inference ~{default="joint" date_inference} \
            ~{true="--date-confidence" false="" date_confidence} \
            --branch-length-inference ~{default="auto" branch_length_inference} \
            ~{true="--clock-filter-iqd" false="" clock_filter_iqd} \
            --divergence-units ~{default="mutations-per-site" divergence_units} \
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

task ancestral_tree {
    input {
        File     refined_tree
        File     aligned_fasta
        String   basename

        String?  inference
        Boolean? keep_ambiguous = false
        Boolean? infer_ambiguous = false
        Boolean? keep_overhangs = false
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
            ~{true="--keep-0verhands" false="" keep_overhangs} \
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
    input {
        File   auspice_config
        File   metadata
        File   refined_tree
        File?   branch_lengths
        File?   nt_muts
        File?   aa_muts
        String basename

        Int?   machine_mem_gb
        String docker = "nextstrain/base"
    }
    command {
        augur export v2 --tree ~{refined_tree} \
            --metadata ~{metadata} \
            --node-data ~{branch_lengths} ~{nt_muts} ~{aa_muts}\
            --auspice-config ~{auspice_config} \
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
