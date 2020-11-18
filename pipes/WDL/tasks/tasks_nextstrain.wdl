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

task fasta_to_ids {
    meta {
        description: "Return the headers only from a fasta file"
    }
    input {
        File sequences_fasta
    }
    String basename = basename(sequences_fasta, ".fasta")
    command {
        cat "~{sequences_fasta}" | grep \> | cut -c 2- > "~{basename}.txt"
    }
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File ids_txt = "~{basename}.txt"
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
        File     sample_metadata_tsv

        Int?     sequences_per_group
        String?  group_by
        File?    include
        File?    exclude

        Boolean  non_nucleotide=true

        Float?   min_date
        Float?   max_date
        Int?     min_length
        File?    priority
        Int?     subsample_seed
        Array[String]?  exclude_where
        Array[String]?  include_where

        String   docker = "nextstrain/base:build-20200629T201240Z"
    }
    parameter_meta {
        sequences_fasta: {
          description: "Set of sequences (unaligned fasta or aligned fasta -- one sequence per genome) or variants (vcf format) to subsample using augur filter.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
        sample_metadata_tsv: {
          description: "Metadata in tab-separated text format. See https://nextstrain-augur.readthedocs.io/en/stable/faq/metadata.html for details.",
          patterns: ["*.txt", "*.tsv"]
        }
    }
    String out_fname = sub(sub(basename(sequences_fasta), ".vcf", ".filtered.vcf"), ".fasta$", ".filtered.fasta")
    command {
        set -e
        augur version > VERSION

        touch wherefile
        VALS="~{write_lines(select_first([exclude_where, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--exclude-where" >> wherefile
            cat $VALS >> wherefile
        fi
        VALS="~{write_lines(select_first([include_where, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--include-where" >> wherefile
            cat $VALS >> wherefile
        fi

        cat wherefile | tr '\n' '\0' | xargs -0 -t augur filter \
            --sequences ~{sequences_fasta} \
            --metadata ~{sample_metadata_tsv} \
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
            --output "~{out_fname}" | tee STDOUT
        #cat ~{sequences_fasta} | grep \> | wc -l > IN_COUNT
        grep "sequences were dropped during filtering" STDOUT | cut -f 1 -d ' ' > DROP_COUNT
        grep "sequences have been written out to" STDOUT | cut -f 1 -d ' ' > OUT_COUNT
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu :   4
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x4"
        preemptible: 1
    }
    output {
        File   filtered_fasta    = out_fname
        String augur_version     = read_string("VERSION")
        Int    sequences_dropped = read_int("DROP_COUNT")
        Int    sequences_out     = read_int("OUT_COUNT")
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
    }
}

task filter_sequences_to_list {
    meta {
        description: "Filter and subsample a sequence set to a specific list of ids in a text file (one id per line)."
    }
    input {
        File          sequences
        Array[File]?  keep_list

        String   docker = "nextstrain/base:build-20200629T201240Z"
    }
    parameter_meta {
        sequences: {
          description: "Set of sequences (unaligned fasta or aligned fasta -- one sequence per genome) or variants (vcf format) to subsample using augur filter.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
        keep_list: {
          description: "List of strain ids.",
          patterns: ["*.txt", "*.tsv"]
        }
    }
    String out_fname = sub(sub(basename(sequences), ".vcf", ".filtered.vcf"), ".fasta$", ".filtered.fasta")
    Int mem_size = ceil(size(sequences, "GB") * 2 + 0.001)
    command {
        set -e
        augur version > VERSION
        KEEP_LISTS="~{sep=' ' select_first([keep_list, []])}"

        if [ -n "$KEEP_LISTS" ]; then
            echo "strain" > keep_list.txt
            cat $KEEP_LISTS >> keep_list.txt
            augur filter \
                --sequences "~{sequences}" \
                --metadata keep_list.txt \
                --output "~{out_fname}" | tee STDOUT
            grep "sequences were dropped during filtering" STDOUT | cut -f 1 -d ' ' > DROP_COUNT
            grep "sequences have been written out to" STDOUT | cut -f 1 -d ' ' > OUT_COUNT
        else
            cp "~{sequences}" "~{out_fname}"
            echo "0" > DROP_COUNT
            echo "-1" > OUT_COUNT
        fi
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: mem_size + " GB"
        cpu :   2
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
    output {
        File   filtered_fasta    = out_fname
        String augur_version     = read_string("VERSION")
        Int    sequences_dropped = read_int("DROP_COUNT")
        Int    sequences_out     = read_int("OUT_COUNT")
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
    }
}

task mafft_one_chr {
    meta {
        description: "Align multiple sequences from FASTA. Only appropriate for closely related (within 99% nucleotide conservation) genomes. See https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html"
    }
    input {
        File     sequences
        File?    ref_fasta
        String   basename
        Boolean  remove_reference = false
        Boolean  keep_length = true
        Boolean  large = false
        Boolean  memsavetree = false

        String   docker = "quay.io/broadinstitute/viral-phylo:2.1.10.0"
        Int      mem_size = 60
        Int      cpus = 32
    }
    command {
        set -e
        touch args.txt

        # boolean options
        echo "~{true='--large' false='' large}" >> args.txt
        echo "~{true='--memsavetree' false='' memsavetree}" >> args.txt
        echo "--auto" >> args.txt

        # if ref_fasta is specified, use "closely related" mode
        # see https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html
        if [ -f "~{ref_fasta}" ]; then
            echo --addfragments >> args.txt
            echo "~{sequences}" >> args.txt
            echo "~{ref_fasta}" >> args.txt
        else
            echo "~{sequences}" >> args.txt
        fi

        # mafft align to reference in "closely related" mode
        cat args.txt | grep . | xargs -d '\n' mafft --thread -1 \
            ~{true='--keeplength --mapout' false='' keep_length} \
            > msa.fasta

        # remove reference sequence
        python3 <<CODE
        import Bio.SeqIO
        seq_it = Bio.SeqIO.parse('msa.fasta', 'fasta')
        print("dumping " + str(seq_it.__next__().id))
        Bio.SeqIO.write(seq_it, 'msa_drop_one.fasta', 'fasta')
        CODE
        REMOVE_REF="~{true='--remove-reference' false='' remove_reference}"
        if [ -n "$REMOVE_REF" -a -f "~{ref_fasta}" ]; then
            mv msa_drop_one.fasta "~{basename}_aligned.fasta"
        else
            mv msa.fasta "~{basename}_aligned.fasta"
        fi

        # profiling and stats
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: mem_size + " GB"
        cpu :   cpus
        disks:  "local-disk 100 HDD"
        preemptible: 0
        dx_instance_type: "mem1_ssd1_v2_x36"
    }
    output {
        File   aligned_sequences = "~{basename}_aligned.fasta"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
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

        String   docker = "nextstrain/base:build-20200629T201240Z"
    }
    command {
        set -e
        augur version > VERSION
        augur align --sequences ~{sequences} \
            --reference-sequence ~{ref_fasta} \
            --output ~{basename}_aligned.fasta \
            ~{true="--fill-gaps" false="" fill_gaps} \
            ~{"--existing-alignment " + existing_alignment} \
            ~{true="--remove-reference" false="" remove_reference} \
            --debug \
            --nthreads auto
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "180 GB"
        cpu :   64
        disks:  "local-disk 750 LOCAL"
        preemptible: 0
        dx_instance_type: "mem3_ssd2_v2_x32"
    }
    output {
        File   aligned_sequences = "~{basename}_aligned.fasta"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task snp_sites {
    input {
        File   msa_fasta
        String docker = "quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0"
    }
    String out_basename = basename(msa_fasta, ".fasta")
    command {
        snp-sites -V > VERSION
        snp-sites -v -c -o ~{out_basename}.vcf ~{msa_fasta}
    }
    runtime {
        docker: docker
        memory: "1 GB"
        cpu :   1
        disks:  "local-disk 50 HDD"
        preemptible: 0
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File   snps_vcf = "~{out_basename}.vcf"
        String snp_sites_version = read_string("VERSION")
    }
}

task augur_mask_sites {
    meta {
        description: "Mask unwanted positions from alignment or SNP table. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/mask.html"
    }
    input {
        File     sequences
        File?    mask_bed

        String   docker = "nextstrain/base:build-20200629T201240Z"
    }
    parameter_meta {
        sequences: {
          description: "Set of alignments (fasta format) or variants (vcf format) to mask.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_fname = sub(sub(basename(sequences), ".vcf", ".masked.vcf"), ".fasta$", ".masked.fasta")
    command {
        set -e
        augur version > VERSION
        BEDFILE=~{select_first([mask_bed, "/dev/null"])}
        if [ -s "$BEDFILE" ]; then
            augur mask --sequences ~{sequences} \
                --mask "$BEDFILE" \
                --output "~{out_fname}"
        else
            cp "~{sequences}" "~{out_fname}"
        fi
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu :   4
        disks:  "local-disk 100 HDD"
        preemptible: 1
        dx_instance_type: "mem1_ssd1_v2_x4"
    }
    output {
        File   masked_sequences = out_fname
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version  = read_string("VERSION")
    }
}

task draft_augur_tree {
    meta {
        description: "Build a tree using iqTree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/tree.html"
    }
    input {
        File     msa_or_vcf

        String   method = "iqtree"
        String   substitution_model = "GTR"
        File?    exclude_sites
        File?    vcf_reference
        String?  tree_builder_args

        Int?     cpus
        String   docker = "nextstrain/base:build-20200629T201240Z"
    }
    parameter_meta {
        msa_or_vcf: {
          description: "Set of alignments (fasta format) or variants (vcf format) to construct a tree from using augur tree (iqTree).",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_basename = basename(basename(basename(msa_or_vcf, '.gz'), '.vcf'), '.fasta')
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur tree --alignment ~{msa_or_vcf} \
            --output ~{out_basename}_~{method}.nwk \
            --method ~{method} \
            --substitution-model ~{default="GTR" substitution_model} \
            ~{"--exclude-sites " + exclude_sites} \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--tree-builder-args " + tree_builder_args} \
            --nthreads auto
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "32 GB"
        cpu:    select_first([cpus, 64])
        disks:  "local-disk 750 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x36"
        preemptible: 0
    }
    output {
        File   aligned_tree = "~{out_basename}_~{method}.nwk"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task refine_augur_tree {
    meta {
        description: "Refine an initial tree using sequence metadata and Treetime. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/refine.html"
    }
    input {
        File     raw_tree
        File     msa_or_vcf
        File     metadata

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

        String   docker = "nextstrain/base:build-20200629T201240Z"
    }
    parameter_meta {
        msa_or_vcf: {
          description: "Set of alignments (fasta format) or variants (vcf format) to use to guide Treetime.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_basename = basename(basename(basename(msa_or_vcf, '.gz'), '.vcf'), '.fasta')
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur refine \
            --tree ~{raw_tree} \
            --alignment ~{msa_or_vcf} \
            --metadata ~{metadata} \
            --output-tree ~{out_basename}_timetree.nwk \
            --output-node-data ~{out_basename}_branch_lengths.json \
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
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "50 GB"
        cpu :   2
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem3_ssd1_v2_x8"
        preemptible: 0
    }
    output {
        File   tree_refined  = "~{out_basename}_timetree.nwk"
        File   branch_lengths = "~{out_basename}_branch_lengths.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task ancestral_traits {
    meta {
        description: "Infer ancestral traits based on a tree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/traits.html"
    }
    input {
        File           tree
        File           metadata
        Array[String]  columns

        Boolean        confidence = true
        File?          weights
        Float?         sampling_bias_correction

        String   docker = "nextstrain/base:build-20200629T201240Z"
    }
    String out_basename = basename(tree, '.nwk')
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur traits \
            --tree ~{tree} \
            --metadata ~{metadata} \
            --columns ~{sep=" " columns} \
            --output-node-data "~{out_basename}_ancestral_traits.json" \
            ~{"--weights " + weights} \
            ~{"--sampling-bias-correction " + sampling_bias_correction} \
            ~{true="--confidence" false="" confidence}
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu :   2
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
    output {
        File   node_data_json = "~{out_basename}_ancestral_traits.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task ancestral_tree {
    meta {
        description: "Infer ancestral sequences based on a tree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/ancestral.html"
    }
    input {
        File     tree
        File     msa_or_vcf

        String   inference = "joint"
        Boolean  keep_ambiguous = false
        Boolean  infer_ambiguous = false
        Boolean  keep_overhangs = false
        File?    vcf_reference
        File?    output_vcf

        String   docker = "nextstrain/base:build-20200629T201240Z"
    }
    parameter_meta {
        msa_or_vcf: {
          description: "Set of alignments (fasta format) or variants (vcf format) to use to guide Treetime.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_basename = basename(basename(basename(msa_or_vcf, '.gz'), '.vcf'), '.fasta')
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur ancestral \
            --tree ~{tree} \
            --alignment ~{msa_or_vcf} \
            --output-node-data ~{out_basename}_nt_muts.json \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--output-vcf " + output_vcf} \
            --output-sequences ~{out_basename}_ancestral_sequences.fasta \
            ~{true="--keep-overhangs" false="" keep_overhangs} \
            --inference ~{default="joint" inference} \
            ~{true="--keep-ambiguous" false="" keep_ambiguous} \
            ~{true="--infer-ambiguous" false="" infer_ambiguous}
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "50 GB"
        cpu :   4
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem3_ssd1_v2_x8"
        preemptible: 0
    }
    output {
        File   nt_muts_json = "~{out_basename}_nt_muts.json"
        File   sequences    = "~{out_basename}_ancestral_sequences.fasta"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task translate_augur_tree {
    meta {
        description: "export augur files to json suitable for auspice visualization. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/translate.html"
    }
    input {
        File   tree
        File   nt_muts
        File   genbank_gb

        File?  genes
        File?  vcf_reference_output
        File?  vcf_reference

        String docker = "nextstrain/base:build-20200629T201240Z"
    }
    String out_basename = basename(tree, '.nwk')
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur translate --tree ~{tree} \
            --ancestral-sequences ~{nt_muts} \
            --reference-sequence ~{genbank_gb} \
            ~{"--vcf-reference-output " + vcf_reference_output} \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--genes " + genes} \
            --output-node-data ~{out_basename}_aa_muts.json
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "2 GB"
        cpu :   1
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
    output {
        File   aa_muts_json = "~{out_basename}_aa_muts.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        String augur_version = read_string("VERSION")
    }
}

task assign_clades_to_nodes {
    meta {
        description: "Assign taxonomic clades to tree nodes based on mutation information"
    }
    input {
        File tree_nwk
        File nt_muts_json
        File aa_muts_json
        File ref_fasta
        File clades_tsv

        String docker = "nextstrain/base:build-20200629T201240Z"
    }
    String out_basename = basename(basename(tree_nwk, ".nwk"), "_timetree")
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur clades \
        --tree ~{tree_nwk} \
        --mutations ~{nt_muts_json} ~{aa_muts_json} \
        --reference ~{ref_fasta} \
        --clades ~{clades_tsv} \
        --output-node-data ~{out_basename}_clades.json
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "2 GB"
        cpu :   1
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
    output {
        File   node_clade_data_json = "~{out_basename}_clades.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        String augur_version      = read_string("VERSION")
    }
}

task augur_import_beast {
    meta {
        description: "Import BEAST tree into files ready for augur export, including a Newick-formatted tree and node data in json format. See both https://nextstrain-augur.readthedocs.io/en/stable/faq/import-beast.html and https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/import.html for details."
    }
    input {
        File    beast_mcc_tree

        Float?  most_recent_tip_date
        String? tip_date_regex
        String? tip_date_format
        String? tip_date_delimiter

        Int?    machine_mem_gb
        String  docker = "nextstrain/base:build-20200629T201240Z"
    }
    String tree_basename = basename(beast_mcc_tree, ".tree")
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur import beast \
            --mcc "~{beast_mcc_tree}" \
            --output-tree "~{tree_basename}.nwk" \
            --output-node-data "~{tree_basename}.json" \
            ~{"--most-recent-tip-date " + most_recent_tip_date} \
            ~{"--tip-date-regex " + tip_date_regex} \
            ~{"--tip-date-format " + tip_date_format} \
            ~{"--tip-date-delimeter " + tip_date_delimiter}
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu :   2
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
    output {
        File   tree_newick    = "~{tree_basename}.nwk"
        File   node_data_json = "~{tree_basename}.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task export_auspice_json {
    meta {
        description: "export augur files to json suitable for auspice visualization. The metadata tsv input is generally required unless the node_data_jsons comprehensively capture all of it. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/export.html"
    }
    input {
        File        auspice_config
        File?       sample_metadata
        File        tree
        Array[File] node_data_jsons

        File?          lat_longs_tsv
        File?          colors_tsv
        Array[String]? geo_resolutions
        Array[String]? color_by_metadata
        File?          description_md
        Array[String]? maintainers
        String?        title

        String docker = "nextstrain/base:build-20200629T201240Z"
    }
    String out_basename = basename(basename(tree, ".nwk"), "_timetree")
    command {
        set -e -o pipefail
        augur version > VERSION
        touch exportargs

        # --node-data
        if [ -n "~{sep=' ' node_data_jsons}" ]; then
            echo "--node-data" >> exportargs
            cat "~{write_lines(node_data_jsons)}" >> exportargs
        fi

        # --geo-resolutions
        VALS="~{write_lines(select_first([geo_resolutions, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--geo-resolutions" >> exportargs;
        fi
        cat $VALS >> exportargs

        # --color-by-metadata
        VALS="~{write_lines(select_first([color_by_metadata, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--color-by-metadata" >> exportargs;
        fi
        cat $VALS >> exportargs

        # --title
        if [ -n "~{title}" ]; then
            echo "--title" >> exportargs
            echo "~{title}" >> exportargs
        fi

        # --maintainers
        VALS="~{write_lines(select_first([maintainers, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--maintainers" >> exportargs;
        fi
        cat $VALS >> exportargs

        (export AUGUR_RECURSION_LIMIT=10000; cat exportargs | tr '\n' '\0' | xargs -0 -t augur export v2 \
            --tree ~{tree} \
            ~{"--metadata " + sample_metadata} \
            --auspice-config ~{auspice_config} \
            ~{"--lat-longs " + lat_longs_tsv} \
            ~{"--colors " + colors_tsv} \
            ~{"--description " + description_md} \
            --output ~{out_basename}_auspice.json)
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu :   2
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 0
    }
    output {
        File   virus_json = "~{out_basename}_auspice.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}
