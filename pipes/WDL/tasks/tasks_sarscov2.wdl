version 1.0

task nextclade_one_sample {
    meta {
        description: "Nextclade classification of one sample. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
    }
    input {
        File   genome_fasta
        File?  root_sequence
        File?  auspice_reference_tree_json
        File?  qc_config_json
        File?  gene_annotations_json
        File?  pcr_primers_csv
    }
    String basename = basename(genome_fasta, ".fasta")
    command {
        set -e
        nextclade.js --version > VERSION
        nextclade.js \
            --input-fasta "~{genome_fasta}" \
            ~{"--input-root-seq " + root_sequence} \
            ~{"--input-tree " + auspice_reference_tree_json} \
            ~{"--input-qc-config " + qc_config_json} \
            ~{"--input-gene-map " + gene_annotations_json} \
            ~{"--input-pcr-primers " + pcr_primers_csv} \
            --output-json "~{basename}".nextclade.json \
            --output-tsv  "~{basename}".nextclade.tsv \
            --output-tree "~{basename}".nextclade.auspice.json
        cp "~{basename}".nextclade.tsv input.tsv
        python3 <<CODE
        with open('input.tsv', 'rt') as inf:
            with open('transposed.tsv', 'wt') as outf:
                for c in zip(*(l.rstrip().split('\t') for l in inf)):
                    outf.write('\t'.join(c)+'\n')
        CODE
    }
    runtime {
        docker: "neherlab/nextclade:0.10.0"
        memory: "3 GB"
        cpu:    2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        String nextclade_version  = read_string("VERSION")
        File   nextclade_json     = "~{basename}.nextclade.json"
        File   auspice_json       = "~{basename}.nextclade.auspice.json"
        File   nextclade_tsv      = "~{basename}.nextclade.tsv"
        String nextclade_clade    = read_map("transposed.tsv")["clade"]
    }
}

task pangolin_one_sample {
    meta {
        description: "Pangolin classification of one SARS-CoV-2 sample."
    }
    input {
        File    genome_fasta
        Int?    min_length
        Float?  max_ambig
        Boolean include_putative = true
    }
    String basename = basename(genome_fasta, ".fasta")
    command {
        set -e
        pangolin -v > VERSION_PANGOLIN
        pangolin -lv > VERSION_LINEAGES
        pangolin -pv > VERSION_PANGOLEARN

        pangolin "~{genome_fasta}" \
            --outfile "~{basename}.pangolin_report.csv" \
            -t "$(nproc)" \
            --include-putative \
            ~{"--min-length " + min_length} \
            ~{"--max-ambig " + max_ambig} \
            ~{true="--include-putative" false="" include_putative} \
            --verbose

        cp "~{basename}.pangolin_report.csv" input.csv
        python3 <<CODE
        with open('input.csv', 'rt') as inf:
            with open('transposed.tsv', 'wt') as outf:
                for c in zip(*(l.rstrip().split(',') for l in inf)):
                    outf.write('\t'.join(c)+'\n')
        CODE
    }
    runtime {
        docker: "staphb/pangolin:2.1.1"
        memory: "3 GB"
        cpu:    2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        String pangolin_version   = read_string("VERSION_PANGOLIN")
        String lineages_version   = read_string("VERSION_LINEAGES")
        String pangolearn_version = read_string("VERSION_PANGOLEARN")
        File   pangolin_csv       = "~{basename}.pangolin_report.csv"
        String pangolin_clade     = read_map("transposed.tsv")["lineage"]
    }
}
