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
        # transpose table
        import codecs
        with codecs.open('input.tsv', 'r', encoding='utf-8') as inf:
            with codecs.open('transposed.tsv', 'w', encoding='utf-8') as outf:
                for c in zip(*(l.rstrip().split('\t') for l in inf)):
                    outf.write('\t'.join(c)+'\n')
        CODE
        grep ^clade transposed.tsv | cut -f 2 | grep -v clade > NEXTCLADE_CLADE
        grep ^aaSubstitutions transposed.tsv | cut -f 2 | grep -v aaSubstitutions > NEXTCLADE_AASUBS
        grep ^aaDeletions transposed.tsv | cut -f 2 | grep -v aaDeletions > NEXTCLADE_AADELS
    }
    runtime {
        docker: "nextstrain/nextclade:0.13.0"
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
        String nextclade_clade    = read_string("NEXTCLADE_CLADE")
        String aa_subs_csv        = read_string("NEXTCLADE_AASUBS")
        String aa_dels_csv        = read_string("NEXTCLADE_AADELS")
    }
}

task nextclade_many_samples {
    meta {
        description: "Nextclade classification of many samples. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
    }
    input {
        Array[File]+ genome_fastas
        File?  root_sequence
        File?  auspice_reference_tree_json
        File?  qc_config_json
        File?  gene_annotations_json
        File?  pcr_primers_csv
        String basename
    }
    command {
        set -e
        nextclade.js --version > VERSION
        cat ~{sep=" " genome_fastas} > genomes.fasta
        nextclade.js \
            --input-fasta genomes.fasta \
            ~{"--input-root-seq " + root_sequence} \
            ~{"--input-tree " + auspice_reference_tree_json} \
            ~{"--input-qc-config " + qc_config_json} \
            ~{"--input-gene-map " + gene_annotations_json} \
            ~{"--input-pcr-primers " + pcr_primers_csv} \
            --output-json "~{basename}".nextclade.json \
            --output-tsv  "~{basename}".nextclade.tsv \
            --output-tree "~{basename}".nextclade.auspice.json
    }
    runtime {
        docker: "nextstrain/nextclade:0.13.0"
        memory: "14 GB"
        cpu:    16
        disks: "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x16"
    }
    output {
        String nextclade_version  = read_string("VERSION")
        File   nextclade_json     = "~{basename}.nextclade.json"
        File   auspice_json       = "~{basename}.nextclade.auspice.json"
        File   nextclade_tsv      = "~{basename}.nextclade.tsv"
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
        String  docker = "staphb/pangolin:2.3.8-pangolearn-2021-04-01"
    }
    String basename = basename(genome_fasta, ".fasta")
    command {
        set -e
        pangolin -v > VERSION_PANGOLIN
        pangolin -pv > VERSION_PANGOLEARN

        pangolin "~{genome_fasta}" \
            --outfile "~{basename}.pangolin_report.csv" \
            ~{"--min-length " + min_length} \
            ~{"--max-ambig " + max_ambig} \
            --alignment \
            --verbose

        cp sequences.aln.fasta "~{basename}.pangolin_msa.fasta"
        cp "~{basename}.pangolin_report.csv" input.csv
        python3 <<CODE
        # transpose table and convert csv to tsv
        with open('input.csv', 'rt') as inf:
            with open('transposed.tsv', 'wt') as outf:
                for c in zip(*(l.rstrip().split(',') for l in inf)):
                    outf.write('\t'.join(c)+'\n')
        CODE
        grep ^lineage transposed.tsv | cut -f 2 | grep -v lineage > PANGOLIN_CLADE
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu:    2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        String pangolin_version   = read_string("VERSION_PANGOLIN")
        String pangolearn_version = read_string("VERSION_PANGOLEARN")
        File   pangolin_csv       = "~{basename}.pangolin_report.csv"
        String pango_lineage      = read_string("PANGOLIN_CLADE")
        File   msa_fasta          = "~{basename}.pangolin_msa.fasta"
        String pangolin_docker    = docker
    }
}


task sequencing_report {
    meta {
        description: "Produce sequencing progress report."
    }
    input {
        File           assembly_stats_tsv
        File?          collab_ids_tsv

        String?        sequencing_lab = "Broad Institute"
        String?        intro_blurb = "The Broad Institute Viral Genomics group, in partnership with the Genomics Platform and Data Sciences Platform, has been engaged in viral sequencing of COVID-19 patients since March 2020."
        String?        max_date
        String?        min_date
        Int?           min_unambig
        String?        voc_list
        String?        voi_list

        String  docker = "quay.io/broadinstitute/sc2-rmd:0.1.10"
    }
    command {
        set -e
        /docker/reports.py \
            "~{assembly_stats_tsv}" \
            ~{'--collab_tsv="' + collab_ids_tsv + '"'} \
            ~{'--sequencing_lab="' + sequencing_lab + '"'} \
            ~{'--intro_blurb="' + intro_blurb + '"'} \
            ~{'--max_date=' + max_date} \
            ~{'--min_date=' + min_date} \
            ~{'--min_unambig=' + min_unambig} \
            ~{'--voc_list=' + voc_list} \
            ~{'--voi_list=' + voi_list}
        zip all_reports.zip *.pdf *.xlsx
    }
    runtime {
        docker: docker
        memory: "2 GB"
        cpu:    2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        Array[File] reports = glob("*.pdf")
        Array[File] sheets = glob("*.xlsx")
        File        all_zip = "all_reports.zip"
    }
}



