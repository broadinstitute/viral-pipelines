version 1.0 

task vadr_tool {
    meta {
        description: "VADR test run"
    }
    parameter_meta {
        genome_fasta: {
            description: "Input fasta file to annotate."
            category: "required"
        }
        library_accession_no: {
            description: "Accesion number to build VADR DB from."
            category: "required"
        }
    }
    input {
        File    genome_fasta
        String  library_accession_no
        String  docker = "quay.io/staphb/vadr:1.5.1"
        Int     minlen = 50
        Int     maxlen = 30000
        Int    mem_size = 4
        Int     disk_size_gb = ceil(2*20) + 5
        Int    cpu = 2
  }
  String out_base = basename(genome_fasta, '.fasta')
  command <<<
    set -ex -o pipefail

    # remove terminal ambiguous nucleotides
    /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
      "~{genome_fasta}" \
      --minlen ~{minlen} \
      --maxlen ~{maxlen} \
      > "~{out_base}.fasta"

    #build VADR database
    v-build.pl "~{library_accession_no}" "~{library_accession_no}"
    
    # run VADR
    v-annotate.pl \
    --mdir /opt/vadr/vadr-models/ \
    {genome_fasta} "~vadr_{genome_fasta}"

    # package everything for output
    tar -C "~{out_base}" -czvf "~{out_base}.vadr.tar.gz" .
    
    # prep alerts into a tsv file for parsing
    cat "~{out_base}/~{out_base}.vadr.alt.list" | cut -f 5 | tail -n +2 \
      > "~{out_base}.vadr.alerts.tsv"
    cat "~{out_base}.vadr.alerts.tsv" | wc -l > NUM_ALERTS
  >>>
  output{ 
    File                 feature_tbl = "~{out_base}/~{out_base}.vadr.pass.tbl"
    Int                  num_alerts  = read_int("NUM_ALERTS")
    File                 alerts_list = "~{out_base}/~{out_base}.vadr.alt.list"
    Array[Array[String]] alerts      = read_tsv("~{out_base}.vadr.alerts.tsv")
    File                 outputs_tgz = "~{out_base}.vadr.tar.gz"
    Boolean              pass        = num_alerts==0
    String               vadr_docker = docker
  }
    runtime {
        docker: docker
        memory:  mem_size + " GB"
        cpu: cpu
        disk: disk_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
   }
    
