version 1.0

task krakenuniq {
  meta {
    description: "Runs Krakenuniq classification"
  }

  input {
    Array[File] reads_unmapped_bam
    File        krakenuniq_db_tar_lz4  # {database.kdb,taxonomy}
    File        krona_taxonomy_db_tgz  # taxonomy.tab

    Int         machine_mem_gb = 320
    String      docker = "quay.io/broadinstitute/viral-classify:2.1.33.0" #skip-global-version-pin
  }

  Int disk_size = 750

  parameter_meta {
    reads_unmapped_bam: {
      description: "Reads to classify. May be unmapped or mapped or both, paired-end or single-end.",
      patterns: ["*.bam"],
      category: "required" }
    krakenuniq_db_tar_lz4: {
      description: "Pre-built Kraken database tarball.",
      patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"],
      category:"required"
    }
    krona_taxonomy_db_tgz: {
      description: "Krona taxonomy database containing a single file: taxonomy.tab, or possibly just a compressed taxonomy.tab",
      patterns: ["*.tab.zst", "*.tab.gz", "*.tab", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"],
      category: "required"
    }
  }

  command <<<
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    DB_DIR=$(mktemp -d --suffix _db)
    mkdir -p $DB_DIR/krakenuniq $DB_DIR/krona

    metagenomics --version | tee VERSION

    # decompress DB to $DB_DIR
    read_utils extract_tarball \
      "~{krakenuniq_db_tar_lz4}" $DB_DIR/krakenuniq \
      --loglevel=DEBUG
    # Support old db tar format
    if [ -d "$DB_DIR/krakenuniq/krakenuniq" ]; then
      mv $DB_DIR/krakenuniq/krakenuniq/* $DB_DIR/krakenuniq
    fi

    # unpack krona taxonomy.tab
    if [[ "~{krona_taxonomy_db_tgz}" == *.tar.* ]]; then
      read_utils extract_tarball \
        "~{krona_taxonomy_db_tgz}" $DB_DIR/krona \
        --loglevel=DEBUG &  # we don't need this until later
    else
      if [[ "~{krona_taxonomy_db_tgz}" == *.zst ]]; then
        cat "~{krona_taxonomy_db_tgz}" | zstd -d > $DB_DIR/krona/taxonomy.tab &
      elif [[ "~{krona_taxonomy_db_tgz}" == *.gz ]]; then
        cat "~{krona_taxonomy_db_tgz}" | pigz -dc > $DB_DIR/krona/taxonomy.tab &
      elif [[ "~{krona_taxonomy_db_tgz}" == *.bz2 ]]; then
        cat "~{krona_taxonomy_db_tgz}" | bzip -dc > $DB_DIR/krona/taxonomy.tab &
      else
        cp "~{krona_taxonomy_db_tgz}" $DB_DIR/krona/taxonomy.tab &
      fi
    fi

    # prep input and output file names
    OUT_READS=fnames_outreads.txt
    OUT_REPORTS=fnames_outreports.txt
    OUT_BASENAME=basenames_reports.txt
    for bam in "~{sep='" "' reads_unmapped_bam}"; do
      echo "$(basename $bam .bam).krakenuniq-reads.txt.gz" >> $OUT_READS
      echo "$(basename $bam .bam)" >> $OUT_BASENAME
      echo "$(basename $bam .bam).krakenuniq-summary_report.txt" >> $OUT_REPORTS
    done

    # execute on all inputs and outputs serially, but with a single
    # database load into ram
    metagenomics krakenuniq \
      $DB_DIR/krakenuniq \
      "~{sep='" "' reads_unmapped_bam}" \
      --outReads $(cat $OUT_READS) \
      --outReport $(cat $OUT_REPORTS) \
      --loglevel=DEBUG

    wait # for krona_taxonomy_db_tgz to download and extract
    # Support old db tar format
    if [ -d $DB_DIR/krona/taxonomy ]; then
      mv $DB_DIR/krona/taxonomy/* $DB_DIR/krona
    fi

    # run single-threaded krona on up to nproc samples at once
    parallel -I ,, \
      "metagenomics krona \
        ,,.krakenuniq-summary_report.txt \
        $DB_DIR/krona \
        ,,.krakenuniq-krona.html \
        --sample_name ,, \
        --noRank --noHits --inputType krakenuniq \
        --loglevel=DEBUG" \
      ::: $(cat $OUT_BASENAME)

    # merge all krona reports
    ktImportKrona -o krakenuniq.krona.combined.html *.krakenuniq-krona.html

    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
  >>>

  output {
    Array[File] krakenuniq_classified_reads = glob("*.krakenuniq-reads.txt.gz")
    Array[File] krakenuniq_summary_reports  = glob("*.krakenuniq-summary_report.txt")
    Array[File] krona_report_html           = glob("*.krakenuniq-krona.html")
    File        krona_report_merged_html    = "krakenuniq.krona.combined.html"

    Int         max_ram_gb                  = ceil(read_float("MEM_BYTES")/1000000000)

    String      viralngs_version            = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    cpu: 32
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem3_ssd1_v2_x48"
    preemptible: 0
  }
}

task build_krakenuniq_db {
  input {
    File     genome_fastas_tarball
    File     taxonomy_db_tarball
    String   db_basename

    Boolean? subsetTaxonomy
    Int?     minimizerLen
    Int?     kmerLen
    Int?     maxDbSize
    Int?     zstd_compression_level

    Int      machine_mem_gb = 240
    String   docker = "quay.io/broadinstitute/viral-classify:2.1.33.0" #skip-global-version-pin
  }

  Int disk_size = 750

  command <<<
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    TAXDB_DIR=$(mktemp -d --suffix _taxdb)
    FASTAS_DIR=$(mktemp -d --suffix fasta)
    DB_DIR="$TMPDIR/~{db_basename}"
    mkdir -p $DB_DIR

    metagenomics --version | tee VERSION

    # decompress input tarballs
    read_utils extract_tarball \
      "~{genome_fastas_tarball}" $FASTAS_DIR \
      --loglevel=DEBUG
    read_utils extract_tarball \
      "~{taxonomy_db_tarball}" $TAXDB_DIR \
      --loglevel=DEBUG

    # build database
    metagenomics krakenuniq_build \
      $DB_DIR --library $FASTAS_DIR --taxonomy $TAXDB_DIR \
      ~{if select_first([subsetTaxonomy, false]) then '--subsetTaxonomy=' else ''} \
      ~{'--minimizerLen=' + minimizerLen} \
      ~{'--kmerLen=' + kmerLen} \
      ~{'--maxDbSize=' + maxDbSize} \
      --clean \
      --loglevel=DEBUG

    # tar it up
    tar -c -C $DB_DIR . | zstd ~{"-" + zstd_compression_level} > "~{db_basename}.tar.zst"
  >>>

  output {
    File   krakenuniq_db    = "~{db_basename}.tar.zst"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    cpu: 32
    dx_instance_type: "mem3_ssd1_v2_x32"
    preemptible: 0
  }
}

task kraken2 {
  meta {
    description: "Runs Kraken2 classification"
  }

  input {
    File   reads_bam
    File   kraken2_db_tgz         # {database.kdb,taxonomy}
    File   krona_taxonomy_db_tgz  # taxonomy.tab
    Float? confidence_threshold = 0.05
    Int?   min_base_qual

    Int    machine_mem_gb = 90
    String docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }

  parameter_meta {
    reads_bam: {
      description: "Reads or contigs to classify. May be unmapped or mapped or both, paired-end or single-end.",
      patterns: ["*.bam", "*.fasta"]
    }
    kraken2_db_tgz: {
      description: "Pre-built Kraken database tarball containing three files: hash.k2d, opts.k2d, and taxo.k2d.",
      patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    krona_taxonomy_db_tgz: {
      description: "Krona taxonomy database containing a single file: taxonomy.tab, or possibly just a compressed taxonomy.tab",
      patterns: ["*.tab.zst", "*.tab.gz", "*.tab", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    confidence_threshold: {
      description: "Kraken2 confidence score threshold (0.0-1.0). See https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#confidence-scoring"
    }
    min_base_qual: {
      description: "Minimum base quality used in classification"
    }
  }

  String out_basename = basename(basename(reads_bam, '.bam'), '.fasta')

  # Disk autoscaling: BAM->FASTQ expansion is ~7-8x, plus kraken2 reads output (~1x input),
  # plus kraken2 database (1x localized tarball + 2x decompressed = 3x), plus overhead for krona and temp files.
  # Minimum 750GB to accommodate typical database sizes.
  # Note: GCP local SSDs must be allocated in pairs (2, 4, 8, 16, 24 × 375GB), so we round to 750GB multiples.
  Int disk_size_auto = ceil((8 * size(reads_bam, "GB") + 3 * size(kraken2_db_tgz, "GB") + 50) / 750.0) * 750
  Int disk_size = if disk_size_auto < 750 then 750 else disk_size_auto

  command <<<
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    DB_DIR=$(mktemp -d --suffix _db)
    mkdir -p $DB_DIR/kraken2 $DB_DIR/krona

    # decompress DB to $DB_DIR
    read_utils extract_tarball \
      "~{kraken2_db_tgz}" $DB_DIR/kraken2 \
      --loglevel=DEBUG
    du -hs $DB_DIR/kraken2

    # unpack krona taxonomy.tab
    if [[ "~{krona_taxonomy_db_tgz}" == *.tar.* ]]; then
      read_utils extract_tarball \
        "~{krona_taxonomy_db_tgz}" $DB_DIR/krona \
        --loglevel=DEBUG &  # we don't need this until later
    else
      if [[ "~{krona_taxonomy_db_tgz}" == *.zst ]]; then
        cat "~{krona_taxonomy_db_tgz}" | zstd -d > $DB_DIR/krona/taxonomy.tab &
      elif [[ "~{krona_taxonomy_db_tgz}" == *.gz ]]; then
        cat "~{krona_taxonomy_db_tgz}" | pigz -dc > $DB_DIR/krona/taxonomy.tab &
      elif [[ "~{krona_taxonomy_db_tgz}" == *.bz2 ]]; then
        cat "~{krona_taxonomy_db_tgz}" | bzip -dc > $DB_DIR/krona/taxonomy.tab &
      else
        cp "~{krona_taxonomy_db_tgz}" $DB_DIR/krona/taxonomy.tab &
      fi
    fi

    metagenomics --version | tee VERSION

    if [[ "~{reads_bam}" == *.bam ]]; then
        metagenomics kraken2 \
          $DB_DIR/kraken2 \
          "~{reads_bam}" \
          --outReads   "~{out_basename}".kraken2.reads.txt \
          --outReports "~{out_basename}".kraken2.report.txt \
          ~{"--confidence " + confidence_threshold} \
          ~{"--min_base_qual " + min_base_qual} \
          --loglevel=DEBUG
    else # fasta input file: call kraken2 directly
        kraken2 \
          --db $DB_DIR/kraken2 \
          "~{reads_bam}" \
          --output "~{out_basename}".kraken2.reads.txt \
          --report "~{out_basename}".kraken2.report.txt \
          ~{"--confidence " + confidence_threshold} \
          ~{"--min_base_qual " + min_base_qual}
    fi

    wait # for krona_taxonomy_db_tgz to download and extract
    pigz "~{out_basename}".kraken2.reads.txt &

    metagenomics krona \
      "~{out_basename}".kraken2.report.txt \
      $DB_DIR/krona \
      "~{out_basename}".kraken2.krona.html \
      --sample_name "~{out_basename}" \
      --noRank --noHits --inputType kraken2 \
      --loglevel=DEBUG

    wait # pigz reads.txt

    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
  >>>

  output {
    File   kraken2_reads_report   = "~{out_basename}.kraken2.reads.txt.gz"
    File   kraken2_summary_report = "~{out_basename}.kraken2.report.txt"
    File   krona_report_html      = "~{out_basename}.kraken2.krona.html"

    Int    max_ram_gb             = ceil(read_float("MEM_BYTES")/1000000000)

    String viralngs_version       = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    cpu: 16
    cpuPlatform: "Intel Ice Lake"
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TESs
    dx_instance_type: "mem3_ssd1_v2_x8"
    preemptible: 3
  }
}

task report_primary_kraken_taxa {
  meta {
    description: "Interprets a kraken (or kraken2 or krakenuniq) summary report file and emits the primary contributing taxa under a focal taxon of interest."
  }
  input {
    File          kraken_summary_report
    String        focal_taxon = "Viruses"

    String        docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }
  String out_basename = basename(kraken_summary_report, '.txt')
  Int disk_size = 50
  Int machine_mem_gb = 2

  command <<<
    set -e
    metagenomics taxlevel_plurality "~{kraken_summary_report}" "~{focal_taxon}" "~{out_basename}.ranked_focal_report.tsv"
    cat "~{out_basename}.ranked_focal_report.tsv" | head -2 | tail +2 > TOPROW
    cut -f 2 TOPROW > NUM_FOCAL
    cut -f 4 TOPROW > PCT_OF_FOCAL
    cut -f 7 TOPROW > NUM_READS
    cut -f 8 TOPROW > TAX_RANK
    cut -f 9 TOPROW > TAX_ID
    cut -f 10 TOPROW > TAX_NAME
  >>>

  output {
    String focal_tax_name = focal_taxon
    File   ranked_focal_report = "~{out_basename}.ranked_focal_report.tsv"
    Int    total_focal_reads = read_int("NUM_FOCAL")
    Float  percent_of_focal = read_float("PCT_OF_FOCAL")
    Int    num_reads = read_int("NUM_READS")
    String tax_rank = read_string("TAX_RANK")
    String tax_id = read_string("TAX_ID")
    String tax_name = read_string("TAX_NAME")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    cpu: 1
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TESs
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}

task filter_refs_to_found_taxa {
  meta {
    description: "Filters a taxid_to_ref_accessions_tsv to the set of taxa found in a focal_report."
  }
  input {
    File          taxid_to_ref_accessions_tsv
    File          focal_report_tsv
    File          taxdump_tgz
    Int           min_read_count = 100

    String        docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }
  String ref_basename = basename(taxid_to_ref_accessions_tsv, '.tsv')
  String hits_basename = basename(focal_report_tsv, '.tsv')
  Int disk_size = 50

  command <<<
    set -e
    mkdir -p taxdump
    read_utils extract_tarball "~{taxdump_tgz}" taxdump
    metagenomics filter_taxids_to_focal_hits "~{taxid_to_ref_accessions_tsv}" "~{focal_report_tsv}" taxdump ~{min_read_count} "~{ref_basename}-~{hits_basename}.tsv"
  >>>

  output {
    File   filtered_taxid_to_ref_accessions_tsv = "~{ref_basename}-~{hits_basename}.tsv"
  }

  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TESs
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 2
  }
}

task build_kraken2_db {
  meta {
    description: "Builds a custom kraken2 database. Outputs tar.zst tarballs of kraken2 database, associated krona taxonomy db, and an ncbi taxdump.tar.gz. Requires live internet access if any standard_libraries are specified or if taxonomy_db_tgz is absent."
  }

  input {
    String        db_basename
    File?         taxonomy_db_tgz
    Array[String] standard_libraries = [
                      "archaea", "bacteria", "plasmid",
                      "viral", "human", "fungi", "protozoa",
                      "UniVec_Core"]
    Array[File]   custom_libraries = []
    Boolean       protein = false

    Int?          kmerLen
    Int?          minimizerLen
    Int?          minimizerSpaces
    Int?          maxDbSize
    Int?          zstd_compression_level

    Int           machine_mem_gb = 100
    String        docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }

  Int disk_size = 750

  parameter_meta {
    db_basename: { description: "A descriptive string used in output filenames. Outputs will be called kraken2-<db_basename>.tar.zst, krona-<db_basename>.tar.zst, and taxdump-<db_basename>.tar.gz" }
    taxonomy_db_tgz: {
       description: "Optional tarball of kraken2 taxonomy database directory. Omitting this input will cause a fresh download from NCBI at the time of build.",
       patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    standard_libraries: {
      description: "A list of 'standard' kraken2 databases to include in this build. Including any values here will cause fresh downloads of data at the time of build. A list of acceptable names is available at https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#custom-databases"
    }
    custom_libraries: {
      description: "A list of 'custom' kraken2 databases to include in this build. Headers must be formatted as described in the kraken2 documentation. These are fastas or tarball collections of such fastas--multiple may be provided here.",
      patterns: ["*.fasta", "*.fasta.zst", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    protein: {
      description: "Build a protein (translated search) database. Default is nucleotide."
    }
    kmerLen: {
      description: "(k) K-mer length in bp/aa (Kraken2 defaults: 35 nt, 15 aa)"
    }
    minimizerLen: {
      description: "(l) Minimizer length in bp/aa (Kraken2 defaults: 31 nt, 12 aa)"
    }
    minimizerSpaces: {
      description: "(s) Number of bp/aa in minimizer that are ignored in comparisons (Kraken2 defaults: 7 nt, 0 aa)"
    }
  }

  command <<<
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    TAXDB_DIR=$(mktemp -d)
    FASTAS_DIR=$(mktemp -d)
    KRONA_DIR=$(mktemp -d)
    DB_DIR=$(mktemp -d)

    metagenomics --version | tee VERSION

    # prep input taxonomy db, if specified
    if [ -n "~{taxonomy_db_tgz}" ]; then
      read_utils extract_tarball \
        "~{taxonomy_db_tgz}" $TAXDB_DIR \
        --loglevel=DEBUG &
      TAX_INPUT_CMD="--tax_db=$TAXDB_DIR"
    else
      TAX_INPUT_CMD=""
    fi

    # prep input custom fastas, if specified
    CUSTOM_INPUT_CMD=""
    if [ -n "~{sep=' ' custom_libraries}" ]; then
      CUSTOM_INPUT_CMD="--custom_libraries "
      for TGZ in ~{sep=' ' custom_libraries}; do
        if [[ ($TGZ == *.tar.*) || ($TGZ == *.tgz) ]]; then
          read_utils extract_tarball \
            $TGZ $FASTAS_DIR \
            --loglevel=DEBUG &
        else
          if [[ $TGZ == *.zst ]]; then
            cat $TGZ | zstd -d > $FASTAS_DIR/$(basename $TGZ .zst) &
          elif [[ $TGZ == *.gz ]]; then
            cat $TGZ | pigz -dc > $FASTAS_DIR/$(basename $TGZ .gz) &
          elif [[ $TGZ == *.bz2 ]]; then
            cat $TGZ | bzip -dc > $FASTAS_DIR/$(basename $TGZ .bz2) &
          else
            CUSTOM_INPUT_CMD="$CUSTOM_INPUT_CMD $TGZ"
          fi
        fi
      done
      wait # wait for all decompressions to finish
      for FASTA in $FASTAS_DIR/*; do
        CUSTOM_INPUT_CMD="$CUSTOM_INPUT_CMD $FASTA"
      done
    fi

    # prep standard libraries, if specified
    STD_INPUT_CMD=""
    if [ -n "~{sep=' ' standard_libraries}" ]; then
      STD_INPUT_CMD="--standard_libraries ~{sep=' ' standard_libraries}"
    fi

    # build kraken2 database
    wait # wait for all decompressions to finish
    metagenomics kraken2_build \
      $DB_DIR \
      $TAX_INPUT_CMD \
      $STD_INPUT_CMD \
      $CUSTOM_INPUT_CMD \
      --taxdump_out "taxdump-~{db_basename}.tar.gz" \
      ~{true='--protein' false='' protein} \
      ~{'--kmerLen=' + kmerLen} \
      ~{'--minimizerLen=' + minimizerLen} \
      ~{'--minimizerSpaces=' + minimizerSpaces} \
      ~{'--maxDbSize=' + maxDbSize} \
      --loglevel=DEBUG
    tar -c -C $DB_DIR . | zstd ~{"-" + zstd_compression_level} > "kraken2-~{db_basename}.tar.zst" &

    # build matching krona db
    metagenomics krona_build \
      $KRONA_DIR --taxdump_tar_gz "taxdump-~{db_basename}.tar.gz"
    cat $KRONA_DIR/taxonomy.tab | zstd -19 > "krona-~{db_basename}-taxonomy.tab.zst"

    wait # tar/zst of kraken2 db
  >>>

  output {
    File   kraken2_db       = "kraken2-~{db_basename}.tar.zst"
    File   taxdump_tgz      = "taxdump-~{db_basename}.tar.gz"
    File   krona_db         = "krona-~{db_basename}-taxonomy.tab.zst"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    cpu: 16
    dx_instance_type: "mem3_ssd1_v2_x16"
    preemptible: 0
  }
}

task blastx {
  meta {
    description: "Runs BLASTx classification"
  }

  input {
    File   contigs_fasta
    File   blast_db_tgz
    File   krona_taxonomy_db_tgz

    Int    machine_mem_gb = 8
    String docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }

  parameter_meta {
    contigs_fasta: {
      description: "Sequences to classify. Use for a small number of longer query sequences (e.g. contigs)",
      patterns: ["*.fasta"] }
    blast_db_tgz: {
      description: "Pre-built BLAST database tarball containing an indexed blast database named 'nr'",
      patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    krona_taxonomy_db_tgz: {
      description: "Krona taxonomy database: a tarball containing a taxonomy.tab file as well as accession to taxid mapping (a kraken-based taxonomy database will not suffice).",
      patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
  }

  String out_basename = basename(contigs_fasta, '.fasta')
  Int    disk_size = 375

  command <<<
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    DB_DIR=$(mktemp -d --suffix _db)
    mkdir -p $DB_DIR/blast $DB_DIR/krona

    # decompress DB to $DB_DIR
    read_utils extract_tarball \
      "~{blast_db_tgz}" $DB_DIR/blast \
      --loglevel=DEBUG

    # unpack krona taxonomy database
    read_utils extract_tarball \
      "~{krona_taxonomy_db_tgz}" $DB_DIR/krona \
      --loglevel=DEBUG &  # we don't need this until later

    blastx -version | tee VERSION

    blastx \
      -query "~{contigs_fasta}" \
      -db $DB_DIR/blast/nr \
      -out "~{out_basename}.blastx.contigs.txt" \
      -outfmt 7 \
      -num_threads $(nproc)

    wait # for krona_taxonomy_db_tgz to download and extract

    ktImportBLAST \
      -i -k \
      -tax $DB_DIR/krona \
      -o "~{out_basename}.blastx.krona.html" \
      "~{out_basename}.blastx.contigs.txt","~{out_basename}"

    pigz "~{out_basename}".blastx.contigs.txt
  >>>

  output {
    File    blast_report       = "~{out_basename}.blastx.contigs.txt.gz"
    File    krona_report_html  = "~{out_basename}.blastx.krona.html"
    String  blastx_version     = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    cpu: 32
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x36"
    preemptible: 1
  }
}

task krona {
  input {
    Array[File]+ reports_txt_gz
    File         krona_taxonomy_db_tgz
    String       out_basename = basename(basename(reports_txt_gz[0], '.gz'), '.txt')

    String?      input_type
    Int?         query_column
    Int?         taxid_column
    Int?         score_column
    Int?         magnitude_column

    Int          machine_mem_gb = 3
    String       docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }

  Int disk_size = 50

  command <<<
    set -ex -o pipefail
    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    DB_DIR=$(mktemp -d --suffix _db)
    mkdir -p $DB_DIR/krona

    metagenomics --version | tee VERSION

    # unpack krona taxonomy.tab
    if [[ "~{krona_taxonomy_db_tgz}" == *.tar.* ]]; then
      read_utils extract_tarball \
        "~{krona_taxonomy_db_tgz}" $DB_DIR/krona \
        --loglevel=DEBUG
    else
      if [[ "~{krona_taxonomy_db_tgz}" == *.zst ]]; then
        cat "~{krona_taxonomy_db_tgz}" | zstd -d > $DB_DIR/krona/taxonomy.tab
      elif [[ "~{krona_taxonomy_db_tgz}" == *.gz ]]; then
        cat "~{krona_taxonomy_db_tgz}" | pigz -dc > $DB_DIR/krona/taxonomy.tab
      elif [[ "~{krona_taxonomy_db_tgz}" == *.bz2 ]]; then
        cat "~{krona_taxonomy_db_tgz}" | bzip -dc > $DB_DIR/krona/taxonomy.tab
      else
        cp "~{krona_taxonomy_db_tgz}" $DB_DIR/krona/taxonomy.tab
      fi
    fi

    metagenomics krona \
      "~{sep='" "' reports_txt_gz}" \
      $DB_DIR/krona \
      "~{out_basename}.html" \
      ~{'--inputType=' + input_type} \
      ~{'--queryColumn=' + query_column} \
      ~{'--taxidColumn=' + taxid_column} \
      ~{'--scoreColumn=' + score_column} \
      ~{'--magnitudeColumn=' + magnitude_column} \
      --noRank --noHits \
      --loglevel=DEBUG
  >>>

  output {
    File   krona_report_html = "~{out_basename}.html"
    String viralngs_version  = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    cpu: 1
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd2_v2_x2"
  }
}

task krona_merge {
  input {
    Array[File] krona_reports
    String      out_basename

    Int         machine_mem_gb = 3
    String      docker = "biocontainers/krona:v2.7.1_cv1"
  }

  Int disk_size = 50

  command <<<
    set -ex -o pipefail
    ktImportKrona | head -2 | tail -1 | cut -f 2-3 -d ' ' | tee VERSION
    ktImportKrona -o "~{out_basename}.html" ~{sep=' ' krona_reports}
  >>>

  output {
    File   krona_report_html = "~{out_basename}.html"
    String krona_version     = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    cpu: 1
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd2_v2_x2"
  }
}

task filter_bam_to_taxa {
  input {
    File           classified_bam
    File           classified_reads_txt_gz
    File           ncbi_taxonomy_db_tgz # nodes.dmp names.dmp
    Array[String]? taxonomic_names
    Array[Int]?    taxonomic_ids
    Int?           minimum_hit_groups
    Boolean        withoutChildren = false
    Boolean        exclude_taxa = false
    String         out_filename_suffix = "filtered"

    Int            machine_mem_gb = 8
    String         docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }

  String out_basename = basename(classified_bam, ".bam") + "." + out_filename_suffix
  Int disk_size = ceil((2 * size(classified_bam, "GB") + 100) / 375.0) * 375

  command <<<
    set -ex -o pipefail
    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi

    # decompress taxonomy DB to CWD
    read_utils extract_tarball \
      "~{ncbi_taxonomy_db_tgz}" . \
      --loglevel=DEBUG
    if [ -d "taxonomy" ]; then mv taxonomy/* .; fi

    touch taxfilterargs
    TAXNAMELIST="~{write_lines(select_first([taxonomic_names, []]))}"
    if [ -n "$(cat $TAXNAMELIST)" ]; then
      echo "--taxNames" >> taxfilterargs
      cat $TAXNAMELIST >> taxfilterargs
      echo "" >> taxfilterargs # cromwell write_lines lacks a final newline, so add one manually
    fi

    TAXIDLIST="~{write_lines(select_first([taxonomic_ids, []]))}"
    if [ -n "$(cat $TAXIDLIST)" ]; then
      echo "--taxIDs" >> taxfilterargs
      cat $TAXIDLIST >> taxfilterargs
      echo "" >> taxfilterargs # cromwell write_lines lacks a final newline, so add one manually
    fi

    echo "taxfilterargs:"
    cat taxfilterargs

    metagenomics --version | tee VERSION

    samtools view -c "~{classified_bam}" | tee classified_taxonomic_filter_read_count_pre &

    cat taxfilterargs | grep . | xargs -d '\n' metagenomics filter_bam_to_taxa \
      "~{classified_bam}" \
      "~{classified_reads_txt_gz}" \
      "~{out_basename}.bam" \
      nodes.dmp \
      names.dmp \
      ~{true='--exclude' false='' exclude_taxa} \
      ~{true='--without-children' false='' withoutChildren} \
      ~{'--minimum_hit_groups=' + minimum_hit_groups} \
      --out_count COUNT \
      --loglevel=DEBUG

    samtools view -c "~{out_basename}.bam" | tee classified_taxonomic_filter_read_count_post
    wait

    cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
    cat /proc/loadavg | cut -f 3 -d ' ' > LOAD_15M
    set +o pipefail
    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi; } > MEM_BYTES
  >>>

  output {
    File    bam_filtered_to_taxa                        = "~{out_basename}.bam"
    Int     reads_matching_taxa                         = read_int("COUNT")
    Int     classified_taxonomic_filter_read_count_pre  = read_int("classified_taxonomic_filter_read_count_pre")
    Int     classified_taxonomic_filter_read_count_post = read_int("classified_taxonomic_filter_read_count_post")

    Int     max_ram_gb                                  = ceil(read_float("MEM_BYTES")/1000000000)
    Int     runtime_sec                                 = ceil(read_float("UPTIME_SEC"))
    Int     cpu_load_15min                              = ceil(read_float("LOAD_15M"))
    String  viralngs_version                            = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    disks: "local-disk ~{disk_size} SSD"
    disk: "~{disk_size} GB" # TES
    cpu: 8
    dx_instance_type: "mem1_ssd1_v2_x8"
    preemptible: 3
  }
}

task kaiju {
  input {
    File   reads_unmapped_bam
    File   kaiju_db_lz4  # <something>.fmi
    File   ncbi_taxonomy_db_tgz # taxonomy/{nodes.dmp, names.dmp}
    File   krona_taxonomy_db_tgz  # taxonomy/taxonomy.tab

    Int    machine_mem_gb = 100
    String docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }

  String   input_basename = basename(reads_unmapped_bam, ".bam")
  Int      disk_size = 375

  command <<<
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    DB_DIR=$(mktemp -d --suffix _db)
    mkdir -p $DB_DIR/kaiju $DB_DIR/krona $DB_DIR/taxonomy

    lz4 -dc "~{kaiju_db_lz4}" > $DB_DIR/kaiju/kaiju.fmi

    read_utils extract_tarball \
      "~{ncbi_taxonomy_db_tgz}" $DB_DIR/taxonomy \
      --loglevel=DEBUG
    # Support old db tar format
    if [ -d "$DB_DIR/taxonomy/taxonomy" ]; then
      mv $DB_DIR/taxonomy/taxonomy/* $DB_DIR/taxonomy
    fi

    read_utils extract_tarball \
      "~{krona_taxonomy_db_tgz}" $DB_DIR/krona \
      --loglevel=DEBUG

    metagenomics --version | tee VERSION

    # classify contigs
    metagenomics kaiju \
      "~{reads_unmapped_bam}" \
      $DB_DIR/kaiju/kaiju.fmi \
      $DB_DIR/taxonomy \
      "~{input_basename}.kaiju.summary_report.txt" \
      --outReads "~{input_basename}.kaiju.reads.txt.gz" \
      --loglevel=DEBUG

    # run krona
    metagenomics krona \
      "~{input_basename}.kaiju.summary_report.txt" \
      $DB_DIR/krona \
      "~{input_basename}.kaiju-krona.html" \
      --inputType kaiju \
      --noRank --noHits \
      --loglevel=DEBUG
  >>>

  output {
    File   kaiju_report      = "~{input_basename}.kaiju-summary_report.txt"
    File   kaiju_reads       = "~{input_basename}.kaiju-reads.txt.gz"
    File   krona_report_html = "~{input_basename}.kaiju-krona.html"
    String viralngs_version  = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    cpu: 16
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem3_ssd1_v2_x16"
  }
}

task kallisto {
  meta {
    description: "Runs kb count classification tool"
  }

  input {
    File     reads_bam
    File     kb_index
    File     t2g
    Int      kmer_size=31

    String   technology
    String   parity
    Boolean  h5ad=false
    Boolean  loom=false
    Boolean  protein=false

    String   docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }

  parameter_meta {
    reads_bam: {
      description: "Reads to classify. Must be un-mapped reads, paired-end or single-end.",
      patterns: ["*.bam", "*.fastq", "*.fastq.gz"]
    }
    kb_index: {
      description: "Pre-built kb index tarball",
      patterns: ["*.idx", "*.index"]
    }
    t2g: {
      description: "Transcript-to-gene mapping file. Two-column TSV with transcript IDs in first column and gene IDs in second column.",
      patterns: ["*.tsv", "*.txt"]
    }
    kmer_size: {
      description: "K-mer size used in kb index. Must match the k-mer size used to build the index. Default is 31."
    }
    technology: {
      description: "Single-cell technology used  {10xv2, 10xv3, dropseq, inDrops, seqwell, smartseq2, bulk}."
    }
    h5ad: {
      description: "Output an h5ad file (requires scanpy). Default is false"
    }
    loom: {
      description: "Output a loom file (requires loompy). Default is false"
    }
    protein: {
      description: "Indicates that the kb index is built from protein sequences. Default is false"
    } 
  }

  String out_basename=basename(basename(reads_bam, '.bam'), '.fasta')

  command <<<
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi

    printf "%s\n" "$(kb -h 2>&1 | grep "kb_python" | tee VERSION)"
    kallisto version | tee -a VERSION
    metagenomics --version | tee -a VERSION

    paste -sd ';' VERSION | sed 's/;/; /g' > VERSION.tmp && mv VERSION.tmp VERSION

    if [[ ~{reads_bam} == *.bam ]]; then
        metagenomics kb \
          "~{reads_bam}" \
          --index "~{kb_index}" \
          --t2g "~{t2g}" \
          --kmer_len ~{kmer_size} \
          --technology ~{technology} \
          --parity ~{parity} \
          ~{if h5ad then "--h5ad" else ""} \
          ~{if loom then "--loom" else ""} \
          ~{true='--protein' false='' protein} \
          --out_dir ~{out_basename}_count \
          --loglevel=DEBUG
    else # we have a single-ended fastq file so just call it directly
        kb count \
          --kallisto kallisto \
          -t `nproc` \
          -k ~{kmer_size} \
          --parity single \
          -i "~{kb_index}" \
          -g "~{t2g}" \
          -o "~{out_basename}_count" \
          -x ~{technology} \
          ~{if h5ad then "--h5ad" else ""} \
          ~{if loom then "--loom" else ""} \
          ~{true='--aa' false='' protein} \
          "~{reads_bam}"
    fi

    # Since we are running this file-by-file we need to add our sample name to the matrix.cells file
    bn=$(basename "~{reads_bam}")
    sample_name=$(echo "$bn" | cut -d'.' -f1)
    echo "$sample_name" > "~{out_basename}_count/matrix.cells"

    tar -c -C "~{out_basename}_count" . | zstd > "~{out_basename}_kb_count.tar.zst"
  >>>

  output {
    File    kb_count_tar  = "~{out_basename}_kb_count.tar.zst"
    String  viralngs_version    = read_string("VERSION")
  }

  runtime {
    docker: "~{docker}"
    memory: "32 GB"
    cpu: 16
    disks: "local-disk 350 LOCAL"
    dx_instance_type: "mem3_ssd1_v2_x16"
    preemptible: 2
  }
}

task build_kallisto_db {
  meta {
    description: "Builds a custom kb index from provided input FASTA file."
  }

  input {
    String      out_basename
    File        reference_sequences
    Boolean     protein=false

    Int        kmer_size=31
    String?     workflow_type

    String      docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }

  parameter_meta {
    out_basename: { description: "A descriptive string used in output index filename. Output will be called <out_basename>.idx" }
    reference_sequences: {
      description: "FASTA file of reference sequences to index.",
      patterns: ["*.fasta", "*.fa", "*.fa.gz", "*.fasta.gz"]
    }
    protein: {
      description: "Generate an index from a FASTA file containing protein sequences. Default is nucleotide."
    }
    kmer_size: {
      description: "K-mer size to use in index. Default is 31."
    }
    workflow_type: {
      description: "Type of index to create {standard, nac, kite, custom}. Default is 'standard'."
    }
  }

  command <<<
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi

    printf "kb_python %s\n" "$(kb -h 2>&1 | grep "kb_python" | cut -d" " -f2)" | tee VERSION
    kallisto version | tee -a VERSION    
    metagenomics --version | tee -a VERSION

    if [[ ~{reference_sequences} == *.gz ]]; then
      gunzip -c "~{reference_sequences}" > reference_sequences.fasta
      REF_FASTA=reference_sequences.fasta
    else
      REF_FASTA="~{reference_sequences}"
    fi

    # build kb database
    metagenomics kb_build \
      ~{true='--protein' false='' protein} \
      --kmer_len=~{kmer_size} \
      ~{if defined(workflow_type) then "--workflow=" + workflow_type else ""} \
      --index="~{out_basename}.idx" \
      $REF_FASTA \
      --loglevel=DEBUG


       tar -c "~{out_basename}.idx" | zstd > "~{out_basename}.idx.tar.zst"
  >>>

  output {
    File        kb_index   = "~{out_basename}.idx.tar.zst"
    String      viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "~{docker}"
    memory: "32 GB"
    disks: "local-disk 750 LOCAL"
    cpu: 16
    dx_instance_type: "mem3_ssd1_v2_x16"
    preemptible: 0
  }
}

task kallisto_extract {
  meta {
    description: "Extracts reads that pseudoalign to a kb index"
  }

  input {
    File            reads_bam
    File            kb_index
    File            t2g
    Array[String]?  target_ids
    File?           h5ad_file
    Int             threads=8
    Int             threshold=1
    Boolean         protein=false

     String          docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }

  parameter_meta {
    reads_bam: {
      description: "Reads used by kb count to psuedoalign. Must be un-mapped or single-end.",
      patterns: ["*.bam", "*.fastq", "*.fastq.gz"]
    }
    kb_index: {
      description: "kb index used to psuedoalign reads.",
      patterns: ["*.index"]
    }
    t2g: {
      description: "Transcript-to-gene mapping file. Two-column TSV with transcript IDs in first column and gene IDs in second column.",
      patterns: ["*.tsv", "*.txt"]
    }
    target_ids: {
      description: "List of target transcript or gene IDs to extract reads for."
    }
    h5ad_file: {
      description: "If no target IDs are provided, an h5ad file may be provided from which to extract target IDs. The h5ad file must contain a 'gene_ids' column in the .var dataframe.",
      patterns: ["*.h5ad"]
    }
    protein: {
      description: "Input sequences contain amino acid sequences."
    }
    threshold: {
      description: "Minimum number of pseudoalignments to a target ID for a read to be extracted. Default is 1."
    }
    threads: {
      description: "Number of threads to use. Default is 8."
    }
  }

  String out_basename = sub(sub(sub(sub(basename(reads_bam), "\\.bam$", ""), "\\.gz$", ""), "\\.fastq$", ""), "\\.fq$", "")

  command <<<
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi

    printf "%s\n" "$(kb -h 2>&1 | grep "kb_python" | tee VERSION)"
    kallisto version | tee -a VERSION
    metagenomics --version | tee -a VERSION

    paste -sd ';' VERSION | sed 's/;/; /g' > VERSION.tmp && mv VERSION.tmp VERSION

    # Determine source of target IDs
    TARGET_IDS="~{sep=' ' select_first([target_ids, []])}"
    if [ -z "$TARGET_IDS" ]; then
      echo "No target IDs provided, will attempt to extract from h5ad!" >&2
      TARGET_SOURCE="--h5ad ~{h5ad_file}"
    else
      TARGET_SOURCE="--targets $TARGET_IDS"
    fi

    metagenomics kb_extract \
      "~{reads_bam}" \
      --index "~{kb_index}" \
      --t2g "~{t2g}" \
      --out_dir "~{out_basename}_extract" \
      ~{if protein then "--protein" else ""} \
      --threshold ~{threshold} \
      $TARGET_SOURCE \
      --loglevel=DEBUG

    tar -c -C "~{out_basename}_extract" . | zstd > "~{out_basename}_kb_extract.tar.zst"
  >>>

  output {
    File    kb_extract_tar    = "~{out_basename}_kb_extract.tar.zst"
    String  viralngs_version  = read_string("VERSION")
  }

  runtime {
    docker: "~{docker}"
    memory: "32 GB"
    cpu: 16
    disks: "local-disk 350 LOCAL"
    dx_instance_type: "mem3_ssd1_v2_x16"
    preemptible: 2
  }
}

task report_primary_kallisto_taxa {
  meta {
    description: "Interprets kb count output file and emits the primary contributing taxa under a focal taxon of interest."
  }
  input {
    File          kb_count_tar
    File          id_to_taxon_map
    String        focal_taxon = "Viruses"

    String        docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }
  String out_basename = sub(basename(kb_count_tar, ".tar.zst"), "_kb_count", "")
  Int disk_size = 200
  Int machine_mem_gb = 16

  command <<<
    set -e

    printf "%s\n" "$(kb -h 2>&1 | grep "kb_python" | tee VERSION)"
    kallisto version | tee -a VERSION
    metagenomics --version | tee -a VERSION

    paste -sd ';' VERSION | sed 's/;/; /g' > VERSION.tmp && mv VERSION.tmp VERSION

    metagenomics kb_top_taxa \
      "~{kb_count_tar}" \
      "~{out_basename}.ranked_focal_report.tsv" \
      --id-to-tax-map "~{id_to_taxon_map}" \
      --target-taxon "~{focal_taxon}"
    cat "~{out_basename}.ranked_focal_report.tsv" | head -2 | tail +2 > TOPROW
    cut -f 2 TOPROW > NUM_FOCAL           # focal_taxon_count
    cut -f 7 TOPROW > PCT_OF_FOCAL        # pct_of_focal
    cut -f 6 TOPROW > NUM_READS           # hit_reads
    cut -f 3 TOPROW > TAX_ID              # palmdb_id
    cut -f 4 TOPROW > TAX_NAME            # hit_id (species name)
    echo "" > TAX_RANK                    # Not provided by kb_top_taxa
  >>>

  output {
    File   ranked_focal_report = "~{out_basename}.ranked_focal_report.tsv"
    String focal_tax_name = focal_taxon
    Int    total_focal_reads = read_int("NUM_FOCAL")
    Float  percent_of_focal = read_float("PCT_OF_FOCAL")
    Int    num_reads = read_int("NUM_READS")
    String tax_rank = read_string("TAX_RANK")
    String tax_id = read_string("TAX_ID")
    String tax_name = read_string("TAX_NAME")
  }

  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    cpu: 16
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TESs
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 2
    maxRetries: 2
  }
}

task kallisto_merge_h5ads {
  meta {
    description: "Merges multiple kb count output tarballs into a single .h5ad file with sample metadata."
  }

  input {
    Array[File]     in_count_tars
    String          out_basename

    String          docker = "quay.io/broadinstitute/viral-ngs:3.0.10-classify"
  }

  parameter_meta {
    in_count_tars: {
      description: "List of kb count output tarballs to merge. Each tarball should contain counts_unfiltered/*.h5ad and matrix.cells.",
      patterns: ["*.tar.zst", "*.tar.gz"]
    }
    out_basename: {
      description: "Basename for output merged .h5ad file. Output will be named <out_basename>.h5ad."
    }

  }

  command <<<
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi

    printf "%s\n" "$(kb -h 2>&1 | grep "kb_python" | tee VERSION)"
    kallisto version | tee -a VERSION
    metagenomics --version | tee -a VERSION

    paste -sd ';' VERSION | sed 's/;/; /g' > VERSION.tmp && mv VERSION.tmp VERSION

    metagenomics kb_merge_h5ads \
      "~{sep='" "' in_count_tars}" \
      --out-h5ad "~{out_basename}.h5ad" \
      --loglevel=DEBUG
   >>>

  output {
    File    kb_merged_h5ad        = "~{out_basename}.h5ad"
    String  viralngs_version      = read_string("VERSION")
  }

  runtime {
    docker: "~{docker}"
    memory: "16 GB"
    cpu: 16
    disks: "local-disk 350 LOCAL"
    dx_instance_type: "mem3_ssd1_v2_x16"
    preemptible: 2
  }
}

task classify_virnucpro {
  meta {
    description: "Runs VirNucPro deep learning viral classification on unmapped reads using GPU acceleration."
  }

  input {
    File      reads_input
    Int       expected_length = 300

    Boolean   parallel=false
    Boolean   persist_model=false
    Boolean   use_gpu=false
    Int       batch_size=256
    Int?      esm_batch_size
    Int?      dnabert_batch_size

    Boolean   resume = false

    String?   accelerator_type
    Int?      accelerator_count
    String?   gpu_type
    Int?      gpu_count
    String?   vm_size
    String    docker = "ghcr.io/broadinstitute/virnucpro-cuda:1.0.9"
  }

  parameter_meta {
    reads_input: {
      description: "Reads to classify. Must be unmapped reads, paired-end or single-end. Accepts BAM or FASTA formats.",
      patterns: ["*.bam", "*.fasta", "*.fa", "*.fasta.gz", "*.fa.gz"]
    }
    expected_length: {
      description: "Expected sequence length in bp. Must be 300 or 500 to match bundled models (300_model.pth, 500_model.pth). VirNucPro silently accepts other values but produces invalid predictions."
    }
    parallel: {
      description: "If true, enables parallel data loading to improve throughput. May increase memory usage."
    }
    persist_model: {
      description: "If true, keeps the model loaded in GPU memory between batches to improve throughput. May increase GPU memory usage."
    }
    use_gpu: {
      description: "If true, uses GPU acceleration for inference. Requires a GPU-enabled runtime."
    }
    batch_size: {
      description: "Number of sequences to process in each batch. Larger batch sizes may improve throughput but increase memory usage."
    }
    esm_batch_size: {
      description: "Batch size specifically for the ESM2 model component. If not provided, defaults to the value of batch_size."
    }
    dnabert_batch_size: {
      description: "Batch size specifically for the DNABERT model component. If not provided, defaults to the value of batch_size."
    }
    resume: {
      description: "If true, enables checkpoint-based resume for preempted or interrupted runs."
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

  String basename = sub(basename(reads_input), "\\.(bam|fasta\\.gz|fa\\.gz|fasta|fa)$", "")

  command <<<
    set -ex -o pipefail

    export TMPDIR=/tmp
    export TEMP=/tmp
    export TMP=/tmp

    /opt/virnucpro_cli.py "~{reads_input}" "~{basename}.virnucpro.tsv" --expected-length ~{expected_length} \
      ~{true='--use-gpu' false='' use_gpu} \
      ~{true='--parallel' false='' parallel} \
      ~{true='--persistent-models' false='' persist_model} \
      --batch-size ~{batch_size} \
      ~{'--esm-batch-size=' + esm_batch_size} \
      ~{'--dnabert-batch-size=' + dnabert_batch_size} \
      --threads $(nproc) \
      ~{true='--resume' false='' resume}

    # Tarball both output files
    tar -czf "~{basename}.virnucpro.tgz" \
      "~{basename}.virnucpro.tsv" \
      "~{basename}.virnucpro_highestscore.csv"
  >>>

  output {
    File virnuc_pro_scores = "~{basename}.virnucpro.tgz"
  }

  # GPU multi-platform support: ALL platform attributes required (GCP: acceleratorType/acceleratorCount, Terra: gpuType/gpuCount, DNAnexus: gpu/dx_instance_type, Azure: vm_size). Missing attributes cause silent CPU fallback.
  runtime {
    docker: docker
    dockerRunOptions: "--tmpfs /dev/shm:rw,nosuid,nodev,size=16g"
    memory: "30 GB"
    cpu: 8
    disks: "local-disk 120 SSD"
    disk: "120 GB"
    gpu: true
    dx_instance_type: "mem1_ssd1_gpu2_x8"
    dx_timeout: "6H"
    acceleratorType: select_first([accelerator_type, "nvidia-tesla-t4"])
    acceleratorCount: select_first([accelerator_count, 1])
    gpuType: select_first([gpu_type, "nvidia-tesla-t4"])
    gpuCount: select_first([gpu_count, 1])
    vm_size: select_first([vm_size, "Standard_NC6"])
    maxRetries: 1
  }
}

task classify_virnucpro_contigs {
  meta {
    description: "Classify contigs as viral or non-viral based on confidence-weighted delta scores from VirNucPro chunk-level predictions."
  }

  input {
    File    virnucpro_scores_tsv

    Float   min_viral_prop    = 0.1
    Float   min_nonviral_prop = 0.1
    Int     min_chunks        = 5
    String  id_col            = "Modified_ID"
    String  id_pattern        = "(NODE_\\d+)"

    String  docker            = "quay.io/broadinstitute/py3-bio:0.1.3"
  }

  parameter_meta {
    virnucpro_scores_tsv: {
      description: "TSV of chunk-level VirNucPro scores with columns: Modified_ID, max_score_0, max_score_1. Also accepts a tar.zst/tar.gz tarball containing one such file.",
      patterns: ["*.tsv", "*.tar.zst", "*.tar.gz"],
      category: "required"
    }
    min_viral_prop: {
      description: "Minimum proportion of confident viral chunks (among non-ambiguous) required for viral call.",
      category: "advanced"
    }
    min_nonviral_prop: {
      description: "Minimum proportion of confident non-viral chunks (among non-ambiguous) required for non-viral call.",
      category: "advanced"
    }
    min_chunks: {
      description: "Minimum number of chunks for high/moderate confidence tiers. Contigs with fewer chunks are downgraded.",
      category: "advanced"
    }
    id_col: {
      description: "Column name containing contig IDs in the input TSV.",
      category: "advanced"
    }
    id_pattern: {
      description: "Regex pattern to extract contig group ID from the ID column. Default extracts NODE_N from SPAdes contig names.",
      category: "advanced"
    }
  }

  String out_filename = basename(virnucpro_scores_tsv, ".tsv") + ".contigs_classified.tsv"
  Int disk_size = 50

  command <<<
    set -e
    pip install pandas zstandard --quiet --no-cache-dir
    python3<<CODE
import re
import sys
import tarfile
import tempfile

import pandas as pd
import zstandard as zstd


def _resolve_file(path):
    if not any(path.endswith(ext) for ext in ('.tar.zst', '.tar.gz', '.tar.bz2', '.tar')):
        return path
    extract_dir = tempfile.mkdtemp()
    if path.endswith('.tar.zst'):
        dctx = zstd.ZstdDecompressor()
        with open(path, 'rb') as fh:
            with dctx.stream_reader(fh) as reader:
                with tarfile.open(fileobj=reader, mode='r|') as tar:
                    tar.extractall(path=extract_dir)
    else:
        with tarfile.open(path, 'r:*') as tar:
            tar.extractall(path=extract_dir)
    import os
    files = [os.path.join(extract_dir, f) for f in os.listdir(extract_dir)
             if os.path.isfile(os.path.join(extract_dir, f))]
    if len(files) != 1:
        print(f"ERROR: expected exactly 1 file in tarball {path}, found {len(files)}", file=sys.stderr)
        sys.exit(1)
    return files[0]

def classify_sequence(group, min_viral_proportion=0.1, min_nonviral_proportion=0.1, min_chunk_count=5):
    group = group.copy()
    group['delta'] = group['max_score_1'] - group['max_score_0']
    group['confidence'] = group['delta'].abs().pow(0.5)

    # Weighted mean delta (confidence-weighted)
    if group['confidence'].sum() > 0:
        weighted_delta = (group['delta'] * group['confidence']).sum() / group['confidence'].sum()
    else:
        weighted_delta = group['delta'].mean()

    # Count high-confidence viral chunks
    confident_viral = (group['max_score_1'] > 0.8) & (group['max_score_0'] < 0.3)
    # Count high-confidence non-viral chunks
    confident_nonviral = (group['max_score_0'] > 0.8) & (group['max_score_1'] < 0.3)
    # Count ambiguous chunks
    ambiguous = (group['max_score_1'] > 0.7) & (group['max_score_0'] > 0.7)

    n_chunks = len(group)
    n_confident_viral = confident_viral.sum()
    n_confident_nonviral = confident_nonviral.sum()
    n_ambiguous = ambiguous.sum()
    # exclude ambiguous chunks from denominator so they don't dilute
    # the proportion of confident viral/nonviral chunks
    n_effective = n_chunks - n_ambiguous
    viral_proportion = n_confident_viral / n_effective if n_effective > 0 else 0
    nonviral_proportion = n_confident_nonviral / n_effective if n_effective > 0 else 0

    # weighted_delta is the primary signal; chunk evidence determines tier
    if weighted_delta > 0.3:
        call = 'Viral'
        if n_confident_viral >= 1 and viral_proportion >= min_viral_proportion:
            tier = 'high_confidence' if weighted_delta > 0.6 else 'moderate_confidence'
        else:
            tier = 'low_confidence'
    elif weighted_delta < -0.3:
        call = 'Non-viral'
        if n_confident_nonviral >= 1 and nonviral_proportion >= min_nonviral_proportion:
            tier = 'high_confidence' if weighted_delta < -0.6 else 'moderate_confidence'
        else:
            tier = 'low_confidence'
    else:
        call = 'Ambiguous'
        tier = 'review'

    # Apply low chunk count penalty
    if n_chunks < min_chunk_count:
        if tier in ['high_confidence', 'moderate_confidence']:
            tier = 'low_confidence'
        elif tier == 'low_confidence':
            tier = 'review'
    return pd.Series({
        'call': call,
        'tier': tier,
        'weighted_delta': round(weighted_delta, 3),
        'n_chunks': n_chunks,
        'n_confident_viral': n_confident_viral,
        'n_confident_nonviral': n_confident_nonviral,
        'n_ambiguous': n_ambiguous,
        'viral_proportion': round(viral_proportion, 3),
        'nonviral_proportion': round(nonviral_proportion, 3)
    })

df = pd.read_csv(_resolve_file("~{virnucpro_scores_tsv}"), sep='\t')

required_cols = ["~{id_col}", 'max_score_0', 'max_score_1']
missing = [c for c in required_cols if c not in df.columns]
if missing:
    print(f"Error: Missing required columns: {missing}", file=sys.stderr)
    sys.exit(1)

df['ID'] = df["~{id_col}"].str.replace(r'_chunk_\d+$', '', regex=True)
df['_group'] = df["~{id_col}"].str.extract("~{id_pattern}")

n_unmatched = df['_group'].isna().sum()
if n_unmatched == len(df):
    print("Error: No valid NODE IDs extracted. Check id_pattern.", file=sys.stderr)
    sys.exit(1)
elif n_unmatched > 0:
    print(f"Warning: {n_unmatched} of {len(df)} rows did not match id_pattern and were excluded.", file=sys.stderr)

for col in ['max_score_0', 'max_score_1']:
    if df[col].isna().any():
        print(f"Warning: {df[col].isna().sum()} NaN values in {col}.", file=sys.stderr)

group_to_id = df.dropna(subset=['_group']).groupby('_group')['ID'].first()
results = df.groupby('_group').apply(
    lambda g: classify_sequence(g, ~{min_viral_prop}, ~{min_nonviral_prop}, ~{min_chunks}),
    include_groups=False
)
results = results.reset_index()
results['ID'] = results['_group'].map(group_to_id)
results = results.drop(columns=['_group'])
results = results[['ID'] + [c for c in results.columns if c != 'ID']]

def natural_sort_key(val):
    return [int(s) if s.isdigit() else s.lower() for s in re.split(r'(\d+)', str(val))]
results = results.sort_values('ID', key=lambda col: col.map(natural_sort_key))

results.to_csv("~{out_filename}", sep='\t', index=False)
CODE
  >>>

  output {
    File contig_classifications = "~{out_filename}"
  }

  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 2
    maxRetries: 2
  }
}

task classify_reads_by_contig {
  meta {
    description: "Classify reads by contig mapping using PAF alignments and VirNucPro contig classifications via DuckDB SQL joins."
  }

  input {
    File    paf_file
    File    contig_classifications

    Int     min_mapq       = 5
    Float   min_identity   = 90.0
    Float   min_query_cov  = 80.0

    String  docker         = "quay.io/broadinstitute/py3-bio:0.1.3"
  }

  parameter_meta {
    paf_file: {
      description: "Gzipped PAF file from minimap2 alignment (with pct_identity and pct_query_cov appended columns).",
      patterns: ["*.paf.gz", "*.paf"],
      category: "required"
    }
    contig_classifications: {
      description: "TSV of contig classifications from classify_virnucpro_contigs (columns: ID, call, tier, weighted_delta, etc.).",
      patterns: ["*.tsv"],
      category: "required"
    }
    min_mapq: {
      description: "Minimum mapping quality threshold for well-mapped reads.",
      category: "advanced"
    }
    min_identity: {
      description: "Minimum percent identity threshold for well-mapped reads (0-100 scale).",
      category: "advanced"
    }
    min_query_cov: {
      description: "Minimum percent query coverage threshold for well-mapped reads (0-100 scale).",
      category: "advanced"
    }
  }

  String out_filename = basename(paf_file, ".paf.gz") + ".reads_classified.tsv"
  Int disk_size = 100

  command <<<
    set -e
    pip install duckdb --quiet --no-cache-dir
    python3<<CODE
    import gzip
    import os
    import sys
    import tempfile

    import duckdb


    OUTPUT_COLUMNS = [
        'read_id', 'read_length', 'contig_id', 'contig_length', 'strand',
        'mapping_quality', 'pct_identity', 'pct_query_cov', 'mapped_well',
        'call', 'tier', 'weighted_delta', 'n_chunks', 'n_confident_viral',
        'n_confident_nonviral', 'n_ambiguous', 'viral_proportion', 'nonviral_proportion',
    ]


    def prepare_paf_file(paf_path):
        """Normalize PAF to fixed 15-column TSV for duckdb ingestion.

        Standard PAF has 12 columns, followed by optional SAM-like tags (key:type:value),
        with pct_identity and pct_query_cov appended as the last two tab-separated fields.
        Tags can appear as a single semicolon-delimited field or as multiple tab-separated
        fields -- either way they sit between column 12 and the two trailing numeric columns.

        Output: 12 standard PAF cols + tags + pct_identity + pct_query_cov
        """
        opener = gzip.open if paf_path.endswith('.gz') else open
        tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
        with opener(paf_path, 'rt') as f:
            for line in f:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 14:
                    continue
                # First 12 are always the standard PAF columns
                row = fields[:12]
                # Last 2 are always pct_identity and pct_query_cov
                pct_identity = fields[-2]
                pct_query_cov = fields[-1]
                # Everything in between is tags
                tag_fields = fields[12:-2]
                tags = '\t'.join(tag_fields) if tag_fields else ''
                row.extend([tags, pct_identity, pct_query_cov])
                tmp.write('\t'.join(row) + '\n')
        tmp.close()
        return tmp.name


    con = duckdb.connect()

    # Step 1: Parse PAF -- normalize variable-width file then read as fixed 15-col TSV
    print(f"Reading PAF file: ~{paf_file}", file=sys.stderr)
    normalized_paf = prepare_paf_file("~{paf_file}")

    con.execute("""
        CREATE TABLE paf AS
        SELECT * FROM read_csv($1,
            delim='\t', header=false, all_varchar=true,
            columns={
                'query_name': 'VARCHAR',
                'query_length': 'VARCHAR',
                'query_start': 'VARCHAR',
                'query_end': 'VARCHAR',
                'strand': 'VARCHAR',
                'target_name': 'VARCHAR',
                'target_length': 'VARCHAR',
                'target_start': 'VARCHAR',
                'target_end': 'VARCHAR',
                'num_matches': 'VARCHAR',
                'alignment_block_length': 'VARCHAR',
                'mapping_quality': 'VARCHAR',
                '_tags': 'VARCHAR',
                'pct_identity': 'VARCHAR',
                'pct_query_cov': 'VARCHAR'
            }
        )
    """, [normalized_paf])

    os.unlink(normalized_paf)

    n_total = con.execute("SELECT count(*) FROM paf").fetchone()[0]
    print(f"  {n_total} alignments loaded.", file=sys.stderr)

    # Cast numeric columns
    con.execute("""
        CREATE TABLE paf_typed AS
        SELECT
            query_name,
            CAST(query_length AS INTEGER) AS query_length,
            CAST(query_start AS INTEGER) AS query_start,
            CAST(query_end AS INTEGER) AS query_end,
            strand,
            target_name,
            CAST(target_length AS INTEGER) AS target_length,
            CAST(target_start AS INTEGER) AS target_start,
            CAST(target_end AS INTEGER) AS target_end,
            CAST(num_matches AS INTEGER) AS num_matches,
            CAST(alignment_block_length AS INTEGER) AS alignment_block_length,
            CAST(mapping_quality AS INTEGER) AS mapping_quality,
            CAST(pct_identity AS DOUBLE) AS pct_identity,
            CAST(pct_query_cov AS DOUBLE) AS pct_query_cov,
            _tags
        FROM paf
    """)
    con.execute("DROP TABLE paf")

    # Auto-detect scale of pct_identity and pct_query_cov (0-1 vs 0-100)
    max_id, max_cov = con.execute(
        "SELECT max(pct_identity), max(pct_query_cov) FROM paf_typed"
    ).fetchone()
    fractional_scale = (max_id is not None and max_id <= 1.0
                        and max_cov is not None and max_cov <= 1.0)
    if fractional_scale:
        print("  Detected fractional scale (0-1) for pct_identity/pct_query_cov.", file=sys.stderr)
        # Normalize thresholds from percentage to fractional
        min_identity = ~{min_identity} / 100.0
        min_query_cov = ~{min_query_cov} / 100.0
    else:
        min_identity = ~{min_identity}
        min_query_cov = ~{min_query_cov}

    # Step 2: Filter secondary/supplementary and flag mapping quality
    con.execute("""
        CREATE TABLE paf_filtered AS
        SELECT
            query_name, query_length, query_start, query_end, strand,
            target_name, target_length, target_start, target_end,
            num_matches, alignment_block_length, mapping_quality,
            pct_identity, pct_query_cov,
            (mapping_quality >= $1 AND pct_identity >= $2 AND pct_query_cov >= $3) AS mapped_well
        FROM paf_typed
        WHERE _tags IS NULL
           OR NOT (_tags LIKE '%tp:A:S%' OR _tags LIKE '%tp:A:I%')
    """, [~{min_mapq}, min_identity, min_query_cov])

    n_secondary = n_total - con.execute("SELECT count(*) FROM paf_filtered").fetchone()[0]
    if n_secondary > 0:
        print(f"Removed {n_secondary} secondary/supplementary alignments.", file=sys.stderr)
    con.execute("DROP TABLE paf_typed")

    n_primary = con.execute("SELECT count(*) FROM paf_filtered").fetchone()[0]
    print(f"  {n_primary} primary alignments retained.", file=sys.stderr)

    # Read classification file
    print(f"Reading classification file: ~{contig_classifications}", file=sys.stderr)
    con.execute("""
        CREATE TABLE clf AS
        SELECT * FROM read_csv($1, delim='\t', header=true, auto_detect=true)
    """, ["~{contig_classifications}"])
    n_clf = con.execute("SELECT count(*) FROM clf").fetchone()[0]
    print(f"  {n_clf} classified contigs loaded.", file=sys.stderr)

    # Step 3: Left join + fill defaults for unclassified
    con.execute("""
        CREATE TABLE merged AS
        SELECT
            p.query_name, p.query_length, p.target_name, p.target_length, p.strand,
            p.mapping_quality, p.pct_identity, p.pct_query_cov, p.mapped_well,
            COALESCE(c.call, 'Unclassified') AS call,
            COALESCE(c.tier, '') AS tier,
            COALESCE(c.weighted_delta, 0) AS weighted_delta,
            COALESCE(c.n_chunks, 0) AS n_chunks,
            COALESCE(c.n_confident_viral, 0) AS n_confident_viral,
            COALESCE(c.n_confident_nonviral, 0) AS n_confident_nonviral,
            COALESCE(c.n_ambiguous, 0) AS n_ambiguous,
            COALESCE(c.viral_proportion, 0) AS viral_proportion,
            COALESCE(c.nonviral_proportion, 0) AS nonviral_proportion
        FROM paf_filtered p
        LEFT JOIN clf c ON p.target_name = c.ID
    """)
    con.execute("DROP TABLE paf_filtered")
    con.execute("DROP TABLE clf")

    # Step 4: Aggregate to one row per read
    # Pick best alignment per read (highest MAPQ, then highest pct_identity).
    # If a read maps to contigs with different classifications, flag as Multi-mapped.
    print("Aggregating to one row per read...", file=sys.stderr)
    con.execute("""
        CREATE TABLE result AS
        WITH ranked AS (
            SELECT *,
                ROW_NUMBER() OVER (
                    PARTITION BY query_name
                    ORDER BY mapping_quality DESC, pct_identity DESC
                ) AS rn
            FROM merged
        ),
        best AS (
            SELECT * FROM ranked WHERE rn = 1
        ),
        call_counts AS (
            SELECT query_name, COUNT(DISTINCT call) AS n_distinct_calls
            FROM merged
            GROUP BY query_name
        )
        SELECT
            b.query_name AS read_id,
            b.query_length AS read_length,
            b.target_name AS contig_id,
            b.target_length AS contig_length,
            b.strand,
            b.mapping_quality,
            b.pct_identity,
            b.pct_query_cov,
            b.mapped_well,
            CASE WHEN cc.n_distinct_calls > 1 THEN 'Multi-mapped' ELSE b.call END AS call,
            CASE WHEN cc.n_distinct_calls > 1 THEN 'review' ELSE b.tier END AS tier,
            b.weighted_delta,
            b.n_chunks,
            b.n_confident_viral,
            b.n_confident_nonviral,
            b.n_ambiguous,
            b.viral_proportion,
            b.nonviral_proportion
        FROM best b
        JOIN call_counts cc ON b.query_name = cc.query_name
    """)
    con.execute("DROP TABLE merged")

    # Stats
    n_reads = con.execute("SELECT count(*) FROM result").fetchone()[0]
    n_well = con.execute("SELECT count(*) FROM result WHERE mapped_well").fetchone()[0]
    n_multi = con.execute("SELECT count(*) FROM result WHERE call = 'Multi-mapped'").fetchone()[0]
    print(f"  {n_reads} reads in output.", file=sys.stderr)
    pct_well = (n_well / n_reads * 100) if n_reads > 0 else 0
    print(f"  {n_well} reads mapped well ({pct_well:.1f}%).", file=sys.stderr)
    if n_multi > 0:
        print(f"  {n_multi} reads flagged as Multi-mapped.", file=sys.stderr)

    # Step 5: Write output as TSV
    con.execute(f"""
        COPY (SELECT {', '.join(OUTPUT_COLUMNS)} FROM result)
        TO $1 (FORMAT CSV, DELIMITER '\t', HEADER TRUE)
    """, ["~{out_filename}"])
    print(f"Output written to: ~{out_filename}", file=sys.stderr)

    con.close()
CODE
  >>>

  output {
    File read_classifications = "~{out_filename}"
  }

  runtime {
    docker: docker
    memory: "4 GB"
    cpu: 2
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem2_ssd1_v2_x4"
    preemptible: 2
    maxRetries: 2
  }
}

task parse_kraken2_reads {
  meta {
    description: "Parse Kraken2 per-read output and annotate each read with NCBI taxonomy (name, kingdom, rank) using a pre-built DuckDB taxonomy database."
  }

  input {
    File    kraken2_reads_output
    File    taxonomy_db
    String  sample_id = sub(basename(kraken2_reads_output), "\\.kraken2\\.reads\\.txt(\\.gz)?$", "")
    Boolean resolve_strains = false

    String  docker = "quay.io/broadinstitute/py3-bio:0.1.3"
  }

  parameter_meta {
    kraken2_reads_output: {
      description: "Kraken2 per-read classification output file (may be gzipped). Format: C/U <read_id> <taxid> <length> <kmer_info>.",
      patterns: ["*.txt.gz", "*.txt"],
      category: "required"
    }
    taxonomy_db: {
      description: "Pre-built DuckDB taxonomy database file from build_taxonomy_db.py.",
      patterns: ["*.duckdb"],
      category: "required"
    }
    sample_id: {
      description: "Sample identifier. Defaults to the filename stem after stripping .kraken2.reads.txt(.gz) extension.",
      category: "common"
    }
    resolve_strains: {
      description: "When true, reclassify 'no rank' nodes whose lineage passes through a species or strain ancestor as 'strain' in TAX_RANK.",
      category: "advanced"
    }
  }

  String out_filename = "~{sample_id}.read_taxonomy.tsv"

  command <<<
    set -e
    pip install duckdb --quiet --no-cache-dir
    python3<<CODE
    import csv
    import gzip
    import sys

    import duckdb


    class DuckDBTaxonomyDatabase:
        """NCBI taxonomy database backed by a pre-built DuckDB file.

        Provides the same interface as TaxonomyDatabase (get_name, get_kingdom)
        but loads much faster by reading pre-computed data from DuckDB.
        """

        def __init__(self, db_path, resolve_strains=False):
            """Initialize taxonomy database from a DuckDB file.

            Args:
                db_path: Path to ncbi_taxonomy.duckdb file
                resolve_strains: If True, resolve 'no rank' nodes below species to 'strain'
            """
            print(f"Loading taxonomy from DuckDB: {db_path}...", file=sys.stderr)

            con = duckdb.connect(str(db_path), read_only=True)
            try:
                rows = con.execute("SELECT taxid, name, kingdom, rank, parent_id FROM taxonomy").fetchall()
                has_parent = True
            except Exception:
                print("Warning: 'rank' column not found in taxonomy table; defaulting ranks to 'unknown'",
                      file=sys.stderr)
                rows_no_rank = con.execute("SELECT taxid, name, kingdom FROM taxonomy").fetchall()
                rows = [(taxid, name, kingdom, 'unknown', None) for taxid, name, kingdom in rows_no_rank]
                has_parent = False
            con.close()

            self.names = {}
            self.kingdoms = {}
            self.ranks = {}
            parents = {}
            for taxid, name, kingdom, rank, parent_id in rows:
                self.names[taxid] = name
                self.kingdoms[taxid] = kingdom
                self.ranks[taxid] = rank
                if parent_id is not None:
                    parents[taxid] = parent_id

            if resolve_strains and has_parent:
                # Resolve "no rank" -> "strain" for nodes below species
                for taxid, rank in list(self.ranks.items()):
                    if rank == 'no rank':
                        current = parents.get(taxid)
                        while current and current != 1:
                            ancestor_rank = self.ranks.get(current, 'unknown')
                            if ancestor_rank in ('species', 'strain'):
                                self.ranks[taxid] = 'strain'
                                break
                            current = parents.get(current)

            print(f"Loaded {len(self.names)} taxonomy entries from DuckDB", file=sys.stderr)

        def get_rank(self, taxid, resolve_strains=False):
            """Get taxonomic rank for taxid."""
            if taxid == 0:
                return 'unclassified'
            return self.ranks.get(taxid, 'unknown')

        def get_name(self, taxid):
            """Get scientific name for taxid."""
            if taxid == 0:
                return 'Unclassified'
            return self.names.get(taxid, f'Unknown (taxid:{taxid})')

        def get_kingdom(self, taxid):
            """Get kingdom/domain for a taxonomy ID."""
            if taxid == 0:
                return 'Unclassified'
            return self.kingdoms.get(taxid, 'Other')


    def parse_kraken2_output(kraken_file, tax_db, output_file, sample_id=None, resolve_strains=False):
        """Parse Kraken2 read classification output and annotate with taxonomy.

        Args:
            kraken_file: Path to Kraken2 output file (may be gzipped)
            tax_db: DuckDBTaxonomyDatabase instance
            output_file: Path to output TSV file
            sample_id: Sample identifier (if None, extracted from filename)
            resolve_strains: If True, resolve 'no rank' nodes below species to 'strain'
        """
        # Extract sample ID from filename if not provided
        if sample_id is None:
            sample_id = kraken_file.split('/')[-1]
            for ext in ['.kraken2', '.output', '.txt', '.gz']:
                if sample_id.endswith(ext):
                    sample_id = sample_id[:-len(ext)]

        print(f"Processing sample: {sample_id}", file=sys.stderr)

        # Handle gzipped files
        if kraken_file.endswith('.gz'):
            f = gzip.open(kraken_file, 'rt')
        else:
            f = open(kraken_file, 'r')

        classified_count = 0
        unclassified_count = 0

        # Collect rows
        rows = []

        try:
            for line in f:
                # Skip empty lines
                line = line.strip()
                if not line:
                    continue

                # Parse Kraken2 output format
                # Format: C/U <read_id> <taxid> <length> <kmer_info>
                parts = line.split('\t')

                if len(parts) < 3:
                    print(f"Warning: Skipping malformed line: {line[:100]}", file=sys.stderr)
                    continue

                classification = parts[0].strip()  # C or U
                read_id = parts[1].strip()
                taxid_str = parts[2].strip()

                # Handle unclassified reads (taxid = 0)
                try:
                    taxid = int(taxid_str)
                except ValueError:
                    print(f"Warning: Invalid taxid '{taxid_str}' for read {read_id}", file=sys.stderr)
                    continue

                if classification == 'U':
                    unclassified_count += 1
                    tax_name = 'Unclassified'
                    kingdom = 'Unclassified'
                    tax_rank = 'unclassified'
                else:
                    classified_count += 1
                    tax_name = tax_db.get_name(taxid)
                    kingdom = tax_db.get_kingdom(taxid)
                    tax_rank = tax_db.get_rank(taxid, resolve_strains=resolve_strains)

                rows.append((sample_id, read_id, taxid, tax_name, kingdom, tax_rank))

        finally:
            f.close()

        # Write output as TSV
        _write_tsv(rows, output_file)

        print(f"\nProcessing complete:", file=sys.stderr)
        print(f"  Classified reads: {classified_count}", file=sys.stderr)
        print(f"  Unclassified reads: {unclassified_count}", file=sys.stderr)
        print(f"  Total reads: {classified_count + unclassified_count}", file=sys.stderr)


    def _write_tsv(rows, output_file):
        """Write rows as TSV."""
        with open(output_file, 'w', newline='') as out_f:
            writer = csv.writer(out_f, delimiter='\t')
            writer.writerow(['SAMPLE_ID', 'READ_ID', 'TAXONOMY_ID', 'TAX_NAME', 'KINGDOM', 'TAX_RANK'])
            for row in rows:
                writer.writerow(row)


    tax_db = DuckDBTaxonomyDatabase("~{taxonomy_db}", resolve_strains=~{true="True" false="False" resolve_strains})
    parse_kraken2_output("~{kraken2_reads_output}", tax_db, "~{out_filename}", "~{sample_id}")
    CODE
  >>>

  output {
    File read_taxonomy = "~{out_filename}"
  }

  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 1
    disks: "local-disk ~{ceil(size(kraken2_reads_output)*3 + size(taxonomy_db) + 20)} HDD"
    disk: "~{ceil(size(kraken2_reads_output)*3 + size(taxonomy_db) + 20)} GB" # TES
    dx_instance_type: "mem2_ssd1_v2_x2"
    preemptible: 2
    maxRetries: 2
  }
}

task summarize_kb_extract_reads {
  meta {
    description: "Summarize kb extract read output with taxonomy annotation. Extracts reads from tarball, joins with taxonomy mapping, and outputs zstd-compressed TSV."
  }

  input {
    File    extract_reads_tar
    File    taxonomy_map_csv
    String  taxonomy_level = "highest"

    String  docker = "quay.io/broadinstitute/py3-bio:0.1.3"
  }

  parameter_meta {
    extract_reads_tar: {
      description: "Tarball containing sample subdirectories, each with <UID>.fastq.gz files. Top-level structure must be <sample_name>/<UID>.fastq.gz (no enclosing directory required).",
      patterns: ["*.tar", "*.tar.gz", "*.tar.bz2", "*.tar.zst"],
      category: "required"
    }
    taxonomy_map_csv: {
      description: "CSV file mapping palmDB_ID to taxonomy lineage. Columns: palmDB_ID,palmDB_ID,tax_level_1,...,tax_level_n,strand",
      patterns: ["*.csv", "*.tsv"],
      category: "required"
    }
    taxonomy_level: {
      description: "Taxonomy level to report: 'highest' (most general) or 'deepest' (most specific). Default 'highest'.",
      category: "advanced"
    }
  }

  String out_filename = "extracted_reads_summary.tsv.zst"
  Int disk_size = 50

  command <<<
    set -e
    pip install zstandard --quiet --no-cache-dir
    python3<<CODE
    import csv
    import gzip
    import os
    import sys
    import tarfile
    import zstandard as zstd

    TAXONOMY_LEVEL = "~{taxonomy_level}"
    OUTPUT_FILE = "~{out_filename}"

    # Derive sample ID from tarball filename, stripping all extensions (e.g. .tar.gz, .tar.zst)
    tar_basename = os.path.basename("~{extract_reads_tar}")
    sample_id = tar_basename.split('.')[0]

    # Load taxonomy map: palmDB_ID -> taxonomy lineage
    print("Loading taxonomy map...", file=sys.stderr)
    taxonomy_map = {}
    with open("~{taxonomy_map_csv}", 'r') as f:
        reader = csv.reader(f)
        header = next(reader, None)  # palmDB_ID,palmDB_ID,tax1,tax2,...,strand
        if header is None:
            print("Error: Empty taxonomy map file", file=sys.stderr)
            sys.exit(1)
        # Determine which columns contain taxonomy (exclude first 2 palmDB_ID duplicates and last strand)
        # Header: [palmDB_ID, palmDB_ID_dup, tax_level_1, ..., tax_level_n, strand]
        tax_cols = header[2:-1]  # All columns between the two IDs and strand

        for row in reader:
            if len(row) < 3:
                continue
            palm_id = row[0]
            # Extract taxonomy values (skip first 2 cols and last col which is strand)
            tax_values = row[2:-1] if len(row) > 2 else []
            taxonomy_map[palm_id] = tax_values

    print(f"  Loaded {len(taxonomy_map)} taxonomy entries", file=sys.stderr)

    # Extract tarball into isolated subdirectory to avoid iterating over unrelated files
    extract_dir = "extract_contents"
    os.makedirs(extract_dir, exist_ok=True)
    print("Extracting tarball...", file=sys.stderr)
    tar_path = "~{extract_reads_tar}"
    if tar_path.endswith('.zst'):
        dctx = zstd.ZstdDecompressor()
        with open(tar_path, 'rb') as zst_f:
            with dctx.stream_reader(zst_f) as reader:
                with tarfile.open(fileobj=reader, mode='r|') as tar:
                    tar.extractall(path=extract_dir)
    else:
        with tarfile.open(tar_path, 'r:*') as tar:
            tar.extractall(path=extract_dir)
    print("  Extraction complete", file=sys.stderr)

    # Open output file with zstd compression
    cctx = zstd.ZstdCompressor()
    with open(OUTPUT_FILE, 'wb') as out_f:
        compressor = cctx.stream_writer(out_f)

        # Write header
        header_line = "SAMPLE_ID\tREAD_ID\tDB_ID\tTAXONOMY_LINEAGE\tTAXONOMY_NAME\tSEQUENCE_LENGTH\n"
        compressor.write(header_line.encode('utf-8'))

        samples_processed = 0
        reads_processed = 0

        for sample_name in os.listdir(extract_dir):
            sample_path = os.path.join(extract_dir, sample_name)
            if not os.path.isdir(sample_path):
                continue

            samples_processed += 1

            # Process all FASTQ files in sample directory
            for filename in os.listdir(sample_path):
                if not filename.endswith('.fastq.gz'):
                    continue

                filepath = os.path.join(sample_path, filename)
                # Directory name is the DB_ID (palmDB_ID); sample ID comes from the tarball name
                db_id = sample_name

                # Look up taxonomy by DB_ID
                tax_values = taxonomy_map.get(db_id, [])

                if tax_values:
                    if TAXONOMY_LEVEL == "deepest":
                        # Get the last non-empty taxonomy level
                        tax_lineage = [t for t in tax_values if t and t.strip()]
                        tax_name = tax_lineage[-1] if tax_lineage else "Unclassified RdRP"
                    else:  # highest
                        # Get the first non-empty taxonomy level
                        tax_lineage = [t for t in tax_values if t and t.strip()]
                        tax_name = tax_lineage[0] if tax_lineage else "Unclassified RdRP"
                    taxonomy_str = ";".join(tax_lineage) if tax_lineage else "Unclassified RdRP"
                else:
                    tax_name = "Unclassified RdRP"
                    taxonomy_str = "Unclassified RdRP"

                # Parse FASTQ to get read IDs and sequence lengths
                with gzip.open(filepath, 'rt') as fq:
                    while True:
                        # Read 4 lines of FASTQ
                        header_line = fq.readline()
                        if not header_line:
                            break
                        seq_line = fq.readline().strip()
                        plus_line = fq.readline()
                        qual_line = fq.readline()

                        if not header_line.startswith('@'):
                            continue

                        read_id = header_line[1:].split()[0]  # Remove @ and take first token
                        seq_length = len(seq_line)

                        # Write output line
                        out_line = f"{sample_id}\t{read_id}\t{db_id}\t{taxonomy_str}\t{tax_name}\t{seq_length}\n"
                        compressor.write(out_line.encode('utf-8'))
                        reads_processed += 1

        compressor.close()

    print(f"Processed {samples_processed} samples, {reads_processed} reads", file=sys.stderr)
    print(f"Output written to: {OUTPUT_FILE}", file=sys.stderr)
    CODE
  >>>

  output {
    File summary_tsv_zst = "~{out_filename}"
  }

  runtime {
    docker: docker
    memory: "4 GB"
    cpu: 2
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem2_ssd1_v2_x4"
    preemptible: 2
    maxRetries: 2
  }
}

task centrifuger {
  meta {
    description: "Runs Centrifuger taxonomic classification on one or more BAM files. Each BAM is converted to FASTQ via picard SamToFastq, classified, and a Kraken-style kreport is generated. The pre-built index is loaded once and reused across all samples."
  }

  input {
    File         centrifuger_db_tgz  # pre-built Centrifuger index tarball (.tar.gz, .tar.lz4, .tar.zst)
    String       db_name             # index prefix (common stem of .1.cfr/.2.cfr/.3.cfr/.4.cfr files)

    Array[File]  reads_bams

    Int          machine_mem_gb = 240
    Int?         cpu
    String       docker = "ghcr.io/broadinstitute/docker-centrifuger:1.0.0"
  }

  Int disk_size = ceil((8 * size(reads_bams, "GB") + 3 * size(centrifuger_db_tgz, "GB") + 400) / 400.0) * 400

  parameter_meta {
    centrifuger_db_tgz: {
      description: "Pre-built Centrifuger index as a compressed tarball (.tar.gz, .tar.lz4, or .tar.zst). The tarball must contain index files sharing a common prefix (db_name). Extracted at runtime into a temporary directory.",
      patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.zst", "*.tar.bz2"],
      category: "required"
    }
    db_name: {
      description: "Centrifuger index prefix (the common filename stem of the .1.cfr, .2.cfr, .3.cfr, .4.cfr files inside the tarball).",
      category: "required"
    }
    reads_bams: {
      description: "Reads in BAM format, one file per sample. Each BAM is converted to FASTQ internally via picard SamToFastq. Paired-end BAMs produce R1/R2; single-end BAMs produce one FASTQ. Sample names are derived from BAM filenames.",
      patterns: ["*.bam"],
      category: "required"
    }
    machine_mem_gb: {
      description: "Memory in GB. Default 240 GB sized for NT-scale centrifuger index.",
      category: "other"
    }
    cpu: {
      description: "Number of CPUs. Default 8.",
      category: "other"
    }
    docker: {
      description: "Docker image with centrifuger and picard.",
      category: "other"
    }
  }

  command <<<
    set -ex -o pipefail

    # Extract centrifuger index tarball
    DB_DIR=$(mktemp -d --suffix _centrifuger_db)
    if [[ "~{centrifuger_db_tgz}" == *.tar.lz4 ]]; then
      lz4 -dc "~{centrifuger_db_tgz}" | tar -x -C "$DB_DIR"
    elif [[ "~{centrifuger_db_tgz}" == *.tar.zst ]]; then
      zstd -dc "~{centrifuger_db_tgz}" | tar -x -C "$DB_DIR"
    elif [[ "~{centrifuger_db_tgz}" == *.tar.bz2 ]]; then
      tar -xjf "~{centrifuger_db_tgz}" -C "$DB_DIR"
    else
      tar -xzf "~{centrifuger_db_tgz}" -C "$DB_DIR"
    fi
    INDEX_PREFIX="$DB_DIR/~{db_name}"

    # Process each BAM: convert to FASTQ, classify, generate kreport
    THREADS="~{select_first([cpu, 8])}"

    for bam in "~{sep='" "' reads_bams}"; do
      SAMPLE="$(basename $bam .bam)"

      # BAM -> FASTQ via picard SamToFastq (per CFGR-02 / D-02)
      picard SamToFastq \
        I="$bam" \
        FASTQ="${SAMPLE}_R1.fq" \
        SECOND_END_FASTQ="${SAMPLE}_R2.fq" \
        READ1_SUFFIX=/1 \
        READ2_SUFFIX=/2 \
        VALIDATION_STRINGENCY=LENIENT

      if [ -s "${SAMPLE}_R2.fq" ]; then
        # Paired-end BAM
        centrifuger \
          -x "$INDEX_PREFIX" \
          -1 "${SAMPLE}_R1.fq" \
          -2 "${SAMPLE}_R2.fq" \
          -t "$THREADS" \
          > "${SAMPLE}.centrifuger.tsv" \
          2> "${SAMPLE}.centrifuger.log"
      else
        # Single-end BAM
        centrifuger \
          -x "$INDEX_PREFIX" \
          -u "${SAMPLE}_R1.fq" \
          -t "$THREADS" \
          > "${SAMPLE}.centrifuger.tsv" \
          2> "${SAMPLE}.centrifuger.log"
      fi

      # Generate Kraken-style report (requires index on disk -- per Phase 4 D-02)
      centrifuger-kreport \
        -x "$INDEX_PREFIX" \
        "${SAMPLE}.centrifuger.tsv" \
        > "${SAMPLE}.centrifuger.kreport"

      # Clean up per-sample FASTQ to free disk
      rm -f "${SAMPLE}_R1.fq" "${SAMPLE}_R2.fq"
    done
  >>>

  output {
    Array[File] classification_tsvs = glob("*.centrifuger.tsv")
    Array[File] kreports            = glob("*.centrifuger.kreport")
    Array[File] centrifuger_logs    = glob("*.centrifuger.log")
  }

  runtime {
    docker:           docker
    memory:           "~{machine_mem_gb} GB"
    cpu:              select_first([cpu, 8])
    disks:            "local-disk ~{disk_size} SSD"
    disk:             "~{disk_size} GB" # TES
    dx_instance_type: "mem3_ssd1_v2_x8"
    preemptible:      2
    maxRetries:       2
  }
}

task join_read_classifications {
  meta {
    description: "Join read-level classifications from Kallisto, Kraken2, VirNucPro, and geNomad into a single ZSTD-compressed Parquet file keyed on SAMPLE_ID + READ_ID. Uses DuckDB FULL OUTER JOINs so reads from any source are included."
  }

  input {
    File?   kallisto_summary          # Kallisto summary parquet (multi-sample, filtered by sample_id)
    File?   kraken2_reads             # Kraken2 reads parquet (filtered by sample_id)
    File?   vnp_reads                 # VirNucPro reads parquet (all rows used)
    File?   genomad_virus_summary     # geNomad virus summary TSV (LEFT JOIN via VNP_CONTIG_ID)
    String  sample_id                 # Required — filters Kallisto/K2 tables, stamps SAMPLE_ID column

    String  docker = "quay.io/broadinstitute/py3-bio:0.1.3"
  }

  parameter_meta {
    kallisto_summary: {
      description: "Kallisto summary TSV for a single sample, or a tar.zst/tar.gz tarball containing one such file. Skipped if not provided or empty.",
      patterns: ["*.tsv", "*.tar.zst", "*.tar.gz"],
      category: "common"
    }
    kraken2_reads: {
      description: "Kraken2 per-read classifications in Parquet or TSV format, or a tar.zst/tar.gz tarball containing one such file. Filtered by sample_id at runtime. Skipped if not provided or empty.",
      patterns: ["*.parquet", "*.tsv", "*.tar.zst", "*.tar.gz"],
      category: "common"
    }
    vnp_reads: {
      description: "VirNucPro per-read classifications in Parquet or TSV format, or a tar.zst/tar.gz tarball containing one such file. All rows used (not filtered by sample_id). Skipped if not provided or empty.",
      patterns: ["*.parquet", "*.tsv", "*.tar.zst", "*.tar.gz"],
      category: "common"
    }
    genomad_virus_summary: {
      description: "geNomad virus summary in TSV format. LEFT JOINed to reads via VNP_CONTIG_ID to annotate viral contigs. Skipped if not provided or empty.",
      patterns: ["*.tsv"],
      category: "common"
    }
    sample_id: {
      description: "Sample identifier used to filter Kallisto and Kraken2 tables and stamp the SAMPLE_ID column on every output row.",
      category: "required"
    }
  }

  String out_filename = "~{sample_id}.read_classifications.parquet"

  command <<<
    set -e
    pip install duckdb zstandard --quiet --no-cache-dir
    python3<<CODE
import os
import sys
import tarfile
import tempfile

import duckdb
import zstandard as zstd


def _file_is_usable(path):
    """Return True if path is provided, exists, and is non-empty."""
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def _resolve_file(path):
    """If path is a tarball, extract it to a temp dir and return the single inner file path.
    Supports .tar.zst, .tar.gz, .tar.bz2, and .tar. Otherwise returns path unchanged."""
    if not any(path.endswith(ext) for ext in ('.tar.zst', '.tar.gz', '.tar.bz2', '.tar')):
        return path
    extract_dir = tempfile.mkdtemp()
    if path.endswith('.tar.zst'):
        dctx = zstd.ZstdDecompressor()
        with open(path, 'rb') as fh:
            with dctx.stream_reader(fh) as reader:
                with tarfile.open(fileobj=reader, mode='r|') as tar:
                    tar.extractall(path=extract_dir)
    else:
        with tarfile.open(path, 'r:*') as tar:
            tar.extractall(path=extract_dir)
    files = [os.path.join(extract_dir, f) for f in os.listdir(extract_dir)
             if os.path.isfile(os.path.join(extract_dir, f))]
    if len(files) != 1:
        print(f"ERROR: expected exactly 1 file in tarball {path}, found {len(files)}", file=sys.stderr)
        sys.exit(1)
    return files[0]


def _duckdb_reader(path):
    """Return the appropriate DuckDB reader function call string for a given file path."""
    if path.endswith('.parquet'):
        return f"read_parquet('{path}')"
    return f"read_csv('{path}', delim='\\t', header=true, auto_detect=true)"


summary      = "~{default='__NONE__' kallisto_summary}"
vnp          = "~{default='__NONE__' vnp_reads}"
kraken2      = "~{default='__NONE__' kraken2_reads}"
virus_summary = "~{default='__NONE__' genomad_virus_summary}"
sample_id    = "~{sample_id}"
output       = "~{out_filename}"

has_kallisto = _file_is_usable(summary)
has_vnp      = _file_is_usable(vnp)
has_kraken2  = _file_is_usable(kraken2)
has_genomad  = _file_is_usable(virus_summary)

if has_kallisto:
    summary = _resolve_file(summary)
if has_vnp:
    vnp = _resolve_file(vnp)
if has_kraken2:
    kraken2 = _resolve_file(kraken2)

if not any([has_kallisto, has_vnp, has_kraken2, has_genomad]):
    print("ERROR: All input files are missing or empty — nothing to join.", file=sys.stderr)
    sys.exit(1)

con = duckdb.connect()

# ── 1. Kallisto summary ──────────────────────────────────────────────
if has_kallisto:
    print(f"Reading Kallisto summary: {summary}", file=sys.stderr)
    con.execute(f"""
        CREATE TABLE kallisto AS
        SELECT
            SAMPLE_ID,
            READ_ID,
            DB_ID            AS KALLISTO_DB_ID,
            TAXONOMY_LINEAGE AS KALLISTO_TAXONOMY_LINEAGE,
            TAXONOMY_NAME    AS KALLISTO_TAXONOMY_NAME,
            SEQUENCE_LENGTH  AS KALLISTO_SEQUENCE_LENGTH
        FROM {_duckdb_reader(summary)}
    """)
    n_kallisto = con.execute("SELECT count(*) FROM kallisto").fetchone()[0]
    print(f"  {n_kallisto:,} Kallisto reads loaded.", file=sys.stderr)
else:
    print(f"  Skipping Kallisto (file missing or empty).", file=sys.stderr)
    con.execute("""
        CREATE TABLE kallisto (
            SAMPLE_ID VARCHAR, READ_ID VARCHAR,
            KALLISTO_DB_ID VARCHAR, KALLISTO_TAXONOMY_LINEAGE VARCHAR,
            KALLISTO_TAXONOMY_NAME VARCHAR, KALLISTO_SEQUENCE_LENGTH BIGINT
        )
    """)

# ── 2. VirNucPro reads ───────────────────────────────────────────────
if has_vnp:
    print(f"Reading VirNucPro reads: {vnp}", file=sys.stderr)
    con.execute(f"""
        CREATE TABLE vnp AS
        SELECT
            read_id           AS READ_ID,
            contig_id         AS VNP_CONTIG_ID,
            contig_length     AS VNP_CONTIG_LENGTH,
            strand            AS VNP_STRAND,
            read_length       AS VNP_READ_LENGTH,
            mapping_quality   AS VNP_MAPPING_QUALITY,
            pct_identity      AS VNP_PCT_IDENTITY,
            pct_query_cov     AS VNP_PCT_QUERY_COV,
            mapped_well       AS VNP_MAPPED_WELL,
            call              AS VNP_CALL,
            tier              AS VNP_TIER,
            weighted_delta    AS VNP_WEIGHTED_DELTA,
            n_chunks          AS VNP_N_CHUNKS,
            n_confident_viral       AS VNP_N_CONFIDENT_VIRAL,
            n_confident_nonviral    AS VNP_N_CONFIDENT_NONVIRAL,
            n_ambiguous             AS VNP_N_AMBIGUOUS,
            viral_proportion        AS VNP_VIRAL_PROPORTION,
            nonviral_proportion     AS VNP_NONVIRAL_PROPORTION
        FROM {_duckdb_reader(vnp)}
    """)
    n_vnp = con.execute("SELECT count(*) FROM vnp").fetchone()[0]
    print(f"  {n_vnp:,} VNP reads loaded.", file=sys.stderr)
else:
    print(f"  Skipping VirNucPro (file missing or empty).", file=sys.stderr)
    con.execute("""
        CREATE TABLE vnp (
            READ_ID VARCHAR, VNP_CONTIG_ID VARCHAR, VNP_CONTIG_LENGTH BIGINT,
            VNP_STRAND VARCHAR, VNP_READ_LENGTH BIGINT, VNP_MAPPING_QUALITY BIGINT,
            VNP_PCT_IDENTITY DOUBLE, VNP_PCT_QUERY_COV DOUBLE, VNP_MAPPED_WELL BOOLEAN,
            VNP_CALL VARCHAR, VNP_TIER BIGINT, VNP_WEIGHTED_DELTA DOUBLE,
            VNP_N_CHUNKS BIGINT, VNP_N_CONFIDENT_VIRAL BIGINT,
            VNP_N_CONFIDENT_NONVIRAL BIGINT, VNP_N_AMBIGUOUS BIGINT,
            VNP_VIRAL_PROPORTION DOUBLE, VNP_NONVIRAL_PROPORTION DOUBLE
        )
    """)

# Kallisto and VNP both use READ_ID with /1 or /2 mate suffix
con.execute("""
    CREATE TABLE joined_kv AS
    SELECT
        COALESCE(k.READ_ID, v.READ_ID) AS READ_ID,
        k.KALLISTO_DB_ID,
        k.KALLISTO_TAXONOMY_LINEAGE,
        k.KALLISTO_TAXONOMY_NAME,
        k.KALLISTO_SEQUENCE_LENGTH,
        v.VNP_CONTIG_ID,
        v.VNP_CONTIG_LENGTH,
        v.VNP_STRAND,
        v.VNP_READ_LENGTH,
        v.VNP_MAPPING_QUALITY,
        v.VNP_PCT_IDENTITY,
        v.VNP_PCT_QUERY_COV,
        v.VNP_MAPPED_WELL,
        v.VNP_CALL,
        v.VNP_TIER,
        v.VNP_WEIGHTED_DELTA,
        v.VNP_N_CHUNKS,
        v.VNP_N_CONFIDENT_VIRAL,
        v.VNP_N_CONFIDENT_NONVIRAL,
        v.VNP_N_AMBIGUOUS,
        v.VNP_VIRAL_PROPORTION,
        v.VNP_NONVIRAL_PROPORTION
    FROM kallisto k
    FULL OUTER JOIN vnp v ON k.READ_ID = v.READ_ID
""")
con.execute("DROP TABLE kallisto")
con.execute("DROP TABLE vnp")
n_kv = con.execute("SELECT count(*) FROM joined_kv").fetchone()[0]
print(f"  {n_kv:,} rows after Kallisto + VNP full outer join.", file=sys.stderr)

# ── 3. Kraken2 reads ─────────────────────────────────────────────────
# Kraken2 READ_IDs lack the /1 or /2 mate suffix, so strip it from
# the Kallisto/VNP side when joining.
if has_kraken2:
    print(f"Reading Kraken2 reads: {kraken2}", file=sys.stderr)
    con.execute(f"""
        CREATE TABLE k2 AS
        SELECT
            READ_ID,
            TAXONOMY_ID AS K2_TAXONOMY_ID,
            TAX_NAME    AS K2_TAX_NAME,
            KINGDOM     AS K2_KINGDOM
        FROM {_duckdb_reader(kraken2)}
        WHERE SAMPLE_ID = $1
    """, [sample_id])
    n_k2 = con.execute("SELECT count(*) FROM k2").fetchone()[0]
    print(f"  {n_k2:,} Kraken2 reads loaded.", file=sys.stderr)
else:
    print(f"  Skipping Kraken2 (file missing or empty).", file=sys.stderr)
    con.execute("""
        CREATE TABLE k2 (
            READ_ID VARCHAR, K2_TAXONOMY_ID VARCHAR,
            K2_TAX_NAME VARCHAR, K2_KINGDOM VARCHAR
        )
    """)

con.execute("""
    CREATE TABLE joined_all AS
    SELECT
        COALESCE(jkv.READ_ID, k.READ_ID) AS READ_ID,
        jkv.KALLISTO_DB_ID,
        jkv.KALLISTO_TAXONOMY_LINEAGE,
        jkv.KALLISTO_TAXONOMY_NAME,
        jkv.KALLISTO_SEQUENCE_LENGTH,
        k.K2_TAXONOMY_ID,
        k.K2_TAX_NAME,
        k.K2_KINGDOM,
        jkv.VNP_CONTIG_ID,
        jkv.VNP_CONTIG_LENGTH,
        jkv.VNP_STRAND,
        jkv.VNP_READ_LENGTH,
        jkv.VNP_MAPPING_QUALITY,
        jkv.VNP_PCT_IDENTITY,
        jkv.VNP_PCT_QUERY_COV,
        jkv.VNP_MAPPED_WELL,
        jkv.VNP_CALL,
        jkv.VNP_TIER,
        jkv.VNP_WEIGHTED_DELTA,
        jkv.VNP_N_CHUNKS,
        jkv.VNP_N_CONFIDENT_VIRAL,
        jkv.VNP_N_CONFIDENT_NONVIRAL,
        jkv.VNP_N_AMBIGUOUS,
        jkv.VNP_VIRAL_PROPORTION,
        jkv.VNP_NONVIRAL_PROPORTION
    FROM joined_kv jkv
    FULL OUTER JOIN k2 k
        ON regexp_replace(jkv.READ_ID, '/[12]$', '') = k.READ_ID
""")
con.execute("DROP TABLE joined_kv")
con.execute("DROP TABLE k2")
n_all = con.execute("SELECT count(*) FROM joined_all").fetchone()[0]
print(f"  {n_all:,} rows after Kraken2 full outer join.", file=sys.stderr)

# ── 4. geNomad virus summary (LEFT JOIN on contig ID) ────────────────
# geNomad is contig-level annotation, not a read source — LEFT JOIN
# enriches reads that mapped to viral contigs via VNP.
if has_genomad:
    print(f"Reading geNomad virus summary: {virus_summary}", file=sys.stderr)
    con.execute("""
        CREATE TABLE genomad AS
        SELECT
            seq_name          AS CONTIG_ID,
            length            AS GENOMAD_LENGTH,
            topology          AS GENOMAD_TOPOLOGY,
            coordinates       AS GENOMAD_COORDINATES,
            n_genes           AS GENOMAD_N_GENES,
            genetic_code      AS GENOMAD_GENETIC_CODE,
            virus_score       AS GENOMAD_VIRUS_SCORE,
            fdr               AS GENOMAD_FDR,
            n_hallmarks       AS GENOMAD_N_HALLMARKS,
            marker_enrichment AS GENOMAD_MARKER_ENRICHMENT,
            taxonomy          AS GENOMAD_TAXONOMY
        FROM read_csv($1, delim='\t', header=true, auto_detect=true)
    """, [virus_summary])
    n_genomad = con.execute("SELECT count(*) FROM genomad").fetchone()[0]
    print(f"  {n_genomad:,} geNomad contigs loaded.", file=sys.stderr)
else:
    print(f"  Skipping geNomad (file missing or empty).", file=sys.stderr)
    con.execute("""
        CREATE TABLE genomad (
            CONTIG_ID VARCHAR, GENOMAD_LENGTH BIGINT, GENOMAD_TOPOLOGY VARCHAR,
            GENOMAD_COORDINATES VARCHAR, GENOMAD_N_GENES BIGINT,
            GENOMAD_GENETIC_CODE BIGINT, GENOMAD_VIRUS_SCORE DOUBLE,
            GENOMAD_FDR DOUBLE, GENOMAD_N_HALLMARKS BIGINT,
            GENOMAD_MARKER_ENRICHMENT DOUBLE, GENOMAD_TAXONOMY VARCHAR
        )
    """)

con.execute("""
    CREATE TABLE result AS
    SELECT
        $1 AS SAMPLE_ID,
        j.*,
        g.GENOMAD_LENGTH,
        g.GENOMAD_TOPOLOGY,
        g.GENOMAD_COORDINATES,
        g.GENOMAD_N_GENES,
        g.GENOMAD_GENETIC_CODE,
        g.GENOMAD_VIRUS_SCORE,
        g.GENOMAD_FDR,
        g.GENOMAD_N_HALLMARKS,
        g.GENOMAD_MARKER_ENRICHMENT,
        g.GENOMAD_TAXONOMY
    FROM joined_all j
    LEFT JOIN genomad g ON j.VNP_CONTIG_ID = g.CONTIG_ID
""", [sample_id])
con.execute("DROP TABLE joined_all")
con.execute("DROP TABLE genomad")

# ── Stats ────────────────────────────────────────────────────────────
n_result = con.execute("SELECT count(*) FROM result").fetchone()[0]
n_with_kallisto = con.execute("SELECT count(*) FROM result WHERE KALLISTO_DB_ID IS NOT NULL").fetchone()[0]
n_with_k2 = con.execute("SELECT count(*) FROM result WHERE K2_TAXONOMY_ID IS NOT NULL").fetchone()[0]
n_with_vnp = con.execute("SELECT count(*) FROM result WHERE VNP_CONTIG_ID IS NOT NULL").fetchone()[0]
n_with_genomad = con.execute("SELECT count(*) FROM result WHERE GENOMAD_VIRUS_SCORE IS NOT NULL").fetchone()[0]
print(f"\n  Final table: {n_result:,} rows", file=sys.stderr)
print(f"  With Kallisto hit: {n_with_kallisto:,} ({n_with_kallisto/n_result*100:.1f}%)", file=sys.stderr)
print(f"  With Kraken2 hit:  {n_with_k2:,} ({n_with_k2/n_result*100:.1f}%)", file=sys.stderr)
print(f"  With VNP hit:      {n_with_vnp:,} ({n_with_vnp/n_result*100:.1f}%)", file=sys.stderr)
print(f"  With geNomad hit:  {n_with_genomad:,} ({n_with_genomad/n_result*100:.1f}%)", file=sys.stderr)

# ── Write output ─────────────────────────────────────────────────────
con.execute("""
    COPY result TO $1 (FORMAT PARQUET, COMPRESSION ZSTD)
""", [output])
print(f"\nOutput written to: {output}", file=sys.stderr)

con.close()
CODE
  >>>

  output {
    File classifications_parquet = "~{out_filename}"
  }

  runtime {
    docker:           docker
    memory:           "16 GB"
    cpu:              1
    disks:            "local-disk ~{ceil((size(kallisto_summary, 'GB') + size(kraken2_reads, 'GB') + size(vnp_reads, 'GB') + size(genomad_virus_summary, 'GB')) * 4 + 10)} HDD"
    disk:             "~{ceil((size(kallisto_summary, 'GB') + size(kraken2_reads, 'GB') + size(vnp_reads, 'GB') + size(genomad_virus_summary, 'GB')) * 4 + 10)} GB"
    dx_instance_type: "mem2_ssd1_v2_x2"
    preemptible:      2
    maxRetries:       2
  }
}
