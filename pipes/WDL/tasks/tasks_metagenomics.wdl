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
    String docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
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

    String        docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
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

    String        docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
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
    String        docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
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

task kallisto {
  meta {
    description: "Runs Kallisto/kb pseudoalignment classification on one BAM or single FASTQ file."
  }

  input {
    File    reads_bam
    File    kallisto_index
    File    t2g
    String? sample_id

    Int     kmer_size = 31
    String  parity = "single"
    String  technology = "bulk"
    Boolean loom = false
    Boolean protein = false

    Int     machine_mem_gb = 32
    Int     cpu = 16
    String  docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
  }

  parameter_meta {
    reads_bam: {
      description: "Reads to classify. May be unaligned BAM or single FASTQ/FASTQ.GZ.",
      patterns: ["*.bam", "*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"],
      category: "required"
    }
    kallisto_index: {
      description: "Pre-built Kallisto index file.",
      patterns: ["*.idx", "*.index"],
      category: "required"
    }
    t2g: {
      description: "Transcript-to-gene mapping file. Two-column TSV with transcript IDs in the first column and collapsed target IDs in the second column.",
      patterns: ["*.tsv", "*.txt"],
      category: "required"
    }
    sample_id: {
      description: "Optional sample identifier to stamp into the long-form counts TSV. Defaults to the input reads basename with common BAM/FASTQ extensions removed.",
      category: "common"
    }
    kmer_size: {
      description: "K-mer size used by the Kallisto index. Must match the value used to build the index. Default 31.",
      category: "advanced"
    }
    parity: {
      description: "Library parity passed to kb count. Common values are single or paired. Default single.",
      category: "advanced"
    }
    technology: {
      description: "Technology preset passed to kb count. Default bulk.",
      category: "advanced"
    }
    loom: {
      description: "Also emit kb-generated Loom output when supported by the runtime image. Default false.",
      category: "advanced"
    }
    protein: {
      description: "Indicates that the Kallisto index was built from protein sequences. Default false.",
      category: "advanced"
    }
    machine_mem_gb: {
      description: "Memory allocation in GB.",
      category: "runtime"
    }
    cpu: {
      description: "Number of CPUs to request and pass to Kallisto/kb.",
      category: "runtime"
    }
    docker: {
      description: "viral-ngs classify-flavored image containing the refactored `metagenomics kallisto` command.",
      category: "runtime"
    }
  }

  String out_basename = sub(sub(sub(sub(sub(basename(reads_bam), "\\.bam$", ""), "\\.fastq\\.gz$", ""), "\\.fq\\.gz$", ""), "\\.fastq$", ""), "\\.fq$", "")
  String output_sample_id = select_first([sample_id, out_basename])
  String count_dir = out_basename + ".kallisto_count"
  Int disk_size = ceil((8 * size(reads_bam, "GB") + 2 * size(kallisto_index, "GB") + size(t2g, "GB") + 100) / 375.0) * 375

  command <<<
    set -ex -o pipefail

    if [ -z "${TMPDIR:-}" ]; then
      TMPDIR=$(pwd)
      export TMPDIR
    fi

    metagenomics --version | tee VERSION

    metagenomics kallisto \
      "~{reads_bam}" \
      --index "~{kallisto_index}" \
      --t2g "~{t2g}" \
      --kmer_len "~{kmer_size}" \
      --technology "~{technology}" \
      --parity "~{parity}" \
      --sample-id "~{output_sample_id}" \
      ~{true="--loom" false="" loom} \
      ~{true="--protein" false="" protein} \
      --out_dir "~{count_dir}" \
      --threads "~{cpu}" \
      --loglevel=DEBUG

    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
  >>>

  output {
    File        kallisto_counts_tsv  = "~{count_dir}/counts.tsv"
    File        kallisto_counts_h5ad = "~{count_dir}/adata.h5ad"
    Array[File] kallisto_loom_files  = glob("~{count_dir}/*.loom")
    Int         max_ram_gb           = ceil(read_float("MEM_BYTES")/1000000000)
    String      viralngs_version     = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem3_ssd1_v2_x16"
    preemptible: 2
    maxRetries: 2
  }
}

task build_kallisto_db {
  meta {
    description: "Builds a Kallisto index from reference sequences."
  }

  input {
    File    reference_sequences
    String? out_basename

    Boolean protein = false
    Int     kmer_size = 31
    String  workflow_type = "standard"

    Int     machine_mem_gb = 32
    Int     cpu = 16
    String  docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
  }

  parameter_meta {
    reference_sequences: {
      description: "FASTA file of reference sequences to index.",
      patterns: ["*.fasta", "*.fa", "*.fna", "*.fasta.gz", "*.fa.gz", "*.fna.gz"],
      category: "required"
    }
    out_basename: {
      description: "Output basename for the index. Defaults to the reference FASTA basename.",
      category: "common"
    }
    protein: {
      description: "Generate an index from amino-acid reference sequences. Default false.",
      category: "advanced"
    }
    kmer_size: {
      description: "K-mer size for the Kallisto index. Default 31.",
      category: "advanced"
    }
    workflow_type: {
      description: "Kallisto/kb workflow type. Valid values are standard, nac, kite, or custom. Default standard.",
      category: "advanced"
    }
    machine_mem_gb: {
      description: "Memory allocation in GB.",
      category: "runtime"
    }
    cpu: {
      description: "Number of CPUs to request and pass to Kallisto/kb.",
      category: "runtime"
    }
    docker: {
      description: "viral-ngs classify-flavored image containing the refactored `metagenomics kallisto_build` command.",
      category: "runtime"
    }
  }

  String reference_basename = sub(sub(sub(sub(basename(reference_sequences), "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", ""), "\\.fna$", "")
  String output_basename = select_first([out_basename, reference_basename])
  Int disk_size = ceil((4 * size(reference_sequences, "GB") + 100) / 375.0) * 375

  command <<<
    set -ex -o pipefail

    if [ -z "${TMPDIR:-}" ]; then
      TMPDIR=$(pwd)
      export TMPDIR
    fi

    metagenomics --version | tee VERSION

    if [[ "~{reference_sequences}" == *.gz ]]; then
      pigz -dc "~{reference_sequences}" > reference_sequences.fasta
      REF_FASTA=reference_sequences.fasta
    else
      REF_FASTA="~{reference_sequences}"
    fi

    metagenomics kallisto_build \
      "$REF_FASTA" \
      --index "~{output_basename}.kallisto.idx" \
      --workflow "~{workflow_type}" \
      --kmer_len "~{kmer_size}" \
      ~{true="--protein" false="" protein} \
      --threads "~{cpu}" \
      --loglevel=DEBUG

    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
  >>>

  output {
    File   kallisto_index  = "~{output_basename}.kallisto.idx"
    Int    max_ram_gb      = ceil(read_float("MEM_BYTES")/1000000000)
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem3_ssd1_v2_x16"
    preemptible: 0
    maxRetries: 1
  }
}

task kallisto_read_summary {
  meta {
    description: "Runs Kallisto extract internally and emits a read-level Kallisto summary TSV."
  }

  input {
    File           reads_bam
    File           kallisto_index
    File           t2g
    String?        sample_id
    Array[String]? target_ids
    File?          h5ad_file
    File?          id_to_taxon_map

    Int            threshold = 1
    String         taxonomy_level = "highest"
    Boolean        protein = false

    Int            machine_mem_gb = 32
    Int            cpu = 16
    String         docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
  }

  parameter_meta {
    reads_bam: {
      description: "Reads to extract from. May be unaligned BAM or single FASTQ/FASTQ.GZ.",
      patterns: ["*.bam", "*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"],
      category: "required"
    }
    kallisto_index: {
      description: "Kallisto index used to pseudoalign reads.",
      patterns: ["*.idx", "*.index"],
      category: "required"
    }
    t2g: {
      description: "Transcript-to-gene mapping file. Two-column TSV with transcript IDs in the first column and collapsed target IDs in the second column.",
      patterns: ["*.tsv", "*.txt"],
      category: "required"
    }
    sample_id: {
      description: "Optional sample identifier to stamp into the Kallisto summary TSV. Defaults to the input reads basename with common BAM/FASTQ extensions removed.",
      category: "common"
    }
    target_ids: {
      description: "Target transcript or collapsed hit IDs to extract. If omitted, h5ad_file is used to discover targets over threshold.",
      category: "common"
    }
    h5ad_file: {
      description: "Kallisto count h5ad file used to discover target IDs when target_ids is omitted.",
      patterns: ["*.h5ad"],
      category: "common"
    }
    id_to_taxon_map: {
      description: "Optional CSV/TSV mapping Kallisto hit IDs to taxonomy columns. When provided, taxonomy lineage and selected taxonomy name are added to the Kallisto summary TSV.",
      patterns: ["*.csv", "*.tsv", "*.csv.gz", "*.tsv.gz"],
      category: "common"
    }
    threshold: {
      description: "Minimum count threshold when extracting target IDs from h5ad_file. Default 1.",
      category: "advanced"
    }
    taxonomy_level: {
      description: "Taxonomy level to report from id_to_taxon_map. Valid values are highest or deepest. Default highest.",
      category: "advanced"
    }
    protein: {
      description: "Indicates that the Kallisto index was built from protein sequences. Default false.",
      category: "advanced"
    }
    machine_mem_gb: {
      description: "Memory allocation in GB.",
      category: "runtime"
    }
    cpu: {
      description: "Number of CPUs to request and pass to Kallisto/kb.",
      category: "runtime"
    }
    docker: {
      description: "viral-ngs classify-flavored image containing the refactored `metagenomics kallisto_extract` command.",
      category: "runtime"
    }
  }

  String out_basename = sub(sub(sub(sub(sub(basename(reads_bam), "\\.bam$", ""), "\\.fastq\\.gz$", ""), "\\.fq\\.gz$", ""), "\\.fastq$", ""), "\\.fq$", "")
  String output_sample_id = select_first([sample_id, out_basename])
  String extract_dir = out_basename + ".kallisto_extract"
  String summary_tsv_filename = output_sample_id + ".kallisto_summary.tsv"
  Int disk_size = ceil((8 * size(reads_bam, "GB") + 2 * size(kallisto_index, "GB") + size(t2g, "GB") + size(h5ad_file, "GB") + size(id_to_taxon_map, "GB") + 100) / 375.0) * 375

  command <<<
    set -ex -o pipefail

    if [ -z "${TMPDIR:-}" ]; then
      TMPDIR=$(pwd)
      export TMPDIR
    fi

    metagenomics --version | tee VERSION

    TARGET_IDS="~{sep=',' select_first([target_ids, []])}"
    H5AD_FILE="~{default='' h5ad_file}"
    ID_TO_TAXON_MAP="~{default='' id_to_taxon_map}"
    if [ -n "$TARGET_IDS" ]; then
      TARGET_ARGS=(--targets "$TARGET_IDS")
    elif [ -n "$H5AD_FILE" ]; then
      TARGET_ARGS=(--h5ad "$H5AD_FILE" --threshold "~{threshold}")
    else
      echo "Either target_ids or h5ad_file must be provided to kallisto_extract." >&2
      exit 1
    fi
    if [ -n "$ID_TO_TAXON_MAP" ]; then
      TAXONOMY_ARGS=(--id-to-tax-map "$ID_TO_TAXON_MAP")
    else
      TAXONOMY_ARGS=()
    fi

    metagenomics kallisto_extract \
      "~{reads_bam}" \
      --index "~{kallisto_index}" \
      --t2g "~{t2g}" \
      --out_dir "~{extract_dir}" \
      ~{true="--protein" false="" protein} \
      --sample-id "~{output_sample_id}" \
      "${TAXONOMY_ARGS[@]}" \
      --taxonomy-level "~{taxonomy_level}" \
      --threads "~{cpu}" \
      "${TARGET_ARGS[@]}" \
      --loglevel=DEBUG

    cp "~{extract_dir}/summary.tsv" "~{summary_tsv_filename}"

    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
  >>>

  output {
    File   summary_tsv          = "~{summary_tsv_filename}"
    Int    max_ram_gb           = ceil(read_float("MEM_BYTES")/1000000000)
    String viralngs_version     = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem3_ssd1_v2_x16"
    preemptible: 2
    maxRetries: 2
  }
}

task report_primary_kallisto_taxa {
  meta {
    description: "Interprets Kallisto counts TSV output and emits the primary contributing hit under a focal taxon."
  }

  input {
    File   kallisto_counts_tsv
    File?  id_to_taxon_map
    String focal_taxon = "Viruses"

    Int    machine_mem_gb = 16
    String docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
  }

  parameter_meta {
    kallisto_counts_tsv: {
      description: "Long-form Kallisto counts TSV emitted by the kallisto task.",
      patterns: ["*.tsv", "*.tsv.gz"],
      category: "required"
    }
    id_to_taxon_map: {
      description: "Optional CSV/TSV mapping Kallisto hit IDs to taxonomy columns.",
      patterns: ["*.csv", "*.tsv", "*.csv.gz", "*.tsv.gz"],
      category: "common"
    }
    focal_taxon: {
      description: "Taxonomic category to summarize. Default Viruses.",
      category: "common"
    }
    machine_mem_gb: {
      description: "Memory allocation in GB.",
      category: "runtime"
    }
    docker: {
      description: "viral-ngs classify-flavored image containing the refactored `metagenomics kallisto_top_taxa` command.",
      category: "runtime"
    }
  }

  String out_basename = sub(sub(basename(kallisto_counts_tsv), "\\.tsv\\.gz$", ""), "\\.tsv$", "")
  Int disk_size = ceil((4 * size(kallisto_counts_tsv, "GB") + size(id_to_taxon_map, "GB") + 100) / 375.0) * 375

  command <<<
    set -ex -o pipefail

    metagenomics --version | tee VERSION

    ID_TO_TAXON_MAP="~{default='' id_to_taxon_map}"
    if [ -n "$ID_TO_TAXON_MAP" ]; then
      TAXONOMY_ARGS=(--id-to-tax-map "$ID_TO_TAXON_MAP")
    else
      TAXONOMY_ARGS=()
    fi

    metagenomics kallisto_top_taxa \
      "~{kallisto_counts_tsv}" \
      "~{out_basename}.ranked_focal_report.tsv" \
      "${TAXONOMY_ARGS[@]}" \
      --target-taxon "~{focal_taxon}" \
      --loglevel=DEBUG

    head -2 "~{out_basename}.ranked_focal_report.tsv" | tail +2 > TOPROW
    cut -f 2 TOPROW > NUM_FOCAL
    cut -f 3 TOPROW > TOP_HIT_ID
    cut -f 4 TOPROW > TOP_HIT_NAME
    cut -f 5 TOPROW > TOP_HIT_LOWEST_TAXA_NAME
    cut -f 6 TOPROW > TOP_HIT_READS
    cut -f 7 TOPROW > PCT_OF_FOCAL
  >>>

  output {
    File   ranked_focal_report     = "~{out_basename}.ranked_focal_report.tsv"
    String focal_tax_name          = focal_taxon
    Int    total_focal_reads       = read_int("NUM_FOCAL")
    String top_hit_id              = read_string("TOP_HIT_ID")
    String top_hit_name            = read_string("TOP_HIT_NAME")
    String top_hit_lowest_tax_name = read_string("TOP_HIT_LOWEST_TAXA_NAME")
    Int    top_hit_reads           = read_int("TOP_HIT_READS")
    Float  percent_of_focal        = read_float("PCT_OF_FOCAL")
    String viralngs_version        = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GB"
    cpu: 16
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 2
    maxRetries: 2
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
    String docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
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
    String       docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
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
    String         docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
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
    String docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
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
