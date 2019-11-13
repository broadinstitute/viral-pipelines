task krakenuniq {
  Array[File] reads_unmapped_bam
  File        krakenuniq_db_tar_lz4  # krakenuniq/{database.kdb,taxonomy}
  File        krona_taxonomy_db_tgz  # taxonomy/taxonomy.tab

  String?     docker="quay.io/broadinstitute/viral-classify"

#  parameter_meta {
#    krakenuniq_db_tar_lz4: "stream"
#    krona_taxonomy_db_tgz: "stream"
#  }

  command {
    set -ex -o pipefail

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    DB_DIR=$(mktemp -d)

    metagenomics.py --version | tee VERSION

    # decompress DB to $DB_DIR
    read_utils.py extract_tarball \
      ${krakenuniq_db_tar_lz4} $DB_DIR \
      --loglevel=DEBUG
    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tgz} . \
      --loglevel=DEBUG &  # we don't need this until later

    # prep input and output file names
    OUT_READS=fnames_outreads.txt
    OUT_REPORTS=fnames_outreports.txt
    OUT_BASENAME=basenames_reports.txt
    for bam in ${sep=' ' reads_unmapped_bam}; do
      echo "$(basename $bam .bam).krakenuniq-reads.txt.gz" >> $OUT_READS
      echo "$(basename $bam .bam).krakenuniq" >> $OUT_BASENAME
      echo "$(basename $bam .bam).krakenuniq-summary_report.txt" >> $OUT_REPORTS
    done

    # execute on all inputs and outputs serially, but with a single
    # database load into ram
    metagenomics.py krakenuniq \
      $DB_DIR/krakenuniq \
      ${sep=' ' reads_unmapped_bam} \
      --outReads `cat $OUT_READS` \
      --outReport `cat $OUT_REPORTS` \
      --loglevel=DEBUG

    wait # for krona_taxonomy_db_tgz to download and extract

    # run single-threaded krona on up to nproc samples at once
    parallel -I ,, \
      "metagenomics.py krona \
        ,,-summary_report.txt \
        taxonomy \
        ,,-krona.html \
        --noRank --noHits --inputType krakenuniq \
        --loglevel=DEBUG" \
      ::: `cat $OUT_BASENAME`
  }

  output {
    Array[File] krakenuniq_classified_reads   = glob("*.krakenuniq-reads.txt.gz")
    Array[File] krakenuniq_summary_reports    = glob("*.krakenuniq-summary_report.txt")
    Array[File] krona_report_html             = glob("*.krakenuniq-krona.html")
    String      viralngs_version              = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "200 GB"
    cpu: 32
    dx_instance_type: "mem3_ssd1_v2_x32"
    preemptible: 0
  }
}

task build_krakenuniq_db {
  File        genome_fastas_tarball
  File        taxonomy_db_tarball
  String      db_basename

  Boolean?    subsetTaxonomy
  Int?        minimizerLen
  Int?        kmerLen
  Int?        maxDbSize
  String?     docker="quay.io/broadinstitute/viral-classify"

  command {
    set -ex -o pipefail

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    TAXDB_DIR=$(mktemp -d)
    FASTAS_DIR=$(mktemp -d)
    DB_DIR="$TMPDIR/${db_basename}"
    mkdir -p $DB_DIR

    metagenomics.py --version | tee VERSION

    # decompress input tarballs
    read_utils.py extract_tarball \
      ${genome_fastas_tarball} $FASTAS_DIR \
      --loglevel=DEBUG
    read_utils.py extract_tarball \
      ${taxonomy_db_tarball} $TAXDB_DIR \
      --loglevel=DEBUG

    # build database
    metagenomics.py krakenuniq_build \
      $DB_DIR --library $FASTAS_DIR --taxonomy $TAXDB_DIR \
      ${true='--subsetTaxonomy=' false='' subsetTaxonomy} \
      ${'--minimizerLen=' + minimizerLen} \
      ${'--kmerLen=' + kmerLen} \
      ${'--maxDbSize=' + maxDbSize} \
      --clean \
      --loglevel=DEBUG

    # tar it up
    tar -c ${db_basename} -C $TMP_DIR | zstd -19 > ${db_basename}.tar.zst
  }

  output {
    File        krakenuniq_db    = "${db_basename}.tar.zst"
    String      viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "200 GB"
    disks: "local-disk 350 HDD"
    cpu: 32
    dx_instance_type: "mem3_ssd1_v2_x32"
    preemptible: 0
  }
}

task kraken {
  Array[File] reads_unmapped_bam
  File        kraken_db_tar_lz4      # kraken/{database.kdb,taxonomy}
  File        krona_taxonomy_db_tgz  # taxonomy/taxonomy.tab

  String?     docker="quay.io/broadinstitute/viral-classify"

#  parameter_meta {
#    kraken_db_tar_lz4:     "stream"
#    krona_taxonomy_db_tgz: "stream"
#  }

  command {
    set -ex -o pipefail

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    DB_DIR=$(mktemp -d)

    # decompress DB to $DB_DIR
    read_utils.py extract_tarball \
      ${kraken_db_tar_lz4} $DB_DIR \
      --loglevel=DEBUG
    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tgz} . \
      --loglevel=DEBUG &  # we don't need this until later

    # prep input and output file names
    OUT_READS=fnames_outreads.txt
    OUT_REPORTS=fnames_outreports.txt
    OUT_BASENAME=basenames_reads.txt
    for bam in ${sep=' ' reads_unmapped_bam}; do
      echo "$(basename $bam .bam).kraken-reads" >> $OUT_BASENAME
      echo "$(basename $bam .bam).kraken-reads.txt.gz" >> $OUT_READS
      echo "$(basename $bam .bam).kraken-summary_report.txt" >> $OUT_REPORTS
    done

    metagenomics.py --version | tee VERSION

    # execute on all inputs and outputs serially, but with a single
    # database load into ram
    metagenomics.py kraken \
      $DB_DIR \
      ${sep=' ' reads_unmapped_bam} \
      --outReads `cat $OUT_READS` \
      --outReport `cat $OUT_REPORTS` \
      --loglevel=DEBUG

    wait # for krona_taxonomy_db_tgz to download and extract

    # run single-threaded krona on up to nproc samples at once
    parallel -I ,, \
      "metagenomics.py krona \
        ,,-summary_report.txt \
        taxonomy \
        ,,-krona.html \
        --noRank --noHits --inputType tsv \
        --loglevel=DEBUG" \
      ::: `cat $OUT_BASENAME`
  }

  output {
    Array[File] kraken_classified_reads = glob("*.kraken-reads.txt.gz")
    Array[File] kraken_summary_reports  = glob("*.kraken-summary_report.txt")
    Array[File] krona_report_html       = glob("*.kraken-krona.html")
    String      viralngs_version        = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "200 GB"
    cpu: 32
    dx_instance_type: "mem3_ssd1_v2_x32"
    preemptible: 0
  }
}

task build_kraken_db {
  File        genome_fastas_tarball
  File        taxonomy_db_tarball
  String      db_basename

  Boolean?    subsetTaxonomy
  Int?        minimizerLen
  Int?        kmerLen
  Int?        maxDbSize
  String?     docker="quay.io/broadinstitute/viral-classify"

  command {
    set -ex -o pipefail

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    TAXDB_DIR=$(mktemp -d)
    FASTAS_DIR=$(mktemp -d)
    DB_DIR="$TMPDIR/${db_basename}"
    mkdir -p $DB_DIR

    metagenomics.py --version | tee VERSION

    # decompress input tarballs
    read_utils.py extract_tarball \
      ${genome_fastas_tarball} $FASTAS_DIR \
      --loglevel=DEBUG
    read_utils.py extract_tarball \
      ${taxonomy_db_tarball} $TAXDB_DIR \
      --loglevel=DEBUG

    # build database
    metagenomics.py kraken_build \
      $DB_DIR --library $FASTAS_DIR --taxonomy $TAXDB_DIR \
      ${true='--subsetTaxonomy=' false='' subsetTaxonomy} \
      ${'--minimizerLen=' + minimizerLen} \
      ${'--kmerLen=' + kmerLen} \
      ${'--maxDbSize=' + maxDbSize} \
      --clean \
      --loglevel=DEBUG

    # tar it up
    tar -c ${db_basename} -C $TMP_DIR | zstd -19 > ${db_basename}.tar.zst
  }

  output {
    File        kraken_db        = "${db_basename}.tar.zst"
    String      viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "200 GB"
    disks: "local-disk 350 HDD"
    cpu: 32
    dx_instance_type: "mem3_ssd1_v2_x32"
    preemptible: 0
  }
}

task krona {
  File     classified_reads_txt_gz
  File     krona_taxonomy_db_tgz

  String?  input_type
  Int?     query_column
  Int?     taxid_column
  Int?     score_column
  Int?     magnitude_column
  String?  docker="quay.io/broadinstitute/viral-classify"

  String  input_basename = basename(classified_reads_txt_gz, ".txt.gz")

#  parameter_meta {
#    krona_taxonomy_db_tgz : "stream"
#  }

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tgz} . \
      --loglevel=DEBUG

    metagenomics.py --version | tee VERSION

    metagenomics.py krona \
      ${classified_reads_txt_gz} \
      taxonomy \
      ${input_basename}.html \
      ${'--inputType=' + input_type} \
      ${'--queryColumn=' + query_column} \
      ${'--taxidColumn=' + taxid_column} \
      ${'--scoreColumn=' + score_column} \
      ${'--magnitudeColumn=' + magnitude_column} \
      --noRank --noHits \
      --loglevel=DEBUG
  }

  output {
    File    krona_report_html  = "${input_basename}.html"
    String  viralngs_version   = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "4 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd2_v2_x2"
  }
}

task filter_bam_to_taxa {
  File           classified_bam
  File           classified_reads_txt_gz
  File           ncbi_taxonomy_db_tgz # nodes.dmp names.dmp
  Array[String]? taxonomic_names
  Array[Int]?    taxonomic_ids
  Boolean?       withoutChildren=false

  String?        docker="quay.io/broadinstitute/viral-classify"

  String         input_basename = basename(classified_bam, ".bam")

#  parameter_meta {
#    ncbi_taxonomy_db_tgz: "stream"
#  }

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    read_utils.py extract_tarball \
      ${ncbi_taxonomy_db_tgz} . \
      --loglevel=DEBUG

    TAX_NAMES="${sep=' ' taxonomic_names}"
    if [ -n "$TAX_NAMES" ]; then TAX_NAMES="--taxNames $TAX_NAMES"; fi

    TAX_IDs="${sep=' ' taxonomic_ids}"
    if [ -n "$TAX_IDs" ]; then TAX_IDs="--taxIDs $TAX_IDs"; fi

    metagenomics.py --version | tee VERSION

    metagenomics.py filter_bam_to_taxa \
      ${classified_bam} \
      ${classified_reads_txt_gz} \
      "${input_basename}_filtered.bam" \
      taxonomy/nodes.dmp \
      taxonomy/names.dmp \
      $TAX_NAMES \
      $TAX_IDs \
      ${true='--without-children' false='' withoutChildren} \
      --loglevel=DEBUG

      samtools view -c ${classified_bam} | tee classified_taxonomic_filter_read_count_pre
      samtools view -c "${input_basename}_filtered.bam" | tee classified_taxonomic_filter_read_count_post
  }

  output {
    File    bam_filtered_to_taxa                        = "${input_basename}_filtered.bam"
    Int     classified_taxonomic_filter_read_count_pre  = read_int("classified_taxonomic_filter_read_count_pre")
    Int     classified_taxonomic_filter_read_count_post = read_int("classified_taxonomic_filter_read_count_post")
    String  viralngs_version                            = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "4 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd2_v2_x2"
  }

}

task kaiju {
  File     reads_unmapped_bam
  File     kaiju_db_lz4  # <something>.fmi
  File     ncbi_taxonomy_db_tgz # taxonomy/{nodes.dmp, names.dmp}
  File     krona_taxonomy_db_tgz  # taxonomy/taxonomy.tab

  String?  docker="quay.io/broadinstitute/viral-classify"

  String   input_basename = basename(reads_unmapped_bam, ".bam")

#  parameter_meta {
#    kaiju_db_lz4            : "stream" # for DNAnexus, until WDL implements the File| type
#    ncbi_taxonomy_db_tgz    : "stream"
#    krona_taxonomy_db_tgz   : "stream"
#  }

  command {
    set -ex -o pipefail

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    DB_DIR=$(mktemp -d)

    lz4 -d ${kaiju_db_lz4} $DB_DIR/kaiju.fmi

    read_utils.py extract_tarball \
      ${ncbi_taxonomy_db_tgz} $DB_DIR \
      --loglevel=DEBUG

    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tgz} . \
      --loglevel=DEBUG

    metagenomics.py --version | tee VERSION

    # classify contigs
    metagenomics.py kaiju \
      ${reads_unmapped_bam} \
      $DB_DIR/kaiju.fmi \
      $DB_DIR/taxonomy \
      ${input_basename}.kaiju.summary_report.txt \
      --outReads ${input_basename}.kaiju.reads.txt.gz \
      --loglevel=DEBUG

    # run krona
    metagenomics.py krona \
      ${input_basename}.kaiju.summary_report.txt \
      taxonomy \
      ${input_basename}.kaiju-krona.html \
      --inputType kaiju \
      --noRank --noHits \
      --loglevel=DEBUG
  }

  output {
    File    kaiju_report       = "${input_basename}.kaiju-summary_report.txt"
    File    kaiju_reads        = "${input_basename}.kaiju-reads.txt.gz"
    File    krona_report_html  = "${input_basename}.kaiju-krona.html"
    String  viralngs_version   = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "100 GB"
    cpu: 16
    dx_instance_type: "mem3_ssd1_v2_x16"
  }
}
