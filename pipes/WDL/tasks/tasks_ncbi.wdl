version 1.0

task download_fasta {
  input {
    String         out_prefix
    Array[String]+ accessions
    String         emailAddress

    String         docker = "quay.io/broadinstitute/viral-phylo:2.5.1.0"
  }

  command {
    ncbi.py --version | tee VERSION
    ncbi.py fetch_fastas \
        ${emailAddress} \
        . \
        ${sep=' ' accessions} \
        --combinedFilePrefix ${out_prefix} \
  }

  output {
    File   sequences_fasta  = "${out_prefix}.fasta"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "7 GB"
    cpu: 2
    dx_instance_type: "mem2_ssd1_v2_x2"
    maxRetries: 2
  }
}

task download_fasta_from_accession_string {
  meta {
    description: "Download GenBank sequences from a space-separated string of accessions (optionally with coverage specifications like 'accession:12.5x') and combine into a single fasta file."
  }

  input {
    String accession_string
    String out_prefix
    String emailAddress

    String docker = "quay.io/broadinstitute/viral-phylo:2.5.1.0"
  }

  command <<<
    set -e
    ncbi.py --version | tee VERSION

    # Parse accession string to extract just the accession IDs (strip coverage values if present)
    python3 <<CODE
    import re
    accession_string = "~{accession_string}"
    # Split by whitespace and extract accession (before colon if present)
    accessions = []
    for pair in accession_string.split():
        if ':' in pair:
            accessions.append(pair.split(':')[0])
        else:
            accessions.append(pair)

    with open('accessions.txt', 'w') as f:
        for acc in accessions:
            f.write(acc + '\n')
    CODE

    # Download sequences using extracted accessions
    ncbi.py fetch_fastas \
        ~{emailAddress} \
        . \
        $(cat accessions.txt | tr '\n' ' ') \
        --combinedFilePrefix ~{out_prefix}
  >>>

  output {
    File   sequences_fasta  = "~{out_prefix}.fasta"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "7 GB"
    cpu: 2
    dx_instance_type: "mem2_ssd1_v2_x2"
    maxRetries: 2
  }
}

task download_annotations {
  input {
    Array[String]+ accessions
    String         emailAddress
    String         combined_out_prefix

    String         docker = "quay.io/broadinstitute/viral-phylo:2.5.1.0"
  }

  command <<<
    set -ex -o pipefail
    ncbi.py --version | tee VERSION
    ncbi.py fetch_feature_tables \
        ~{emailAddress} \
        ./ \
        ~{sep=' ' accessions} \
        --loglevel DEBUG
    mkdir -p combined
    ncbi.py fetch_fastas \
        ~{emailAddress} \
        ./ \
        ~{sep=' ' accessions} \
        --combinedFilePrefix "combined/~{combined_out_prefix}" \
        --forceOverwrite \
        --loglevel DEBUG
  >>>

  output {
    File        combined_fasta   = "combined/~{combined_out_prefix}.fasta"
    Array[File] genomes_fasta    = glob("*.fasta")
    Array[File] features_tbl     = glob("*.tbl")
    String      viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "7 GB"
    cpu: 2
    dx_instance_type: "mem2_ssd1_v2_x2"
    maxRetries: 2
  }
}

task download_ref_genomes_from_tsv {
  input {
    File      ref_genomes_tsv    # [tax_id, isolate_prefix, taxname, colon_delim_accession_list]
    String    emailAddress

    String    docker = "quay.io/broadinstitute/viral-phylo:2.5.1.0"
  }

  command <<<
    set -ex -o pipefail
    ncbi.py --version | tee VERSION
    mkdir -p combined

    python3<<CODE
    import csv
    import phylo.genbank
    with open("~{ref_genomes_tsv}", 'rt') as inf:
      reader = csv.DictReader(inf, delimiter='\t',
        fieldnames=['tax_id', 'isolate_prefix', 'taxname', 'accessions']) # backwards support for headerless tsvs
      # for the future: batch all the downloads in a single call and re-organize output files afterwards
      for ref_genome in reader:
        if ref_genome['tax_id'] != 'tax_id': # skip header
          accessions = ref_genome['accessions'].split(':')
          phylo.genbank.fetch_fastas_from_genbank(
            accessionList=accessions,
            destinationDir=".",
            emailAddress="~{emailAddress}",
            forceOverwrite=True,
            combinedFilePrefix="combined/" + '-'.join(accessions),
            removeSeparateFiles=False,
            chunkSize=500)
    CODE
  >>>

  output {
    Array[File] ref_genomes_fastas  = glob("combined/*.fasta")
    Int         num_references      = length(ref_genomes_fastas)
  }

  runtime {
    docker: docker
    memory: "7 GB"
    cpu: 2
    dx_instance_type: "mem2_ssd1_v2_x2"
    maxRetries: 2
  }
}

task sequencing_platform_from_bam {
  input {
    File    bam

    String  docker = "quay.io/broadinstitute/viral-core:2.5.9"
  }

  command <<<
    set -ex -o pipefail
    samtools view -H "~{bam}" | grep '^@RG' | grep -o 'PL:[^[:space:]]*' | cut -d':' -f2 | sort | uniq > BAM_PLATFORMS
    if [ $(wc -l < BAM_PLATFORMS) -ne 1 ]; then
      echo "Other: hybrid" > GENBANK_SEQ_TECH
    elif grep -qi 'ILLUMINA' BAM_PLATFORMS; then
      echo "Illumina" > GENBANK_SEQ_TECH
    elif grep -qi 'ONT' BAM_PLATFORMS; then
      echo "Oxford Nanopore" > GENBANK_SEQ_TECH
    elif grep -qi 'PACBIO' BAM_PLATFORMS; then
      echo "PacBio" > GENBANK_SEQ_TECH
    elif grep -qi 'IONTORRENT' BAM_PLATFORMS; then
      echo "IonTorrent" > GENBANK_SEQ_TECH
    elif grep -qi 'SOLID' BAM_PLATFORMS; then
      echo "SOLiD" > GENBANK_SEQ_TECH
    elif grep -qi 'ULTIMA' BAM_PLATFORMS; then
      echo "Ultima" > GENBANK_SEQ_TECH
    elif grep -qi 'ELEMENT' BAM_PLATFORMS; then
      echo "Element" > GENBANK_SEQ_TECH
    elif grep -qi 'CAPILLARY' BAM_PLATFORMS; then
      echo "Sanger" > GENBANK_SEQ_TECH
    else
      echo "Haven't seen this one before!"
      exit 1
    fi
  >>>

  output {
    String  genbank_sequencing_technology  = read_string("GENBANK_SEQ_TECH")
  }

  runtime {
    docker: docker
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task align_and_annot_transfer_single {
  meta {
    description: "Given a reference genome annotation in TBL format (e.g. from Genbank or RefSeq) and new genome not in Genbank, produce new annotation files (TBL format with appropriate coordinate conversions) for the new genome. Resulting output can be fed to tbl2asn for Genbank submission."
  }

  input {
    File         genome_fasta
    Array[File]+ reference_fastas
    Array[File]+ reference_feature_tables

    String       out_basename = basename(genome_fasta, '.fasta')
    Int          machine_mem_gb = 30
    String       docker = "quay.io/broadinstitute/viral-phylo:2.5.1.0"
  }

  parameter_meta {
    genome_fasta: {
      description: "New genome, all segments/chromosomes in one fasta file. Must contain the same number of sequences as reference_fasta",
      patterns: ["*.fasta"]
    }
    reference_fastas: {
      description: "Reference genome, each segment/chromosome in a separate fasta file, in the exact same count and order as the segments/chromosomes described in genome_fasta. Headers must be Genbank accessions.",
      patterns: ["*.fasta"]
    }
    reference_feature_tables: {
      description: "NCBI Genbank feature table, each segment/chromosome in a separate TBL file, in the exact same count and order as the segments/chromosomes described in genome_fasta and reference_fastas. Accession numbers in the TBL files must correspond exactly to those in reference_fasta.",
      patterns: ["*.tbl"]
    }
  }

  command <<<
    set -e
    ncbi.py --version | tee VERSION
    mkdir -p out
    ncbi.py tbl_transfer_multichr \
        "~{genome_fasta}" \
        out \
        --ref_fastas ~{sep=' ' reference_fastas} \
        --ref_tbls ~{sep=' ' reference_feature_tables} \
        --oob_clip \
        --loglevel DEBUG
    cat out/*.tbl > "~{out_basename}.tbl"
  >>>

  output {
    File         feature_tbl           = "~{out_basename}.tbl"
    #Array[File]+ genome_per_chr_tbls   = glob("out/*.tbl")
    #Array[File]+ genome_per_chr_fastas = glob("out/*.fasta")
    String       viralngs_version      = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    cpu: 8
    dx_instance_type: "mem2_ssd1_v2_x4"
    preemptible: 1
    maxRetries: 2
  }
}

task structured_comments {
  input {
    File   assembly_stats_tsv

    File?  filter_to_ids

    String docker = "quay.io/broadinstitute/viral-core:2.5.9"
  }
  String out_base = basename(assembly_stats_tsv, '.txt')
  command <<<
    set -e

    python3 << CODE
    import util.file

    samples_to_filter_to = set()
    if "~{default='' filter_to_ids}":
        with open("~{default='' filter_to_ids}", 'rt') as inf:
            samples_to_filter_to = set(line.strip() for line in inf)

    out_headers = ('SeqID', 'StructuredCommentPrefix', 'Assembly Method', 'Coverage', 'Sequencing Technology', 'StructuredCommentSuffix')
    with open("~{out_base}.cmt", 'wt') as outf:
        outf.write('\t'.join(out_headers)+'\n')

        for row in util.file.read_tabfile_dict("~{assembly_stats_tsv}"):
            outrow = dict((h, row.get(h, '')) for h in out_headers)

            if samples_to_filter_to:
              if row['SeqID'] not in samples_to_filter_to:
                  continue

            if outrow['Coverage']:
              outrow['Coverage'] = "{}x".format(round(float(outrow['Coverage'])))
            outrow['StructuredCommentPrefix'] = 'Assembly-Data'
            outrow['StructuredCommentSuffix'] = 'Assembly-Data'
            outf.write('\t'.join(outrow[h] for h in out_headers)+'\n')
    CODE
  >>>
  output {
    File   structured_comment_table = "~{out_base}.cmt"
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task structured_comments_from_aligned_bam {
  input {
    File    aligned_bam
    String  assembly_method
    String  assembly_method_version

    String  out_basename = basename(aligned_bam, '.bam')
    Boolean is_genome_assembly = true
    Boolean sanitize_ids = true
    String  docker = "quay.io/broadinstitute/viral-core:2.5.9"
  }
  # see https://www.ncbi.nlm.nih.gov/genbank/structuredcomment/
  command <<<
    set -e -o pipefail

    cp "~{aligned_bam}" aligned.bam
    reports.py coverage_only aligned.bam coverage.txt
    cat coverage.txt

    samtools view -H "~{aligned_bam}" | grep '^@SQ' | grep -o 'SN:[^[:space:]]*' | cut -d':' -f2 | tee SEQ_IDS

    samtools view -H "~{aligned_bam}" | grep '^@RG' | grep -o 'PL:[^[:space:]]*' | cut -d':' -f2 | sort | uniq | tee BAM_PLATFORMS
    if [ $(wc -l < BAM_PLATFORMS) -ne 1 ]; then
      echo "Other: hybrid" > GENBANK_SEQ_TECH
    elif grep -qi 'ILLUMINA' BAM_PLATFORMS; then
      echo "Illumina" > GENBANK_SEQ_TECH
    elif grep -qi 'ONT' BAM_PLATFORMS; then
      echo "Oxford Nanopore" > GENBANK_SEQ_TECH
    elif grep -qi 'PACBIO' BAM_PLATFORMS; then
      echo "PacBio" > GENBANK_SEQ_TECH
    elif grep -qi 'IONTORRENT' BAM_PLATFORMS; then
      echo "IonTorrent" > GENBANK_SEQ_TECH
    elif grep -qi 'SOLID' BAM_PLATFORMS; then
      echo "SOLiD" > GENBANK_SEQ_TECH
    elif grep -qi 'ULTIMA' BAM_PLATFORMS; then
      echo "Ultima" > GENBANK_SEQ_TECH
    elif grep -qi 'ELEMENT' BAM_PLATFORMS; then
      echo "Element" > GENBANK_SEQ_TECH
    elif grep -qi 'CAPILLARY' BAM_PLATFORMS; then
      echo "Sanger" > GENBANK_SEQ_TECH
    else
      echo "Haven't seen this one before!"
      exit 1
    fi
    cat GENBANK_SEQ_TECH

    python3 << CODE
    import Bio.SeqIO
    import csv
    import re

    # get list of sequence IDs from BAM header
    with open("SEQ_IDS") as inf:
        seqids = set(line.strip() for line in inf)
        if ~{true="True" false="False" sanitize_ids}:
          seqids = set(re.sub(r'[^0-9A-Za-z!_-]', '-', x) for x in seqids)

    # get sequencing technology from BAM header    
    with open("GENBANK_SEQ_TECH", "rt") as inf:
        sequencing_tech = inf.read().strip()

    # this header has to be in this specific order -- don't reorder the columns!
    out_headers = ('SeqID', 'StructuredCommentPrefix', 'Assembly Method', "~{true='Genome ' false='' is_genome_assembly}Coverage", 'Sequencing Technology', 'StructuredCommentSuffix')
    with open("~{out_basename}.cmt", 'wt') as outf:
      outf.write('\t'.join(out_headers)+'\n')

      with open("coverage.txt", "rt") as inf:
        for row in csv.DictReader(inf, delimiter='\t'):
          if row.get('sample') and row.get('aln2self_cov_median') and row['sample'] in seqids:
            outrow = {
              'SeqID': row['sample'],
              'Assembly Method': "~{assembly_method} v. ~{assembly_method_version}",  # note: the <tool name> v. <version name> format is required by NCBI, don't remove the " v. "
              'Sequencing Technology': sequencing_tech,
              "~{true='Genome ' false='' is_genome_assembly}Coverage": "{}x".format(round(float(row['aln2self_cov_median']))),
              'StructuredCommentPrefix': "~{true='Genome-' false='' is_genome_assembly}Assembly-Data",
              'StructuredCommentSuffix': "~{true='Genome-' false='' is_genome_assembly}Assembly-Data",
            }
            outf.write('\t'.join(outrow[h] for h in out_headers)+'\n')
    CODE
  >>>
  output {
    File   structured_comment_file = "~{out_basename}.cmt"
  }
  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task prefix_fasta_header {
  input {
    File   genome_fasta
    String prefix
    String out_basename = basename(genome_fasta, ".fasta")
  }
  command <<<
    set -e
    python3 <<CODE
    with open('~{genome_fasta}', 'rt') as inf:
      with open('~{out_basename}.fasta', 'wt') as outf:
        for line in inf:
          if line.startswith('>'):
            line = ">{}{}\n".format('~{prefix}', line.rstrip()[1:])
          outf.write(line)
    CODE
  >>>
  output {
    File renamed_fasta = "~{out_basename}.fasta"
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task rename_fasta_header {
  input {
    File   genome_fasta
    String new_name

    String out_basename = basename(genome_fasta, ".fasta")

    String docker = "quay.io/broadinstitute/viral-core:2.5.9"
  }
  command {
    set -e
    file_utils.py rename_fasta_sequences \
      "~{genome_fasta}" "~{out_basename}.fasta" "~{new_name}"
  }
  output {
    File renamed_fasta = "~{out_basename}.fasta"
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task gisaid_meta_prep {
  input {
    File    source_modifier_table
    File    structured_comments
    String  out_name
    String  continent = "North America"
    Boolean strict = true
    String? username
    String  submitting_lab_name
    String? fasta_filename

    String  address_map = '{}'
    String  authors_map = '{}'
  }
  command <<<
    python3 << CODE
    import os.path
    import csv
    import json

    strict = ~{true="True" false="False" strict}

    # institutional mappings
    address_map = json.loads('~{address_map}')
    authors_map = json.loads('~{authors_map}')
    assert "~{submitting_lab_name}" in address_map, f"error: institution '~{submitting_lab_name}' not found in address_map"

    # lookup table files to dicts
    sample_to_cmt = {}
    with open('~{structured_comments}', 'rt') as inf:
      for row in csv.DictReader(inf, delimiter='\t'):
        sample_to_cmt[row['SeqID']] = row

    out_headers = ('submitter', 'fn', 'covv_virus_name', 'covv_type', 'covv_passage', 'covv_collection_date', 'covv_location', 'covv_add_location', 'covv_host', 'covv_add_host_info', 'covv_sampling_strategy', 'covv_gender', 'covv_patient_age', 'covv_patient_status', 'covv_specimen', 'covv_outbreak', 'covv_last_vaccinated', 'covv_treatment', 'covv_seq_technology', 'covv_assembly_method', 'covv_coverage', 'covv_orig_lab', 'covv_orig_lab_addr', 'covv_provider_sample_id', 'covv_subm_lab', 'covv_subm_lab_addr', 'covv_subm_sample_id', 'covv_authors', 'covv_comment', 'comment_type')

    with open('~{out_name}', 'w', newline='') as outf:
      writer = csv.DictWriter(outf, out_headers, dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
      writer.writeheader()

      with open('~{source_modifier_table}', 'rt') as inf:
        for row in csv.DictReader(inf, delimiter='\t'):

          isolation_source = row['isolation_source'].lower()
          #covv_specimen
          if strict:
            valid_isolation_sources = ('clinical', 'environmental')
            assert isolation_source in valid_isolation_sources, f"Metadata error: 'isolation_source' not one of: {valid_isolation_sources}\n{row}"
            assert row['host'] == 'Homo sapiens' or isolation_source == 'environmental', f"Metadata error: 'host' must be 'Homo sapiens' if 'isolation_source' is not 'Environmental'\n{row}"
            assert row['organism']         == 'Severe acute respiratory syndrome coronavirus 2', f"'organism' != 'Severe acute respiratory syndrome coronavirus 2'\n{row}"
            assert row['db_xref']          == 'taxon:2697049', f"Metadata error: 'db_xref' != 'taxon:2697049'\n{row}"

            collected_by = row['collected_by']
            assert collected_by in address_map, f"error: institution '{collected_by}' not found in address_map"
            assert collected_by in authors_map, f"error: institution '{collected_by}' not found in authors_map"

          # PHA4GE/INSDC controlled vocabulary for source information
          # from "Vocabulary" tab of this sheet:
          #   https://github.com/pha4ge/SARS-CoV-2-Contextual-Data-Specification/blob/master/PHA4GE%20SARS-CoV-2%20Contextual%20Data%20Template.xlsx
          gisaid_specimen_source = "unknown"
          if isolation_source == 'clinical':
            gisaid_specimen_source = row.get("body_product",row.get("anatomical_material",row.get("anatomical_part","missing")))
          if isolation_source == 'environmental':
            gisaid_specimen_source = row.get("environmental_material",row.get("environmental_site","missing"))

          writer.writerow({
            'covv_virus_name'     : 'hCoV-19/' +row['Sequence_ID'],
            'covv_collection_date': row['collection_date'],
            'covv_location'       : '~{continent} / ' + row['country'].replace(':',' /'),

            'covv_type'           : 'betacoronavirus',
            'covv_passage'        : 'Original',
            'covv_host'           : 'Human' if isolation_source == 'clinical' else isolation_source.replace("environmental","Environment"),
            'covv_add_host_info'  : 'unknown',
            'covv_gender'         : 'unknown',
            'covv_patient_age'    : 'unknown',
            'covv_patient_status' : 'unknown',
            'covv_specimen'       : gisaid_specimen_source.capitalize(), # capitalization of the first word seems to be the norm for GISAID

            'covv_assembly_method': sample_to_cmt[row['Sequence_ID']]['Assembly Method'],
            'covv_coverage'       : sample_to_cmt[row['Sequence_ID']]['Coverage'],
            'covv_seq_technology' : sample_to_cmt[row['Sequence_ID']]['Sequencing Technology'],

            'covv_orig_lab'       : row['collected_by'],
            'covv_subm_lab'       : "~{submitting_lab_name}",
            'covv_authors'        : authors_map.get(row['collected_by'], 'REQUIRED'),
            'covv_orig_lab_addr'  : address_map.get(row['collected_by'], 'REQUIRED'),
            'covv_subm_lab_addr'  : address_map.get("~{submitting_lab_name}", 'REQUIRED'),

            'submitter'           : "~{default='REQUIRED' username}",
            'fn'                  : "~{default='REQUIRED' fasta_filename}",

            'covv_sampling_strategy'  : row.get('note',''),
          })

    CODE
  >>>
  output {
    File meta_csv = "~{out_name}"
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task lookup_table_by_filename {
  input {
    String id
    File   mapping_tsv
    Int    return_col = 2

    String docker = "ubuntu"
  }
  command {
    set -e -o pipefail
    grep ^"~{id}" ~{mapping_tsv} | cut -f ~{return_col} > OUTVAL
  }
  output {
    String value = read_string("OUTVAL")
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task sra_meta_prep {
  meta {
    description: "Prepare tables for submission to NCBI's SRA database. This only works on bam files produced by illumina.py illumina_demux --append_run_id in viral-core."
  }
  input {
    Array[File] cleaned_bam_filepaths
    File        biosample_map
    Array[File] library_metadata
    String      platform
    String      instrument_model
    String      title
    Boolean     paired

    String      out_name = "sra_metadata.tsv"
    String      docker="quay.io/broadinstitute/viral-core:2.5.9"
  }
  Int disk_size = 100
  parameter_meta {
    cleaned_bam_filepaths: {
      description: "Unaligned bam files containing cleaned (submittable) reads.",
      localization_optional: true,
      stream: true,
      patterns: ["*.bam"]
    }
    biosample_map: {
      description: "Tab text file with a header and at least two columns named accession and sample_name. 'accession' maps to the BioSample accession number. Any samples without an accession will be omitted from output. 'sample_name' maps to the internal lab sample name used in filenames, samplesheets, and library_metadata files.",
      patterns: ["*.txt", "*.tsv"]
    }
    library_metadata: {
      description: "Tab text file with a header and at least six columns (sample, library_id_per_sample, library_strategy, library_source, library_selection, design_description). See 3rd tab of https://www.ncbi.nlm.nih.gov/core/assets/sra/files/SRA_metadata_acc_example.xlsx for controlled vocabulary and term definition.",
      patterns: ["*.txt", "*.tsv"]
    }
    platform: {
      description: "Sequencing platform (one of _LS454, ABI_SOLID, BGISEQ, CAPILLARY, COMPLETE_GENOMICS, HELICOS, ILLUMINA, ION_TORRENT, OXFORD_NANOPORE, PACBIO_SMRT)."
    }
    instrument_model: {
      description: "Sequencing instrument model (examples for platform=ILLUMINA: HiSeq X Five, HiSeq X Ten, Illumina Genome Analyzer, Illumina Genome Analyzer II, Illumina Genome Analyzer IIx, Illumina HiScanSQ, Illumina HiSeq 1000, Illumina HiSeq 1500, Illumina HiSeq 2000, Illumina HiSeq 2500, Illumina HiSeq 3000, Illumina HiSeq 4000, Illumina iSeq 100, Illumina NovaSeq 6000, Illumina MiniSeq, Illumina MiSeq, NextSeq 500, NextSeq 550)."
    }
    title: {
      description: "Descriptive sentence of the form <method> of <organism>, e.g. Metagenomic RNA-seq of SARS-CoV-2."
    }
  }
  command <<<
    python3 << CODE
    import os.path
    import csv
    import util.file

    # WDL arrays to python arrays
    bam_uris = list(x for x in '~{sep="*" cleaned_bam_filepaths}'.split('*') if x)
    library_metadata = list(x for x in '~{sep="*" library_metadata}'.split('*') if x)

    # lookup table files to dicts
    lib_to_bams = {}
    sample_to_biosample = {}
    for bam in bam_uris:
      # filename must be <libraryname>.<flowcell>.<lane>.cleaned.bam or <libraryname>.<flowcell>.<lane>.bam
      bam_base = os.path.basename(bam)
      bam_parts = bam_base.split('.')
      assert bam_parts[-1] == 'bam', "filename does not end in .bam -- {}".format(bam) 
      bam_parts = bam_parts[:-1]
      if bam_parts[-1] == 'cleaned':
        bam_parts = bam_parts[:-1]
      assert len(bam_parts) >= 3, "filename does not conform to <libraryname>.<flowcell>.<lane>.cleaned.bam -- {}".format(bam_base)
      lib = '.'.join(bam_parts[:-2]) # drop flowcell and lane
      lib_to_bams.setdefault(lib, [])
      lib_to_bams[lib].append(bam_base)
      print("debug: registering lib={} bam={}".format(lib, bam_base))
    with open('~{biosample_map}', 'rt') as inf:
      for row in csv.DictReader(inf, delimiter='\t'):
        sample_to_biosample[row['sample_name']] = row['accession']

    # set up SRA metadata table
    outrows = []
    out_headers = ['biosample_accession', 'library_ID', 'title', 'library_strategy', 'library_source', 'library_selection', 'library_layout', 'platform', 'instrument_model', 'design_description', 'filetype', 'assembly', 'filename']

    # iterate through library_metadata entries and produce an output row for each entry
    libs_written = set()
    for libfile in library_metadata:
      with open(libfile, 'rt') as inf:
        for row in csv.DictReader(inf, delimiter='\t'):
          lib = util.file.string_to_file_name("{}.l{}".format(row['sample'], row['library_id_per_sample']))
          biosample = sample_to_biosample.get(row['sample'],'')
          bams = lib_to_bams.get(lib,[])
          print("debug: sample={} lib={} biosample={}, bams={}".format(row['sample'], lib, biosample, bams))
          if biosample and bams and lib not in libs_written:
            libs_written.add(lib)
            outrows.append({
              'biosample_accession': sample_to_biosample[row['sample']],
              'library_ID': lib,
              'title': "~{title}",
              'library_strategy': row.get('library_strategy',''),
              'library_source': row.get('library_source',''),
              'library_selection': row.get('library_selection',''),
              'library_layout': '~{true="paired" false="single" paired}',
              'platform': '~{platform}',
              'instrument_model': '~{instrument_model}',
              'design_description': row.get('design_description',''),
              'filetype': 'bam',
              'assembly': 'unaligned',
              'files': lib_to_bams[lib],
            })
    assert outrows, "failed to prepare any metadata -- output is empty!"

    # find library with the most files and add col headers
    n_cols = max(len(row['files']) for row in outrows)
    for i in range(n_cols-1):
      out_headers.append('filename{}'.format(i+2))

    # write output file
    with open('~{out_name}', 'wt') as outf:
      outf.write('\t'.join(out_headers)+'\n')
      for row in outrows:
        row['filename'] = row['files'][0]
        for i in range(len(row['files'])):
          row['filename{}'.format(i+1)] = row['files'][i]
        outf.write('\t'.join(row.get(h,'') for h in out_headers)+'\n')
    CODE
  >>>
  output {
    File sra_metadata     = "~{out_name}"
    File cleaned_bam_uris = write_lines(cleaned_bam_filepaths)
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task biosample_to_table {
  meta {
    description: "Reformats a BioSample registration attributes table (attributes.tsv) for ingest into a Terra table."
  }
  input {
    File        biosample_attributes_tsv
    Array[File] raw_bam_filepaths
    File        demux_meta_json

    String  sample_table_name  = "sample"
    String  docker = "python:slim"
  }
  String  sanitized_id_col = "entity:~{sample_table_name}_id"
  String base = basename(basename(biosample_attributes_tsv, ".txt"), ".tsv")
  parameter_meta {
    raw_bam_filepaths: {
      description: "Unaligned bam files containing raw reads.",
      localization_optional: true,
      stream: true,
      patterns: ["*.bam"]
    }
  }
  command <<<
    set -ex -o pipefail
    python3 << CODE
    import os.path
    import csv
    import json

    # load demux metadata
    with open("~{demux_meta_json}", 'rt') as inf:
      demux_meta_by_file = json.load(inf)

    # load list of bams surviving filters
    bam_fnames = list(os.path.basename(x) for x in '~{sep="*" raw_bam_filepaths}'.split('*'))
    bam_fnames = list(x[:-len('.bam')] if x.endswith('.bam') else x for x in bam_fnames)
    bam_fnames = list(x[:-len('.cleaned')] if x.endswith('.cleaned') else x for x in bam_fnames)
    print("bam basenames ({}): {}".format(len(bam_fnames), bam_fnames))
    sample_to_sanitized = {demux_meta_by_file.get(x, {}).get('sample_original'): demux_meta_by_file.get(x, {}).get('sample') for x in bam_fnames}
    if None in sample_to_sanitized:
      del sample_to_sanitized[None]
    sample_names_seen = sample_to_sanitized.keys()
    print("samples seen ({}): {}".format(len(sample_names_seen), sorted(sample_names_seen)))

    # load biosample metadata
    biosample_attributes = []
    biosample_headers = ['biosample_accession']
    with open('~{biosample_attributes_tsv}', 'rt') as inf:
      for row in csv.DictReader(inf, delimiter='\t'):
        if row['sample_name'] in sample_names_seen and row['message'] == "Successfully loaded":
          row['biosample_accession'] = row.get('accession')
          row = dict({k:v for k,v in row.items() if v.strip().lower() not in ('missing', 'na', 'not applicable', 'not collected', '')})
          for k,v in row.items():
            if v and (k not in biosample_headers) and k not in ('message', 'accession'):
              biosample_headers.append(k)
          row['biosample_json'] = json.dumps({k: v for k,v in row.items() if k in biosample_headers})
          biosample_attributes.append(row)
    biosample_headers.append('biosample_json')

    print("biosample headers ({}): {}".format(len(biosample_headers), biosample_headers))
    print("biosample output rows ({})".format(len(biosample_attributes)))
    samples_seen_without_biosample = set(sample_names_seen) - set(row['sample_name'] for row in biosample_attributes)
    print("samples seen in bams without biosample entries ({}): {}".format(len(samples_seen_without_biosample), sorted(samples_seen_without_biosample)))

    # write reformatted table
    with open('~{base}.entities.tsv', 'w', newline='') as outf:
      writer = csv.DictWriter(outf, delimiter='\t', fieldnames=["~{sanitized_id_col}"]+biosample_headers, dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
      writer.writeheader()
      for row in biosample_attributes:
        outrow = {h: row.get(h, '') for h in biosample_headers}
        outrow["~{sanitized_id_col}"] = sample_to_sanitized[row['sample_name']]
        writer.writerow(outrow)
    CODE
  >>>
  output {
    File sample_meta_tsv = "~{base}.entities.tsv"
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task biosample_to_genbank {
  meta {
    description: "Prepares two input metadata files for Genbank submission based on a BioSample registration attributes table (attributes.tsv) since all of the necessary values are there. This produces both a Genbank Source Modifier Table and a BioSample ID map file that can be fed into the prepare_genbank task."
  }
  input {
    File    biosample_attributes
    Int     taxid
    Int     num_segments = 1
    String  biosample_col_for_fasta_headers = "sample_name"

    File?   filter_to_ids
    String? filter_to_accession
    Map[String,String] src_to_attr_map = {}
    String?  organism_name_override
    String?  sequence_id_override
    String?  isolate_prefix_override
    File?    source_overrides_json # optional set of key-value pairs to override source metadata

    Boolean sanitize_seq_ids = true

    String  out_basename = basename(basename(biosample_attributes, ".txt"), ".tsv")
    String  docker = "python:slim"
  }
  command <<<
    set -e
    python3<<CODE
    import csv
    import json
    import re
    import os.path

    header_key_map = {
        'Sequence_ID':'~{biosample_col_for_fasta_headers}',
        'BioProject':'bioproject_accession',
        'BioSample':'accession',
    }
    with open("~{write_json(src_to_attr_map)}", 'rt') as inf:
        more_header_key_map = json.load(inf)
    for k,v in more_header_key_map.items():
        header_key_map[k] = v
    print("header_key_map: {}".format(header_key_map))

    out_headers_total = ['Sequence_ID', 'isolate', 'collection_date', 'geo_loc_name', 'collected_by', 'isolation_source', 'organism', 'host', 'note', 'db_xref', 'BioProject', 'BioSample']
    samples_to_filter_to = set()
    if "~{default='' filter_to_ids}":
        with open("~{filter_to_ids}", 'rt') as inf:
            samples_to_filter_to = set(line.strip() for line in inf)
            print("filtering to samples: {}".format(samples_to_filter_to))
    only_accession = "~{default='' filter_to_accession}" if "~{default='' filter_to_accession}" else None
    if only_accession:
        print("filtering to biosample: {}".format(only_accession))

    # read entire tsv -> biosample_attributes, filtered to only the entries we keep
    with open("~{biosample_attributes}", 'rt') as inf_biosample:
      biosample_attributes_reader = csv.DictReader(inf_biosample, delimiter='\t')
      in_headers = biosample_attributes_reader.fieldnames
      if 'accession' not in in_headers:
        assert 'biosample_accession' in in_headers, "no accession column found in ~{biosample_attributes}"
        header_key_map['BioSample'] = 'biosample_accession'
      biosample_attributes = list(row for row in biosample_attributes_reader
        if row.get('message', 'Success').startswith('Success')
        and (not only_accession or row[header_key_map['BioSample']] == only_accession)
        and (not samples_to_filter_to or row[header_key_map['Sequence_ID']] in samples_to_filter_to))
      print("filtered to {} samples".format(len(biosample_attributes)))

    # clean out some na values
    for row in biosample_attributes:
      for k in ('strain', 'isolate', 'host', 'geo_loc_name', 'collection_date', 'serotype', 'food_origin'):
        if row.get(k, '').lower().strip() in ('na', 'not applicable', 'missing', 'not collected', 'not available', ''):
            row[k] = ''

    # override organism_name if provided (this allows us to submit Genbank assemblies for
    # specific species even though the metagenomic BioSample may have been registered with a different
    # species or none at all)
    if "~{default='' organism_name_override}":
        for row in biosample_attributes:
            row['organism'] = "~{default='' organism_name_override}"

    # handle special submission types: flu, sc2, noro, dengue
    special_bugs = ('Influenza A virus', 'Influenza B virus', 'Influenza C virus',
                    'Severe acute respiratory syndrome coronavirus 2',
                    'Norovirus', 'Dengue virus')
    for special in special_bugs:
        # sanitize organism name if it's a special one
        for row in biosample_attributes:
          if row['organism'].startswith(special):
              row['organism'] = special

        # enforce that special submissions are all the same special thing
        if any(row['organism'] == special for row in biosample_attributes):
          print("special organism found " + special)
          assert all(row['organism'] == special for row in biosample_attributes), "if any samples are {}, all samples must be {}".format(special, special)
          ### Influenza-specific requirements
          if special.startswith('Influenza'):
            print("special organism is Influenza A/B/C")
            # simplify isolate name
            if row.get('strain'):
              header_key_map['isolate'] = 'strain'
            if 'serotype' not in out_headers_total:
              out_headers_total.append('serotype')
            for row in biosample_attributes:
              # simplify isolate name more if it still looks structured with metadata (not allowed for flu submissions)
              if len(row[header_key_map.get('isolate', 'isolate')].split('/')) >= 2:
                row[header_key_map.get('isolate', 'isolate')] = row[header_key_map.get('isolate', 'isolate')].split('/')[-2]
              # populate serotype from name parsing
              match = re.search(r'\(([^()]+)\)+$', row['sample_name'])
              if match:
                  row['serotype'] = match.group(1)
                  print("found serotype {}". format(row['serotype']))
              # populate host field from name parsing if empty, override milk
              if not row.get('host','').strip():
                match = re.search(r'[^/]+/([^/]+)/[^/]+/[^/]+/[^/]+', row['sample_name'])
                if match:
                    row['host'] = match.group(1)
                    if row['host'] == 'bovine_milk':
                      row['host'] = 'Cattle'
                    assert 'host' in out_headers_total
              # override geo_loc_name if food_origin exists
              if row.get('food_origin','').strip():
                  print("overriding geo_loc_name with food_origin")
                  row['geo_loc_name'] = row['food_origin']

    # prepare output headers
    out_headers = list(h for h in out_headers_total if header_key_map.get(h,h) in in_headers)
    if 'db_xref' not in out_headers:
        out_headers.append('db_xref')
    if 'note' not in out_headers:
        out_headers.append('note')
    if 'serotype' in out_headers_total and 'serotype' not in out_headers:
        out_headers.append('serotype')
    if 'strain' not in out_headers and any(row['organism'].startswith('Influenza') for row in biosample_attributes):
        out_headers.append('strain')

    # manual overrides if provided
    if os.path.isfile("~{source_overrides_json}"):
      with open("~{source_overrides_json}", 'rt') as inf:
        for k,v in json.load(inf).items():
          print("overriding {} with {}".format(k,v))
          if k not in out_headers:
            out_headers.append(k)
          for row in biosample_attributes:
            row[k] = v

    # write output source modifier table
    with open("~{out_basename}.src", 'wt') as outf_smt:
      outf_smt.write('\t'.join(out_headers)+'\n')

      with open("~{out_basename}.sample_ids.txt", 'wt') as outf_ids:
          for row in biosample_attributes:
            # populate output row as a dict
            outrow = dict((h, row.get(header_key_map.get(h,h), '')) for h in out_headers)

            ### Taxon-specific genome naming rules go here
            if outrow['organism'].startswith('Influenza'):
              ### Influenza-specific isolate naming was handled above already and is required to be metadata-free
              # FTP submission pathway requires us to do the naming via the strain field
              type = outrow['organism'].split()[1] # A, B, C, D, etc
              state = outrow['geo_loc_name'].split(':')[1].strip() if ':' in outrow['geo_loc_name'] else outrow['geo_loc_name']
              year = outrow['collection_date'].split('-')[0]
              outrow['strain'] = '/'.join([type, state, outrow['isolate'], year])
              print("new strain name: {}".format(outrow['strain']))
            elif outrow['organism'].startswith('Special taxon with special naming rules'):
              ### Example special case here
              pass
            else:
              ### Most viral taxa can follow this generic model containing structured metadata
              # <short organism name>/<host lowercase>/Country/ST-Institution-LabID/Year
              # e.g. RSV-A/human/USA/MA-Broad-1234/2020
              # e.g. SARS-CoV-2/human/USA/MA-Broad-1234/2020
              # assume that the current isolate name has the structure:
              # <something wrong to be overwritten>/<host>/<country>/<state>-<institution>-<labid>/<year>
              print("old isolate name: {}".format(outrow['isolate']))
              year = outrow['collection_date'].split('-')[0]
              country = outrow['geo_loc_name'].split(':')[0]
              host = outrow['host'].lower()

              # host name mapping to the informal names used by NCBI in sequence titles for certain species
              host_to_informal_name_map = {
                'homo sapiens': 'human',
                'sus scrofa domesticus': 'swine',
                'sus scrofa': 'swine'
              }
              host = host_to_informal_name_map.get(host, host)

              if len(outrow['isolate'].split('/')) >= 2:
                state_inst_labid = outrow['isolate'].split('/')[-2]
              else:
                state_inst_labid = outrow['Sequence_ID']
              name_prefix = outrow['organism'].split('(')[0].split('/')[0].strip()
              if name_prefix.lower().endswith('refseq'):
                name_prefix = name_prefix[:-6].strip()
              if "~{default='' isolate_prefix_override}".strip():
                name_prefix = "~{default='' isolate_prefix_override}".strip()
              outrow['isolate'] = '/'.join([name_prefix, host, country, state_inst_labid, year])
              print("new isolate name: {}".format(outrow['isolate']))

            # some fields are not allowed to be empty
            if not outrow.get('geo_loc_name'):
                outrow['geo_loc_name'] = 'missing'
            if not outrow.get('host'):
                outrow['host'] = 'not applicable'

            # custom db_xref/taxon
            outrow['db_xref'] = "taxon:{}".format(~{taxid})

            # load the purpose of sequencing (or if not, the purpose of sampling) in the note field
            outrow['note'] = row.get('purpose_of_sequencing', row.get('purpose_of_sampling', ''))

            # overwrite sequence ID if requested
            if "~{default='' sequence_id_override}":
                outrow['Sequence_ID'] = "~{default='' sequence_id_override}"

            # sanitize sequence IDs to match fasta headers
            if "~{sanitize_seq_ids}".lower() == 'true':
                outrow['Sequence_ID'] = re.sub(r'[^0-9A-Za-z!_-]', '-', outrow['Sequence_ID'])

            # write isolate name for this sample (it will overwrite and only keep the last one if there are many)
            with open("isolate_name.str", 'wt') as outf_isolate:
                print("writing isolate name: {}".format(outrow['isolate']))
                outf_isolate.write(outrow['isolate']+'\n')

            # write entry for this sample
            sample_name = outrow['Sequence_ID']
            if ~{num_segments}>1:
                for i in range(~{num_segments}):
                    outrow['Sequence_ID'] = "{}-{}".format(sample_name, i+1)
                    outf_smt.write('\t'.join(outrow[h] for h in out_headers)+'\n')
                    outf_ids.write(outrow['Sequence_ID']+'\n')
            else:
                outf_smt.write('\t'.join(outrow[h] for h in out_headers)+'\n')
                outf_ids.write(outrow['Sequence_ID']+'\n')
    CODE
  >>>
  output {
    File   genbank_source_modifier_table = "~{out_basename}.src"
    File   sample_ids                    = "~{out_basename}.sample_ids.txt"
    String isolate_name                  = read_string("isolate_name.str")
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task generate_author_sbt_file {
  meta {
    description: "Generate an NCBI-compatible author sbt file for submission of sequence data to GenBank. Accepts an author string, a defaults yaml file, and a jinja2-format template. Output is comparable to what is generated by http://www.ncbi.nlm.nih.gov/WebSub/template.cgi"
  }

  input {
    String? author_list
    File    j2_template
    File?   defaults_yaml
    String  out_base = "authors"

    String  docker = "quay.io/broadinstitute/py3-bio:0.1.2"
  }

  parameter_meta {
    author_list: {
      description: "A string containing a space-delimited list with of author surnames separated by first name and (optional) middle initial. Ex. 'Lastname,Firstname, Last-hypenated,First,M., Last,F.'"
    }
    j2_template: {
      description: "an sbt file (optionally) with Jinja2 variables to be filled in based on values present in author_sbt_defaults_yaml, if provided. If no yaml is provided, this file is passed through verbatim. Example: gs://pathogen-public-dbs/other-related/author_template.sbt.j2"
    }
    defaults_yaml: {
      description: "A YAML file with default values to use for the submitter, submitter affiliation, and author affiliation. Optionally including authors at the start and end of the author_list. Example: gs://pathogen-public-dbs/other-related/default_sbt_values.yaml",
      patterns: ["*.yaml","*.yml"]
    }
    out_base: {
      description: "prefix to use for the generated *.sbt output file"
    }
  }
  
  command <<<
    set -e

    # blank yaml file to be used if the optional input is not specified
    touch blank.yml

    python3 << CODE
    # generates an sbt file of the format returned by:
    # http://www.ncbi.nlm.nih.gov/WebSub/template.cgi
    import re
    import shutil
    # external dependencies
    import yaml # pyyaml
    from jinja2 import Template #jinja2

    def render_sbt(author_string, defaults_yaml=None, sbt_out_path="authors.sbt", j2_template="author_template.sbt.j2"):
        # simple version for only initials: #author_re=re.compile(r"\s?(?P<lastname>[\w\'\-\ ]+),(?P<initials>(?:[A-Z]\.){1,3})")
        author_re=re.compile(r"\s?(?P<lastname>[\w\'\-\ ]+),((?P<first>\w[\w\'\-\ ]+\.?),?|(?P<initials>(?:[A-Z]\.)+))(?P<initials_ext>(?:[A-Z]\.)*)")

        authors=[]
        defaults_data_last_authors=[]
        defaults_data = {}

        authors_affil = None
        submitter     = None
        bioproject    = None
        title         = None
        citation      = None

        if defaults_yaml is not None:
            with open(defaults_yaml) as defaults_yaml:
                defaults_data = yaml.load(defaults_yaml, Loader=yaml.FullLoader)

                if defaults_data is not None:
                    submitter     = defaults_data.get("submitter")
                    bioproject    = defaults_data.get("bioproject")
                    title         = defaults_data.get("title")
                    citation      = defaults_data.get("citation")
                    authors_affil = defaults_data.get("authors_affil")
                    
                    defaults_data_authors = defaults_data.get("authors_start",[])
                    for author in defaults_data_authors:
                        authors.extend(author)

                    defaults_data_last_authors = defaults_data.get("authors_last",[])
                    for author in defaults_data_last_authors:
                        last_authors.append(author)
        
        for author_match in author_re.finditer(author_string):
            author = {}
            lastname=author_match.group("lastname")
            initials=[]
            if author_match.group("initials"):
                initials.extend(list(filter(None,author_match.group("initials").split("."))))
            if author_match.group("initials_ext"):
                initials.extend(list(filter(None,author_match.group("initials_ext").split("."))))

            first=""
            if author_match.group("first"):
                first=author_match.group("first")
            else:
                first=initials[0]+"."
            author["last"]     = author_match.group("lastname")
            author["first"]    = first
            author["initials"]   = ".".join(initials[1:]) if not author_match.group("first") else ".".join(initials)
            author["initials"]   = author["initials"]+"." if len(author["initials"])>0 else author["initials"]
            
            if author not in authors: # could use less exact match
                authors.append(author)

        for author in defaults_data_last_authors:
            if author not in authors:
                authors.append(author)

        jinja_rendering_kwargs={}
        if authors_affil is not None:
            jinja_rendering_kwargs["authors_affil"]=authors_affil
        if title is not None:
            jinja_rendering_kwargs["title"]=title
        if submitter is not None:
            jinja_rendering_kwargs["submitter"]=submitter
        if citation is not None:
            jinja_rendering_kwargs["citation"]=citation
        if bioproject is not None:
            jinja_rendering_kwargs["bioproject"]=bioproject

        if len(authors) >= 1 or len(jinja_rendering_kwargs) >= 1:
            with open(j2_template) as sbt_template:
                template = Template(sbt_template.read())

            rendered = template.render( authors=authors, 
                                        **jinja_rendering_kwargs)
        
            #print(rendered)
            with open(sbt_out_path,"w") as sbt_out:
                sbt_out.write(rendered)
        else:
            # if no authors were specified, simply copy the template to the output
            shutil.copyfile(j2_template, sbt_out_path)

    render_sbt("~{author_list}", defaults_yaml="~{default='blank.yml' defaults_yaml}", sbt_out_path="~{out_base}.sbt", j2_template="~{j2_template}")
    CODE
  >>>
  output {
    File sbt_file = "~{out_base}.sbt"
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task table2asn {
  meta {
    description: "This task runs NCBI's table2asn, the artist formerly known as tbl2asn"
  }

  input {
    File         assembly_fasta
    File         annotations_tbl
    File         authors_sbt
    File?        source_modifier_table
    File?        structured_comment_file
    String?      comment
    String       organism
    String       mol_type = "cRNA"
    Int          genetic_code = 1

    String       out_basename = basename(assembly_fasta, ".fasta")
    Int          machine_mem_gb = 8
    String       docker = "quay.io/broadinstitute/viral-phylo:2.5.1.0"  # this could be a simpler docker image, we don't use anything beyond table2asn itself
  }
  Int disk_size = 50

  parameter_meta {
    assembly_fasta: {
      description: "Assembled genome. All chromosomes/segments in one file.",
      patterns: ["*.fasta"]
    }
    annotations_tbl: {
      description: "Gene annotations in TBL format, one per fasta file. Filename basenames must match the assemblies_fasta basenames. These files are typically output from the ncbi.annot_transfer task.",
      patterns: ["*.tbl"]
    }
    authors_sbt: {
      description: "A genbank submission template file (SBT) with the author list, created at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/",
      patterns: ["*.sbt"]
    }
    source_modifier_table: {
      description: "A tab-delimited text file containing requisite metadata for Genbank (a 'source modifier table'). https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html",
      patterns: ["*.src"]
    }
    structured_comment_file: {
      patterns: ["*.cmt"]
    }
    organism: {
      description: "The scientific name for the organism being submitted. This is typically the species name and should match the name given by the NCBI Taxonomy database. For more info, see: https://www.ncbi.nlm.nih.gov/Sequin/sequin.hlp.html#Organism"
    }
    mol_type: {
      description: "The type of molecule being described. Any value allowed by the INSDC controlled vocabulary may be used here. Valid values are described at http://www.insdc.org/controlled-vocabulary-moltype-qualifier"
    }
    comment: {
      description: "Optional comments that can be displayed in the COMMENT section of the Genbank record. This may include any disclaimers about assembly quality or notes about pre-publication availability or requests to discuss pre-publication use with authors."
    }
  }

  command <<<
    set -ex -o pipefail
    table2asn -version | cut -f 2 -d ' ' > TABLE2ASN_VERSION
    cp "~{assembly_fasta}" "~{out_basename}.fsa" # input fasta must be in CWD so output files end up next to it
    touch "~{out_basename}.val"  # this file isn't produced if no errors/warnings

    table2asn \
      -t "~{authors_sbt}" \
      -i "~{out_basename}.fsa" \
      -f "~{annotations_tbl}" \
      ~{'-w "' + structured_comment_file + '"'} \
      -j '[gcode=~{genetic_code}][moltype=~{mol_type}][organism=~{organism}]' \
      ~{'-src-file "' + source_modifier_table + '"'} \
      ~{'-y "' + comment + '"'} \
      -a s -V vb

    set +x
    echo "table2asn complete"
    cat "~{out_basename}.val" | { grep -vi '^Info:' || test $? = 1; } | tee "~{out_basename}.val.no_info"
    cat "~{out_basename}.val.no_info" | { grep -vi '^Warning: valid' || test $? = 1; } | tee "~{out_basename}.val.errors_only"
  >>>

  output {
    File          genbank_submission_sqn   = "~{out_basename}.sqn"
    File          genbank_preview_file     = "~{out_basename}.gbf"
    File          genbank_validation_file  = "~{out_basename}.val"
    Array[String] table2asn_errors         = read_lines("~{out_basename}.val")
    String        table2asn_version        = read_string("TABLE2ASN_VERSION")
    Boolean       table2asn_passing        = length(read_lines("~{out_basename}.val.errors_only")) == 0
  }

  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
  }
}

task package_special_genbank_ftp_submission {
  meta {
    description: "Prepares a zip and xml file for FTP-based NCBI Genbank submission according to instructions at https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/genbank/SARS-CoV-2/."
  }
  input {
    File   sequences_fasta
    File   structured_comment_table
    File   source_modifier_table
    File   author_template_sbt
    String submission_name
    String submission_uid
    String spuid_namespace
    String account_name
    String wizard="BankIt_SARSCoV2_api"

    String  docker = "quay.io/broadinstitute/viral-baseimage:0.3.0"
  }
  command <<<
    set -e

    # make the submission zip file
    cp "~{sequences_fasta}" sequence.fsa
    cp "~{structured_comment_table}" comment.cmt
    cp "~{source_modifier_table}" source.src
    cp "~{author_template_sbt}" template.sbt
    zip "~{submission_uid}.zip" sequence.fsa comment.cmt source.src template.sbt

    # make the submission xml file
    SUB_NAME="~{submission_name}"
    ACCT_NAME="~{account_name}"
    SPUID="~{submission_uid}"
    cat << EOF > submission.xml
    <?xml version="1.0"?>
    <Submission>
      <Description>
        <Comment>$SUB_NAME</Comment>
        <Organization type="center" role="owner">
          <Name>$ACCT_NAME</Name>
        </Organization>
      </Description>
      <Action>
        <AddFiles target_db="GenBank">
          <File file_path="$SPUID.zip">
            <DataType>genbank-submission-package</DataType>
          </File>
          <Attribute name="wizard">~{wizard}</Attribute>
          <Identifier>
            <SPUID spuid_namespace="~{spuid_namespace}">$SPUID</SPUID>
          </Identifier>
        </AddFiles>
      </Action>
    </Submission>
    EOF

    # make the (empty) ready file
    touch submit.ready
  >>>
  output {
    File submission_zip = "~{submission_uid}.zip"
    File submission_xml = "submission.xml"
    File submit_ready   = "submit.ready"
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task genbank_special_taxa {
  meta {
    description: "this task tells you if you have a special taxon that NCBI treats differently"
  }

  input {
    Int     taxid
    File    taxdump_tgz
    File    vadr_by_taxid_tsv # "gs://pathogen-public-dbs/viral-references/annotation/vadr/vadr-by-taxid.tsv"
    String  docker = "quay.io/broadinstitute/viral-classify:2.5.1.0"
  }

  command <<<
    set -e

    # unpack the taxdump tarball
    mkdir -p taxdump
    read_utils.py extract_tarball "~{taxdump_tgz}" taxdump

    python3 << CODE
    import csv
    import metagenomics
    import tarfile
    import urllib.request
    import json
    import re
    import sys
    taxid = ~{taxid}

    # load taxdb and retrieve full hierarchy leading to this taxid
    taxdb = metagenomics.TaxonomyDb(tax_dir="taxdump", load_nodes=True, load_names=True, load_gis=False)
    ancestors = taxdb.get_ordered_ancestors(taxid)
    this_and_ancestors = [taxid] + ancestors

    # Genbank prohibits normal submissions for SC2, Flu A/B/C, Noro, and Dengue
    table2asn_prohibited = {
      11320: "Influenza A virus",
      11520: "Influenza B virus",
      11552: "Influenza C virus",
      2697049: "Severe acute respiratory syndrome coronavirus 2",
      11983: "Norovirus",
      3052464: "Dengue virus"
    }
    prohibited = any(node in table2asn_prohibited for node in this_and_ancestors)
    with open("table2asn_allowed.boolean", "wt") as outf:
      outf.write("false" if prohibited else "true")
    with open("genbank_submission_mechanism.str", "wt") as outf:
      if any(node == 11320 for node in this_and_ancestors):
        outf.write("InfluenzaA")
      elif any(node == 11520 for node in this_and_ancestors):
        outf.write("InfluenzaB")
      elif any(node == 11552 for node in this_and_ancestors):
        outf.write("InfluenzaC")
      elif any(node == 2697049 for node in this_and_ancestors):
        outf.write("SARS-CoV-2")
      elif any(node == 11983 for node in this_and_ancestors):
        outf.write("Norovirus")
      elif any(node == 3052464 for node in this_and_ancestors):
        outf.write("Dengue")
      else:
        outf.write("table2asn")

    # Genbank requires certain fields to be populated for Flu A and Dengue
    if any(node == 11320 for node in this_and_ancestors):
      # flu A needs subtype specified in the serotype column
      with open("genbank_source_overrides.json", "wt") as outf:
        # This taxid matches a specific named strain
        match = re.search(r'\(([^()]+)\)+$', taxdb.names[taxid])
        if match:
          subtype = match.group(1)
          print("found Influenza A subtype {}". format(subtype))
          json.dump({'serotype':subtype}, outf)
        else:
          # Influenza A has a number of taxids at the subtype level
          match = re.search(r'(\w+)\s+subtype$', taxdb.names[taxid])
          if match:
            subtype = match.group(1)
            print("found Influenza A subtype {}". format(subtype))
            json.dump({'serotype':subtype}, outf)
          else:
            print("failed to find Influenza A subtype from taxid {}, taxname {}".format(taxid, taxdb.names[taxid]))
            sys.exit(1)
    elif any(node == 3052464 for node in this_and_ancestors):
      # dengue needs serotype specified in the genotype column
      with open("genbank_source_overrides.json", "wt") as outf:
        match = re.search(r'(\d)$', taxdb.names[taxid])
        if match:
          serotype = match.group(1)
          print("found Dengue serotype {}". format(serotype))
          json.dump({'genotype':serotype}, outf)
        else:
          print("failed to find Dengue serotype from taxid {}, taxname {}".format(taxid, taxdb.names[taxid]))
          sys.exit(1)

    # VADR is an annotation tool that supports SC2, Flu A/B/C/D, Noro, Dengue, RSV A/B, MPXV, etc
    # https://github.com/ncbi/vadr/wiki/Available-VADR-model-files
    # Note this table includes some taxa that are subtaxa of others and in those cases, it is ORDERED from
    # more specific to less specific (e.g. noro before calici, dengue before flavi, sc2 before corona)
    # so use the *first* hit in the table.
    # tsv header: tax_id, taxon_name, min_seq_len, max_seq_len, vadr_min_ram_gb, vadr_opts, vadr_model_tar_url
    with open("~{vadr_by_taxid_tsv}", 'rt') as inf:
      vadr_supported = list(row for row in csv.DictReader(inf, delimiter='\t'))

    out_vadr_supported = False
    out_vadr_cli_options = ""
    out_vadr_model_tar_url = ""
    out_min_genome_length = 0
    out_max_genome_length = 1000000
    out_vadr_taxid = 0
    out_vadr_min_ram_gb = 8
    out_vadr_model_tar_subdir = ""
    for row in vadr_supported:
      if any(node == int(row['tax_id']) for node in this_and_ancestors):
        out_vadr_taxid = int(row['tax_id'])
        out_vadr_supported = True
        out_vadr_cli_options = row['vadr_opts']
        out_vadr_model_tar_url = row['vadr_model_tar_url']
        if row['min_seq_len']:
          out_min_genome_length = int(row['min_seq_len'])
        if row['max_seq_len']:
          out_max_genome_length = int(row['max_seq_len'])
        if row['vadr_min_ram_gb']:
          out_vadr_min_ram_gb = int(row['vadr_min_ram_gb'])
        if row['vadr_model_tar_subdir']:
          out_vadr_model_tar_subdir = row['vadr_model_tar_subdir']
        break
    with open("vadr_supported.boolean", "wt") as outf:
      outf.write("true" if out_vadr_supported else "false")
    with open("vadr_cli_options.string", "wt") as outf:
      outf.write(out_vadr_cli_options)
    with open("min_genome_length.int", "wt") as outf:
      outf.write(str(out_min_genome_length))
    with open("max_genome_length.int", "wt") as outf:
      outf.write(str(out_max_genome_length))
    with open("vadr_taxid.int", "wt") as outf:
      outf.write(str(out_vadr_taxid))
    with open("vadr_min_ram_gb.int", "wt") as outf:
      outf.write(str(out_vadr_min_ram_gb))
    with open("vadr_model_tar_subdir.str", "wt") as outf:
      outf.write(out_vadr_model_tar_subdir)

    if out_vadr_model_tar_url:
      urllib.request.urlretrieve(out_vadr_model_tar_url, "vadr_model-~{taxid}.tar.gz")

    CODE
  >>>

  output {
    Boolean  table2asn_allowed     = read_boolean("table2asn_allowed.boolean")
    String   genbank_submission_mechanism = read_string("genbank_submission_mechanism.str")
    File?    genbank_source_overrides_json = "genbank_source_overrides.json"
    Boolean  vadr_supported        = read_boolean("vadr_supported.boolean")
    String   vadr_cli_options      = read_string("vadr_cli_options.string")
    File?    vadr_model_tar        = "vadr_model-~{taxid}.tar.gz"
    String   vadr_model_tar_subdir = read_string("vadr_model_tar_subdir.str")
    Int      vadr_taxid            = read_int("vadr_taxid.int")
    Int      vadr_min_ram_gb       = read_int("vadr_min_ram_gb.int")
    Int      min_genome_length     = read_int("min_genome_length.int")
    Int      max_genome_length     = read_int("max_genome_length.int")
  }

  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task vadr {
  meta {
    description: "Runs NCBI's Viral Annotation DefineR for annotation and QC."
  }
  input {
    File    genome_fasta
    String? vadr_opts
    Int?    minlen
    Int?    maxlen
    File?   vadr_model_tar
    String? vadr_model_tar_subdir

    String out_basename = basename(genome_fasta, '.fasta')
    String docker = "mirror.gcr.io/staphb/vadr:1.6.4"
    Int    mem_size = 16  # the RSV model in particular seems to consume 15GB RAM
    Int    cpus = 4
  }
  command <<<
    set -e
    # unpack custom VADR models or use default
    if [ -n "~{vadr_model_tar}" ]; then
      mkdir -p vadr-untar
      tar -C vadr-untar -xzvf "~{vadr_model_tar}"
      # just link the model subdirectory, not the outer wrapper
      ln -s vadr-untar/*/ vadr-models
    else
      # use default (distributed with docker image) models
      ln -s /opt/vadr/vadr-models vadr-models
    fi
    if [ -n "~{default='' vadr_model_tar_subdir}" ]; then
      VADR_MODEL_DIR="vadr-models/~{default='' vadr_model_tar_subdir}"
    else
      VADR_MODEL_DIR="vadr-models"
    fi

    # remove terminal ambiguous nucleotides
    /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
      "~{genome_fasta}" \
      ~{'--minlen ' + minlen} \
      ~{'--maxlen ' + maxlen} \
      > "~{out_basename}.fasta"

    # run VADR
    v-annotate.pl \
      ~{default='' vadr_opts} \
      --split --cpu `nproc` \
      --mdir "$VADR_MODEL_DIR" \
      "~{out_basename}.fasta" \
      "~{out_basename}"

    # package everything for output
    tar -C "~{out_basename}" -czvf "~{out_basename}.vadr.tar.gz" .

    # get the gene annotations (feature table)
    # the vadr.fail.tbl sometimes contains junk / invalid content that needs to be filtered out
    cat "~{out_basename}/~{out_basename}.vadr.pass.tbl" \
        "~{out_basename}/~{out_basename}.vadr.fail.tbl" \
       | sed '/Additional note/,$d' \
        > "~{out_basename}.tbl"

    # prep alerts into a tsv file for parsing
    cat "~{out_basename}/~{out_basename}.vadr.alt.list" | cut -f 5 | tail -n +2 \
      > "~{out_basename}.vadr.alerts.tsv"
    cat "~{out_basename}.vadr.alerts.tsv" | wc -l > NUM_ALERTS

    # record peak memory usage
    set +o pipefail
    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } | tee MEM_BYTES
  >>>
  output {
    File                 feature_tbl = "~{out_basename}.tbl"
    Int                  num_alerts  = read_int("NUM_ALERTS")
    File                 alerts_list = "~{out_basename}/~{out_basename}.vadr.alt.list"
    Array[String]        alerts      = read_lines("~{out_basename}.vadr.alerts.tsv")
    File                 outputs_tgz = "~{out_basename}.vadr.tar.gz"
    Boolean              pass        = num_alerts==0
    String               vadr_docker = docker
    Int                  max_ram_gb  = ceil(read_float("MEM_BYTES")/1000000000)
  }
  runtime {
    docker: docker
    memory: mem_size + " GB"
    cpu: cpus
    dx_instance_type: "mem2_ssd1_v2_x4"
    maxRetries: 2
  }
}

task package_genbank_submissions {
  input {
    Array[String]   genbank_file_manifest
    Array[File]     genbank_submit_files

    String          spuid_base
    String          spuid_namespace
    String          account_name
    File            authors_sbt

    String          docker = "python:slim"
    Int             mem_size = 2
    Int             cpus = 1
  }
  command <<<
    set -e

    # bring everything into CWD -- at this point we are sure all files are uniquely named / no collisions
    cp "~{sep='" "' genbank_submit_files}" .
    cp "~{authors_sbt}" template.sbt

    python3 << CODE
    import json
    import zipfile
    import csv

    # initialize counters and file lists for each submission category
    counts_by_type = {}
    files_by_type = {}
    groups_total = list(m+"_"+v for m in ('table2asn', 'InfluenzaA', 'InfluenzaB', 'InfluenzaC', 'SARS-CoV-2', 'Norovirus', 'Dengue') for v in ('clean', 'warnings'))
    for group in groups_total:
      counts_by_type[group] = 0
      files_by_type[group] = []

    # read manifest from genbank_single
    for genome in json.loads('[~{sep="," genbank_file_manifest}]'):
      group = genome['submission_type'] + '_' + ('clean' if genome['validation_passing'] else 'warnings')
      counts_by_type[group] += 1
      files_by_type[group].extend(genome['files'])

    # function for writing submission.xml files
    def write_submission_xml(spuid, spuid_namespace, account_name, wizard=None):
      fname = spuid + ".xml"
      submission_name = spuid
      with open(fname, 'wt') as outf:
        outf.write('<?xml version="1.0"?>\n')
        outf.write('<Submission>\n')
        outf.write('  <Description>\n')
        outf.write('    <Comment>{}</Comment>\n'.format(submission_name))
        outf.write('    <Organization type="center" role="owner">\n')
        outf.write('      <Name>{}</Name>\n'.format(account_name))
        outf.write('    </Organization>\n')
        outf.write('  </Description>\n')
        outf.write('  <Action>\n')
        outf.write('    <AddFiles target_db="GenBank">\n')
        outf.write('      <File file_path="{}.zip">\n'.format(spuid))
        outf.write('        <DataType>genbank-submission-package</DataType>\n')
        outf.write('      </File>\n')
        if wizard:
          outf.write('      <Attribute name="wizard">{}</Attribute>\n'.format(wizard))
        outf.write('      <Identifier>\n')
        outf.write('        <SPUID spuid_namespace="{}">{}</SPUID>\n'.format(spuid_namespace, spuid))
        outf.write('      </Identifier>\n')
        outf.write('    </AddFiles>\n')
        outf.write('  </Action>\n')
        outf.write('</Submission>\n')

    # write and bundle outputs
    for group in groups_total:
      with open("count_"+group+".int", 'wt') as outf:
        outf.write(str(counts_by_type[group]))
      if counts_by_type[group]:
        if group.startswith('table2asn'):
          # bundle the sqn files together and that's it
          with zipfile.ZipFile("~{spuid_base}_" + group + ".zip", 'w') as zf:
            for fname in files_by_type[group]:
              zf.write(fname)
        else:
          # for special submissions, merge each individual file type
          with open("sequences.fsa", 'wt') as out_fsa:
            with open("comment.cmt", 'wt') as out_cmt:
              with open("source.src", 'wt') as out_src:
                src_header = ""
                cmt_header = ""
                for fname in files_by_type[group]:
                  if fname.endswith('.fsa') or fname.endswith('.fasta'):
                    with open(fname) as inf:
                      out_fsa.write(inf.read())
                  elif fname.endswith('.cmt'):
                    with open(fname) as inf:
                      header = inf.readline()
                      if cmt_header:
                        assert header == cmt_header
                      else:
                        cmt_header = header
                        out_cmt.write(header)
                      out_cmt.write(inf.read())
                  elif fname.endswith('.src'):
                    with open(fname) as inf:
                      header = inf.readline()
                      if src_header:
                        assert header == src_header
                      else:
                        src_header = header
                        out_src.write(header)
                      out_src.write(inf.read())
          # make FTP zips and xmls
          if group.startswith('SARS-CoV-2'):
            with zipfile.ZipFile("~{spuid_base}_FTP_" + group + ".zip", 'w') as zf:
              for fname in ('sequences.fsa', 'comment.cmt', 'source.src', 'template.sbt'):
                zf.write(fname)
            wizard = "BankIt_SARSCoV2_api"
            write_submission_xml("~{spuid_base}_FTP_" + group, "~{spuid_namespace}", "~{account_name}", wizard=wizard)
          elif group.startswith('Influenza'):
            with zipfile.ZipFile("~{spuid_base}_FTP_" + group + ".zip", 'w') as zf:
              for fname in ('sequences.fsa', 'comment.cmt', 'source.src', 'template.sbt'):
                zf.write(fname)
            wizard = "BankIt_influenza_api"
            write_submission_xml("~{spuid_base}_FTP_" + group, "~{spuid_namespace}", "~{account_name}", wizard=wizard)
            # make a different src file for web after making FTP zip
            # because for some reason, you can't use the same source.src file for web and FTP for flu
            # strain column is required for FTP and prohibited for web
            with open('source.src', 'rt') as inf:
              reader = csv.DictReader(inf, delimiter='\t')
              out_header = list(h for h in reader.fieldnames if h != 'strain')
              rows = []
              for row in reader:
                del row['strain']
                rows.append(row)
            with open('source.src', 'wt') as outf:
              writer = csv.DictWriter(outf, delimiter='\t', fieldnames=out_header, dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
              writer.writeheader()
              for row in rows:
                writer.writerow(row)

          # make web portal zips
          with zipfile.ZipFile("~{spuid_base}_" + group + ".zip", 'w') as zf:
            for fname in ('sequences.fsa', 'comment.cmt', 'source.src', 'template.sbt'):
              zf.write(fname)
    CODE
  >>>

  output {
    Array[File] ftp_submission_files  = glob("~{spuid_base}_FTP_*_clean.*")
    File? submit_sqns_clean_zip       = "~{spuid_base}_table2asn_clean.zip"
    File? submit_sqns_warnings_zip    = "~{spuid_base}_table2asn_warnings.zip"
    Int   num_sqns_clean              = read_int("count_table2asn_clean.int")
    Int   num_sqns_warnings           = read_int("count_table2asn_warnings.int")
    File? submit_fluA_clean_zip       = "~{spuid_base}_InfluenzaA_clean.zip"
    File? submit_fluA_warnings_zip    = "~{spuid_base}_InfluenzaA_warnings.zip"
    Int   num_fluA_clean              = read_int("count_InfluenzaA_clean.int")
    Int   num_fluA_warnings           = read_int("count_InfluenzaA_warnings.int")
    File? submit_fluB_clean_zip       = "~{spuid_base}_InfluenzaB_clean.zip"
    File? submit_fluB_warnings_zip    = "~{spuid_base}_InfluenzaB_warnings.zip"
    Int   num_fluB_clean              = read_int("count_InfluenzaB_clean.int")
    Int   num_fluB_warnings           = read_int("count_InfluenzaB_warnings.int")
    File? submit_fluC_clean_zip       = "~{spuid_base}_InfluenzaC_clean.zip"
    File? submit_fluC_warnings_zip    = "~{spuid_base}_InfluenzaC_warnings.zip"
    Int   num_fluC_clean              = read_int("count_InfluenzaC_clean.int")
    Int   num_fluC_warnings           = read_int("count_InfluenzaC_warnings.int")
    File? submit_sc2_clean_zip        = "~{spuid_base}_SARS-CoV-2_clean.zip"
    File? submit_sc2_warnings_zip     = "~{spuid_base}_SARS-CoV-2_warnings.zip"
    Int   num_sc2_clean               = read_int("count_SARS-CoV-2_clean.int")
    Int   num_sc2_warnings            = read_int("count_SARS-CoV-2_warnings.int")
    File? submit_noro_clean_zip       = "~{spuid_base}_Norovirus_clean.zip"
    File? submit_noro_warnings_zip    = "~{spuid_base}_Norovirus_warnings.zip"
    Int   num_noro_clean              = read_int("count_Norovirus_clean.int")
    Int   num_noro_warnings           = read_int("count_Norovirus_warnings.int")
    File? submit_dengue_clean_zip     = "~{spuid_base}_Dengue_clean.zip"
    File? submit_dengue_warnings_zip  = "~{spuid_base}_Dengue_warnings.zip"
    Int   num_dengue_clean            = read_int("count_Dengue_clean.int")
    Int   num_dengue_warnings         = read_int("count_Dengue_warnings.int")
  }

  runtime {
    docker: docker
    memory: mem_size + " GB"
    cpu: cpus
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 1
  }
}
