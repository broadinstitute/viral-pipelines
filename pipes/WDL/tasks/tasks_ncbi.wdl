version 1.0

task download_fasta {
  input {
    String         out_prefix
    Array[String]+ accessions
    String         emailAddress

    String         docker = "quay.io/broadinstitute/viral-phylo:2.1.19.1"
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
  }
}

task download_annotations {
  input {
    Array[String]+ accessions
    String         emailAddress
    String         combined_out_prefix

    String         docker = "quay.io/broadinstitute/viral-phylo:2.1.19.1"
  }

  command {
    set -ex -o pipefail
    ncbi.py --version | tee VERSION
    ncbi.py fetch_feature_tables \
        ${emailAddress} \
        ./ \
        ${sep=' ' accessions} \
        --loglevel DEBUG
    ncbi.py fetch_fastas \
        ${emailAddress} \
        ./ \
        ${sep=' ' accessions} \
        --combinedFilePrefix "${combined_out_prefix}" \
        --loglevel DEBUG
  }

  output {
    File        combined_fasta   = "${combined_out_prefix}.fasta"
    Array[File] genomes_fasta    = glob("*.fasta")
    Array[File] features_tbl     = glob("*.tbl")
    String      viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "7 GB"
    cpu: 2
    dx_instance_type: "mem2_ssd1_v2_x2"
  }
}

task annot_transfer {
  meta {
    description: "Given a reference genome annotation in TBL format (e.g. from Genbank or RefSeq) and a multiple alignment of that reference to other genomes, produce new annotation files (TBL format with appropriate coordinate conversions) for each sequence in the multiple alignment. Resulting output can be fed to tbl2asn for Genbank submission."
  }

  input {
    File         multi_aln_fasta
    File         reference_fasta
    Array[File]+ reference_feature_table

    String       docker = "quay.io/broadinstitute/viral-phylo:2.1.19.1"
  }

  parameter_meta {
    multi_aln_fasta: {
      description: "multiple alignment of sample sequences against a reference genome -- for a single chromosome",
      patterns: ["*.fasta"]
    }
    reference_fasta: {
      description: "Reference genome, all segments/chromosomes in one fasta file. Headers must be Genbank accessions.",
      patterns: ["*.fasta"]
    }
    reference_feature_table: {
      description: "NCBI Genbank feature tables, one file for each segment/chromosome described in reference_fasta.",
      patterns: ["*.tbl"]
    }
  }

  command {
    set -e
    ncbi.py --version | tee VERSION
    ncbi.py tbl_transfer_prealigned \
        ${multi_aln_fasta} \
        ${reference_fasta} \
        ${sep=' ' reference_feature_table} \
        . \
        --oob_clip \
        --loglevel DEBUG
  }

  output {
    Array[File] transferred_feature_tables = glob("*.tbl")
    String      viralngs_version           = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
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

    String       docker = "quay.io/broadinstitute/viral-phylo:2.1.19.1"
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

  command {
    set -e
    ncbi.py --version | tee VERSION
    mkdir -p out
    ncbi.py tbl_transfer_multichr \
        "${genome_fasta}" \
        out \
        --ref_fastas ${sep=' ' reference_fastas} \
        --ref_tbls ${sep=' ' reference_feature_tables} \
        --oob_clip \
        --loglevel DEBUG
  }

  output {
    Array[File]+ genome_per_chr_tbls   = glob("out/*.tbl")
    Array[File]+ genome_per_chr_fastas = glob("out/*.fasta")
    String       viralngs_version      = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "15 GB"
    cpu: 4
    dx_instance_type: "mem2_ssd1_v2_x4"
    preemptible: 1
  }
}

task structured_comments {
  input {
    File   assembly_stats_tsv

    File?  filter_to_ids

    String docker = "quay.io/broadinstitute/viral-core:2.1.32"
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
  }
}

task rename_fasta_header {
  input {
    File   genome_fasta
    String new_name

    String out_basename = basename(genome_fasta, ".fasta")

    String docker = "quay.io/broadinstitute/viral-core:2.1.32"
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
    String? submitting_lab_name
    String? submitting_lab_addr
    String? originating_lab_addr
    String? authors
    String? fasta_filename
  }
  command <<<
    python3 << CODE
    import os.path
    import csv

    strict = ~{true="True" false="False" strict}

    # lookup table files to dicts
    sample_to_cmt = {}
    with open('~{structured_comments}', 'rt') as inf:
      for row in csv.DictReader(inf, delimiter='\t'):
        sample_to_cmt[row['SeqID']] = row

    out_headers = ('submitter', 'fn', 'covv_virus_name', 'covv_type', 'covv_passage', 'covv_collection_date', 'covv_location', 'covv_add_location', 'covv_host', 'covv_add_host_info', 'covv_gender', 'covv_patient_age', 'covv_patient_status', 'covv_specimen', 'covv_outbreak', 'covv_last_vaccinated', 'covv_treatment', 'covv_seq_technology', 'covv_assembly_method', 'covv_coverage', 'covv_orig_lab', 'covv_orig_lab_addr', 'covv_provider_sample_id', 'covv_subm_lab', 'covv_subm_lab_addr', 'covv_subm_sample_id', 'covv_authors', 'covv_comment', 'comment_type')

    with open('~{out_name}', 'wt') as outf:
      writer = csv.DictWriter(outf, out_headers, dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
      writer.writeheader()

      with open('~{source_modifier_table}', 'rt') as inf:
        for row in csv.DictReader(inf, delimiter='\t'):

          #covv_specimen
          if strict:
            valid_isolation_sources = ('Clinical', 'Environmental')
            assert row['isolation_source'] in valid_isolation_sources, f"Metadata error: 'isolation_source' not one of: {valid_isolation_sources}\n{row}"
            assert row['host'] == 'Homo sapiens' or row['isolation_source'] == 'Environmental', f"Metadata error: 'host' must be 'Homo sapiens' if 'isolation_source' is not 'Environmental'\n{row}"
            assert row['organism']         == 'Severe acute respiratory syndrome coronavirus 2', f"'organism' != 'Severe acute respiratory syndrome coronavirus 2'\n{row}"
            assert row['db_xref']          == 'taxon:2697049', f"Metadata error: 'db_xref' != 'taxon:2697049'\n{row}"

          # PHA4GE/INSDC controlled vocabulary for source information
          # from "Vocabulary" tab of this sheet:
          #   https://github.com/pha4ge/SARS-CoV-2-Contextual-Data-Specification/blob/master/PHA4GE%20SARS-CoV-2%20Contextual%20Data%20Template.xlsx
          gisaid_specimen_source = "unknown"
          if row['isolation_source'] == 'Clinical':
            gisaid_specimen_source = row.get("body_product",row.get("anatomical_material",row.get("anatomical_part","missing")))
          if row['isolation_source'] == 'Environmental':
            gisaid_specimen_source = row.get("environmental_material",row.get("environmental_site","missing"))

          writer.writerow({
            'covv_virus_name'     : 'hCoV-19/' +row['Sequence_ID'],
            'covv_collection_date': row['collection_date'],
            'covv_location'       : '~{continent} / ' + row['country'].replace(':',' /'),

            'covv_type'           : 'betacoronavirus',
            'covv_passage'        : 'Original',
            'covv_host'           : 'Human' if row['isolation_source'] == 'Clinical' else row['isolation_source'].replace("Environmental","Environment"),
            'covv_gender'         : 'unknown',
            'covv_patient_age'    : 'unknown',
            'covv_patient_status' : 'unknown',
            'covv_specimen'       : gisaid_specimen_source,

            'covv_assembly_method': sample_to_cmt[row['Sequence_ID']]['Assembly Method'],
            'covv_coverage'       : sample_to_cmt[row['Sequence_ID']]['Coverage'],
            'covv_seq_technology' : sample_to_cmt[row['Sequence_ID']]['Sequencing Technology'],

            'covv_orig_lab'       : row['collected_by'],
            'covv_subm_lab'       : "~{default='REQUIRED' submitting_lab_name}",
            'covv_authors'        : "~{default='REQUIRED' authors}",
            'covv_orig_lab_addr'  : "~{default='REQUIRED' originating_lab_addr}",
            'covv_subm_lab_addr'  : "~{default='REQUIRED' submitting_lab_addr}",

            'submitter'           : "~{default='REQUIRED' username}",
            'fn'                  : "~{default='REQUIRED' fasta_filename}",

            'covv_add_host_info'  : row.get('note',''),
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
    String      docker="quay.io/broadinstitute/viral-core:2.1.32"
  }
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
      # filename must be <libraryname>.<flowcell>.<lane>.cleaned.bam
      assert bam.endswith('.cleaned.bam'), "filename does not end in .cleaned.bam: {}".format(bam)
      bam_base = os.path.basename(bam)
      bam_parts = bam_base.split('.')
      assert len(bam_parts) >= 5, "filename does not conform to <libraryname>.<flowcell>.<lane>.cleaned.bam -- {}".format(bam_base)
      lib = '.'.join(bam_parts[:-4])
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
    disks: "local-disk 100 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}

task biosample_to_genbank {
  meta {
    description: "Prepares two input metadata files for Genbank submission based on a BioSample registration attributes table (attributes.tsv) since all of the necessary values are there. This produces both a Genbank Source Modifier Table and a BioSample ID map file that can be fed into the prepare_genbank task."
  }
  input {
    File    biosample_attributes
    Int     num_segments = 1
    Int     taxid

    File?   filter_to_ids

    Boolean s_dropout_note = true
    String  docker = "quay.io/broadinstitute/viral-phylo:2.1.19.1"
  }
  String base = basename(biosample_attributes, ".txt")
  command {
    set -ex -o pipefail
    ncbi.py --version | tee VERSION
    ncbi.py biosample_to_genbank \
        "${biosample_attributes}" \
        ${num_segments} \
        ${taxid} \
        "${base}".genbank.src \
        "${base}".biosample.map.txt \
        ${'--filter_to_samples ' + filter_to_ids} \
        --biosample_in_smt \
        --iso_dates \
        ~{true="--sgtf_override" false="" s_dropout_note} \
        --loglevel DEBUG
    cut -f 1 "${base}.genbank.src" | tail +2 > "${base}.sample_ids.txt"
  }
  output {
    File genbank_source_modifier_table = "${base}.genbank.src"
    File biosample_map                 = "${base}.biosample.map.txt"
    File sample_ids                    = "${base}.sample_ids.txt"
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}

task generate_author_sbt_file {
  meta {
    description: "Generate an NCBI-compatible author sbt file for submission of sequence data to GenBank. Accepts an author string, a defaults yaml file, and a jinja2-format template. Output is comparable to what is generated by http://www.ncbi.nlm.nih.gov/WebSub/template.cgi"
  }

  input {
    String? author_list
    File    j2_template
    File    defaults_yaml
    String? out_base = "authors"

    String  docker = "quay.io/broadinstitute/py3-bio:0.1.2"
  }

  parameter_meta {
    author_list: {
      description: "A string containing a space-delimited list with of author surnames separated by first name and (optional) middle initial. Ex. 'Lastname,Firstname, Last-hypenated,First,M., Last,F.'"
    }
    j2_template: {
      description: "A jinja2-format template for the sbt file expected by NCBI. Example: gs://pathogen-public-dbs/other-related/author_template.sbt.j2"
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

    python3 << CODE
    # generates an sbt file of the format returned by:
    # http://www.ncbi.nlm.nih.gov/WebSub/template.cgi
    import re
    # external dependencies
    import yaml # pyyaml
    from jinja2 import Template #jinja2

    def render_sbt(author_string, defaults_yaml=None, sbt_out_path="authors.sbt", j2_template="author_template.sbt.j2"):
        # simple version for only initials: #author_re=re.compile(r"\s?(?P<lastname>[\w\'\-\ ]+),(?P<initials>(?:[A-Z]\.){1,3})")
        author_re=re.compile(r"\s?(?P<lastname>[\w\'\-\ ]+),((?P<first>\w[\w\'\-\ ]+\.?),?|(?P<initials>(?:[A-Z]\.)+))(?P<initials_ext>(?:[A-Z]\.)*)")

        defaults_data = {}
        if defaults_yaml is not None:
            with open(defaults_yaml) as defaults_yaml:
                defaults_data = yaml.load(defaults_yaml, Loader=yaml.FullLoader)

        authors=[]
        submitter     = defaults_data.get("submitter")
        bioproject    = defaults_data.get("bioproject")
        title         = defaults_data.get("title")
        citation      = defaults_data.get("citation")
        authors_affil = defaults_data.get("authors_affil")
        
        authors.extend(defaults_data.get("authors_start",[]))
        
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

        for author in defaults_data.get("authors_last",[]):
            if author not in authors:
                authors.append(author)

        with open(j2_template) as sbt_template:
            template = Template(sbt_template.read())
        rendered = template.render( authors=authors, 
                                    authors_affil=authors_affil, 
                                    title=title, 
                                    submitter=submitter, 
                                    citation=citation, 
                                    bioproject=bioproject)
        
        #print(rendered)
        with open(sbt_out_path,"w") as sbt_out:
            sbt_out.write(rendered)

    render_sbt("~{author_list}", defaults_yaml="~{defaults_yaml}", sbt_out_path="~{out_base}.sbt", j2_template="~{j2_template}")
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
  }
}

task prepare_genbank {
  meta {
    description: "this task runs NCBI's tbl2asn"
  }

  input {
    Array[File]+ assemblies_fasta
    Array[File]  annotations_tbl
    File         authors_sbt
    File?        biosampleMap
    File?        genbankSourceTable
    File?        coverage_table
    String?      sequencingTech
    String?      comment
    String?      organism
    String?      molType
    String?      assembly_method
    String?      assembly_method_version

    Int?         machine_mem_gb
    String       docker = "quay.io/broadinstitute/viral-phylo:2.1.19.1"
  }

  parameter_meta {
    assemblies_fasta: {
      description: "Assembled genomes. One chromosome/segment per fasta file.",
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
    biosampleMap: {
      description: "A two column tab text file mapping sample IDs (first column) to NCBI BioSample accession numbers (second column). These typically take the format 'SAMN****' and are obtained by registering your samples first at https://submit.ncbi.nlm.nih.gov/",
      patterns: ["*.txt", "*.tsv"]
    }
    genbankSourceTable: {
      description: "A tab-delimited text file containing requisite metadata for Genbank (a 'source modifier table'). https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html",
      patterns: ["*.txt", "*.tsv"]
    }
    coverage_table: {
      description: "A two column tab text file mapping sample IDs (first column) to average sequencing coverage (second column, floating point number).",
      patterns: ["*.txt", "*.tsv"]
    }
    sequencingTech: {
      description: "The type of sequencer used to generate reads. NCBI has a controlled vocabulary for this value which can be found here: https://submit.ncbi.nlm.nih.gov/structcomment/nongenomes/"
    }
    organism: {
      description: "The scientific name for the organism being submitted. This is typically the species name and should match the name given by the NCBI Taxonomy database. For more info, see: https://www.ncbi.nlm.nih.gov/Sequin/sequin.hlp.html#Organism"
    }
    molType: {
      description: "The type of molecule being described. Any value allowed by the INSDC controlled vocabulary may be used here. Valid values are described at http://www.insdc.org/controlled-vocabulary-moltype-qualifier"
    }
    assembly_method: {
      description: "Very short description of the software approach used to assemble the genome. We typically provide a github link here. If this is specified, assembly_method_version should also be specified."
    }
    assembly_method_version: {
      description: "The version of the software used. If this is specified, assembly_method should also be specified."
    }
    comment: {
      description: "Optional comments that can be displayed in the COMMENT section of the Genbank record. This may include any disclaimers about assembly quality or notes about pre-publication availability or requests to discuss pre-publication use with authors."
    }

  }

  command {
    set -ex -o pipefail
    ncbi.py --version | tee VERSION
    cp ${sep=' ' annotations_tbl} .

    touch special_args
    if [ -n "${comment}" ]; then
      echo "--comment" >> special_args
      echo "${comment}" >> special_args
    fi
    if [ -n "${sequencingTech}" ]; then
      echo "--sequencing_tech" >> special_args
      echo "${sequencingTech}" >> special_args
    fi
    if [ -n "${organism}" ]; then
      echo "--organism" >> special_args
      echo "${organism}" >> special_args
    fi
    if [ -n "${molType}" ]; then
      echo "--mol_type" >> special_args
      echo "${molType}" >> special_args
    fi
    if [ -n "${assembly_method}" -a -n "${assembly_method_version}" ]; then
      echo "--assembly_method" >> special_args
      echo "${assembly_method}" >> special_args
      echo "--assembly_method_version" >> special_args
      echo "${assembly_method_version}" >> special_args
    fi
    if [ -n "${coverage_table}" ]; then
      echo -e "sample\taln2self_cov_median" > coverage_table.txt
      cat ${coverage_table} >> coverage_table.txt
      echo "--coverage_table" >> special_args
      echo coverage_table.txt >> special_args
    fi

    cat special_args | xargs -d '\n' ncbi.py prep_genbank_files \
        ${authors_sbt} \
        ${sep=' ' assemblies_fasta} \
        . \
        ${'--biosample_map ' + biosampleMap} \
        ${'--master_source_table ' + genbankSourceTable} \
        --loglevel DEBUG
    zip sequins_only.zip *.sqn
    zip all_files.zip *.sqn *.cmt *.gbf *.src *.fsa *.val
    mv errorsummary.val errorsummary.val.txt # to keep it separate from the glob
  }

  output {
    File        submission_zip           = "sequins_only.zip"
    File        archive_zip              = "all_files.zip"
    Array[File] sequin_files             = glob("*.sqn")
    Array[File] structured_comment_files = glob("*.cmt")
    Array[File] genbank_preview_files    = glob("*.gbf")
    Array[File] source_table_files       = glob("*.src")
    Array[File] fasta_per_chr_files      = glob("*.fsa")
    Array[File] validation_files         = glob("*.val")
    File        errorSummary             = "errorsummary.val.txt"
    String      viralngs_version         = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 3]) + " GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}

task package_genbank_ftp_submission {
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

    String  docker = "quay.io/broadinstitute/viral-baseimage:0.1.20"
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
          <Attribute name="wizard">BankIt_SARSCoV2_api</Attribute>
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
  }
}

task vadr {
  meta {
    description: "Runs NCBI's Viral Annotation DefineR for annotation and QC. Defaults here are for SARS-CoV-2 (see https://github.com/ncbi/vadr/wiki/Coronavirus-annotation), but VADR itself is applicable to a larger number of viral taxa (change the vadr_opts accordingly)."
  }
  input {
    File   genome_fasta
    String vadr_opts = "--glsearch -s -r --nomisc --mkey sarscov2 --alt_fail lowscore,fstukcnf,insertnn,deletinn"

    String docker = "staphb/vadr:1.2"
    Int    minlen = 50
    Int    maxlen = 30000
  }
  String out_base = basename(genome_fasta, '.fasta')
  command <<<
    set -e

    # remove terminal ambiguous nucleotides
    /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
      "~{genome_fasta}" \
      --minlen ~{minlen} \
      --maxlen ~{maxlen} \
      > "~{out_base}.fasta"

    # run VADR
    v-annotate.pl \
      ~{vadr_opts} \
      --mdir /opt/vadr/vadr-models/ \
      "~{out_base}.fasta" \
      "~{out_base}"

    # package everything for output
    tar -C "~{out_base}" -czvf "~{out_base}.vadr.tar.gz" .

    # prep alerts into a tsv file for parsing
    cat "~{out_base}/~{out_base}.vadr.alt.list" | cut -f 2 | tail -n +2 \
      > "~{out_base}.vadr.alerts.tsv"
    cat "~{out_base}.vadr.alerts.tsv" | wc -l > NUM_ALERTS
  >>>
  output {
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
    memory: "2 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}

