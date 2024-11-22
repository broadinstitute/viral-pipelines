version 1.0

task download_fasta {
  input {
    String         out_prefix
    Array[String]+ accessions
    String         emailAddress

    String         docker = "quay.io/broadinstitute/viral-phylo:2.3.6.0"
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

task download_annotations {
  input {
    Array[String]+ accessions
    String         emailAddress
    String         combined_out_prefix

    String         docker = "quay.io/broadinstitute/viral-phylo:2.3.6.0"
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

task download_viral_taxon_genomes {
  # Download viral genomes for a given taxid from NCBI using the NCBI Datasets tool
  #   see:
  #     https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/virus/
  #     https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/

  input {
    String  taxid_or_organism_name
    String? date_released_after
    String? date_released_before
    String? geo_location_country
    String? geo_location_us_state

    String? segment_num_or_name
    
    String? host_taxid_or_organism_name     = "Homo sapiens"
    Int?    max_num_ambig_bases
    Int?    min_sequence_length_end_to_end

    String?  pango_lineage
    String?  genotype_or_subtype_name_query

    Boolean refseq_only                     = false
    Boolean exclude_lab_passaged            = true
    Boolean exclude_vaccine_strains         = true
    Boolean exclude_environmental_specimens = true

    String? ncbi_api_key

    String metadata_columns = "accession,length,isolate-collection-date,geo-location,geo-region,geo-state,biosample-acc,sra-accs,completeness,isolate-lineage-source,release-date,host-common-name,virus-common-name,virus-name,virus-intraspecific-strain,virus-tax-id,is-vaccine-strain,is-lab-host,purpose-of-sequencing,submitter-names,submitter-affiliation"

    String data_package_basename = taxid_or_organism_name+"genome_data_package"


    # ncbi datasets is in the viral-core image but could also use:
    #   quay.io/biocontainers/ncbi-datasets-cli:14.26.0
    String docker = "quay.io/biocontainers/ncbi-datasets-cli:14.26.0"
  }

  meta {
        description: "This uses the NCBI datasets CLI tool to download a set of genomes for a given taxonomic ID or organism name."
        volatile: true
  }

  parameter_meta {
    taxid_or_organism_name: {
      description: "The taxonomic ID or taxnomic name for which genomes should be downloaded. This must match a value in NCBI Taxonomy.",
    }
    date_released_after: {
      description: "Only include sequences released after this date ('YYYY-MM-DD' format).",
    }
    ncbi_api_key: {
      description: "NCBI API key; request volume limited to 5 rps without an API key, 10+ with an API key. To obtain a key, see: https://support.nlm.nih.gov/kbArticle/?pn=KA-05317"
    }
    metadata_columns: {
      description: "Columns to include in the metadata output. For valid column names, see the 'Table Field Mnemonic' values in the Virus data package report schema: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/virus/"
    }
    segment_num_or_name: {
        description: "limit to sequences of a specified segment number (ex. '4','6' for Influenza A) or name ('L' for Lassa mammarenavirus); this finds sequences where each has a '/segment' field in the GenBank record and the value of '/segment' matches exactly with the desired segment number or name. Values of '/segment' are entered at sequence submission and are not curated or canonicalized. NB: not all sequences have '/segment' metadata, and the values do not currently conform to a controlled vocabulary. Specifying a segment may unintentionally exclude sequences if they lack metadata or the '/segment' value differs from the requested segment number or name."
    }
    refseq_only: {
      description: "limit to only refseq sequences"
    }
    genotype_or_subtype_name_query: {
        description: "limit to sequences with a specified genotype or subtype name; wildcards allowed. (Ex. 'H*N1')"
    }
    host_taxid_or_organism_name: {
      description: "Limit to virus genomes isolated from a specified host species (organism name, or its numeric taxid)"
    }
    geo_location_country: {
      description: "Limit to virus genomes isolated from the specified country or countries"
    }
    geo_location_state: {
      description: "Two letter abbreviation of the state of the virus specifime collection (if United States)"
    }
  }

  command <<<
    export NCBI_API_KEY="~{ncbi_api_key}"

    

    


    
    function get_tax_id_for_organism_name_or_taxid_input {
        requested_taxid_or_organism_name="$1"

        tax_id_to_download=""
        # check if the input is a numeric taxid or an organism name
        int_validation_re='^[0-9]+$'
        if ! [[ $requested_taxid_or_organism_name =~ $int_validation_re ]] ; then
            echo "Warning: taxid_or_organism_name '"${requested_taxid_or_organism_name}"' does not look like a numeric NCBI taxonomy ID; treating as organism name and attempting to look up a taxid via entrez" >&2; exit 1

            # look up taxid from organism name using esearch, efetch, and xtract
            # currently disabled since NCBI tools are not available in the viral-core docker image
            # taxid="$(esearch -db taxonomy -query '~{taxid_or_organism_name}' | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId)"
            
            # use esearch to query NCBI Taxonomy for the organism name specified
            # then use jq to parse out the result count and ID, 
            # and pass to read as a pipe-delimited string,
            # which in turn stores the values in bash variables
            IFS='|' read esearch_count esearch_taxid < <(curl --silent 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&retmax=1&retmode=json&term="'${requested_taxid_or_organism_name}'"&field=name~{"&api_key=" + ncbi_api_key}' | jq -r '. | [.esearchresult] | map(.count?,.idlist?[0]) | join("|")')
            # check if the esearch results include only one result and that the taxid is defined
            if [ ! -z $esearch_taxid ] && [ $esearch_count -eq 1 ]; then 
                echo "taxid is present as one value"
                tax_id_to_download="$esearch_taxid"
            else 
                echo "Error: could not find a single taxid for organism name '~{taxid_or_organism_name}'" >&2; exit 1
            fi
        else
          # treat as numeric taxid
          tax_id_to_download="$requested_taxid_or_organism_name"
        fi
        echo $tax_id_to_download
    }

    tax_id_to_download=$(get_tax_id_for_organism_name_or_taxid_input "~{taxid_or_organism_name}")

    # use the NCBI Virus API to obtain a list of accessions for the given taxid
    # this is a workaround for the lack of filtering options in NCBI datasets CLI
    #accession_list_download_url='https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?fq={!tag=SeqType_s}SeqType_s:("Nucleotide")&fq=VirusLineageId_ss:(2697049)&fq={!tag=CollectionDate_dr}CollectionDate_dr:[2024-11-01T00:00:00.00Z TO 2024-11-15T23:59:59.00Z]&fq={!tag=LabHost_s} NOT LabHost_s:*&fq={!tag=VacStrain_s} NOT VacStrain_s:*&fq={!tag=EnvSample_s} NOT EnvSample_s:*&fq=HostLineageId_ss:(9606)&fq={!tag=SLen_i}SLen_i:([26000 TO 3000000])&fq={!tag=QualNum_i}QualNum_i:([0 TO 1400])&cmd=download&sort=SourceDB_s desc,CreateDate_dt desc,id asc&dlfmt=acc&fl=AccVer_s'
    #https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?fq=%7B!tag=SeqType_s%7DSeqType_s:(%22Nucleotide%22)&fq=VirusLineageId_ss:(2697049)&fq=%7B!tag=CollectionDate_dr%7DCollectionDate_dr:%5B2024-11-01T00:00:00.00Z%20TO%202024-11-15T23:59:59.00Z%5D&fq=%7B!tag=LabHost_s%7D%20NOT%20LabHost_s:*&fq=%7B!tag=VacStrain_s%7D%20NOT%20VacStrain_s:*&fq=%7B!tag=EnvSample_s%7D%20NOT%20EnvSample_s:*&fq=HostLineageId_ss:(9606)&fq=%7B!tag=SLen_i%7DSLen_i:(%5B26000%20TO%203000000%5D)&fq=%7B!tag=QualNum_i%7DQualNum_i:(%5B0%20TO%201400%5D)&cmd=download&sort=SourceDB_s%20desc,CreateDate_dt%20desc,id%20asc&dlfmt=acc&fl=AccVer_s


    # backup of working before variables added
    # accession_list_download_query=(
    #     'q=*:*'
    #     'fq={!tag=SeqType_s}SeqType_s:("Nucleotide")'
    #     'fq=VirusLineageId_ss:('${tax_id_to_download}')'
    #     'fq={!tag=CollectionDate_dr}CollectionDate_dr:[2024-11-01T00:00:00.00Z TO 2024-11-15T23:59:59.00Z]'
    #     'fq={!tag=LabHost_s} NOT LabHost_s:*'
    #     'fq={!tag=VacStrain_s} NOT VacStrain_s:*'
    #     'fq={!tag=EnvSample_s} NOT EnvSample_s:*'
    #     'fq=HostLineageId_ss:(9606)'
    #     'fq={!tag=SLen_i}SLen_i:([26000 TO 3000000])'
    #     'fq={!tag=QualNum_i}QualNum_i:([0 TO 1400])'
    #     'cmd=download'
    #     'sort=SourceDB_s desc,CreateDate_dt desc,id asc'
    #     'dlfmt=acc'
    #     'fl=AccVer_s'
    # )

    

    # ToDo: add parameters to array iff they are defined WDL inputs

    accession_list_download_query=(
        'q=*:*'
        'fq={!tag=SeqType_s}SeqType_s:("Nucleotide")'
        'cmd=download'
        'sort=SourceDB_s desc,CreateDate_dt desc,id asc'
        'dlfmt=acc'
        'fl=AccVer_s'
    )

    accession_list_download_query+=('fq=VirusLineageId_ss:('${tax_id_to_download}')')
    
    # limit by date only if a lower or upper bound is specified
    if ~{if defined(date_released_before) || defined(date_released_after) then "true" else "false"}; then {
        # use 1900-01-01 as the start date if no start date is specified
        after_date_to_format=~{if defined(date_released_after) then "-d'~{date_released_before}'" else "1900-01-01"}
        date_released_after_formatted="$(date ${after_date_to_format} +%Y-%m-%d)"
        # use the current date as the end date if no end date is specified
        before_date_to_format=~{if defined(date_released_before) then "-d'~{date_released_before}'" else ""}
        date_released_before_formatted="$(date ${before_date_to_format} +%Y-%m-%d)"

        # add the date range filter to the query
        accession_list_download_query+=('fq={!tag=CollectionDate_dr}CollectionDate_dr:['${date_released_after_formatted}' TO '${date_released_before_formatted}']')
    }
    
    # exclude lab passaged, vaccine strains, or environmental specimens, if specified
    if ~{if exclude_lab_passaged then "true" else "false"}; then {
        accession_list_download_query+=('fq={!tag=LabHost_s} NOT LabHost_s:*')
    }
    if ~{if exclude_vaccine_strains then "true" else "false"}; then {
        accession_list_download_query+=('fq={!tag=VacStrain_s} NOT VacStrain_s:*')    
    }
    if ~{if exclude_environmental_specimens then "true" else "false"}; then {
        accession_list_download_query+=('fq={!tag=EnvSample_s} NOT EnvSample_s:*')
    }

    # limit to virus genomes isolated from a specified host species if specified
    if ~{if defined(host_taxid_or_organism_name) then "true" else "false"}; then {
        tax_id_of_host=$(get_tax_id_for_organism_name_or_taxid_input "~{host_taxid_or_organism_name}")
        accession_list_download_query+=('fq=HostLineageId_ss:('tax_id_of_host')')
    }

    # limit to sequences above a minimum length, if specified
    if ~{if defined(min_sequence_length_end_to_end) then "true" else "false"}; then {
        accession_list_download_query+=('fq={!tag=SLen_i}SLen_i:([~{min_sequence_length_end_to_end} TO 3000000])')
    }


    # limit to sequences below a maximum number of ambiguous bases, if specified
    if ~{if defined(max_num_ambig_bases) then "true" else "false"}; then {
        accession_list_download_query+=('fq={!tag=QualNum_i}QualNum_i:([0 TO ~{max_num_ambig_bases}])')
    }
    
    # limit to sequences of a specified pangolin lineage, if specified
    if ~{if defined(pango_lineage) then "true" else "false"}; then {
        accession_list_download_query+=('fq={!tag=Lineage_s}Lineage_s:("~{pango_lineage}")')
    }

    # limit to refseq sequences only, if specified
    if ~{if refseq_only then "true" else "false"}; then {
        accession_list_download_query+=('fq={!tag=SourceDB_s}SourceDB_s:("RefSeq")')
    }

    # limit to sequences with a specified genotype or subtype name, if specified
    if ~{if defined(genotype_or_subtype_name_query) then "true" else "false"}; then {
        accession_list_download_query+=('fq={!edismax qf=Serotype_s}~{genotype_or_subtype_name_query}')
    }

    # limit to sequences of a specified segment number or name, if specified
    # this finds sequences where each has a '/segment' field in the GenBank record
    # and the value of '/segment' matches exactly with the desired segment number or name.
    # it falls back to 'genome' if no segment is specified
    if ~{if defined(segment_num_or_name) then "true" else "false"}; then {
        accession_list_download_query+=('fq={!tag=Segment_s}Segment_s:("~{segment_num_or_name}")')
    } else {
        accession_list_download_query+=('fq={!tag=Segment_s}Segment_s:("genome")')
    }

    # template
    # if ~{if defined() then "true" else "false"}; then {
    #     accession_list_download_query+=('fq=')
    # }

    # join the url parameters above with the particular url encoding expected by NCBI VirusVariation 'API'
    IFS='' accession_list_download_query_str="$(printf '%s&' ${accession_list_download_query[@]} | jq -sRr @uri | sed 's/%26/\&/g' | sed 's/%3D/=/g' | sed 's/%21/!/g' | sed 's/%3A/:/g' | sed 's/%2A/\*/g' | sed 's/%28/\(/g' | sed 's/%29/\)/g' | sed 's/%2C/,/g')"
    accession_list_download_url='https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?'"${accession_list_download_query_str}"
    echo ${accession_list_download_url}

    curl --silent --retry 3 -L "${accession_list_download_url}" > acc_list.txt

    # check if the accession list download failed or is empty
    if [ ! -s acc_list.txt ]; then
        echo "Error: accession list download failed or is empty" >&2; exit 1
    fi


    #WORKING curl --silent --http1.1 --retry 5 -L 'https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?q=*:*&fq=%7B!tag=SeqType_s%7DSeqType_s:(%22Nucleotide%22)&fq=VirusLineageId_ss:(2697049)&fq=%7B!tag=CollectionDate_dr%7DCollectionDate_dr:%5B2024-11-01T00:00:00.00Z%20TO%202024-11-15T23:59:59.00Z%5D&fq=%7B!tag=LabHost_s%7D%20NOT%20LabHost_s:*&fq=%7B!tag=VacStrain_s%7D%20NOT%20VacStrain_s:*&fq=%7B!tag=EnvSample_s%7D%20NOT%20EnvSample_s:*&fq=HostLineageId_ss:(9606)&fq=%7B!tag=SLen_i%7DSLen_i:(%5B26000%20TO%203000000%5D)&fq=%7B!tag=QualNum_i%7DQualNum_i:(%5B0%20TO%201400%5D)&cmd=download&sort=SourceDB_s%20desc,CreateDate_dt%20desc,id%20asc&dlfmt=acc&fl=AccVer_s' > acc_list.txt
    
    # set the number of workers for downloading data to nproc (-1 core for OS overhead), 
    # with an upper bound of 30 cores as specified by NCBI datasets
    NUM_WORKERS=$(nproc --ignore 1)
    NUM_WORKERS=$((NUM_WORKERS<31 ? NUM_WORKERS : 30))

    # until NCBI datasets offers filters for collection date and segment, we can obtain a list of accessions
    # from the NCBI Virus API and then pass that to NCBI datasets for download via:
    #   datasets download virus genome accession --inputfile acc_list.txt
    # (NCBI datasets CLI filters still apply)
    datasets download virus genome accession --inputfile acc_list.txt \
        --no-progressbar \
        --dehydrated \
        --include genome,biosample \
        --max-workers $NUM_WORKERS \
        ~{"--geo-location" + geo_location_country} \
        ~{"--usa-state" + geo_location_state} \
        --filename ~{data_package_basename}.zip


    # format yyyy-mm-dd to mm/dd/yyyy $(date -d"~{date_released_after}" +%m/%d/%Y)
    # Note: the datasets documentation indictes support for ISO8601 
    #       dates however this does not yet appear to be the case
    #date_released_after_formatted="$(date -d'~{date_released_after}' +%m/%d/%Y)"
    #date_released_before_formatted="$(date -d'~{date_released_after}' +%m/%d/%Y)"

    # download a "dehydrated" data package with metadata but no sequence data
    # datasets download virus genome taxon ~{taxid_or_organism_name} \
    #   --no-progressbar \
    #   --dehydrated \
    #   ~{if defined(date_released_after) then "--released-after" else ""} ${date_released_after_formatted} \
    #   ~{"--refseq" + refseq_only} \
    #   ~{"--host" + host_species_name} \
    #   ~{"--geo-location" + geo_location_country} \
    #   ~{"--usa-state" + geo_location_state} \
    #   --filename ~{data_package_basename}.zip \
    #   --include genome,biosample
      
    # can also use (after obtaining accessions from NCBI Virus API w/ segment filtering):
    #   datasets download virus genome accession --inputfile accessions.txt
    # accession list download example:
    #   https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?fq=%7B!tag=SeqType_s%7DSeqType_s:(%22Nucleotide%22)&fq=VirusLineageId_ss:(102793)&fq=%7B!tag=Segment_s%7DSegment_s:(%226%22)&cmd=download&sort=SourceDB_s%20desc,CreateDate_dt%20desc,id%20asc&dlfmt=acc&fl=AccVer_s

    # decompress the data package
    unzip ~{data_package_basename}.zip -d ~{data_package_basename}

    # download the sequence data described within the data package
    datasets rehydrate \
      --gzip \
      --max-workers $NUM_WORKERS \
      --directory "~{data_package_basename}/"

    # rename the fasta file containing the sequence data
    mv "~{data_package_basename}/ncbi_dataset/data/genomic.fna" \
       "~{data_package_basename}/ncbi_dataset/data/~{data_package_basename}.fasta"

    # output sequence metadata in tsv format, with the desired columns specified via ~{metadata_columns}
    dataformat tsv virus-genome \
      --inputfile "~{data_package_basename}/ncbi_dataset/data/data_report.jsonl" \
      --fields "~{metadata_columns}" > "~{data_package_basename}_metadata.tsv"

    #~{data_package_basename}/ncbi_dataset/data/data_report.jsonl
  >>>

  output {
    File viral_genomes_fasta                = "~{data_package_basename}/ncbi_dataset/data/~{data_package_basename}.fasta"
    File viral_genomes_metadata_report_json = "~{data_package_basename}/ncbi_dataset/data/data_report.jsonl"
    File viral_genomes_metadata_tsv         = "~{data_package_basename}_metadata.tsv"
  }

  runtime {
    docker: docker
    memory: "7 GB"
    cpu: 2
    dx_instance_type: "mem2_ssd1_v2_x2"
    maxRetries: 2
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

    String       docker = "quay.io/broadinstitute/viral-phylo:2.3.6.0"
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

    String       docker = "quay.io/broadinstitute/viral-phylo:2.3.6.0"
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
    maxRetries: 2
  }
}

task structured_comments {
  input {
    File   assembly_stats_tsv

    File?  filter_to_ids

    String docker = "quay.io/broadinstitute/viral-core:2.4.0"
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

    String docker = "quay.io/broadinstitute/viral-core:2.4.0"
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
    String      docker="quay.io/broadinstitute/viral-core:2.4.0"
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
    Array[File] cleaned_bam_filepaths
    File        demux_meta_json

    String  sample_table_name  = "sample"
    String  docker = "python:slim"
  }
  String  sanitized_id_col = "entity:~{sample_table_name}_id"
  String base = basename(basename(biosample_attributes_tsv, ".txt"), ".tsv")
  parameter_meta {
    cleaned_bam_filepaths: {
      description: "Unaligned bam files containing cleaned (submittable) reads.",
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
    bam_fnames = list(os.path.basename(x) for x in '~{sep="*" cleaned_bam_filepaths}'.split('*'))
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
          biosample_attributes.append(row)
          for k,v in row.items():
            if v.strip().lower() in ('missing', 'na', 'not applicable', 'not collected', ''):
              v = None
            if v and (k not in biosample_headers) and k not in ('message', 'accession'):
              biosample_headers.append(k)
    print("biosample headers ({}): {}".format(len(biosample_headers), biosample_headers))
    print("biosample output rows ({})".format(len(biosample_attributes)))
    samples_seen_without_biosample = set(sample_names_seen) - set(row['sample_name'] for row in biosample_attributes)
    print("samples seen in bams without biosample entries ({}): {}".format(len(samples_seen_without_biosample), sorted(samples_seen_without_biosample)))

    # write reformatted table
    with open('~{base}.entities.tsv', 'w', newline='') as outf:
      writer = csv.DictWriter(outf, delimiter='\t', fieldnames=["~{sanitized_id_col}"]+biosample_headers, dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
      writer.writeheader()
      for row in biosample_attributes:
        outrow = {h: row[h] for h in biosample_headers}
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
    Int     num_segments = 1
    Int     taxid

    File?   filter_to_ids

    Boolean s_dropout_note = true
    String  docker = "quay.io/broadinstitute/viral-phylo:2.3.6.0"
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
    String? out_base = "authors"

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
    String       docker = "quay.io/broadinstitute/viral-phylo:2.3.6.0"
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
    maxRetries: 2
  }
}

task package_sc2_genbank_ftp_submission {
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

    String  docker = "quay.io/broadinstitute/viral-baseimage:0.2.0"
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
    maxRetries: 2
  }
}

task vadr {
  meta {
    description: "Runs NCBI's Viral Annotation DefineR for annotation and QC. Defaults here are for SARS-CoV-2 (see https://github.com/ncbi/vadr/wiki/Coronavirus-annotation), but VADR itself is applicable to a larger number of viral taxa (change the vadr_opts accordingly)."
  }
  input {
    File   genome_fasta
    String vadr_opts = "--glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn"

    String docker = "quay.io/staphb/vadr:1.6.3"
    Int    minlen = 50
    Int    maxlen = 30000
    Int    mem_size = 4
    Int    cpus = 2
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
    cat "~{out_base}/~{out_base}.vadr.alt.list" | cut -f 5 | tail -n +2 \
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
    memory: mem_size + " GB"
    cpu: cpus
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task sequence_rename_by_species {
  meta {
    description: "Rename sequences based on species-specific naming conventions for many viral taxa."
  }
  input {
    String sample_id
    String organism_name
    File   biosample_attributes
    String taxid
    File   taxdump_tgz

    String docker = "quay.io/broadinstitute/viral-classify:2.2.5"
  }
  command <<<
    set -e
    mkdir -p taxdump
    read_utils.py extract_tarball "~{taxdump_tgz}" taxdump
    python3 << CODE
    import metagenomics
    taxdb = metagenomics.TaxonomyDb(tax_dir='taxdump', load_nodes=True, load_gis=False)
    taxid = int('~{taxid}')
    ancestors = taxdb.get_ordered_ancestors(taxid)


    if any(node == 3052310 for node in [taxid] + ancestors):
      # LASV
      pass
    elif any(node == 186538 for node in [taxid] + ancestors):
      # ZEBOV
      pass
    elif any(node == 11250 for node in [taxid] + ancestors):
      # RSV -- no real convention! Some coalescence around this:
      # <type>/<host lowercase>/Country/ST-Institution-LabID/Year
      # e.g. RSV-A/human/USA/MA-Broad-1234/2020
      pass
    elif any(node == 2697049 for node in [taxid] + ancestors):
      # SARS-CoV-2
      # SARS-CoV-2/<host lowercase>/Country/ST-Institution-LabID/Year
      # e.g. SARS-CoV-2/human/USA/MA-Broad-1234/2020
      pass
    elif any((node == 11320 or node == 11520) for node in [taxid] + ancestors):
      # Flu A or B
      # <type>/<hostname if not human>/<geoloc>/seqUID/year
      # e.g. A/Massachusetts/Broad_MGH-1234/2001 or A/chicken/Hokkaido/TU25-3/2022 or B/Rhode Island/RISHL-1234/2024
      pass
    elif any(node == 12059 for node in [taxid] + ancestors):
      # Enterovirus (including rhinos)
      pass
    else:
      # everything else
      pass

    CODE
  >>>
  output {
    String assembly_name_genbank = read_string("assembly_name_genbank")
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}
