version 1.0

task md5sum {
  input {
    File in_file
  }
  command {
    md5sum ${in_file} | cut -f 1 | tee MD5
  }
  output {
    String md5 = read_string("MD5")
  }
  runtime {
    docker: "ubuntu"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 100 HDD"
    dx_instance_type: "mem1_ssd2_v2_x2"
  }
}

task download_fasta {
  input {
    String         out_prefix
    Array[String]+ accessions
    String         emailAddress

    String         docker="quay.io/broadinstitute/viral-phylo:2.1.12.0"
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
    docker: "${docker}"
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

    String         docker="quay.io/broadinstitute/viral-phylo:2.1.12.0"
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
    docker: "${docker}"
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

    String  docker="quay.io/broadinstitute/viral-phylo:2.1.12.0"
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
    docker: "${docker}"
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

    String  docker="quay.io/broadinstitute/viral-phylo:2.1.12.0"
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
    docker: "${docker}"
    memory: "15 GB"
    cpu: 4
    dx_instance_type: "mem2_ssd1_v2_x4"
    preemptible: 1
  }
}

task biosample_to_genbank {
  meta {
    description: "Prepares two input metadata files for Genbank submission based on a BioSample registration attributes table (attributes.tsv) since all of the necessary values are there. This produces both a Genbank Source Modifier Table and a BioSample ID map file that can be fed into the prepare_genbank task."
  }
  input {
    File  biosample_attributes
    Int   num_segments=1
    Int   taxid

    String  docker="quay.io/broadinstitute/viral-phylo:2.1.12.0"
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
        --loglevel DEBUG
  }
  output {
    File genbank_source_modifier_table = "${base}.genbank.src"
    File biosample_map                 = "${base}.biosample.map.txt"
  }
  runtime {
    docker: "${docker}"
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
    String       docker="quay.io/broadinstitute/viral-phylo:2.1.12.0"
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
