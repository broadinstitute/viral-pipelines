version 1.0

task nextclade_one_sample {
    meta {
        description: "Nextclade classification of one sample. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
    }
    input {
        File  genome_fasta
        File? root_sequence
        File? auspice_reference_tree_json
        File? qc_config_json
        File? gene_annotations_json
        File? pcr_primers_csv
        String docker = "neherlab/nextclade:0.14.4"
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
        docker: docker
        memory: "3 GB"
        cpu:    2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        String nextclade_version = read_string("VERSION")
        File   nextclade_json    = "~{basename}.nextclade.json"
        File   auspice_json      = "~{basename}.nextclade.auspice.json"
        File   nextclade_tsv     = "~{basename}.nextclade.tsv"
        String nextclade_clade   = read_string("NEXTCLADE_CLADE")
        String aa_subs_csv       = read_string("NEXTCLADE_AASUBS")
        String aa_dels_csv       = read_string("NEXTCLADE_AADELS")
    }
}

task nextclade_many_samples {
    meta {
        description: "Nextclade classification of many samples. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
    }
    input {
        Array[File]+ genome_fastas
        File?        root_sequence
        File?        auspice_reference_tree_json
        File?        qc_config_json
        File?        gene_annotations_json
        File?        pcr_primers_csv
        String       basename
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
        docker: "nextstrain/nextclade:0.14.4"
        memory: "14 GB"
        cpu:    16
        disks: "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x16"
    }
    output {
        String nextclade_version = read_string("VERSION")
        File   nextclade_json    = "~{basename}.nextclade.json"
        File   auspice_json      = "~{basename}.nextclade.auspice.json"
        File   nextclade_tsv     = "~{basename}.nextclade.tsv"
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
        Boolean inference_usher=true
        String  docker = "staphb/pangolin:3.0.5-pangolearn-2021-06-05"
    }
    String basename = basename(genome_fasta, ".fasta")
    command <<<
        date | tee DATE
        conda list -n pangolin | grep "usher" | awk -F ' +' '{print$1, $2}' | tee VERSION_PANGO_USHER
        set -ex
        pangolin -v | tee VERSION_PANGOLIN
        pangolin -pv | tee VERSION_PANGOLEARN

        pangolin "~{genome_fasta}" \
            --outfile "~{basename}.pangolin_report.csv" \
            ~{"--min-length " + min_length} \
            ~{"--max-ambig " + max_ambig} \
            ~{true='--usher' false='' inference_usher} \
            --alignment \
            --verbose

        cp sequences.aln.fasta "~{basename}.pangolin_msa.fasta"
        python3 <<CODE
        import csv
        #grab output values by column header
        with open("~{basename}.pangolin_report.csv", 'rt') as csv_file:
            for line in csv.DictReader(csv_file):
                with open("VERSION", 'wt') as outf:
                    pangolin_version=line["pangolin_version"]
                    version=line["version"]
                    outf.write(f"pangolin {pangolin_version}; {version}")
                with open("PANGO_LINEAGE", 'wt') as outf:
                    outf.write(line["lineage"])
                with open("PANGOLIN_CONFLICTS", 'wt') as outf:
                    outf.write(line["conflict"])
                with open("PANGOLIN_NOTES", 'wt') as outf:
                    outf.write(line["note"])
                break
        CODE
        grep ^lineage transposed.tsv | cut -f 2 | grep -v lineage > PANGOLIN_CLADE
    >>>
    runtime {
        docker: docker
        memory: "3 GB"
        cpu:    2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        String     date                   = read_string("DATE")
        String     version                = read_string("VERSION")
        String     pango_lineage          = read_string("PANGO_LINEAGE")
        String     pangolin_conflicts     = read_string("PANGOLIN_CONFLICTS")
        String     pangolin_notes         = read_string("PANGOLIN_NOTES")
        String     pangolin_usher_version = read_string("VERSION_PANGO_USHER")
        String     pangolin_version       = read_string("VERSION_PANGOLIN")
        String     pangolearn_version     = read_string("VERSION_PANGOLEARN")
        String     pangolin_docker        = docker
        File       pango_lineage_report   = "${basename}.pangolin_report.csv"
        File       msa_fasta              = "~{basename}.pangolin_msa.fasta"
    }
}


task sequencing_report {
    meta {
        description: "Produce sequencing progress report."
    }
    input {
        File    assembly_stats_tsv
        File?   collab_ids_tsv

        String? sequencing_lab = "Broad Institute"
        String? intro_blurb = "The Broad Institute Viral Genomics group, in partnership with the Genomics Platform and Data Sciences Platform, has been engaged in viral sequencing of COVID-19 patients since March 2020."
        String? max_date
        String? min_date
        Int?    min_unambig
        String? voc_list
        String? voi_list

        String  docker = "quay.io/broadinstitute/sc2-rmd:0.1.18"
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
        zip all_reports.zip *.pdf *.xlsx *.tsv
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
        Array[File] sheets  = glob("*.xlsx")
        File        all_zip = "all_reports.zip"
        File        all_tsv = glob("report-*-everything-*.tsv")[0]
    }
}


task sc2_meta_final {
    meta {
        description: "Final wash and cleanup of assembly metadata prior to delivery."
    }
    input {
        File          assembly_stats_tsv
        File?         collab_ids_tsv
        File?         genome_status_json

        String        collab_ids_idcol = 'external_id'
        Array[String] collab_ids_addcols = ['collaborator_id']

        String?       max_date
        String?       min_date
        Int?          min_unambig=24000
        Boolean       drop_file_cols=false

        File?         filter_to_ids

        String        docker = "quay.io/broadinstitute/py3-bio:0.1.2"
    }
    String out_basename = basename(basename(assembly_stats_tsv, '.txt'), '.tsv')
    command <<<
        set -e
        python3<<CODE
        import json
        import os
        import datetime
        import epiweeks
        import pandas as pd
        import numpy as np

        # set inputs
        collab_idcol = '~{collab_ids_idcol}'
        collab_addcols = list(x for x in '~{sep="*" collab_ids_addcols}'.split('*') if x)
        assemblies_tsv = "~{assembly_stats_tsv}"
        collab_tsv = ~{default='None' '"' + collab_ids_tsv + '"'}
        min_unambig = ~{default='0' min_unambig}
        min_date = ~{default='None' '"' + min_date + '"'}
        max_date = ~{default='None' '"' + max_date + '"'}
        drop_file_cols = ~{true='True' false='False' drop_file_cols}
        filter_to_ids = ~{default='None' '"' + filter_to_ids + '"'}
        genome_status_json = '~{default="" genome_status_json}'
        if genome_status_json:
          with open(genome_status_json, 'rt') as inf:
            genome_status = json.load(inf)
        else:
          genome_status = {}

        # read input files
        df_assemblies = pd.read_csv(assemblies_tsv, sep='\t').dropna(how='all')
        if collab_tsv and os.path.isfile(collab_tsv) and os.path.getsize(collab_tsv):
            collab_ids = pd.read_csv(collab_tsv, sep='\t').dropna(how='all').dropna(how='all', axis='columns')
            if 'collection_date' in collab_ids.columns:
                collab_ids.drop(columns=['collection_date'], inplace=True)
            if collab_addcols:
                collab_ids = collab_ids[[collab_idcol] + collab_addcols]
            if collab_idcol != 'sample':
                collab_ids = collab_ids.rename(columns={collab_idcol: 'sample'})
        else:
            collab_ids = pd.DataFrame(columns = ['sample'] + collab_addcols)
        if filter_to_ids:
            with open(filter_to_ids, 'rt') as inf:
                keep_list = set(x.strip() for x in inf)

        # remove columns with File URIs if requested
        if drop_file_cols:
            cols_unwanted = [
                'assembly_fasta','coverage_plot','aligned_bam','replicate_discordant_vcf',
                'variants_from_ref_vcf','nextclade_tsv','nextclade_json',
                'pangolin_csv','vadr_tgz','vadr_alerts',
            ]
            cols_unwanted = list(c for c in cols_unwanted if c in df_assemblies.columns)
            df_assemblies.drop(columns=cols_unwanted, inplace=True)

        # filter to IDs if requested
        if filter_to_ids:
            df_assemblies = df_assemblies[df_assemblies['sample'].isin(keep_list)]

        # format dates properly and remove all rows with missing or bad dates
        df_assemblies = df_assemblies.loc[~df_assemblies['run_date'].isna()]
        df_assemblies.loc[:,'run_date'].replace('missing', pd.NA, inplace=True)
        df_assemblies.loc[:,'collection_date'].replace('missing', pd.NA, inplace=True)
        df_assemblies = df_assemblies.astype({'collection_date':'datetime64[D]','run_date':'datetime64[D]'})

        # fix vadr_num_alerts
        df_assemblies = df_assemblies.astype({'vadr_num_alerts':'Int64'})

        # subset by date range
        if min_date:
            df_assemblies = df_assemblies.loc[np.datetime64(min_date) <= df_assemblies['run_date']]
        if max_date:
            df_assemblies = df_assemblies.loc[df_assemblies['run_date'] <= np.datetime64(max_date)]

        # fix missing data in purpose_of_sequencing
        df_assemblies.loc[:,'purpose_of_sequencing'] = df_assemblies.loc[:,'purpose_of_sequencing'].fillna('Missing').replace('', 'Missing')

        # derived column: genome_status
        df_assemblies.loc[:,'genome_status'] = list(
                'failed_sequencing' if df_assemblies.loc[id, 'assembly_length_unambiguous'] < min_unambig
                else genome_status[df_assemblies.loc[id, 'sample']] if df_assemblies.loc[id, 'sample'] in genome_status
                else 'failed_annotation' if df_assemblies.loc[id, 'vadr_num_alerts'] > 0
                else 'submittable'
                for id in df_assemblies.index)

        # derived columns: geo_country, geo_state, geo_locality, geo_state_abbr
        df_assemblies.loc[:,'geo_country'] = list(g.split(': ')[0] if not pd.isna(g) else '' for g in df_assemblies.loc[:,'geo_loc_name'])
        df_assemblies.loc[:,'geo_state'] = list(g.split(': ')[1].split(', ')[0] if not pd.isna(g) else '' for g in df_assemblies.loc[:,'geo_loc_name'])
        df_assemblies.loc[:,'geo_locality'] = list(g.split(': ')[1].split(', ')[1] if not pd.isna(g) and ', ' in g else '' for g in df_assemblies.loc[:,'geo_loc_name'])
        df_assemblies.loc[:,'geo_state_abbr'] = list(s.split('/')[1].split('-')[0] if '/' in s and '-' in s.split('/')[1] else '' for s in df_assemblies.loc[:,'sample'])

        # derived columns: collection_epiweek, run_epiweek
        df_assemblies.loc[:,'collection_epiweek'] = list(epiweeks.Week.fromdate(x) if not pd.isna(x) else x for x in df_assemblies.loc[:,'collection_date'])
        df_assemblies.loc[:,'run_epiweek'] = list(epiweeks.Week.fromdate(x) if not pd.isna(x) else x for x in df_assemblies.loc[:,'run_date'])
        df_assemblies.loc[:,'collection_epiweek_end'] = list(x.enddate().strftime('%Y-%m-%d') if not pd.isna(x) else '' for x in df_assemblies.loc[:,'collection_epiweek'])
        df_assemblies.loc[:,'run_epiweek_end'] = list(x.enddate().strftime('%Y-%m-%d') if not pd.isna(x) else '' for x in df_assemblies.loc[:,'run_epiweek'])

        # derived column: sample_age_at_runtime
        df_assemblies.loc[:,'sample_age_at_runtime'] = list(
            (x.days if not pd.isna(x) else pd.NA)
            for x in (df_assemblies.loc[:,'run_date'] - df_assemblies.loc[:,'collection_date']))

        # join column: collaborator_id
        df_assemblies = df_assemblies.merge(collab_ids, on='sample', how='left', validate='one_to_one')

        # write final output
        df_assemblies.to_csv("~{out_basename}.final.tsv", sep='\t', index=False)
        CODE
    >>>
    runtime {
        docker: docker
        memory: "2 GB"
        cpu:    2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File meta_tsv = "~{out_basename}.final.tsv"
    }
}


task crsp_meta_etl {
    meta {
        description: "Initial cleanup of CRSP sample metadata."
    }
    input {
        File          sample_meta_crsp
        String        salt
        String        bioproject

        String        country = 'USA'
        String        collected_by = 'Broad Institute Clinical Research Sequencing Platform'
        String        ontology_map_states = '{"AL": "Alabama", "AK": "Alaska", "AZ": "Arizona", "AR": "Arkansas", "CA": "California", "CO": "Colorado", "CT": "Connecticut", "DE": "Delaware", "DC": "District of Columbia", "FL": "Florida", "GA": "Georgia", "HI": "Hawaii", "ID": "Idaho", "IL": "Illinois", "IN": "Indiana", "IA": "Iowa", "KS": "Kansas", "KY": "Kentucky", "LA": "Louisiana", "ME": "Maine", "MD": "Maryland", "MA": "Massachusetts", "MI": "Michigan", "MN": "Minnesota", "MS": "Mississippi", "MO": "Missouri", "MT": "Montana", "NE": "Nebraska", "NV": "Nevada", "NH": "New Hampshire", "NJ": "New Jersey", "NM": "New Mexico", "NY": "New York", "NC": "North Carolina", "ND": "North Dakota", "OH": "Ohio", "OK": "Oklahoma", "OR": "Oregon", "PA": "Pennsylvania", "RI": "Rhode Island", "SC": "South Carolina", "SD": "South Dakota", "TN": "Tennessee", "TX": "Texas", "UT": "Utah", "VT": "Vermont", "VA": "Virginia", "WA": "Washington", "WV": "West Virginia", "WI": "Wisconsin", "WY": "Wyoming"}'
        String        ontology_map_body_part = '{"AN SWAB": "Anterior Nares", "AN Swab": "Anterior Nares", "Swab": "Upper respiratory tract", "Viral": "Upper respiratory tract", "Null": "Anterior Nares", "NP Swab": "Nasopharynx (NP)"}'

        String        docker = "quay.io/broadinstitute/py3-bio:0.1.2"
    }
    String out_basename = basename(basename(basename(sample_meta_crsp, '.txt'), '.tsv'), '.metadata')
    command <<<
        set -e
        python3<<CODE
        import base64
        import datetime
        import hashlib
        import json
        import pandas as pd

        # load some inputs
        salt = '~{salt}'.strip()
        ontology_map_states = json.loads('~{ontology_map_states}')
        ontology_map_body_part = json.loads('~{ontology_map_body_part}')

        # read input files
        sample_meta = pd.read_csv('~{sample_meta_crsp}', sep='\t')

        # clean collection_date
        sample_meta = sample_meta.astype({'collection_date':'datetime64[D]'})
        sample_meta.loc[:,'collection_year'] = list(d.year for d in sample_meta.loc[:,'collection_date'])

        # validation checks
        assert sample_meta.geo_loc_name.isna().sum() == 0, "error: some samples missing geo_loc_name"
        assert sample_meta.collection_date.isna().sum() == 0, "error: some samples missing collection_date"
        assert sample_meta.collected_by.isna().sum() == 0, "error: some samples missing collected_by"
        assert all(sample_meta.collected_by == '~{collected_by}'), "error: not all samples collected by same lab"
        assert sample_meta.anatomical_part.isna().sum() == 0, "error: some samples missing anatomical_part"
        assert sample_meta.hl7_message_id.isna().sum() == 0, "error: some samples missing hl7_message_id"
        assert sample_meta.internal_id.isna().sum() == 0, "error: some samples missing internal_id"

        # clean geoloc
        sample_meta.loc[:,'geo_loc_state_abbr'] = sample_meta.loc[:,'geo_loc_name']
        sample_meta = sample_meta.assign(geo_loc_country = '~{country}')
        sample_meta.loc[:,'geo_loc_name'] = ["~{country}: " + ontology_map_states[x] for x in sample_meta.loc[:,'geo_loc_state_abbr']]

        # translate specimen type
        sample_meta.loc[:,'anatomical_part'] = [ontology_map_body_part[x] for x in sample_meta.loc[:,'anatomical_part']]
        sample_meta = sample_meta.assign(body_product = 'Mucus')

        # construct IDs
        ''' The hl7_message_id is a 10-byte b32encoded unique sample
            identifier provided by CRSP to link and report diagnostic
            results to Public Health and also to submitting labs and
            patients. These IDs may be known to patients and providers.
            We want to report genomic results to Public Health, but not
            to patients or providers. We will use hl7_message_id to
            derive a new hashed ID that we will use for public data release.
            1: No member of the public should be able to infer the
            hl7_message_id from the public ID.
            2: The patient must not be able to derive the public ID from
            information available to them (e.g. hl7_message_id).
            3: The hashes should not collide across the range of
            hl7_message_id inputs.
        '''
        sample_meta.loc[:,'hl7_hashed'] = [
            base64.b32encode(hashlib.pbkdf2_hmac('sha256', id.encode('UTF-8'), salt.encode('UTF-8'), 20000, dklen=10)).decode('UTF-8')
            for id in sample_meta.loc[:,'hl7_message_id']
        ]
        sample_meta['host_subject_id'] = [
            f'CDCBI-CRSP_{id}' for id
            in sample_meta['hl7_hashed']]
        sample_meta['sample_name'] = [
            f'{country}/{state}-{id}/{year}'
            for country, state, id, year
            in zip(sample_meta['geo_loc_country'], sample_meta['geo_loc_state_abbr'], sample_meta['host_subject_id'], sample_meta['collection_year'])]
        sample_meta['isolate'] = [f'SARS-CoV-2/Human/{id}'
            for id in sample_meta['sample_name']]

        # prep biosample submission table
        biosample = sample_meta[['sample_name', 'isolate', 'collected_by', 'collection_date', 'geo_loc_name', 'host_subject_id', 'anatomical_part', 'body_product']]
        biosample = biosample.assign(
            bioproject_accession = '~{bioproject}',
            attribute_package = 'Pathogen.cl',
            organism = 'Severe acute respiratory syndrome coronavirus 2',
            isolation_source = 'Clinical',
            lat_lon = 'missing',
            host = 'Homo sapiens',
            host_disease = 'COVID-19',
            purpose_of_sampling = 'Diagnostic Testing',
            purpose_of_sequencing = 'Baseline surveillance (random sampling)'
        )
        biosample = biosample.reindex(columns= biosample.columns.tolist() + [
            'host_health_state',' host_disease_outcome', 'host_age', 'host_sex',
            'anatomical_material', 'collection_device', 'collection_method',
            'passage_history', 'lab_host', 'passage_method',
            'culture_collection', 'host_speciman_voucher',
            'environmental_material', 'environmental_site',
            'description'
        ])

        # create ID map
        collab_ids = sample_meta[['sample_name','hl7_message_id','internal_id']]
        collab_ids = collab_ids.rename(columns={'sample_name': 'external_id'})

        # write final output
        biosample.to_csv("biosample_meta_submit-~{out_basename}.tsv", sep='\t', index=False)
        collab_ids.to_csv("collab_ids-~{out_basename}.tsv", sep='\t', index=False)
        CODE
    >>>
    runtime {
        docker: docker
        memory: "2 GB"
        cpu:    2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File          biosample_submit_tsv = "biosample_meta_submit-~{out_basename}.tsv"
        File          collab_ids_tsv       = "collab_ids-~{out_basename}.tsv"

        String        collab_ids_idcols    = 'external_id'
        Array[String] collab_ids_addcols   = ['hl7_message_id','internal_id']
    }
}

task gisaid_uploader {
  input {
    File    gisaid_sequences_fasta
    File    gisaid_meta_csv
    File    cli_auth_token
  }
  command {
    set -e
    cp "~{cli_auth_token}" gisaid_uploader.authtoken
    gisaid_uploader CoV upload \
        --fasta "~{gisaid_sequences_fasta}" \
        --csv "~{gisaid_meta_csv}" \
        --failedout failed_metadata.csv
  }
  output {
    File  failed_metadata = "failed_metadata.csv"
  }
  runtime {
    docker: "quay.io/broadinstitute/gisaid-cli:1.0"
    memory: "2 GB"
    cpu: 2
    disks: "local-disk 100 HDD"
  }
}