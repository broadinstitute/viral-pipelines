version 1.1

#DX_SKIP_WORKFLOW

import "../tasks/tasks_terra.wdl" as terra

workflow terra_tsv_to_table {
    meta {
        description: "Merge per-entity metadata tsv files and upload to a Terra data table: insert-or-update on existing rows/columns. Inputs are reconciled by column NAME, not column position, so tsv files carrying the same columns in different orders merge correctly, and the entity id column is forced to column 1 where Terra expects to find it."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File?]+ tsv_files
        Array[String] preferred_col_order = []
        String?       entity_table_name
    }

    parameter_meta {
        tsv_files: {
            description: "Terra entity tsv files to merge and upload. Each must have a header row, and all must describe the same table. The entity id column may be either the 'entity:<table>_id' form Terra expects on upload or the prefix-less '<table>_id' form Terra emits on download -- if NO input carries the prefixed form, supply entity_table_name so the table can be identified. Column order may differ between files, as may the subset of columns present. Nulls are dropped."
        }
        preferred_col_order: {
            description: "Optional canonical column order for the merged tsv, e.g. the assembly_header literal from assemble_denovo_metagenomic.wdl. Only reorders columns, never creates them. Cosmetic: Terra matches columns by name, so only column 1 affects the import."
        }
        entity_table_name: {
            description: "Terra table name (e.g. 'assembly'), needed only when no input file carries an 'entity:<table>_id' column -- which is the case when every input came from terra_table_to_tsv / download_entities_tsv, since Terra strips the prefix on the way out."
        }
    }

    call terra.check_terra_env

    call terra.merge_entities_tsvs {
        input:
            input_tsvs          = select_all(tsv_files),
            preferred_col_order = preferred_col_order,
            entity_table_name   = entity_table_name,
            out_basename        = "terra_upload"
    }

    call terra.upload_entities_tsv {
        input:
            tsv_file         = merge_entities_tsvs.out_tsv,
            workspace_name   = check_terra_env.workspace_name,
            terra_project    = check_terra_env.workspace_namespace
    }

    output {
        Array[String] upload_response = upload_entities_tsv.stdout

        File   merged_tsv           = merge_entities_tsvs.out_tsv
        File   merge_log            = merge_entities_tsvs.merge_log
        String entity_table         = merge_entities_tsvs.entity_table
        String entity_id_col        = merge_entities_tsvs.entity_id_col
        Int    num_rows_uploaded    = merge_entities_tsvs.num_rows
        Int    num_cols_uploaded    = merge_entities_tsvs.num_cols
        Int    num_input_col_orders = merge_entities_tsvs.num_input_col_orders
    }
}
