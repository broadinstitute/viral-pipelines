version 1.0

workflow CreateEntericsQCViz {

    parameter_meta {
        sample_ids: {description: "selected rows of data from data table which will be used for plotting"}
        input_table_name: {description: "name of the Terra data table - root entity type - from where input data is selected"}
        workspace_name: {description: "name of Terra workspace where data lives"}
        workspace_project: {description: "name of Terra project associated with Terra workspace"}
        grouping_column_name: {description: "name of column to be used for grouping/coloring - ex. gambit_predicted_taxon (organism)"}
        output_filename: {description: "name of output file containing visualizations"}
        custom_est_coverage_thresholds: {description: ""}
        custom_contig_thresholds: {description: ""}
        custom_assembly_thresholds: {description: ""}
        custom_mean_q_thresholds: {description: ""}
    }

    input {

        Array[String]   sample_ids
        String          input_table_name
        String          workspace_name
        String          workspace_project

        String?         grouping_column_name
        String?         output_filename

        String?         custom_est_coverage_thresholds
        String?         custom_contig_thresholds
        String?         custom_assembly_thresholds
        String?         custom_mean_q_thresholds
    }

    call create_viz {
        input:
            sample_ids                      =   sample_ids,
            workspace_name                  =   workspace_name,
            workspace_project               =   workspace_project,
            input_table_name                =   input_table_name,
            grouping_column_name            =   grouping_column_name,
            output_filename                 =   output_filename,
            custom_est_coverage_thresholds  =   custom_est_coverage_thresholds,
            custom_contig_thresholds        =   custom_contig_thresholds,
            custom_assembly_thresholds      =   custom_assembly_thresholds,
            custom_mean_q_thresholds        =   custom_mean_q_thresholds
    }

    output {
        File    vizualization_html     =   create_viz.html
    }
}

task create_viz {
    input {
        Array[String]   sample_ids
        String          workspace_name
        String          workspace_project
        String          input_table_name

        String?          grouping_column_name            =   "gambit_predicted_taxon"
        String           output_filename                 =   "QC_visualizations.pdf"

        String?         custom_est_coverage_thresholds
        String?         custom_contig_thresholds
        String?         custom_assembly_thresholds
        String?         custom_mean_q_thresholds

        String          docker                          =   "us-central1-docker.pkg.dev/pgs-automation/enterics-visualizations/create_visualization_html:v1"       
    }

    command {
        python3 /scripts/create_enterics_visualizations_html.py -s ~{sep=' ' sample_ids} \
                                                  -dt ~{input_table_name} \
                                                  -w ~{workspace_name} \
                                                  -bp ~{workspace_project} \
                                                  ~{"-g" + grouping_column_name} \
                                                  ~{"-o" + output_filename} \
                                                  ~{"-ect" + custom_est_coverage_thresholds} \
                                                  ~{"-cnt" + custom_contig_thresholds} \
                                                  ~{"-at" + custom_assembly_thresholds} \
                                                  ~{"-mqt" + custom_mean_q_thresholds} \
    }

    runtime {
        docker: docker
    }

    output {
        File html = "~{output_filename}"
    }
}
