version 1.0

workflow CreateEntericsQCViz {

    meta {
        allowNestedInputs: true
    }

    call create_viz

    output {
        File    visualization_html     =   create_viz.html
    }
}

task create_viz {
    input {
        Array[String]   sample_ids
        String          workspace_name
        String          workspace_project
        String          input_table_name

        String          grouping_column_name            =   "gambit_predicted_taxon"
        String          output_filename                 =   "QC_visualizations.html"

        File?           thresholds_file

        String          docker                           =   "us-central1-docker.pkg.dev/pgs-automation/enterics-visualizations/create_visualization_html:v4"       
    }

    parameter_meta {
        sample_ids: {description: "selected rows of data from data table which will be used for plotting"}
        input_table_name: {description: "name of the Terra data table - root entity type - from where input data is selected"}
        workspace_name: {description: "name of Terra workspace where data lives"}
        workspace_project: {description: "name of Terra project associated with Terra workspace"}
        grouping_column_name: {description: "name of column to be used for grouping/coloring - ex. gambit_predicted_taxon (organism)"}
        output_filename: {description: "name of output file containing visualizations"}
        thresholds_file: {description: "JSON file containing custom thresholds"}
    }

    command {
        python3 /scripts/create_enterics_visualizations_html.py -s "~{sep='" "' sample_ids}" \
                                                                -dt "~{input_table_name}" \
                                                                -w "~{workspace_name}" \
                                                                -bp "~{workspace_project}" \
                                                                ~{'-g "' + grouping_column_name + '"'} \
                                                                ~{'-o "' + output_filename + '"'} \
                                                                ~{'-t "' + thresholds_file + '"'}
    }

    runtime {
        docker: docker
    }

    output {
        File html = output_filename
    }
}
