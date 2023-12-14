version 1.0

workflow CreateEntericsQCViz {

    parameter_meta {
        sample_ids: {description: "selected rows of data from data table which will be used for plotting"}
        input_table_name: {description: "name of the Terra data table - root entity type - from where input data is selected"}
        workspace_name: {description: "name of Terra workspace where data lives"}
        workspace_project: {description: "name of Terra project associated with Terra workspace"}
        grouping_column_name: {description: "name of column to be used for grouping/coloring - ex. gambit_predicted_taxon (organism)"}
        output_filename: {description: "name of output file containing visualizations"}
    }

    input {

        Array[String]   sample_ids
        String          input_table_name
        String          workspace_name
        String          workspace_project

        String?         grouping_column_name
        String?         output_filename
    }

    call create_viz {
        input:
            sample_ids              =   sample_ids,
            workspace_name          =   workspace_name,
            workspace_project       =   workspace_project,
            input_table_name        =   input_table_name,
            grouping_column_name    =   grouping_column_name,
            output_filename         =   output_filename
    }

    output {
        File    viz_pdf     =   create_viz.vizualizations
    }
}

task create_viz {
    input {
        Array[String]   sample_ids
        String          workspace_name
        String          workspace_project
        String          input_table_name

        String          grouping_column_name    = "gambit_predicted_taxon"
        String          output_filename         = "QC_vizualizations.pdf"

        String  docker                          =   "broadinstitute/horsefish:pgs_visualizations_dev"        
    }

    command {
        python3 /scripts/create_visualizations.py -s ~{sep=' ' sample_ids} \
                                                  -t ~{input_table_name} \
                                                  -w ~{workspace_name} \
                                                  -p ~{workspace_project} \
                                                  ~{"-g" + grouping_column_name} \
                                                  ~{"-o" + output_filename}
    }

    runtime {
        docker: docker
    }

    output {
        File vizualizations = "~{output_filename}"
    }
}