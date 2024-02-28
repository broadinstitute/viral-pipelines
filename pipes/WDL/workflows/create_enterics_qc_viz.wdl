version 1.0

workflow CreateEntericsQCViz {

    meta {
        allowNestedInputs: true
    }

    call create_viz {
    }

    output {
        File    visualization_html     =   create_viz.html
    }
}

task create_viz {
    input {
        Array[String]   sample_ids
        Array[String]   gambit_predicted_taxon
        Array[String]   est_coverage_clean
        Array[String]   number_contigs
        Array[String]   assembly_length

        File?           thresholds_file

        String          output_filename                  =   "QC_visualizations.html"

        String          docker                           =   "us-central1-docker.pkg.dev/pgs-automation/enterics-visualizations/create_visualization_html:v5"       
    }

    parameter_meta {
        sample_ids: {description: "selected samples of data"}
        gambit_predicted_taxon: {description: "predicted gambits for selected samples"}
        est_coverage_clean: {description: "estimated coverage metrics for selected samples"}
        number_contigs: {description: "number contigs metrics for selected samples"}
        assembly_length: {description: "assembly length metrics for selected samples"}
        output_filename: {description: "name of output file containing visualizations"}
        thresholds_file: {description: "JSON file containing custom thresholds"}
    }

    command {

        python3 /scripts/create_enterics_visualizations_html.py -s "~{sep='" "' sample_ids}" \
                                                                -g "~{sep='" "' gambit_predicted_taxon}" \
                                                                -ecc ~{sep=' ' est_coverage_clean} \
                                                                -nc ~{sep=' ' number_contigs} \
                                                                -al ~{sep=' ' assembly_length} \
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
