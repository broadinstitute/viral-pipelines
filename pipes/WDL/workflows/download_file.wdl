version 1.0

#DX_SKIP_WORKFLOW

import "../tasks/tasks_utils.wdl" as terra

workflow download_file {
    meta {
        description: "Downloads an http[s] file. Helpful if this is not natively supported by the WDL execution backend for File inputs."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call terra.download_from_url

    output {
        File output_file                      = select_first([download_from_url.downloaded_response_file, download_from_url.passthrough_url])

        # one or the other will be returned, depending on the download method
        # an http[s] url will be downloaded to a file and available via downloaded_response_file
        # other urls (i.e. localizable paths like 'gs://*', 'drs://') will be available via passthrough_url
        File?   downloaded_response_file      = download_from_url.downloaded_response_file
        String? passthrough_url               = download_from_url.passthrough_url
        
        # optional fields only returned in the case of a downloaded file
        File?   downloaded_response_headers   = download_from_url.downloaded_response_headers
        String? md5_sum_of_response_file      = download_from_url.md5_sum_of_response_file
        Int?    file_size_bytes               = download_from_url.file_size_bytes

        # boolean flag to indicate if the download task passed through the input url instead of downloading the file
        Boolean passed_through_input_url_instead_of_downloading  = download_from_url.passed_through_input_url_instead_of_downloading
    }
}
