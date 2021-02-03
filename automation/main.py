"""Google Cloud Function to launch a Terra workflow."""

import os
from typing import Any, Dict
from utils import prepare_and_launch


def launch_workflow(data: Dict[Any, Any], context: Any):
    """Entry point for execution via a Cloud Function.

      This Cloud Function reads configuration from environment variables and the triggering event.

      This example workflow uses entities from a data table so that workflow parameter values do
      not need to be hardcoded here in this script.

      Environment variables:
        WORKSPACE_NAMESPACE: The project id of the Terra billing project in which the workspace resides.
        WORKSPACE_NAME: The name of the workspace in which the workflow resides
        METHOD_NAMESPACE: The namespace of the workflow method.
        METHOD_NAME: The name of the workflow method.
        SECRET_PATH: The 'Resource ID' of the service account key stored in Secret Manager. Or, if
          testing locally, the filepath to the JSON key for the service account.
        TRIGGER_PARAMETER_NAME: The name of the workflow parameter to receive the path to the triggering file.
          Defaults to `MyWorkflowName.aCloudStorageFilePath`.

      Args:
        event:  The dictionary with data specific to this type of event.
                The `data` field contains a description of the event in
                the Cloud Storage `object` format described here:
                https://cloud.google.com/storage/docs/json_api/v1/objects#resource
        context: Metadata of triggering event.
      Returns:
        None; the side effect is the execution of a parameter-parallel Terra workflow.
    """

    # extract file information from the triggering PubSub message.
    file_name = data.get('name')

    # only process .tar.gz files uploaded to the flowcells/broad-walkup/ directory
    if not (file_name.endswith(".tar.gz") and "flowcells/broad-walkup/" in file_name):
        print(f"Invalid file for Broad viral sequencing ({file_name}). Workflow not triggered.")
        return

    bucket_name = data.get('bucket')
    file_path = f"gs://{bucket_name}/{file_name}"
    print(f"Processing input file: {file_name}; full path: {file_path}")
    # default to the parameter name from the example workflow.
    workflow_parameters = {
        os.getenv("TRIGGER_PARAMETER_NAME"): f"\"{file_path}\"",
        }

    prepare_and_launch(
        workspace_namespace=os.getenv("WORKSPACE_NAMESPACE"),
        workspace_name=os.getenv("WORKSPACE_NAME"),
        method_namespace=os.getenv("METHOD_NAMESPACE"),
        method_name=os.getenv("METHOD_NAME"),
        secret_path=os.getenv("SECRET_PATH"),
        workflow_parameters=workflow_parameters
        )


if __name__ == "__main__":
    """Entry point of manual execution for testing purposes."""

    # This example parameter is a world-readable file.
    # gs://genomics-public-data/platinum-genomes/other/platinum_genomes_sample_info.csv

    # launch_workflow(data={"bucket": "viral-sequencing-automation",
    #                       "name": "flowcells/broad-walkup/pizza.txt"},
    #                 context=None)

    file_name = "flowcells/broad-walkup/pizza.txt"
    file_path = "gs://viral-sequencing-automation/flowcells/broad-walkup/pizza.txt"
    workflow_parameters = {
      "HelloWorldPlus.input_file": f"\"{file_path}\""
    }

    print(f"input file: {file_name}; full path: {file_path}")

    prepare_and_launch(
        workspace_namespace="pathogen-genomic-surveillance",
        workspace_name="Viral-Genomics-Lab-private_sushmac_dfe385",
        method_namespace="pathogen-genomic-surveillance",
        method_name="hello-world-plus",
        secret_path="projects/115398405846/secrets/full_viral_workflow_sa_json_key/versions/1",
        workflow_parameters=workflow_parameters
        )
