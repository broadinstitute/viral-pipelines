# -*- coding: utf-8 -*-
"""Cloud Function Utilities."""

import os
import json
import requests

from firecloud.errors import FireCloudServerError

from google.cloud.secretmanager_v1 import SecretManagerServiceClient
from oauth2client.client import GoogleCredentials
from oauth2client.service_account import ServiceAccountCredentials


def prepare_and_launch(workspace_namespace, workspace_name, method_namespace, method_name,
                       secret_path, workflow_parameters):
    """Launch a Terra workflow programmatically.

      Arguments `workflow_parameters` and `entity_set_name` are used to specify the parameters to the
      workflow. The Terra UI is used to control *all other aspects* of the workflow configuration, such as:
      * which version of the method is to be used
      * whether call caching is enabled
      * whether to delete intermediate outputs
      * the type of the root entity
      * how the columns of the root entity map onto the workflow parameters

      Args:
        workspace_namespace: The project id of the Terra billing project in which the workspace resides.
        workspace_name: The name of the workspace in which the workflow resides.
        method_namespace: The namespace of the workflow method.
        method_name: The name of the workflow method.
        secret_path: The 'Resource ID' of the service account key stored in Secret Manager. Or, if
          testing locally, the filepath to the json key for the service account.
        workflow_parameters: A dictionary of key/value pairs to merge with the inputs from the workflow configuration,
          overwriting any duplicate keys.
        entity_set_name: The name of the entity set to be used for all other workflow parameters. Defaults to
          the most recently created entity set of the root entity type.
      Returns:
        None; the side effect is the execution of a parameter-parallel Terra workflow.
    """
    # Get access token and and add to headers for requests.
    headers = {"Authorization": "bearer " + get_access_token(secret_path=secret_path)}

    # Get the workflow config.
    workflow_config_response = get_workflow_method_config(
        workspace_namespace=workspace_namespace,
        workspace_name=workspace_name,
        method_namespace=method_namespace,
        method_name=method_name,
        headers=headers
    )
    check_fapi_response(response=workflow_config_response, success_code=200)
    response_string = json.dumps(workflow_config_response.json()).replace("\n", " ")
    print(f"Current default workflow configuration:  {response_string}")

    # Merge the workflow inputs configured in the Terra UI with the items from 'workflow_parameters'
    # in-place, overwriting existing keys.
    workflow_config_json = workflow_config_response.json()
    workflow_config_json["inputs"].update(workflow_parameters)

    # Update the workflow configuration.
    updated_workflow_response = update_workflow_method_config(
        workspace_namespace=workspace_namespace,
        workspace_name=workspace_name,
        method_namespace=method_namespace,
        method_name=method_name,
        body=workflow_config_json,
        headers=headers
    )
    check_fapi_response(response=updated_workflow_response, success_code=200)

    # Launch the workflow.
    create_submisson_response = create_submission(
        workspace_namespace=workspace_namespace,
        workspace_name=workspace_name,
        method_namespace=method_namespace,
        method_name=method_name,
        headers=headers
    )
    response_string = json.dumps(create_submisson_response.json()).replace("\n", " ")
    print(f"Submission response:  {response_string}")

    check_fapi_response(response=create_submisson_response, success_code=201)
    submission_id = create_submisson_response.json()["submissionId"]
    print(f"Successfully created submission: submissionId = {submission_id}.")


def get_access_token(secret_path):
    """Conditionally create an access token with the minimum necessary scopes.

    When run as a Cloud Function, a service account's JSON credentials file from Secret Manager
    is used to populate the access token. When run locally, either the service account JSON key file or the
    application default credential is used to populate the access token.

    Args:
      secret_path: The 'Resource ID' of the service account key stored in Secret Manager. Or, if
        testing locally, the filepath to the JSON key for the service account.
    Returns:
      An access token.
    """
    scopes = ["https://www.googleapis.com/auth/userinfo.profile", "https://www.googleapis.com/auth/userinfo.email"]

    if not secret_path:  # Running locally as a user.
        credentials = GoogleCredentials.get_application_default()
        credentials = credentials.create_scoped(scopes)
    elif os.path.isfile(secret_path):  # Running locally as the service account.
        credentials = ServiceAccountCredentials.from_json_keyfile_name(secret_path, scopes=scopes)
    else: # Running inside the Cloud Function.
        # Retrieve the secret from the secret manager API.
        client = SecretManagerServiceClient()
        response = client.access_secret_version(secret_path)
        service_account_key = response.payload.data.decode("utf-8")
        json_acct_info = json.loads(service_account_key)
        credentials = ServiceAccountCredentials.from_json_keyfile_dict(json_acct_info, scopes=scopes)

    return credentials.get_access_token().access_token

# Helper methods for Firecloud API calls. We cannot use the methods in firecloud.api because
# we need to explicitly pass a header for auth.


def check_fapi_response(response, success_code):
    if response.status_code != success_code:
        print(response.content)
        raise FireCloudServerError(response.status_code, response.content)


def get_workflow_method_config(workspace_namespace, workspace_name, method_namespace, method_name, headers):
    uri = f"https://api.firecloud.org/api/workspaces/{workspace_namespace}/{workspace_name}/method_configs/{method_namespace}/{method_name}"
    return requests.get(uri, headers=headers)


def update_workflow_method_config(workspace_namespace, workspace_name, method_namespace, method_name, body, headers):
    uri = f"https://api.firecloud.org/api/workspaces/{workspace_namespace}/{workspace_name}/method_configs/{method_namespace}/{method_name}"
    return requests.post(uri, json=body, headers=headers)


def get_entities(workspace_namespace, workspace_name, etype, headers):
    uri = f"https://api.firecloud.org/api/workspaces/{workspace_namespace}/{workspace_name}/entities/{etype}"
    return requests.get(uri, headers=headers)


def create_submission(workspace_namespace, workspace_name, method_namespace, method_name, headers,
                      entity=None, etype=None, expression=None, use_callcache=True):
    uri = f"https://api.firecloud.org/api/workspaces/{workspace_namespace}/{workspace_name}/submissions"
    body = {
        "methodConfigurationNamespace": method_namespace,
        "methodConfigurationName": method_name,
        "useCallCache": use_callcache
    }
    if etype:
        body["entityType"] = etype
    if entity:
        body["entityName"] = entity
    if expression:
        body["expression"] = expression
    return requests.post(uri, json=body, headers=headers)
