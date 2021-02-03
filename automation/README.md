# launch_workflow Cloud Function

The scripts in this directory can be used to deploy a Google Cloud Function (CF) that launches a Terra workflow. The Cloud Function is triggered when a file is created or modified in a specific Google Cloud Storage bucket.

Note: this example is specifically written assuming:

* the triggering file is an input parameter to the Terra workflow
* the most recently created entity set should be used, or override this default via `ENTITY_SET_NAME`
* the Terra UI is used to control *all other aspects* of the workflow configuration, such as:
  * which version of the method is to be used
  * whether call caching is enabled
  * whether to delete intermediate outputs
  * the type of the root entity
  * how the columns of the root entity map onto the workflow parameters

**New to Cloud Functions and Secret Manager?** Don't start here. Instead see:

* [Python Quickstart | Cloud Functions Documentation | Google Cloud](http://cloud/functions/docs/quickstart-python)
* [Quickstart | Secret Manager Documentation | Google Cloud](http://cloud/secret-manager/docs/quickstart)

**New to Terra groups and permissions?** Don't start here. Instead see:

* [Pet accounts and proxy groups – Terra Support](https://support.terra.bio/hc/en-us/articles/360031023592-Understanding-and-setting-up-a-proxy-group)
* [Accessing advanced GCP features in Terra – Terra Support](https://support.terra.bio/hc/en-us/articles/360051229072)

# Setup

## Configure the Terra Workflow

1. Create or identify the Terra workflow you want to trigger via a Cloud Function.
    * For a simple example, see [example.wdl](./example.wdl).
1. Upload an entity set of all but one of the parameters for your workflow. The one parameter absent from the entity set should be the one that is the path to the Cloud Storage file that will trigger the workflow.
    * Create this entity set even if the set would only contain one entity (one row of parameter values).
    * For a simple example, see [example.tsv](./example.tsv).
1. **IMPORTANT** Make sure you have used the Terra User Interface to run your workflow successfully at least once using the entity set.

## Configure the Cloud and Terra execution environments

1. Create or identify the Google Cloud Platform project where this Cloud Function will live / be billed.
    1. Enable the [Cloud Functions API](https://console.developers.google.com/apis/library/cloudfunctions.googleapis.com) for that project.
    1. Enable the [Secret Manager API](https://console.cloud.google.com/marketplace/product/google/secretmanager.googleapis.com) for that project.
    1. Enable the [Cloud Build API](https://console.developers.google.com/apis/library/cloudbuild.googleapis.com) for this project.
1. Create or identify the bucket in that project that will be used for files to trigger the Cloud Function.
1. Create a new [service account](https://console.cloud.google.com/iam-admin/serviceaccounts) within that project. It will look like `<your service account name>@<your project id>.iam.gserviceaccount.com`.
1. Download a JSON key for the newly created service account and store it in a safe place.
1. Upload the JSON key to [Secret Manager](https://console.cloud.google.com/security/secret-manager).
1. For the newly created secret, in column "Actions", choose "Copy Resource ID". You'll need this value for the deployment step and it will look like: `projects/<your project number>/secrets/<your secret name>/versions/1`.
1. Give the service account the needed permissions to items within the project.
    1. Grant role `Secret Manager Secret Accessor` to the service account.
    1. Grant role `Storage Object Viewer` to the service account on the bucket that will trigger the Cloud Function.
1. [Register the service account](https://github.com/broadinstitute/terra-tools/tree/master/scripts/register_service_account) in Terra.
1. [Create a new Terra group](https://app.terra.bio/#groups) to hold only the service account. It will look like `
<your group name>@firecloud.org`.
1. Give the Terra Group the needed permissions to items within the project.
    1. Add the Terra group to the billing project (a.k.a. workspace namespace).
    1. Grant role `Storage Object Viewer` to the Terra group on the bucket that will trigger the Cloud Function. This is because, presumably, the workflow will need to read the contents of the file that triggered the workflow.
    1. If your workflow uses a private Google Container Repository image, also grant role `Storage Object Viewer` to the Terra group on the bucket that holds the image (e.g., `artifacts.<your project id>.appspot.com`).
1. Give the Terra Group the needed permissions within Terra.
    1. Share the Terra workspace with the Terra group, granting `Writer` and `Can compute` access.
    1. If the workspace has an Authorization Domain, ensure that the Terra group is a member of the Authorization Domain.
1. If the workflow method is private, share the workflow method. Note that for the Broad Methods repository, sharing with groups is not yet supported. Instead share the method with the service account.

# Deploy

## Deploy the Cloud Function

1. Update the `config.yaml` file in this directory with the Terra workspace and workflow information you want to use. These values will be used as environment variable inputs to the deployment script.
2. Run the deployment script. Usage:
```
./deploy_cloud_function.sh PROJECT BUCKET SERVICE_ACCOUNT
```
e.g.
```
./deploy_cloud_function.sh launch-workflow-test launch-workflow-test-input \
    launch-workflow-sa@launch-workflow-test.iam.gserviceaccount.com
```

# Test

## Test the Cloud Function trigger

Test the cloud function by placing a file in the bucket to trigger the workflow. You can do this using the [cloud console](https://console.cloud.google.com/storage/browser/) or with the following script:
```
./test_cloud_function.sh BUCKET
```
e.g.
```
./test_cloud_function.sh launch-workflow-test-input
```

## [Optional] Run the Cloud Function code locally.

This is useful for rapid debugging of permissions or other issues.

1. Create a [virtual environment](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment) for Python 3.
```
# Tip: run this outside of your git clone, such as under ${HOME}/my-venvs
cd path/to/where/I/keep/my/venvs
python3 -m venv test-cf-env
```
2. Activate the [virtual environment](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment).
```
source test-cf-env/bin/activate
```
3. Install the dependencies.
```
# Run this from directory scripts/launch_workflow_cf in your clone.
pip3 install -r requirements.txt
```
4. Invoke the workflow manually.
```
# Run this from directory scripts/launch_workflow_cf in your clone.

export WORKSPACE_NAMESPACE="morgan-fieldeng"
export WORKSPACE_NAME="launch-workflow-test"
export METHOD_NAMESPACE="morgan-fieldeng"
export METHOD_NAME="hello-world-plus-prefixes"
export TRIGGER_PARAMETER_NAME="MyWorkflowName.aCloudStorageFilePath"
export SECRET_PATH="/local/filepath/to/service-account-key.json"

python3 main.py
```
