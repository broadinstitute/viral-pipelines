#!/bin/bash
# requires $MODULE_VERSIONS to be set to point to a text file with equal-sign-separated values
set -e

# here is the older script that actually used sed to replace version strings
#while IFS='=' read module version; do
#	OLD_TAG=$module
#	NEW_TAG="$module:$version"
#	echo Replacing $OLD_TAG with $NEW_TAG in all task WDL files
#	sed -i -- "s|$OLD_TAG|$NEW_TAG|g" pipes/WDL/tasks/*.wdl
#done < $MODULE_VERSIONS

# this is the newer script that simply validates existing version strings
while IFS='=' read module version; do
	OLD_TAG=$module
	NEW_TAG="$module:$version"
	echo Validating that $OLD_TAG is $NEW_TAG in all task WDL files
    ! grep $OLD_TAG pipes/WDL/tasks/*.wdl | grep -v $NEW_TAG
done < $MODULE_VERSIONS
