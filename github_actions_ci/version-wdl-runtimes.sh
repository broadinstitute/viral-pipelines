#!/usr/bin/env bash

# use sed to replace version strings of docker images based on versions defined in txt file
#
# skip this replacement for any version string line with the comment "#skip-global-version-pin"
#
# requires $MODULE_VERSIONS to be set to point to a text file with equal-sign-separated values
# export MODULE_VERSIONS="./requirements-modules.txt" && ./github_actions_ci/version-wdl-runtimes.sh

function absolute_path() {
    local SOURCE="$1"
    while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
        DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
        if [[ "$OSTYPE" == "darwin"* ]]; then
            SOURCE="$(readlink "$SOURCE")"
        else
            SOURCE="$(readlink -f "$SOURCE")"
        fi
        [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
    done
    echo "$SOURCE"
}
SOURCE="${BASH_SOURCE[0]}"
SCRIPT=$(absolute_path "$SOURCE")
SCRIPT_DIRNAME="$(dirname "$SOURCE")"
SCRIPTPATH="$(cd -P "$(echo $SCRIPT_DIRNAME)" &> /dev/null && pwd)"
SCRIPT="$SCRIPTPATH/$(basename "$SCRIPT")"
REPO_PATH="$(realpath "$SCRIPTPATH/../")"

# set MODULE_VERSIONS with default value of "$(realpath "$SCRIPTPATH/../requirements-modules.txt")", assumed to be located one level above this script
MODULE_VERSIONS="${MODULE_VERSIONS:-"${REPO_PATH}/requirements-modules.txt"}"
echo $MODULE_VERSIONS

printf "Updating docker image tags in WDL files with those in ${MODULE_VERSIONS}\n"
printf "Performing replacements in ${REPO_PATH}/pipes/WDL/tasks/*.wdl files\n\n"


# check if sed is GNU or BSD (macOS) version
# macOS/BSD sed does not support '--version', and if/when it does it will not emit 'GNU sed'
if ! (sed --version | grep --quiet 'GNU sed') &> /dev/null; then
  echo "using sed with macOS/BSD syntax"
  # if on macOS, use BSD sed's -i '' (empty backup suffix) for true in-place editing
  SED_OPTS=(-i '')
else
  echo "using sed with GNU/Linux syntax"
  # GNU sed: -i alone means no backup file; keep -E for extended regex
  SED_OPTS=(-i -s -e)
fi

while IFS='=' read module version; do
  OLD_TAG=$module
  NEW_TAG="$module:$version"
  NEW_TAG_BOLD="$module:$(tput bold)$version$(tput sgr0)"
  printf "Replacing: %-14s \n with tag: %-14s \n\n" "$OLD_TAG" "$NEW_TAG_BOLD"

  sed "${SED_OPTS[@]}" "/^\(.*\)[[:space:]]*#skip-global-version-pin[[:space:]]*$/!s|$OLD_TAG[^\"']*|$NEW_TAG|g" "${REPO_PATH}/pipes/WDL/tasks/"*.wdl
done < $MODULE_VERSIONS

printf "Replacements skipped for lines marked with '#skip-global-version-pin' \n\n"