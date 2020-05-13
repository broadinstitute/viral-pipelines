#!/bin/bash
set -e -o pipefail

echo "version: 1.2"
echo "workflows:"
for WDL in $*; do
	echo " - name: $(basename $WDL .wdl)"
	echo "   subclass: WDL"
	echo "   primaryDescriptorPath: $WDL"
	echo "   testParameterFiles:"
	echo "    - empty.json"
done
