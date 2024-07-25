#!/bin/bash

echo "pip installing test-related packages (Sphinx, etc.)"
pip install -r docs/requirements.txt

# list what was installed
pip freeze
