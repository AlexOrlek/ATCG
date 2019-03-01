#!/bin/bash
set -e
set -u
set -o pipefail

inputfolder=${1}
find ${inputfolder} -maxdepth 1 -mindepth 1 -type f -name "*.gb" -exec cat {} \;
