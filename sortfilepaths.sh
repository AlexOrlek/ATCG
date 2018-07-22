#!/bin/bash
set -e
set -u
set -o pipefail

sort -k1,1V -o ${1} ${1}  #-k command must precede -o command
