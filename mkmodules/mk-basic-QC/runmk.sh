#!/bin/bash

find -L . \
  -type f \
  -name "*.tsv" \
| sed "s#.tsv#.basicQC.pdf#" \
| xargs mk
