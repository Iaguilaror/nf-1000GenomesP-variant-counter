#!/bin/bash

find -L . \
  -type f \
  -name "*.counts.tsv" \
| sed "s#chr.*_GRCh38.*#allchrom_counts.tsv#" \
| sort -u \
| xargs mk
