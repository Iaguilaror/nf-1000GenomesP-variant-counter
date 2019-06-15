#!/usr/bin/env bash

## Find files with .vcf extension
find -L . \
	-type f \
	-name '*.stats' \
| sed "s#\\.subsampled_.*#\\.allstats\\.tsv#" \
| sort -u \
| xargs mk
