#!/usr/bin/env bash

## Find files with .vcf extension
find -L . \
	-type f \
	-name '*.allstats.tsv' \
| sed 's#.tsv#.tagged.tsv#' \
| xargs mk
