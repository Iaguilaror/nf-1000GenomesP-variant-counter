#!/usr/bin/env bash

## Find files with .vcf extension
find -L . \
	-type f \
	-name '*.stats.tmp' \
| sed 's#.stats.tmp#.stats#' \
| xargs mk
