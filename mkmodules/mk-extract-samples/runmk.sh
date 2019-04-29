#!/usr/bin/env bash

## Find files with .vcf extension
find -L . \
	-type f \
	-name '*.vcf.gz' \
	! -name '*.sample_*.vcf' \
| sed 's#.vcf.gz#.EXTRACT_SAMPLES#' \
| xargs mk
