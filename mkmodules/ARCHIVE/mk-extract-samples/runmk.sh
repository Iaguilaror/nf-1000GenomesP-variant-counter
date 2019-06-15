#!/usr/bin/env bash

## Find files with .vcf extension
find -L . \
	-type f \
	-name '*.vcf.tmp' \
| sed 's#.vcf.tmp#.vcf#' \
| xargs mk
