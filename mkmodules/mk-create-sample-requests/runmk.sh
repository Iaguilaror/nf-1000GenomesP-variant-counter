#!/bin/bash

find -L . \
  -type f \
  -name "*.vcf.gz" \
| sed "s#.vcf.gz#.REQUEST_SAMPLES#" \
| xargs mk

