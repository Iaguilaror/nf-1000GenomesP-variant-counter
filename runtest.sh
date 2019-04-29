echo -e "======\n Testing NF execution \n======" \
&& rm -rf test/results/ \
&& nextflow run calculate-vcf-stats.nf \
	--vcf_dir test/data/ \
  --AN_cutoff 154 \
	--output_dir test/results \
	-resume \
	-with-report test/results/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag test/results/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
