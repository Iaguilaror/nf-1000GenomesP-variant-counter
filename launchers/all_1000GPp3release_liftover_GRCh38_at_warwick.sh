nextflow run calculate-vcf-stats.nf \
	--vcf_dir real-data/GRCh38-liftover-all/ \
	--output_dir real-data/GRCh38-liftover-all/results \
	-resume \
	-with-report real-data/GRCh38-liftover-all/results/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag real-data/GRCh38-liftover-all/results/`date +%Y%m%d_%H%M%S`.DAG.html
