#!/usr/bin/env nextflow

/*================================================================
The MORETT LAB presents...

  The VCF stats calculator

- A variants per sample counter tool

==================================================================
Version: 0.0.1
Project repository: https://bitbucket.org/morett_lab/nf-vcf-table-description/
==================================================================
Authors:

- Bioinformatics Design
 Israel Aguilar-Ordonez (iaguilaror@gmail)

- Bioinformatics Development
 Israel Aguilar-Ordonez (iaguilaror@gmail)

- Nextflow Port
 Israel Aguilar-Ordonez (iaguilaror@gmail)

=============================
Pipeline Processes In Brief:

Pre-processing:

Core-processing:
  _001_create_sample_requests
  _002_count_variants
  _003_concatenate_stats

Post-processing:
  _pos1_tag_samples <<-- also plots basic QC

================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
  The VCF stats calculator
  - A variants per sample counter tool
  v${version}
  ==========================================

	Usage:

  nextflow run calculate-vcf-stats.nf --vcf_dir <path to input 1> [--output_dir path to results]

    --vcf_dir    <- directory containing gz compressed vcf files;
		    vcf must be in .vcf.gz format, and tabix indexed;
        vcf must have sample level GT data;
    --output_dir     <- directory where results, intermediate and log files will be stored;
				default: same level dir where --vcf_dir resides
	  --help           <- Show Pipeline Information
	  --version        <- Show Pipeline version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "countVariants"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.vcf_dir = false  //if no input path is provided, value is false to provoke the error during the parameter validation block
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "VCF variants per sample counter Pipeline v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at FEB 2019
*/
nextflow_required_version = '19.01'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Extend Align required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  Pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  Pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*//////////////////////////////
  INPUT PARAMETER VALIDATION BLOCK
  TODO (iaguilar) check the extension of input queries; see getExtension() at https://www.nextflow.io/docs/latest/script.html#check-file-attributes
*/

/* Check if vcf file was provided and AN_cutoff
    if they were not provided, they keep the 'false' value assigned in the parameter initiation block above
    and this test fails
*/
if ( !params.vcf_dir ) {
  log.error " Please provide the --vcf_dir file \n\n" +
  " For more information, execute: nextflow run calculate-vcf-stats.nf --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --vcf
*/
params.output_dir = file(params.vcf_dir).getParent()

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable (pipeline_name) defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*//////////////////////////////
  LOG RUN INFORMATION
*/
log.info"""
==========================================
The VCF stats calculator
- A variants per sample counter tool
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['VCF dir']			= params.vcf_dir
pipelinesummary['Results Dir']		= results_dir
pipelinesummary['Intermediate Dir']		= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

/*//////////////////////////////
  PIPELINE START
*/

/*
	READ INPUTS
*/
/* Define function to pair vcf and tbi inputs */
def get_chr_prefix = { file -> file.toString().tokenize('.')[0,1].join(".")}

/* Load vcf file into channel */
Channel
  .fromPath("${params.vcf_dir}/*.vcf.gz")
  .into{ vcf_inputs; vcf_inputs_2 }

/* _001_create_sample_requests */
/* Read mkfile module files */
Channel
  .fromPath("${workflow.projectDir}/mkmodules/mk-create-sample-requests/*")
  .toList()
  .set{ mkfiles_001 }

process _001_create_sample_requests {

	input:
	file vcf from vcf_inputs
	file mk_files from mkfiles_001

	output:
	file "*" into results_001_create_sample_requests mode flatten

	"""
	bash runmk.sh
	"""
}

/* _002_count_variants */
/* gather all vcf_inputs_2 to send to nextmodule */
vcf_inputs_2
	.toList()
	.set{ all_vcf_inputs }

/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/mk-count-variants/*")
	.toList()
	.set{ mkfiles_002 }

process _002_count_variants {

	publishDir "${intermediates_dir}/_002_count_variants/",mode:"symlink"

	input:
	file sample_request from results_001_create_sample_requests
	file all_vcf from all_vcf_inputs
	file mk_files from mkfiles_002

	output:
	file "*.stats" into results_002_count_variants

	"""
	bash runmk.sh
	"""

}

/* _003_concatenate_stats */
/* Gather all files before concatneation */
results_002_count_variants
  .toList()
  .set{ inputs_for_003 }


/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/mk-concatenate-stats/*")
	.toList()
	.set{ mkfiles_003 }

process _003_concatenate_stats {

	publishDir "${intermediates_dir}/_003_concatenate_stats/",mode:"symlink"

	input:
	file tsvfiles from inputs_for_003
	file mk_files from mkfiles_003

	output:
	file "*.allstats.tsv" into results_003_concatenate_stats mode flatten

	"""
	bash runmk.sh
	"""

}

/* _pos1_tag_samples */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/mk-tag-samples/*")
	.toList()
	.set{ mkfiles_pos1 }

process _pos1_tag_samples {

	publishDir "${results_dir}/_pos1_tag_samples/",mode:"copy"

	input:
	file tsvfiles from results_003_concatenate_stats
	file mk_files from mkfiles_pos1

	output:
	file "*.tagged.tsv*" into results_pos1_tag_samples

	"""
	export METADATA="${params.metadata}"
	bash runmk.sh
	"""

}
