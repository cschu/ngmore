#!/usr/bin/env nextflow

params.outdir = "outputs"

/*Channel
	.fromFilePairs("/g/scb/bork/schudoma/ngmore/data/*_R{1,2}.fastq.gz") { file -> file.name.replaceAll(/_R[12].fastq.gz$/, '') }
	.set { reads_ch }*/
Channel
	.fromPath("/g/scb/bork/schudoma/ngmore/data/*_R{1,2}.fastq.gz") 
	.map { file -> 
		def sampleId = file.name.replaceAll(/_R[12].fastq.gz$/, '')
		return tuple(sampleId, file)
	}
	.groupTuple()
	.set { reads_ch }

bwa_db = file("/g/scb/bork/schudoma/ngmore/ref/ecoli/*.gz")


process bbduk_preprocess {
	
	conda "bioconda::bbmap"

	publishDir "$params.outdir/$sampleId"

	input:
	set sampleId, file(reads) from reads_ch

	output:
	//stdout result
	//file '*.trim.fastq.gz' into preprocessed_ch
	tuple val(sampleId), file("${sampleId}*.trim.fastq.gz") into preprocessed_ch

	script:
	"""
	# echo ${reads[1]}
	bbduk.sh in1=${sampleId}_R1.fastq.gz in2=${sampleId}_R2.fastq.gz out1=${sampleId}_R1.trim.fastq.gz out2=${sampleId}_R2.trim.fastq.gz
	# bbduk.sh -h
	# touch ${sampleId}_R1.trim.fastq.gz
	"""
}


process bwa_mem {

	//conda "anaconda::openssl bioconda::bwa bioconda::samtools"
	conda "conda/bwa_samtools.yml"

	publishDir "$params.outdir/$sampleId"

	input:
	set sampleId, file(reads) from preprocessed_ch

	output:
	stdout result
	file("${sampleId}.sorted.bam") into aligned_sorted_ch

	script:
	"""
	#echo $sampleId
	bwa mem -a ${bwa_db[0]} $reads | samtools sort -o ${sampleId}.sorted.bam 
	"""

}
result.view { it }


/*#process splitLetters {
#
#    output:
#    file 'chunk_*' into letters
#
#    """
#    printf '${params.str}' | split -b 6 - chunk_
#    """
#}
#
#
#process convertToUpper {
#
#    input:
#    file x from letters.flatten()
#
#    output:
#    stdout result
#
#    """
#    cat $x | tr '[a-z]' '[A-Z]'
#    """
#}
#
#result.view { it.trim() }
#
#
#
#rule qc_bbduk:
#	message:
#			"Preprocessing read data with bbduk..."
#				input:
#						get_sample_files
#							output:
#									r1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
#											r2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
#												log:
#														join(QC_LOGDIR, "{sample}", "{sample}.qc_bbduk.log")
#															threads:
#																	8
#																		params:
#																				cmd = CMD_CALL + "bbduk.sh",
#																						adapters = config["resources"]["bb_adapters"],
#																								bbduk_params = config["params"]["bbduk"]
#																									shell:
#																											TIME_V + " {params.cmd}" + \
#																													" -Xmx30g t={threads} in1={input[0]} in2={input[1]} out1={output.r1} out2={output.r2}" + \
#																															" ref={params.adapters}" + \
#																																	" {params.bbduk_params}" + \
#																																			" &> {log}"
#
#*/
