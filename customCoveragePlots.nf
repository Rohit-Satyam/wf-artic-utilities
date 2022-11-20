nextflow.enable.dsl=2
params.bams = "wf-articresults/*.bam"
params.fqpassDir= "$projectDir/fastq_pass"
params.fasta="wf-articresults/*.fasta"
params.outdir="results/"
params.prefix="summary"
params.maxlen=700
params.minlen=200
params.cpus = 20
params.help = false

// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step 7: BAM Coverage
=============================================
Usage:
	nextflow run 07_coverage.nf --bams ${params.bams} --fqpassDir ${params.fqpassDir} --fasta ${params.fasta} --outdir ${params.outdir} --prefix ${params.prefix}
Input:
	* --bams: Path of .bam files from computing per-base coverage. Defult [${params.bams}]
	* --fqpassDir: Path of fastq_pass directory. Default [${params.fqpassDir}]
	* --fasta: Path of consensus fasta file. Default [${params.fasta}]
	* --outdir: Output directory to save results. Default [${params.outdir}]
	* --prefix: Prefix for fastcat output. Default [${params.prefix}]
	* --maxlen: Maximum length of reads. Default [${params.maxlen}]
	* --minlen: Minimum length of reads. Default [${params.minlen}]
	* --cpus: Number of cores. Default [${params.cpus}]
"""

exit 0
}


process MOSDEPTH{
	publishDir "$params.outdir/mosdepth", mode: 'copy'
	cpus params.cpus
	input:
	tuple val(sid), path(bam)
	output:
	path "*"
	path "*per-base.bed.gz", emit: per_base_cov

	script:
	"""
	mosdepth -t ${task.cpus} -x -b 1 ${sid} ${bam.toRealPath()}
	"""
}

process GETNBED{
	publishDir "$params.outdir", mode: 'copy'
	cpus params.cpus
	input:
	tuple val(sid), path(fasta)
	output:
	path "*.bed", emit: getn_bed

	script:
	"""
	python $projectDir/bin/getmask.py ${fasta} > ${params.prefix}.maskedRegions.bed
	"""


}

process FASTCAT{
	publishDir "$params.outdir", mode: 'copy'
	cpus params.cpus
	output:
	path "*.stats", emit: fastcat_stats
	path "*.txt"

	script:
	"""
	fastcat -r ${params.prefix}.stats -x ${params.fqpassDir} -a ${params.minlen} -b ${params.maxlen} -f ${params.prefix}.txt > /dev/null
	"""
}

process ASSEMBLYQC{
	publishDir "$params.outdir", mode: 'copy'
	cpus params.cpus
	input:
	path(bed)
	path(stats)
	output:
	path "*.html"
	path "temp.cov.csv", emit:barcode_cov_file

	script:
	"""
	Rscript $projectDir/bin/get_coverage.R -b ${bed} -fs ${stats} -o ${params.prefix}
	"""
}

process HORIZONTALCOVERAGE{
	publishDir "$params.outdir", mode: 'copy'
	cpus params.cpus
	input:
	file(list_of_bed)
	path(cov)
	output:
	path "*"

	script:
	"""
	Rscript $projectDir/bin/horizontalCoverageplt.R -l ${list_of_bed} --t ${cov} -o ${params.prefix}_horizontalcov -w 6000 -hi 1000
	"""
}

mosdepth_ch = Channel.fromPath(params.bams).map { file -> tuple(file.simpleName, file) }
getnsbed_ch = Channel.fromPath(params.fasta).map { file -> tuple(file.simpleName, file) }


workflow {
  MOSDEPTH(mosdepth_ch)
	GETNBED(getnsbed_ch)
	FASTCAT()
	ASSEMBLYQC(GETNBED.out.getn_bed,FASTCAT.out.fastcat_stats)
	HORIZONTALCOVERAGE(MOSDEPTH.out.per_base_cov.collect(),ASSEMBLYQC.out.barcode_cov_file)
	//MOSDEPTH.out.per_base_cov.flatten().view()
}
