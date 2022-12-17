nextflow.enable.dsl = 2
params.vcf = "wf-articresults/*.pass.named.vcf.gz"
params.bam = "wf-articresults/*.primertrimmed.rg.sorted.bam"
params.reference = "$projectDir/resources/Sars_cov_2.ASM985889v3.dna.toplevel.fa"
params.cov = 20
params.output = "results/GISAID_revised"
params.help = false

// Help Section

if( params.help ) {
log.info """
wf-artic-utilities@KAUST: Revise GISAID Returned assemblies post manual fintering of variants.

For filtering low quality variants and variants in homopolymer regions
=============================================
Usage:
	nextflow run reviseGisaid.nf --vcf "${params.vcf}" --bam "${params.bam}" --reference ${params.reference} --cov ${params.cov} --output ${params.output}
Input:
	* --vcf: Path to .vcf.gz files from wf-artic results. Defult [${params.vcf}]
	* --bam: Path to indexed .bam files from wf-artic results. Default [${params.bam}]
	* --reference: Path to reference fasta. Default [${params.reference}]
	* --cov: Depth of Coverage Cutoff. Default[${params.cov}X]
	* --output: Path for dumping results. Default[${params.output}]
"""

exit 0
}


process ASSEMBLY {
publishDir "${params.output}", mode: 'copy'

//conda (params.enable_conda ? "bioconda::seqkit bioconda::tabix bioconda::vafator" : null)
input:
tuple val(sid) , path(vcf), path(bam)

output:
path("*.fasta")

script:
"""
zcat ${vcf} | bgzip -c > ${sid}.gz	
tabix -p vcf ${sid}.gz

## Finding the regions of low coverage
	covtobed -x ${params.cov} ${bam.toRealPath()}  > ${sid}_lowcoverage.bed

## Building the Assembly, masking low coverage regions
	bcftools consensus -f ${params.reference}  --include 'FILTER="PASS"' -m ${sid}_lowcoverage.bed --output ${sid}.fasta ${sid}.gz

## Rename Header
	seqtk seq -C ${sid}.fasta > ${sid}.temp

## For some reason seqtk append 1 to fasta header
	seqtk rename ${sid}.temp "${sid} MN908947.3" |seqkit replace -p ".+" -r "${sid} MN908947.3"> ${sid}.fasta
	rm ${sid}.temp
"""
}



workflow{
vcf_ch = Channel.fromPath(params.vcf).map { file -> tuple(file.simpleName, file) }
bam_ch = Channel.fromPath(params.bam).map { file -> tuple(file.simpleName, file) }
all_files_ch = vcf_ch.combine(bam_ch, by: 0)
ASSEMBLY(all_files_ch)
}
