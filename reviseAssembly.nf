nextflow.enable.dsl = 2
params.vcf = "wf-articresults/*.pass.named.vcf.gz"
params.bam = "wf-articresults/*.primertrimmed.rg.sorted.bam"
params.reference = "$projectDir/resources/Sars_cov_2.ASM985889v3.dna.toplevel.fa"
params.cov = 20
params.output = "results/new_assembly"
params.gtf = "$projectDir/resources/Sars_cov_2.ASM985889v3.101.gff3"
params.homopolymer = "$projectDir/resources/homopolymer.csv"
params.help = false
params.skipHomopolymerMasking = true

// Help Section

if( params.help ) {
log.info """
wf-artic-utilities@KAUST: Revise Wf-Artic Assembly

For filtering low quality variants and variants in homopolymer regions
=============================================
Usage:
	nextflow run reviseAssembly.nf --vcf "${params.vcf}" --bam "${params.bam}" --reference ${params.reference} --cov ${params.cov} --output ${params.output}\
    --gtf ${params.gtf} --homopolymer ${params.homopolymer}
Input:
	* --vcf: Path to .vcf.gz files from wf-artic results. Defult [${params.vcf}]
	* --bam: Path to indexed .bam files from wf-artic results. Default [${params.bam}]
	* --reference: Path to reference fasta. Default [${params.reference}]
	* --cov: Depth of Coverage Cutoff. Default[${params.cov}X]
	* --output: Path for dumping results. Default[${params.output}]
  	* --gtf: Path to GTF/GFF file. Default[${params.gtf}]
  	* --homopolymer: Path to Homopolymer regions (>=4bp) BED file. Default[${params.homopolymer}]
	* --skipHomopolymerMasking: If true, variants falling within homopolymer regions will not be filtered. Default [${params.skipHomopolymerMasking}]
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
path("${sid}.vaf.annot.vcf*")

script:
if ("${params.skipHomopolymerMasking}")
"""

## Decompose MNVs

bcftools sort ${vcf.toRealPath()} |  bcftools norm \
--multiallelics -any --check-ref e --fasta-ref ${params.reference} --old-rec-tag OLD_CLUMPED --atomize - | \
bcftools norm --rm-dup exact --output-type z -o ${sid}.normalized.vcf.gz -

## Calculating the Variant Allele Frequency (VAF)
  vafator --input-vcf ${sid}.normalized.vcf.gz --output-vcf ${sid}.vaf.vcf --bam vafator ${bam.toRealPath()} --mapping-quality 0 --base-call-quality 0
  bgzip -c ${sid}.vaf.vcf > ${sid}.vaf.vcf.gz
  tabix -p vcf ${sid}.vaf.vcf.gz

## Bad variant calls Tagging
  bcftools view -Ob ${sid}.vaf.vcf.gz  | \
  bcftools filter --exclude 'INFO/vafator_af < 0.5 || INFO/vafator_dp < 100 && INFO/vafator_ac < 50 ' \
  --soft-filter POOR_CALLS --output-type v - | bcftools norm -d both - > ${sid}.vaf.annot.vcf

  bgzip -c  ${sid}.vaf.annot.vcf >  ${sid}.vaf.annot.vcf.gz
  tabix -p vcf ${sid}.vaf.annot.vcf.gz

## Finding the regions of low coverage
  covtobed -x ${params.cov} ${bam.toRealPath()}  > ${sid}_lowcoverage.bed

## Building the Assembly, masking low coverage regions
  bcftools consensus -f ${params.reference}  --include 'FILTER="PASS"' -m ${sid}_lowcoverage.bed --output ${sid}.fasta ${sid}.vaf.annot.vcf.gz

## Rename Header
	seqtk seq -C ${sid}.fasta > ${sid}.temp

## For some reason seqtk append 1 to fasta header
	seqtk rename ${sid}.temp "${sid} MN908947.3" |seqkit replace -p ".+" -r "${sid} MN908947.3"> ${sid}.fasta
	rm ${sid}.temp

"""
else

"""
bcftools sort ${vcf.toRealPath()} |  bcftools norm \
--multiallelics -any --check-ref e --fasta-ref ${params.reference} --old-rec-tag OLD_CLUMPED --atomize - | \
bcftools norm --rm-dup exact --output-type z -o ${sid}.normalized.vcf.gz -

## Calculating the Variant Allele Frequency (VAF)
  vafator --input-vcf ${sid}.normalized.vcf.gz --output-vcf ${sid}.vaf.vcf --bam vafator ${bam.toRealPath()} --mapping-quality 0 --base-call-quality 0

## Filtering the Variants from Homopolymer regions (regions with repeates longer than 3bp)
  bgzip -c ${sid}.vaf.vcf > ${sid}.vaf.vcf.gz
  tabix -p vcf ${sid}.vaf.vcf.gz
  tabix -R ${params.homopolymer} ${sid}.vaf.vcf.gz -h > ${sid}.vaf.sanatized.vcf

## Bad variant calls Tagging
  bcftools view -Ob ${sid}.vaf.sanatized.vcf  | \
  bcftools filter --exclude 'INFO/vafator_af < 0.5 || INFO/vafator_dp < 100 && INFO/vafator_ac < 50 ' \
  --soft-filter POOR_CALLS --output-type v - | bcftools norm -d both - > ${sid}.vaf.annot.vcf

  bgzip -c  ${sid}.vaf.annot.vcf >  ${sid}.vaf.annot.vcf.gz
  tabix -p vcf ${sid}.vaf.annot.vcf.gz

## Finding the regions of low coverage
  covtobed -x ${params.cov} ${bam.toRealPath()}  > ${sid}_lowcoverage.bed

## Building the Assembly, masking low coverage regions
  bcftools consensus -f ${params.reference}  --include 'FILTER="PASS"' -m ${sid}_lowcoverage.bed --output ${sid}.fasta ${sid}.vaf.annot.vcf.gz

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
