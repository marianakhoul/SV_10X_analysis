"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml R/3.6.2-foss-2019b-fh1
ml Python/3.7.4-foss-2019b-fh1
ml SAMtools/1.10-GCCcore-8.3.0

snakemake -s combineSvabaTitan.snakefile --latency-wait 60 --restart-times 2 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 50 -np
"""

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

import glob
import re
def getLRFullPath(base, filename):
  return glob.glob(''.join([base, filename]))
  
def getTITANpath(base, id, ext):
  return glob.glob(''.join([base, "optimalClusterSolution/", id, "_cluster*", ext]))


CHRS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X']


rule all:
  input: 
  	expand("results/LongRangerSomaticSV/{tumor}/{tumor}.LR.somatic.sv.txt", tumor=config["pairings"]),
  	expand("results/LongRangerSomaticSV/{tumor}/{tumor}.LR.germline.sv.txt", tumor=config["pairings"]),
  	"results/panelOfNormalsSV/PanelOfNormalsSV.txt",
	"results/panelOfNormalsSV/PoNBlacklistBins.txt",
	expand("results/barcodeRescue/{tumor}.bxOverlap.vcf", tumor=config["pairings"]),
  	expand("results/combineSvabaTitan/{tumor}/{tumor}.svabaTitan.sv.txt", tumor=config["pairings"]),
  	expand("results/combineSvabaTitan/{tumor}/{tumor}.svabaTitan.cn.txt", tumor=config["pairings"]),
  	expand("results/combineSvabaTitan/{tumor}/{tumor}.svabaTitan.sv.bedpe", tumor=config["pairings"]),
 	expand("results/combineSvabaTitan/{tumor}/{tumor}.svabaTitan.sv.annotPoN.bedpe", tumor=config["pairings"]),
	expand("results/combineSvabaTitan/{tumor}/{tumor}.svabaTitan.sv.PoNToolFilter.bedpe", tumor=config["pairings"]),
  	#expand("results/plotSvabaTitan/{tumor}/{tumor}_CNA-SV-BX_{type}_chr{chr}.{format}", tumor=config["pairings"], type=config["plot_type"], chr=CHRS, format=config["plot_format"]),
   	#expand("results/plotCircos/{tumor}/{tumor}_Circos.pdf", tumor=config["pairings"])
 		
rule getLongRangerSomaticSV:
	input:
		tumSVFile=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], "large_sv_calls.bedpe"),
		normSVFile=lambda wildcards: getLRFullPath(config["samples"][config["pairings"][wildcards.tumor]], "large_sv_calls.bedpe"),
		tumDelFile=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], "dels.vcf.gz"),
		normDelFile=lambda wildcards: getLRFullPath(config["samples"][config["pairings"][wildcards.tumor]], "dels.vcf.gz")
			
	output:
		outputSVFile="results/LongRangerSomaticSV/{tumor}/{tumor}.LR.somatic.sv.txt",
		outputNormSVFile="results/LongRangerSomaticSV/{tumor}/{tumor}.LR.germline.sv.txt",
	params:
		getLRscript=config["getLRsomaticSV_script"],		
		tenXfuncs=config["tenX_funcs"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		chrs=config["chrs"]
	log:
		"logs/LongRangerSomaticSV/{tumor}.log"
	shell:
		"Rscript {params.getLRscript} --id {wildcards.tumor} --tenX_funcs {params.tenXfuncs} --tumLargeSVFile {input.tumSVFile} --normLargeSVFile {input.normSVFile} --tumDeletionFile {input.tumDelFile} --normDeletionFile {input.normDelFile} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --chrs \"{params.chrs}\" --outDir results/LongRangerSomaticSV/{wildcards.tumor}/ --outputSVFile {output.outputSVFile} --outputNormSVFile {output.outputNormSVFile} > {log} 2> {log}"
		
rule buildPoN:
	input:
		svabaDir=expand("./results/svaba/{tumor}", tumor=config["pairings"]),
		lrDir=expand("./results/LongRangerSomaticSV/{tumor}", tumor=config["pairings"])
	output:
		outputPoNFile="results/panelOfNormalsSV/PanelOfNormalsSV.txt",
		outputBlackListFile="results/panelOfNormalsSV/PoNBlacklistBins.txt"
	params:
		buildPoNscript=config["buildPoN_script"],
		blackListBinWidth=config["PoN_blackListBinWidth"],
		svabafuncs=config["svaba_funcs"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		chrs=config["chrs"]
	log:
		"logs/panelOfNormalsSV/panelOfNormalsSV.log"
	shell:
		"Rscript {params.buildPoNscript} --SVABAdir {input.svabaDir} --LRdir {input.lrDir} --svaba_funcs {params.svabafuncs} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --chrs \"{params.chrs}\" --outputPoNFile {output.outputPoNFile} --outputBlackListFile {output.outputBlackListFile} > {log} 2> {log}"

rule barcodeRescue:
	input:
		tumBam=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], config["bamFileName"]),
		unfiltVCF="results/svaba/{tumor}/{tumor}.svaba.unfiltered.somatic.sv.vcf",
		bps="results/svaba/{tumor}/{tumor}.bps.txt.gz"
	output:
		"results/barcodeRescue/{tumor}.bxOverlap.vcf"
	params:
		bxRescueScript=config["bxRescue_script"],
		id="{tumor}",
		tenXfuncs=config["tenX_funcs"],
		svabafuncs=config["svaba_funcs"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		chrs=config["chrs"],
		minMapQ=config["bxRescue_minMapQ"],
		minLength=config["bxRescue_minLength"],
		windowSize=config["bxRescue_windowSize"],
		minRead=config["bxRescue_minReadOverlapSupport"]	
	log:
		"logs/barcodeRescue/{tumor}.bxOverlap.log"
	shell:
		"Rscript {params.bxRescueScript} --tenX_funcs {params.tenXfuncs} --svaba_funcs {params.svabafuncs} --id {params.id} --tumBam {input.tumBam} --vcf {input.unfiltVCF} --bps {input.bps} --chrs \"{params.chrs}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --minMapQ {params.minMapQ} --minLength {params.minLength} --windowSize {params.windowSize} --minReadOverlapSupport {params.minRead} --outFile {output} > {log} 2> {log}"

rule combineSvabaTitan:
	input:
		LRsummaryFile=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], "summary.csv"),
		svabaVCF="results/barcodeRescue/{tumor}.bxOverlap.vcf",
		titanBinFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.cna.txt"),
		titanSegFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.seg.noSNPs.txt"),
		LRsvFile="results/LongRangerSomaticSV/{tumor}/{tumor}.LR.somatic.sv.txt"
	output:
		outputSVFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.txt",
		outputBedpeFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.bedpe",
		outputCNFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.cn.txt"
	params:
		combineSVCNscript=config["combineSVCN_script"],
		normID=lambda wildcards: config["pairings"][wildcards.tumor],
		tenXfuncs=config["tenX_funcs"],
		svabafuncs=config["svaba_funcs"],
		#manualSVfile=config["manualSVFile"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		chrs=config["chrs"],
		minMapQ=config["bxRescue_minMapQ"],
		minLength=config["bxRescue_minLength"],
		windowSize=config["bxRescue_windowSize"],
		minRead=config["bxRescue_minReadOverlapSupport"]	
	log:
		"logs/combineSvabaGrocsvsTitan/{tumor}.log"
	shell:
		"Rscript {params.combineSVCNscript} --tumID {wildcards.tumor} --normID {params.normID} --tenX_funcs {params.tenXfuncs} --svaba_funcs {params.svabafuncs} --svabaVCF {input.svabaVCF} --manualSVFile {params.manualSVfile} --titanBinFile {input.titanBinFile} --titanSegFile {input.titanSegFile} --LRsummaryFile {input.LRsummaryFile} --LRsvFile {input.LRsvFile} --grocsvsFile {input.grocFile} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --chrs \"{params.chrs}\" --outDir results/combineSvabaGrocsvsTitan/{wildcards.tumor}/ --outputSVFile {output.outputSVFile} --outputCNFile {output.outputCNFile} --outputBedpeFile {output.outputBedpeFile} > {log} 2> {log}"
rule annotatePoNSV:
	input:
		svFile="results/combineSvabaTitan/{tumor}/{tumor}.svabaTitan.sv.txt",
		PoNFile="results/panelOfNormalsSV/PanelOfNormalsSV.txt",
		blackListFile="results/panelOfNormalsSV/PoNBlacklistBins.txt"
	output:
		outputSVAnnotFile="results/combineSvabaTitan/{tumor}/{tumor}.svabaTitan.sv.annotPoN.bedpe",
	params:
		annotScript=config["annotPoNSV_script"],
		svabafuncs=config["svaba_funcs"],
	log:
		"logs/combineSvabaTitan/{tumor}.annotPoNSV.log"
	shell:
		"Rscript {params.annotScript} --id {wildcards.tumor} --svaba_funcs {params.svabafuncs} --svFile {input.svFile} --PoNFile {input.PoNFile} --blackListFile {input.blackListFile} --outputSVAnnotFile {output.outputSVAnnotFile} 2> {log} > {log}"

rule filterSVs:
	input:
		svFile="results/combineSvabaTitan/{tumor}/{tumor}.svabaTitan.sv.annotPoN.bedpe"
	output:
		outputSVFiltFile="results/combineSvabaTitan/{tumor}/{tumor}.svabaTitan.sv.PoNToolFilter.bedpe",
		outputSummaryFile="results/combineSvabaTitan/{tumor}/{tumor}.svabaTitan.sv.PoNToolFilter.summary.txt"
	params:
		filterScript=config["filterSVs_script"],
		minFreqPoNSVBkptOverlap=config["PoN_minFreqSVbkpts"],
		# minFreqPoNCNVBkptOverlap=config["PoN_minFreqCNV"],
		minFreqPoNBlackList=config["PoN_minFreqBlackList"]
	log:
		"logs/combineSvabaTitan/{tumor}.filterSVs.log"
	shell:
		"Rscript {params.filterScript} --id {wildcards.tumor} --svFile {input.svFile} --minFreqPoNSVBkptOverlap {params.minFreqPoNSVBkptOverlap} --minFreqPoNBlackList {params.minFreqPoNBlackList} --outputSVFile {output.outputSVFiltFile} --outputSummary {output.outputSummaryFile} 2> {log} > {log}"


