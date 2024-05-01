#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"])

## Convert all columns to strings
samples = samples.astype(str)

## Set sample names
samples['sn'] = samples[config['fileNamesFrom']].apply('_'.join, axis=1)

## Conditional output - check that mergeBy is valid
outfiles = list()

if config['mergeBy'] == config['fileNamesFrom'] or config['mergeBy'] == '':
	"None"
else:
	## Merge according to mergeBy parameter, define merge name (mn)
	samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)

	## Build dictionary of merged BAM files
	mergeSamples = samples.groupby('mn')['sn'].apply(list).to_dict()

	## Define rule all outputs
	outfiles.append([expand("output/mergeAlign/{mergeName}_{ext}", mergeName=key, ext=['sorted.bam', 'sorted.bam.bai', 'stats.txt']) for key in mergeSamples])
	outfiles.append([expand("output/mergeSignal/unstranded/{mergeName}.bw", mergeName=key) for key in mergeSamples])
	if config['stranded'] == True:
		outfiles.append([expand("output/mergeSignal/stranded/{mergeName}_{dir}.bw", mergeName=key, dir=['fwd', 'rev']) for key in mergeSamples])

## Define actions on success
onsuccess:
	## Success message
	print("RNApipe completed successfully! Wahoo!")

##### Define rules #####
rule all:
	input:
		outfiles

rule mergeAlign:
	input:
		bams = lambda wildcards: ["output/align/{sampleName}_sorted.bam".format(sampleName=value) for value in mergeSamples[wildcards.mergeName]],
		bais = lambda wildcards: ["output/align/{sampleName}_sorted.bam.bai".format(sampleName=value) for value in mergeSamples[wildcards.mergeName]]
	output:
		bam = "output/mergeAlign/{mergeName}_sorted.bam",
		bai = "output/mergeAlign/{mergeName}_sorted.bam.bai",
		stats = "output/mergeAlign/{mergeName}_stats.txt"
	log:
		err = 'output/logs/mergeAlign_{mergeName}.err',
		out = 'output/logs/mergeAlign_{mergeName}.out'
	benchmark: 
		'output/benchmarks/mergeAlign_{mergeName}.tsv'
	params:
		version = config['samtoolsVers']
	shell:
		"""
		module load samtools/{params.version};
		samtools merge {output.bam} {input.bams} 1>> {log.out} 2>> {log.err};
		samtools flagstat {output.bam} > {output.stats} 2>> {log.err};
		samtools index {output.bam} 1>> {log.out} 2>> {log.err}
		"""

rule mergeSignal:
	input:
		bam = rules.mergeAlign.output.bam
	output:
		"output/mergeSignal/unstranded/{mergeName}.bw"
	log:
		err = 'output/logs/mergeSignal_{mergeName}.err',
		out = 'output/logs/mergeSignal_{mergeName}.out'
	benchmark: 
		'output/benchmarks/mergeSignal_{mergeName}.tsv'
	params:
		version = config['deeptoolsVers']
	shell:
		"""
		module load deeptools/{params.version};
		bamCoverage -b {input.bam} -o {output} 1> {log.out} 2> {log.err}
		"""

rule mergeStrandedSignal:
	input:
		bam = rules.mergeAlign.output.bam
	output:
		fwd = "output/mergeSignal/stranded/{mergeName}_fwd.bw",
		rev = "output/mergeSignal/stranded/{mergeName}_rev.bw"
	log:
		err = 'output/logs/mergeStrandedSignal_{mergeName}.err',
		out = 'output/logs/mergeStrandedSignal_{mergeName}.out'
	benchmark: 
		'output/benchmarks/mergeStrandedSignal_{mergeName}.tsv'
	params:
		version = config['deeptoolsVers']
	shell:
		"""
		module load deeptools/{params.version};
		bamCoverage --filterRNAstrand forward -b {input.bam} -o {output.fwd} 1> {log.out} 2> {log.err};
		bamCoverage --filterRNAstrand reverse -b {input.bam} -o {output.rev} 1>> {log.out} 2>> {log.err}
		"""