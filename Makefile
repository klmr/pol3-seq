#
# Recipes for the result files of the Pol III ChIP-seq analysis for repeat elements
#

# Applications used

mapper = bowtie
mkindex = bowtie-build
bsub = scripts/bsub -K
format_repeat_annotation = src/gff-from-repeats

# Filenames of data sources and result targets

genome = Mus_musculus.GRCm38.75
reference = data/${genome}.dna.primary_assembly.fa
sine_reference = data/repbase_sine_all.fasta
index_prefix = $(notdir $(basename ${reference}))
sine_index_prefix = $(notdir $(basename ${sine_reference}))
index_path = data/${mapper}
map_path = results/${mapper}
sine_map_path = results/${mapper}/sines
bigwig_path = ${map_path}
coverage_path = ${map_path}/coverage
trna_coverage_path = ${coverage_path}/trna
script_path = scripts
report_path = results/report
index = ${index_path}/${index_prefix}
sine_index = ${index_path}/${sine_index_prefix}
data_files = $(shell cat data/files-all.txt)
data_base = $(patsubst %/,%,$(dir $(word 1,${data_files})))
mapped_reads = $(addprefix ${map_path}/,$(patsubst %.fq.gz,%.bam,$(notdir ${data_files})))
mapped_sines = $(addprefix ${sine_map_path}/,$(notdir ${mapped_reads}))
genomesize = ${reference}.fai
repeat_annotation = data/${genome}.repeats.gff
trna_annotation = data/${genome}.repeats.sine.trna.gff
line_annotation = data/${genome}.repeats.line.gff
repeat_annotation_repeatmasker = data/combined_repeats.out.gz
bigwig = $(patsubst %.bam,%.bw,${mapped_reads})
coverage = $(addprefix ${coverage_path}/,$(patsubst %.bam,%.counts,$(notdir ${mapped_reads})))
trna_coverage = $(addprefix ${trna_coverage_path}/, $(patsubst %.bam,%.counts,$(notdir ${mapped_reads})))

result_paths = $(sort ${map_path} ${bigwig_path} ${coverage_path} ${trna_coverage_path} ${report_path} ${sine_map_path})

# Other parameters

memlimit = 64000

# Rules to build result files

.PHONY: index
index: ${index}.1.ebwt

# This is inconvenient, since the index actually consists of multiple files,
# which are numbered.
${index}.1.ebwt: ${reference} ${index_path}
	${bsub} -M 16000 -R 'rusage[mem=16000]' "${mkindex} --offrate 3 $< ${index}"

.PHONY: sine_index
sine_index: ${sine_index}.1.ebwt

${sine_index}.1.ebwt: ${sine_reference} ${index_path}
	${bsub} -M 8000 -R 'rusage[mem=8000]' "${mkindex} --offrate 1 $< ${sine_index}"

${index_path}:
	mkdir -p ${index_path}

.PHONY: repeat_annotation
repeat_annotation: ${repeat_annotation}

${repeat_annotation}: ${repeat_annotation_repeatmasker}
	${bsub} "gunzip -c $< | ${format_repeat_annotation} > $@"

${trna_annotation}: ${repeat_annotation}
	fgrep 'SINE/tRNA' $< > $@

${line_annotation}: ${repeat_annotation}
	fgrep 'LINE' $< > $@

.PHONY: genomesize
genomesize: ${genomesize}

${genomesize}: ${reference}
	${bsub} "samtools faidx $<"

.PHONY: mapped-reads
mapped-reads: ${mapped_reads}

${map_path}/%.bam: ${data_base}/%.fq.gz ${index}.1.ebwt ${map_path}
	${bsub} -M ${memlimit} -n 32 -R 'rusage[mem=${memlimit}]' \
		"./scripts/${mapper} ${index} $< $@"

.PHONY: mapped-sines
mapped-sines: ${mapped_sines}

${sine_map_path}/%.bam: ${data_base}/%.fq.gz ${sine_index}.1.ebwt ${sine_map_path}
	${bsub} -M 8000 -n 16 -R 'rusage[mem=8000]' \
		"./scripts/${mapper} ${sine_index} $< $@"

.PHONY: bigwig
bigwig: ${bigwig}

${bigwig_path}/%.bw: ${map_path}/%.bam ${genomesize} ${bigwig_path}
	${bsub} "./scripts/bigwig $< ${genomesize} $@"

.PHONY: coverage
coverage: ${coverage}

${coverage_path}/%.counts: ${map_path}/%.bam ${repeat_annotation} ${coverage_path}
	${bsub} "bedtools coverage -abam $< -b ${repeat_annotation} > $@"

.PHONY: trna-coverage
trna-coverage: ${trna_coverage}

${trna_coverage_path}/%.counts: ${map_path}/%.bam ${trna_annotation} ${trna_coverage_path}
	${bsub} "bedtools coverage -abam $< -b ${trna_annotation} > $@"

# Reports

${report_path}/%.html: ${report_path}/%.md
	$(eval options = c('use_xhtml', 'mathjax', 'highlight_code', 'smartypants'))
	Rscript --vanilla -e "options(markdown.HTML.options = NULL); markdown::markdownToHTML('$<', '$@', options = ${options})"

${report_path}/%.md: ${script_path}/%.rmd ${report_path}
	Rscript --vanilla -e "knitr::knit('$<', '$@')"

${result_paths}:
	mkdir -p $@
