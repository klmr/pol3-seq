#
# Recipes for the result files of the Pol III ChIP-seq analysis for repeat elements
#

# Applications used

mapper = bowtie
mkindex = bowtie-build
bsub = scripts/bsub -K
format_repeat_annotation = src/gff-from-repeats

#define folder =
#	$(shell [[ -d $1 ]] || echo $1)
#endef

folder = $(if $(realpath $1),,$1)

# Filenames of data sources and result targets

genome = Mus_musculus.GRCm38.75
reference = data/${genome}.dna.primary_assembly.fa
nc_reference = data/${genome}.ncrna.fa
sines_reference = data/repbase_sine_all.fasta
index_prefix = $(notdir $(basename ${reference}))
sines_index_prefix = $(notdir $(basename ${sines_reference}))
index_path = data/${mapper}
map_path = results/${mapper}
sines_map_path = results/${mapper}/sines
bigwig_path = ${map_path}
sines_bigwig_path = ${sines_map_path}
coverage_path = ${map_path}/coverage
sines_coverage_path = ${sines_map_path}/express
trna_coverage_path = ${coverage_path}/trna
script_path = scripts
report_path = results/report
index = ${index_path}/${index_prefix}
sines_index = ${index_path}/${sines_index_prefix}
data_files = $(shell cat data/files-all.txt)
data_base = $(patsubst %/,%,$(dir $(word 1,${data_files})))
mapped_reads = $(addprefix ${map_path}/,$(patsubst %.fq.gz,%.bam,$(notdir ${data_files})))
sines_mapped = $(addprefix ${sines_map_path}/,$(notdir ${mapped_reads}))
genomesize = ${reference}.fai
sines_size = ${sines_reference}.fai
all_annotation = data/${genome}.gtf
repeat_annotation = data/${genome}.repeats.gff
trna_annotation = data/${genome}.repeats.sine.trna.gff
line_annotation = data/${genome}.repeats.line.gff
repeat_annotation_repeatmasker = data/combined_repeats.out.gz
bigwig = $(patsubst %.bam,%.bw,${mapped_reads})
sines_bigwig = $(patsubst %.bam,%.bw,${sines_mapped})
coverage = $(addprefix ${coverage_path}/,$(patsubst %.bam,%.counts,$(notdir ${mapped_reads})))
sines_coverage = $(addprefix ${sines_coverage_path}/,$(patsubst %.bam,%.counts,$(notdir ${sines_mapped})))
trna_coverage = $(addprefix ${trna_coverage_path}/, $(patsubst %.bam,%.counts,$(notdir ${mapped_reads})))

result_paths = $(sort \
	${map_path} \
	${bigwig_path} \
	${coverage_path} \
	${trna_coverage_path} \
	${report_path} \
	${sines_map_path} \
	${sines_bigwig_path} \
	${sines_coverage_path} \
)

# Other parameters

memlimit = 64000

# Rules to download data files

${all_annotation}:
	curl -o $@.gz 'ftp://ftp.ensembl.org/pub/release-79/gtf/mus_musculus/Mus_musculus.GRCm38.79.gtf.gz'
	gunzip $@.gz

${reference}:
	curl -o $@.gz 'ftp://ftp.ensembl.org/pub/release-79/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz'
	gunzip $@.gz

${nc_reference}:
	curl -o $@.gz 'ftp://ftp.ensembl.org/pub/release-79/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz'
	gunzip $@.gz

# Rules to build result files

data/${genome}.genes.bed: ${all_annotation}
	awk -vOFS='\t' '!/^#/ && ($$3=="gene") {print $$1, $$4, $$5, $$2}' $< > $@

.PHONY: index
index: ${index}.1.ebwt

# This is inconvenient, since the index actually consists of multiple files,
# which are numbered.
${index}.1.ebwt: ${reference} ${index_path}
	${bsub} -M 16000 -R 'rusage[mem=16000]' "${mkindex} --offrate 3 $< ${index}"

.PHONY: sines_index
sines_index: ${sines_index}.1.ebwt

${sines_index}.1.ebwt: ${sines_reference} ${index_path}
	${bsub} -M 8000 -R 'rusage[mem=8000]' "${mkindex} --offrate 1 $< ${sines_index}"

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

.PHONY: sines-size
sines-size: ${sines_size}

${sines_size}: ${sines_reference}
	${bsub} "samtools faidx $<"

.PHONY: mapped-reads
mapped-reads: ${mapped_reads}

${map_path}/%.bam: ${data_base}/%.fq.gz ${index}.1.ebwt ${map_path}
	${bsub} -M ${memlimit} -n 32 -R 'select[gpfs]' -R 'rusage[mem=${memlimit}]' \
		"./scripts/${mapper} --best ${index} $< $@"

.PHONY: sines-mapped
sines-mapped: ${sines_mapped}

${sines_map_path}/%.bam: ${data_base}/%.fq.gz ${sines_index}.1.ebwt ${sines_map_path}
	${bsub} -M 16000 -n 16 -R 'select[gpfs]' -R 'rusage[mem=16000]' \
		"./scripts/${mapper} --all ${sines_index} $< $@"

.PHONY: bigwig
bigwig: ${bigwig}

${bigwig_path}/%.bw: ${map_path}/%.bam ${genomesize} ${bigwig_path}
	${bsub} "./scripts/bigwig $< ${genomesize} $@"

.PHONY: sines-bigwig
sines-bigwig: ${sines_bigwig}

${sines_bigwig_path}/%.bw: ${sines_map_path}/%.bam ${sines_size} ${sines_bigwig_path}
	${bsub} "./scripts/bigwig $< ${sines_size} $@"

.PHONY: coverage
coverage: ${coverage}

${coverage_path}/%.counts: ${map_path}/%.bam ${repeat_annotation} ${coverage_path}
	${bsub} "bedtools coverage -abam $< -b ${repeat_annotation} > $@"

.PHONY: sines-coverage
sines-coverage: ${sines_coverage}

${sines_coverage_path}/%.counts: ${sines_map_path}/%.bam ${sines_reference} ${sines_coverage_path}
	${bsub} "samtools sort -n $< $(<:.bam=.sorted); \
		mkdir -p $@; \
		express --output-dir $@ ${sines_reference} $(<:.bam=.sorted.bam) > $@/log"

.PHONY: trna-coverage
trna-coverage: ${trna_coverage}

${trna_coverage_path}/%.counts: ${map_path}/%.bam ${trna_annotation} ${trna_coverage_path}
	${bsub} "bedtools coverage -abam $< -b ${trna_annotation} > $@"

# Reports

${report_path}/%.html: ${script_path}/%.rmd ${report_path}
	Rscript -e "knitr::knit2html('$<', '$@')"

${report_path}/%.md: ${script_path}/%.rmd ${report_path}
	Rscript -e "knitr::knit('$<', '$@')"

${result_paths}:
	mkdir -p $@
