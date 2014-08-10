#
# Recipes for the result files of the Pol III ChIP-seq analysis for repeat elements
#

# Applications used

mapper = bowtie
mkindex = bowtie-build
bsub = ./scripts/bsub
format_repeat_annotation = src/gff-from-repeats

# Filenames of data sources and result targets

genome = Mus_musculus.GRCm38.75
reference = data/${genome}.dna.primary_assembly.fa
index_prefix = $(notdir $(basename ${reference}))
index_path = data/${mapper}
map_path = results/${mapper}
bigwig_path = ${map_path}
coverage_path = ${map_path}/coverage
index = ${index_path}/${index_prefix}
data_files = $(shell cat data/files-all.txt)
data_base = $(patsubst %/,%,$(dir $(word 1,${data_files})))
mapped_reads = $(addprefix ${map_path}/,$(patsubst %.fq.gz,%.bam,$(notdir ${data_files})))
genomesize = ${reference}.fai
annotation = data/${genome}.repeats.gff
repeat_annotation = data/${genome}.repeats.gff
repeat_annotation_repeatmasker = data/combined_repeats.out.gz
bigwig = $(patsubst %.bam,%.bw,${mapped_reads})
coverage = $(addprefix ${coverage_path}/,$(patsubst %.bam,%.counts,$(notdir ${mapped_reads})))

result_paths = $(sort ${map_path} ${bigwig_path} ${coverage_path})

# Other parameters

memlimit = 64000

# Rules to build result files

.PHONY: index
index: ${index}.1.ebwt

# This is inconvenient, since the index actually consists of multiple files,
# which are numbered.
${index}.1.ebwt: ${reference} ${index_path}
	${bsub} -M 16000 -R 'rusage[mem=16000]' "${mkindex} --offrate 3 $< ${index}"

${index_path}:
	mkdir -p ${index_path}

.PHONY: repeat_annotation
repeat_annotation: ${repeat_annotation}

${repeat_annotation}: ${repeat_annotation_repeatmasker}
	${bsub} "gunzip -c $< | ${format_repeat_annotation} > $@"

data/${genome}.repeats.sine.trna.gff: ${repeat_annotation}
	fgrep 'SINE/tRNA' $< > $@

.PHONY: genomesize
genomesize: ${genomesize}

${genomesize}: ${reference}
	${bsub} "samtools faidx $<"

.PHONY: mapped-reads
mapped-reads: ${mapped_reads}

${map_path}/%.bam: ${data_base}/%.fq.gz ${index}.1.ebwt ${map_path}
	${bsub} -M ${memlimit} -n 32 -R 'rusage[mem=${memlimit}]' \
		"./scripts/${mapper} ${index} $< $@"

.PHONY: bigwig
bigwig: ${bigwig}

${bigwig_path}/%.bw: ${map_path}/%.bam ${genomesize} ${bigwig_path}
	${bsub} "./scripts/bigwig $< ${genomesize} $@"

.PHONY: coverage
coverage: ${coverage}

${coverage_path}/%.counts: ${map_path}/%.bam ${coverage_path}
	${bsub} -n 32 "bedtools coverage -abam $< -b ${annotation} > $@"

${result_paths}:
	mkdir -p $@
