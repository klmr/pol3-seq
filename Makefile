mapper = bowtie
mkindex = bowtie-build
bsub = ./scripts/bsub

define library_for =
$(addprefix ${data_base},$(patsubst %.bam,%.fq.gz,$(notdir $1)))
endef

define bam_for =
$(addprefix results/${mapper},$(notdir $(addsuffix .bam,$(basename $1))))
endef

genome = Mus_musculus.GRCm38.75
reference = data/${genome}.dna.primary_assembly.fa
index_prefix = $(notdir $(basename ${reference}))
index_path = data/${mapper}
index = ${index_path}/${index_prefix}
data_files = $(shell cat data/files-all.txt)
data_base = $(dir $(word 1,${data_files}))
mapped_reads = $(addprefix results/${mapper}/,$(patsubst %.fq.gz,%.bam, $(notdir ${data_files})))
genomesize = data/${reference}.fai
annotation = data/${genome}.repeats.gff
memlimit = 64000

.PHONY: index
index: ${index}.1.ebwt

# This is inconvenient, since the index actually consists of multiple files,
# which are numbered. Unfortunately I don’t know how many, so I just refer to
# the first file (*.1.ebwt) here.
${index}.1.ebwt: ${reference} ${index_path}
	${bsub} -M 16000 -R 'rusage[mem=16000]' "${mkindex} --offrate 3 $< ${index}"

${index_path}:
	mkdir -p ${index_path}

.PHONY: genomesize
genomesize: ${genomesize}

${genomesize}: ${reference}
	${bsub} "samtools faidx $<"

.PHONY: mapped-reads
mapped-reads: ${mapped_reads}

${mapped_reads}: ${index}.1.ebwt results/${mapper}
	${bsub} -M ${memlimit} -n 32 -R 'rusage[mem=${memlimit}]' "./scripts/${mapper} ${index} $(call library_for,$@) $@"

.PHONY: bigwig
bigwig: ${bigwig}

${bigwig}: ${mapped_reads} ${genomesize}
	echo "./scripts/bigwig $(call bam_for,$@) ${genomesize} $@"

.PHONY: coverage
coverage: ${coverage}

${coverage}: ${mapped_reads} results/${mapper}/coverage
	echo "bedtools coverage -abam $(call bam_for,$@) -b ${annotation} > $@"

results/${mapper}:
	mkdir -p $@

results/${mapper}/coverage:
	mkdir -p $@
