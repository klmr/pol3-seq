#
# Recipes for the result files of the Pol III ChIP-seq analysis for repeat elements
#

# Helpers

keep = $(foreach i,$2,$(if $(findstring $1,$i),$i))

# Applications used

mapper = bowtie
mkindex = bowtie-build
bsub = scripts/bsub -K

folder = $(if $(realpath $1),,$1)

include paths.make

# Other parameters

memlimit = 64000
nthreads = 24

include repeats.make

# Rules to download data files

${all_annotation}:
	mkdir -p data
	curl -o $@.gz 'ftp://ftp.ensembl.org/pub/release-79/gtf/mus_musculus/Mus_musculus.GRCm38.79.gtf.gz'
	gunzip $@.gz

${reference}:
	mkdir -p data
	curl -o $@.gz 'ftp://ftp.ensembl.org/pub/release-79/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz'
	gunzip $@.gz

${trna_data}:
	mkdir -p data
	curl -o $@ http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz

${trna_annotation}: ${trna_data}
	tar xfz $< mm10-tRNAs.bed && mv mm10-tRNAs.bed $@
	./scripts/fix-reference $@

data/%.fa: data/%.bed ${reference}
	bedtools getfasta -name -fi ${reference} -bed $< -fo $@

${trna_prefix}.extended.bed: ${trna_annotation} ${genomesize}
	bedtools slop -i $< -g ${genomesize} -b 100 > $@

${trna_prefix}.flanking.bed: ${trna_annotation} ${genomesize}
	bedtools flank -i $< -g ${genomesize} -b 100 > $@

${trna_prefix}.bare.salmon_index: ${trna_reference}
	salmon index --transcripts $< --index $@

${trna_prefix}.%.salmon_index: ${trna_prefix}.%.fa
	salmon index --transcripts $< --index $@

.PRECIOUS: $(addprefix ${trna_prefix}.,$(addsuffix .fa,flanking extended))

lib=$(shell grep "$1" meta/library-files.txt)

define salmon=
	${bsub} "$$SHELL -c 'salmon quant --index $2 --libType U -r <(gunzip -c $(call lib,$1)) -o $3'"
endef

results/salmon/trna-bare/%: ${trna_prefix}.bare.salmon_index
	$(call salmon,$*,$<,$@)

results/salmon/trna-flanking/%: ${trna_prefix}.flanking.salmon_index
	$(call salmon,$*,$<,$@)

results/salmon/trna-extended/%: ${trna_prefix}.extended.salmon_index
	$(call salmon,$*,$<,$@)

# Rules to build result files

data/${genome}.genes.bed: ${all_annotation}
	awk -vOFS='\t' '!/^#/ && ($$3 == "gene") {print $$1, $$4 - 1, $$5, $$2}' $< > $@

data/${genome}.repeats.bed: ${repeat_annotation}
	awk -vOFS='\t' '!/^#/ && ($$3 != "Simple_repeat") {print $$1, $$4 - 1, $$5, $$3}' $< > $@

.PHONY: index
index: ${index}.1.ebwt

# This is inconvenient, since the index actually consists of multiple files,
# which are numbered.
${index}.1.ebwt: ${reference}
	mkdir -p ${index_path}
	${bsub} -M 16000 -R 'rusage[mem=16000]' "${mkindex} --offrate 3 $< ${index}"

.PHONY: sines_index
sines_index: ${sines_index}.1.ebwt

${sines_index}.1.ebwt: ${sines_reference}
	mkdir -p ${index_path}
	${bsub} -M 8000 -R 'rusage[mem=8000]' "${mkindex} --offrate 1 $< ${sines_index}"

#${trna_annotation}: ${repeat_annotation}
#	fgrep 'SINE/tRNA' $< > $@

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

#${coverage_path}/%.counts: ${map_path}/%.bam ${repeat_annotation} ${coverage_path}
#	${bsub} "bedtools coverage -abam $< -b ${repeat_annotation} > $@"

.PHONY: repeat-coverage
repeat-coverage: ${repeat_coverage}

.PHONY: gene-coverage
gene-coverage: ${gene_coverage}

.PHONY: trna-coverage
trna-coverage: ${trna_coverage}

${coverage_path}/%_repeats.counts: ${map_path}/%.bam data/${genome}.repeats.bed
	${bsub} -M 8000 -R 'rusage[mem=8000]' "bedtools coverage -abam $< -b data/${genome}.repeats.bed > $@"

${coverage_path}/%_genes.counts: ${map_path}/%.bam data/${genome}.genes.bed
	${bsub} -M 8000 -R 'rusage[mem=8000]' "bedtools coverage -abam $< -b data/${genome}.genes.bed > $@"

${coverage_path}/%_trnas.counts: ${map_path}/%.bam data/${genome}.trnas.bed
	${bsub} -M 8000 -R 'rusage[mem=8000]' "bedtools coverage -abam $< -b data/${genome}.trnas.bed > $@"

.PHONY: sines-coverage
sines-coverage: ${sines_coverage}

${sines_coverage_path}/%.counts: ${sines_map_path}/%.bam ${sines_reference} ${sines_coverage_path}
	${bsub} "samtools sort -n $< $(<:.bam=.sorted); \
		mkdir -p $@; \
		express --output-dir $@ ${sines_reference} $(<:.bam=.sorted.bam) > $@/log"

#.PHONY: trna-coverage
#trna-coverage: ${trna_coverage}

#${trna_coverage_path}/%.counts: ${map_path}/%.bam ${trna_annotation} ${trna_coverage_path}
#	${bsub} "bedtools coverage -abam $< -b ${trna_annotation} > $@"

# Reports

${report_path}/%.html: ${script_path}/%.rmd ${report_path}
	Rscript -e "knitr::knit2html('$<', '$@')"

.DELETE_ON_ERROR:

include make-help.make
