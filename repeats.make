.PHONY: repeat-annotation
repeat-annotation: ${repeat_annotation}

rm_engine = crossmatch
rm_threads = ${nthreads}
rm_memlimit = ${memlimit}

rm_repeat_annotation = ${repeat_annotation:%.repeats.gtf=%.RepeatMasker.raw}/${reference:data/%=%}.out

.INTERMEDIATE: ${rm_repeat_annotation}
${rm_repeat_annotation}: ${reference}
	mkdir -p ${@D}
	${bsub} -n${rm_threads} -M${rm_memlimit} -R'rusage[mem=${rm_memlimit}]' \
		RepeatMasker -e ${rm_engine} -pa ${rm_threads} -nolow \
			-species ${species} -dir ${@D} $<

.PRECIOUS: ${repeat_annotation}
${repeat_annotation}: ${rm_repeat_annotation}
	${bsub} "./scripts/repeatmasker-to-gtf $< $@"

.PRECIOUS: ${sine_annotation}
${sine_annotation}: ${repeat_annotation}
	fgrep 'Derives_from="SINE/' $< > $@

.PRECIOUS: ${sine_reference}
${sine_reference}: ${sine_annotation}
	${bsub} "./scripts/gtf-to-fasta $< ${reference} $@"

repeat_reference = data/${genome}.repeats.fa
.PRECIOUS: ${repeat_reference}
${repeat_reference}: ${repeat_annotation}
	${bsub} "./scripts/gtf-to-fasta $< ${reference} $@"

# Conventional read mapping

sine_index = ${sine_reference:.fa=.salmon_index}

sine_quant = $(addprefix results/salmon/sine-bare/,$(patsubst %.fq.gz,%,$(notdir ${data_files})))

sine_rna_quant = $(addprefix results/salmon/sine-rna-bare/,$(patsubst %.fq.gz,%,$(subst p1,,$(call keep,p1,$(notdir ${rna_data_files})))))

.PHONY: sine-index
sine-index: ${sine_index}

${sine_index}: ${sine_reference}
	salmon index --transcripts '$<' --index '$@'

.PHONY: sine-quant
sine-quant: ${sine_quant}

.PHONY: sine-rna-quant
sine-rna-quant: ${sine_rna_quant}

results/salmon/sine-bare/%: ~/nfs/data/trna/bianca/chip/%.fq.gz ${sine_index}
	${bsub} -n8 -R'span[hosts=1]' -M12000 -R'select[mem>12000] rusage[mem=12000]' \
		"$$SHELL -c 'salmon quant --index $(lastword $^) --libType U -r <(gunzip -c $<) -o $@'"

# Really --libType ISF?
results/salmon/sine-rna-bare/%: ~/nfs/data/trna/bianca/rna/%p1.fq.gz ${sine_index}
	${bsub} -n8 -R'span[hosts=1]' -M12000 -R'select[mem>12000] rusage[mem=12000]' \
		"$$SHELL -c 'salmon quant --index $(lastword $^) --libType A \
		-1 <(gunzip -c $<) -2 <(gunzip -c $(subst p1,p2,$<)) -o $@'"

.PHONY: sine-mapped
sine-mapped: ${sine_mapped}

${sine_map_path}/%.bam: ${data_base}/%.fq.gz ${sine_index}
	${bsub} -M ${memlimit} -n 32 -R 'select[gpfs]' -R 'rusage[mem=${memlimit}]' \
		"./scripts/${happer} --best ${index} $< $@"

# For comparison: RNA-seq against whole genome + feature counts on TEs.
#
# See sine-coverage/ subdirectory

# FIXME:

#SINE reference can only be used for flank-less alignment. Produce two additional references:
#1. SINE gene body including flanks
#2. SINE flanks, without gene body; annotate as pseudo exons in a transcript (how?!)

# vim: ft=make
