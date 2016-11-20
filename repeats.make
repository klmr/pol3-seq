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

sine_index = ${sine_reference:.fa=.salmon_index}

sine_quant = $(addprefix results/salmon/sine-bare/,$(patsubst %.fq.gz,%,$(notdir ${data_files})))

.PHONY: sine-index
sine-index: ${sine_index}

${sine_index}: ${sine_reference}
	salmon index --transcripts '$<' --index '$@'

.PHONY: sine-quant
sine-quant: ${sine_quant}

results/salmon/sine-bare/%: ~/nfs/data/trna/bianca/chip/%.fq.gz ${sine_index}
	${bsub} -n8 -R'span[hosts=1]' -M12000 -R'select[mem>12000] rusage[mem=12000]' \
		"$$SHELL -c 'salmon quant --index $(lastword $^) --libType U -r <(gunzip -c $<) -o $@'"

# vim: ft=make
