.PHONY: repeat-annotation
repeat-annotation: ${repeat_annotation}

rm_engine = crossmatch
rm_threads = ${nthreads}
rm_memlimit = ${memlimit}

rm_repeat_annotation = ${repeat_annotation:%.repeats.gtf=%.RepeatMasker.raw}/${reference:data/%=%}.out

${rm_repeat_annotation}: ${reference}
	mkdir -p ${@D}
	${bsub} -n${rm_threads} -M${rm_memlimit} -R'rusage[mem=${rm_memlimit}]' \
		RepeatMasker -e ${rm_engine} -pa ${rm_threads} -nolow \
			-species ${species} -dir ${@D} $<

${repeat_annotation}: ${rm_repeat_annotation}
	${bsub} "./scripts/repeatmasker-to-gtf $< $@"

# vim: ft=make
