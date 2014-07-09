mapper=bowtie2
mkindex=bowtie2-build
bsub=./scripts/bsub

reference=data/Mus_musculus.GRCm38.75.dna.primary_assembly.fa
index_prefix=$(notdir $(basename ${reference}))
index_path=data/${mapper}
index=${index_path}/${index_prefix}

index: ${index}.1.bt2

# This is inconvenient, since the index actually consists of multiple files,
# which are numbered. Unfortunately, I don’t know how many.
${index}.1.bt2: ${reference} ${index_path}
	${bsub} -M 16000 -R 'rusage[mem=16000]' "${mkindex} $< ${index}"

${index_path}:
	mkdir -p ${index_path}

mapped-reads: ${mapped_reads}

${mapped_reads}: ${index}.1.bt2
	bsub -M 16000 -R 'rusage[mem=16000]' "./scripts/${mapper}"
