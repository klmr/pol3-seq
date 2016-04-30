# Filenames of data sources and result targets

species = mouse
genome = Mus_musculus.GRCm38.75
reference = data/${genome}.dna.primary_assembly.fa
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
data_files := $(shell cat meta/library-files.txt)
data_base = $(patsubst %/,%,$(dir $(word 1,${data_files})))
mapped_reads = $(addprefix ${map_path}/,$(patsubst %.fq.gz,%.bam,$(notdir ${data_files})))
sines_mapped = $(addprefix ${sines_map_path}/,$(notdir ${mapped_reads}))
genomesize = ${reference}.fai
sines_size = ${sines_reference}.fai
all_annotation = data/${genome}.gtf
repeat_annotation = data/${genome}.repeats.gtf
sine_annotation = data/${genome}.sines.gtf
sine_reference = data/${genome}.sines.fa
trna_annotation = data/${genome}.repeats.sine.trna.gff
line_annotation = data/${genome}.repeats.line.gff
repeat_annotation_repeatmasker = data/combined_repeats.out.gz
bigwig = $(patsubst %.bam,%.bw,${mapped_reads})
sines_bigwig = $(patsubst %.bam,%.bw,${sines_mapped})
coverage = $(addprefix ${coverage_path}/,$(patsubst %.bam,%.counts,$(notdir ${mapped_reads})))
sines_coverage = $(addprefix ${sines_coverage_path}/,$(patsubst %.bam,%.counts,$(notdir ${sines_mapped})))
#trna_coverage = $(addprefix ${trna_coverage_path}/, $(patsubst %.bam,%.counts,$(notdir ${mapped_reads})))
trna_data = data/mm10-tRNAs.tar.gz
trna_prefix = data/${genome}.trna
trna_annotation = ${trna_prefix}.bed
trna_reference = ${trna_prefix}.fa

repeat_coverage = $(addprefix ${coverage_path}/,$(patsubst %.bam,%_repeats.counts,$(notdir ${mapped_reads})))
gene_coverage = $(addprefix ${coverage_path}/,$(patsubst %.bam,%_genes.counts,$(notdir ${mapped_reads})))
trna_coverage = $(addprefix ${coverage_path}/,$(patsubst %.bam,%_trnas.counts,$(notdir ${mapped_reads})))

# vim: ft=make
