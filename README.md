# Genome-wide coverage analysis of pol III ChIP marks

## Usage

Look at the [`Makefile`][makefile] for possible targets. To get specific feature
coverage, create the `*.count` files in a given coverage output folder.

For example, to investigate tRNA-related SINE repeat coverage, try

```bash
output=results/bowtie/coverage/trna
targets=($(< data/files-all.txt xargs basename --multiple \
    | sed "s#\(.*\)\.fq\.gz\$#$output/\1.counts#"
))
make -j ${targets[@]}
```

Most targets also have `.PHONY` shortcut target names, so the above could also
be written as

```bash
make -j trna-coverage
```

[makefile]: blob/master/Makefile
