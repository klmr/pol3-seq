files = file.path(dir('results/salmon/sine-bare', full.names = TRUE), 'quant.sf')

f = modules::import('klmr/functional')
modules::import('klmr/functional/lambda')
library(dplyr)
library(tidyr)

coldata = tibble(File = files) %>%
    extract(File, c('ChIP', 'Tissue', 'Stage'),
            'do\\d+_([[:alpha:]]+)_([[:alpha:]]+)_(?:[^_]*)_[Mm]mus?BL6([^_]+)',
            remove = FALSE) %>%
    filter(Tissue == 'liver') %>%
    mutate(ChIP = ifelse(ChIP == 'Input', 'Input', 'PolIII')) %>%
    mutate(Stage = toupper(Stage)) %>%
    filter(Stage %in% c('E15.5', 'P22'))

row_list = function (df)
    lapply(seq_len(nrow(df)), i -> df[i, ])

data = coldata$File %>%
    lapply(readr::read_tsv) %>%
    Map(x ~ f -> mutate(x, ChIP = f$ChIP, Tissue = f$Tissue, Stage = f$Stage),
        ., row_list(coldata)) %>%
    bind_rows() %>%
    group_by(Tissue, Stage, ChIP) %>%
    mutate(ID = seq_len(n())) %>%
    group_by(ID, Name, Tissue, Stage, ChIP) %>%
    summarize(TPM = mean(TPM)) %>%
    spread(ChIP, TPM) %>%
    group_by(Tissue, Stage) %>%
    mutate(MinPosInput = min(Input[Input != 0]),
           Value = PolIII / (Input + MinPosInput))

library(ggplot2)

theme_set(theme_minimal())

summary = data %>%
    group_by(Tissue, Stage, Name) %>%
    summarize(Value = sum(PolIII, na.rm = TRUE) / sum(Input)) %>%
    mutate(Name = reorder(Name, -Value))

ggplot(summary, aes(Name, Value)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    facet_wrap(~ Stage) +
    labs(x = 'SINE class', y = 'Ratio signal / input')
