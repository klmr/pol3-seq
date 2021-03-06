f = modules::import('klmr/functional')
modules::import('klmr/functional/lambda')
modules::import_package('dplyr', attach = TRUE)
modules::import_package('tidyr', attach = TRUE)

chip_files = file.path(dir('results/salmon/sine-bare', full.names = TRUE), 'quant.sf')

chip_coldata = tibble(File = chip_files) %>%
    extract(File, c('ChIP', 'Tissue', 'Stage'),
            'do\\d+_([[:alpha:]]+)_([[:alpha:]]+)_(?:[^_]*)_[Mm]mus?BL6([^_]+)',
            remove = FALSE) %>%
    filter(Tissue == 'liver') %>%
    mutate(ChIP = ifelse(ChIP == 'Input', 'Input', 'PolIII')) %>%
    mutate(Stage = toupper(Stage)) %>%
    filter(Stage %in% c('E15.5', 'P22'))

row_list = function (df)
    lapply(seq_len(nrow(df)), i -> df[i, ])

chip_data = chip_coldata$File %>%
    lapply(readr::read_tsv) %>%
    Map(x ~ f -> mutate(x, ChIP = f$ChIP, Tissue = f$Tissue, Stage = f$Stage),
        ., row_list(chip_coldata)) %>%
    bind_rows() %>%
    group_by(Tissue, Stage, ChIP) %>%
    mutate(ID = seq_len(n())) %>%
    group_by(ID, Name, Tissue, Stage, ChIP) %>%
    summarize(TPM = mean(TPM)) %>%
    spread(ChIP, TPM) %>%
    group_by(Tissue, Stage) %>%
    mutate(MinPosInput = min(Input[Input != 0]),
           Value = PolIII / (Input + MinPosInput))

rna_files = file.path(dir('results/salmon/sine-rna-bare', full.names = TRUE), 'quant.sf')

rna_coldata = tibble(File = rna_files) %>%
    extract(File, c('Tissue', 'Stage'), 'RNAseq_([[:alpha:]]+)_mmuBL6([^_]+)',
            remove = FALSE) %>%
    mutate(Stage = toupper(ifelse(Stage == 'Greece', 'e15.5', Stage))) %>%
    filter(Tissue == 'liver', Stage %in% c('E15.5', 'P22'))

rna_data = rna_coldata$File %>%
    lapply(readr::read_tsv) %>%
    Map(x ~ f -> mutate(x, Tissue = f$Tissue, Stage = f$Stage), .,
        row_list(rna_coldata)) %>%
    bind_rows() %>%
    group_by(Tissue, Stage) %>%
    mutate(ID = seq_len(n())) %>%
    group_by(ID, Name, Tissue, Stage) %>%
    summarize(Value = mean(TPM))

modules::import_package('ggplot2', attach = TRUE)

theme_set(theme_minimal() +
          theme(panel.grid.minor = element_blank()))

chip_data %>%
    filter(Stage == 'E15.5') %>%
    group_by(Name) %>%
    summarize(Copies = n()) %>%
    mutate(Name = reorder(Name, -Copies)) %>%
    ggplot(aes(Name, Copies)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    theme(panel.grid.major.y = element_blank())

ggsave('results/salmon/plots/sine-element-numbers-barchart.pdf')

chip_summary = chip_data %>%
    group_by(Tissue, Stage, Name) %>%
    summarize(Value = sum(PolIII, na.rm = TRUE) / sum(Input)) %>%
    mutate(Name = reorder(Name, -Value))

ggplot(chip_summary, aes(Name, Value, fill = Stage)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    labs(x = 'SINE class', y = 'Signal/input') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major.x = element_blank()) +
    scale_fill_manual(values = c('skyblue', 'orange'))

ggsave('results/salmon/plots/chip-sine-expression-barchart.pdf',
       width = 12, height = 7)

rna_summary = rna_data %>%
    group_by(Tissue, Stage, Name) %>%
    summarize(Value = sum(Value)) %>%
    mutate(Name = reorder(Name, -Value))

ggplot(rna_summary, aes(Name, Value, fill = Stage)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    labs(x = 'SINE class', y = 'Transcripts per million') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major.x = element_blank()) +
    scale_fill_manual(values = c('skyblue', 'orange'))

ggsave('results/salmon/plots/rna-sine-expression-barchart.pdf',
       width = 12, height = 7)

summary = inner_join(chip_summary, rna_summary, by = c('Name', 'Tissue', 'Stage'))

models = summary %>%
    split(.$Stage) %>%
    purrr::map(~ lm(Value.y ~ Value.x, data = .)) %>%
    purrr::map(base::summary)

plot_correlation = function (stage) {
    data = summary %>% filter(Stage == stage)
    ggplot(data, aes(Value.x, Value.y, color = Stage)) +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    scale_y_log10() +
    labs(x = 'Pol III ChIP­seq', y = 'RNA­seq') +
    annotate('text', x = min(data$Value.x), y = 1000, hjust = 0,
             label = sprintf('R^2 == %0.2f', models[[stage]]$r.squared),
             parse = TRUE) +
    scale_color_manual(values = c('skyblue', 'orange'),
                       limits = c('E15.5', 'P22'),
                       guide = FALSE)
}

plot_correlation('E15.5')
ggsave('results/salmon/plots/sine-rna-chip-correlation-E15.5-scatter.pdf', width = 5, height = 5)
plot_correlation('P22')
ggsave('results/salmon/plots/sine-rna-chip-correlation-P22-scatter.pdf', width = 5, height = 5)
