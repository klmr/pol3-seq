base = import('ebits/base')
io = import('ebits/io')
import_package('dplyr', attach = TRUE)
import_package('ggplot2', attach = TRUE)
pdf.options(encoding = "CP1250")

data = io$read_table('results/salmon/merged-quant.tsv', header = TRUE)

# Replicate variability:

libraries = names(select(data, -Name, -Feature))
region_types = unique(data$Feature)

rep_scatter = lapply(region_types,
                     f -> ggplot(filter(data, Feature == f)) +
                         do.call(aes_string, as.list(libraries)) +
                         geom_point() +
                         scale_x_log10() +
                         scale_y_log10() +
                         ggtitle(f)) %>%
    setNames(region_types)

for (f in region_types) {
    pdf(sprintf('results/salmon/plots/replicate-%s.pdf', f))
    plot(rep_scatter[[f]])
    dev.off()
}

# Correlate different feature region types

region_contrasts = as.data.frame(combn(region_types, 2))

rep1 = libraries[1]
feature_scatter = lapply(region_contrasts,
       c -> {
           x = filter(data, Feature == c[1])[[rep1]]
           y = filter(data, Feature == c[2])[[rep1]]
           rsq = format(cor(x, y, method = 'spearman') ** 2, digits = 2)
           ggplot(data.frame(x = x, y = y)) +
               aes(x, y) +
               geom_point() +
               scale_x_log10(c[1]) +
               scale_y_log10(c[2]) +
               geom_smooth(method = lm, se = FALSE) +
               ggtitle(paste0(paste(c, collapse = '–'), ' (r^2 = ', rsq, ')'))
       }) %>% setNames(sapply(region_contrasts, paste, collapse = '-'))

for (c in names(feature_scatter)) {
    pdf(sprintf('results/salmon/plots/contrast-%s.pdf', c))
    plot(feature_scatter[[c]])
    dev.off()
}
