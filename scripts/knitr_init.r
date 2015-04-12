library(knitr)

# Adjust knitr paths to my project folder structure
root_dir = normalizePath(file.path(getwd(), '..'))
base_dir = file.path(root_dir, 'results', 'report')
opts_knit$set(root.dir = root_dir,
              base.dir = base_dir)

opts_chunk$set(dev = c('png', 'pdf'),
               cache = TRUE)
