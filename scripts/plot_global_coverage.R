library(tidyverse)

coverage_df <- read_tsv(snakemake@input[[1]],
                        col_names = c("coverage", "count"),
                        col_types = "dd")

maxdepth <- as.integer(snakemake@wildcards[["maxdepth"]])

plot_df <- coverage_df %>% filter(coverage <= maxdepth)
plot_df <- bind_rows(plot_df,
                     coverage_df %>%
                         filter(coverage > maxdepth) %>%
                         summarise(coverage = maxdepth+1,
                                   count = sum(count)))

ggplot(plot_df, aes(x = coverage, y = count)) +
    geom_col(width = 1) +
    labs(x = "Coverage", y = "Number of bases") +
    scale_y_continuous(labels = function(x) format(x / 1e9, big.mark = ","))

ggsave(snakemake@output[[1]], width = 5, height = 4)
