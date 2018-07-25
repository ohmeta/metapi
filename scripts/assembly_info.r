#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(argparser)
library(here)

parse_asm <- function(path_f)
{
    return
    read_delim(path_f, delim = '\t') %>%
      arrange(scaf_L50) %>%
      select(
        filename, n_scaffolds, scaf_bp,
        scaf_N50, scaf_L50,
        scaf_N90, scaf_L90,
        scaf_max, scaf_n_gt50K, scaf_pct_gt50K,
        gc_avg, gc_std)
}

average_asm <- function(asm_df)
{
    return
    asm_df %>%
        select(
          n_scaffolds, scaf_bp,
          scaf_N50, scaf_L50,
          scaf_N90, scaf_L90, scaf_max,
          scaf_n_gt50K, scaf_pct_gt50K,
          gc_avg, gc_std) %>%
        summarise(
            n_scaffolds_average = mean(n_scaffolds),
            scaf_bp_average = mean(scaf_bp),
            scaf_N50_average = mean(scaf_N50),
            scaf_L50_average = mean(scaf_L50),
            scaf_N90_average = mean(scaf_N90),
            scaf_L90_average = mean(scaf_L90),
            scaf_max_average = mean(scaf_max),
            scaf_n_gt50K_average = mean(scaf_n_gt50K),
            scaf_pct_gt50K_average = mean(scaf_pct_gt50K),
            gc_avg_average = mean(gc_avg),
            gc_std_average = mean(gc_std)) %>%
        gather(key, value) %>%
        mutate(value_human = value / 1000)
}

asm_boxplot <- function(df, title)
{
    p <-
    df %>%
        gather(key, value, -filename) %>%
        mutate(key = factor(
            key,
            levels = c(
                "n_scaffolds", "scaf_bp",
                "scaf_N50", "scaf_L50",
                "scaf_N90", "scaf_L90",
                "scaf_max", "scaf_n_gt50K", "scaf_pct_gt50K",
                "gc_avg", "gc_std"))) %>%
        ggplot(., aes(key, value)) +
        geom_boxplot(aes(fill = key), outlier.size = 0.5) +
        geom_jitter(size = 1, width = 0.25) +
        facet_wrap(~ key, scales = "free") +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            legend.title = element_blank()) +
        ggtitle(title)
    return(p)
}

parser <- arg_parser("plot assembly statistics") %>%
  add_argument("--assembly_info", help="assembly statistics info table") %>%
  add_argument("--pdf", help="assembly statistics plot", default="assembly_statistics.pdf")

args <- parse_args(parser)
asm_df <- parse_asm(args$assembly_info)
average_asm_df <- average_asm(asm_df)
plot <- asm_boxplot(asm_df, "8 soil and 2 wood samples megahit assembly statistics")
ggsave(args$pdf, plot, width = 10, height = 10)
