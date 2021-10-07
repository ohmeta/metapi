#!/usr/bin/env Rscript

library(ggplot2)


dada2_stats_barplot <- function(df, stack=FALSE, pretty=FALSE)
{
    df <- df %>% dplyr::arrange(`non-chimeric`)
    df_l <- df %>% 
        dplyr::select("sample-id", "input", "filtered", "denoised", "non-chimeric") %>%
        tidyr::pivot_longer(!"sample-id", names_to="step", values_to="count") %>%
        dplyr::mutate(step=factor(step,
                                  levels=c("input", "filtered", "denoised", "non-chimeric")),
                      `sample-id`=factor(`sample-id`,
                                         levels=df$`sample-id`))

    position = position_dodge(0.8)
    if (stack) { position = "stack" }

    if (pretty) {
        p <- 
            ggpubr::ggbarplot(df_l, x="sample-id", y="count",
                              fill="step", color="step", x.text.angle=90,
                              stat="identity", position=position)
    } else {
        p <-
            ggplot(df_l, aes(x=`sample-id`, y=count)) +
            geom_bar(aes(color=step, fill=step),
                     stat="identity", position=position, width=0.7) +
        theme_classic() +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5,
                                       size=12, color="black"),
              axis.text.y=element_text(size=12, color="black"))
        }

    print(p)
    return(p)
}