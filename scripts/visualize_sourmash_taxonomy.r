library(here)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(metacoder)
library(RColorBrewer)

classify_df <-
    read_csv(here("assay/09.annotation/00.taxonomy/01.bins_annotation/sourmash/sourmash_lca_2913_bins_k31_s10000.classify.csv")) %>%
    group_by(ID) %>%
    mutate(bin_id = str_split(ID, "/")[[1]][12]) %>%
    select(-X1, -ID)

classify_df_summary <-
    classify_df %>%
    group_by(status) %>%
    summarise(count = n())

classify_df_ra <- classify_df %>% filter(str_detect(bin_id, "RA_"))
classify_df_rah <- classify_df %>% filter(!str_detect(bin_id, "RA_"))
classify_df_tooth <- classify_df %>% filter(str_detect(bin_id, "tooth"))
classify_df_saliva <- classify_df %>% filter(str_detect(bin_id, "saliva"))

superkingdom_count <- classify_df %>% group_by(superkingdom) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
phylum_count <- classify_df %>% group_by(phylum) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
class_count <- classify_df %>% group_by(class) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
order_count <- classify_df %>% group_by(order) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
family_count <- classify_df %>% group_by(family) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
genus_count <- classify_df %>% group_by(genus) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
species_count <- classify_df %>% group_by(species) %>% summarise(count = n()) %>% mutate(per = count / sum(count))

superkingdom_count_ra <- classify_df_ra %>% group_by(superkingdom) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
phylum_count_ra <- classify_df_ra %>% group_by(phylum) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
class_count_ra <- classify_df_ra %>% group_by(class) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
order_count_ra <- classify_df_ra %>% group_by(order) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
family_count_ra <- classify_df_ra %>% group_by(family) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
genus_count_ra <- classify_df_ra %>% group_by(genus) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
species_count_ra <- classify_df_ra %>% group_by(species) %>% summarise(count = n()) %>% mutate(per = count / sum(count))

superkingdom_count_rah <- classify_df_rah %>% group_by(superkingdom) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
phylum_count_rah <- classify_df_rah %>% group_by(phylum) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
class_count_rah <- classify_df_rah %>% group_by(class) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
order_count_rah <- classify_df_rah %>% group_by(order) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
family_count_rah <- classify_df_rah %>% group_by(family) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
genus_count_rah <- classify_df_rah %>% group_by(genus) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
species_count_rah <- classify_df_rah %>% group_by(species) %>% summarise(count = n()) %>% mutate(per = count / sum(count))

superkingdom_count_tooth <- classify_df_tooth %>% group_by(superkingdom) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
phylum_count_tooth <- classify_df_tooth %>% group_by(phylum) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
class_count_tooth <- classify_df_tooth %>% group_by(class) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
order_count_tooth <- classify_df_tooth %>% group_by(order) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
family_count_tooth <- classify_df_tooth %>% group_by(family) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
genus_count_tooth <- classify_df_tooth %>% group_by(genus) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
species_count_tooth <- classify_df_tooth %>% group_by(species) %>% summarise(count = n()) %>% mutate(per = count / sum(count))

superkingdom_count_saliva <- classify_df_saliva %>% group_by(superkingdom) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
phylum_count_saliva <- classify_df_saliva %>% group_by(phylum) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
class_count_saliva <- classify_df_saliva %>% group_by(class) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
order_count_saliva <- classify_df_saliva %>% group_by(order) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
family_count_saliva <- classify_df_saliva %>% group_by(family) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
genus_count_saliva <- classify_df_saliva %>% group_by(genus) %>% summarise(count = n()) %>% mutate(per = count / sum(count))
species_count_saliva <- classify_df_saliva %>% group_by(species) %>% summarise(count = n()) %>% mutate(per = count / sum(count))

# sourmash lca classify bin_sig genban_sig
# bar plot
color_values <- brewer.pal(12, "Set3")
genus_count_t12 <-
    ggplot(genus_count %>% filter(!is.na(genus)) %>% arrange(desc(count)) %>% head(12), aes(x = "", y = count, fill = genus)) +
    geom_bar(width = 1, stat = "identity") +
    #geom_text(aes(y = count, label = count), size = 3) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = color_values) +
    labs(title = "oral 2913 bins classfication on genus level(number top 12)") +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank())
genus_count_t12
ggsave(here("assay/09.annotation/00.taxonomy/01.bins_annotation/sourmash/",
            "sourmash_lca_2913_bins_k31_s10000.classify.genus_count_t12.pdf"),
       genus_count_t12)

species_count_t12 <-
    ggplot(species_count %>% filter(!is.na(species)) %>% arrange(desc(count)) %>% head(12), aes(x = "", y = count, fill = species)) +
    geom_bar(width = 1, stat = "identity") +
    #geom_text(aes(y = count, label = count), size = 3) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = color_values) +
    labs(title = "oral 2913 bins classfication on species level(number top 12)") +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank())
species_count_t12
ggsave(here("assay/09.annotation/00.taxonomy/01.bins_annotation/sourmash/",
            "sourmash_lca_2913_bins_k31_s10000.classify.species_count_t12.pdf"),
       species_count_t12)

species_count_ra_ <-
    species_count_ra %>%
    filter(!is.na(species)) %>%
    arrange(desc(count)) %>%
    head(12) %>%
    mutate(case = "RA")
species_count_rah_ <-
    species_count_rah %>%
    filter(!is.na(species)) %>%
    arrange(desc(count)) %>%
    head(12) %>%
    mutate(case = "RAH")
species_count_ra_rah <-
    species_count_ra_ %>%
    bind_rows(species_count_rah_)

color_values <- brewer.pal(17, "Set3")
colors <- colorRampPalette(c("#377EB8", "#E41A1C", "#4DAF4A",
                            "#984EA3", "#FF7F00", "#FFFF33",
                            "#A65628","#F781BF", "#999999",
                            "#1B9E77", "#D95F02","#7570B3",
                            "#E7298A", "#66A61E", "#E6AB02",
                            "#A6761D", "#666666", "#66C2A5",
                            "#FC8D62", "#8DA0CB", "#E78AC3",
                            "#A6D854", "#FFD92F", "#E5C494",
                            "#B3B3B3", "#1A1A1A", "#67001F"))(17)
species_count_ra_rah_p <-
    ggplot(species_count_ra_rah, aes(x = case, y = count, fill = species)) +
    geom_bar(width = 1,  stat = "identity", position = position_fill()) +
    #scale_fill_brewer(palette = c("Set3"))
    scale_fill_manual(values = colors)
species_count_ra_rah_p
ggsave(here("assay/09.annotation/00.taxonomy/01.bins_annotation/sourmash/",
            "sourmash_lca_2913_bins_k31_s10000.classify.species_count_ra_rah.pdf"),
       species_count_ra_rah_p)

colors <- colorRampPalette(c("#377EB8", "#E41A1C", "#4DAF4A",
                             "#984EA3", "#FF7F00", "#FFFF33",
                             "#A65628","#F781BF", "#999999",
                             "#1B9E77", "#D95F02","#7570B3",
                             "#E7298A", "#66A61E", "#E6AB02",
                             "#A6761D", "#666666", "#66C2A5",
                             "#FC8D62", "#8DA0CB", "#E78AC3",
                             "#A6D854", "#FFD92F", "#E5C494",
                             "#B3B3B3", "#1A1A1A", "#67001F"))(22)
species_count_tooth_ <-
    species_count_tooth %>%
    filter(!is.na(species)) %>%
    arrange(desc(count)) %>%
    head(12) %>%
    mutate(case = "tooth")
species_count_saliva_ <-
    species_count_saliva %>%
    filter(!is.na(species)) %>%
    arrange(desc(count)) %>%
    head(12) %>%
    mutate(case = "saliva")
species_count_tooth_saliva <-
    species_count_tooth_ %>%
    bind_rows(species_count_saliva_)

species_count_tooth_saliva_p <-
    ggplot(species_count_tooth_saliva, aes(x = case, y = count, fill = species)) +
    geom_bar(width = 1,  stat = "identity", position = position_fill()) +
    #scale_fill_brewer(palette = c("Set3"))
    scale_fill_manual(values = colors)
species_count_tooth_saliva_p
ggsave(here("assay/09.annotation/00.taxonomy/01.bins_annotation/sourmash/",
            "sourmash_lca_2913_bins_k31_s10000.classify.species_count_tooth_saliva.pdf"),
       species_count_tooth_saliva_p)

# venn data
# http://bioinformatics.psb.ugent.be/webtools/Venn/
# write.table(classify_df_tooth$species, here("assay/09.annotation/00.taxonomy/01.bins_annotation/sourmash/species_tooth.venn.txt"), quote = FALSE, row.names = FALSE)
# write.table(classify_df_saliva$species, here("assay/09.annotation/00.taxonomy/01.bins_annotation/sourmash/species_saliva.venn.txt"), quote = FALSE, row.names = FALSE)
# write.table(classify_df_ra$species, here("assay/09.annotation/00.taxonomy/01.bins_annotation/sourmash/species_ra.venn.txt"), quote = FALSE, row.names = FALSE)
# write.table(classify_df_rah$species, here("assay/09.annotation/00.taxonomy/01.bins_annotation/sourmash/species_rah.venn.txt"), quote = FALSE, row.names = FALSE)

# metacoder analysis
classify_metacoder_df <-
    read_csv(here("assay/09.annotation/00.taxonomy/01.bins_annotation/sourmash/sourmash_lca_2913_bins_k31_s10000.classify.metacoder.csv"))
classify_tax <- parse_tax_data(classify_metacoder_df, class_cols = "lineage", class_sep = ";")
heat_tree(classify_tax, node_label = taxon_names, node_size = n_obs, node_color = n_obs)
