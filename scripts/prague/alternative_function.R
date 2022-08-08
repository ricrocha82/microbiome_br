# For this function to work, it needs as inputs the taxa (species or strain level) of interest and a dataframe containing
# columns corresponding to the taxonomy, the size of the genome im megabases (mb), the 
# GC content percentage (gc_percent) and the level of the assembly, while the lines will correspond
# to the species names as in the NCBI genome database

## The function
genomic_traits <- function(taxa, df){
  
  Assembly_level <- c("Complete Genome","Chromosome","Scaffold","Contig")
  number_of_entries <- df %>% count(taxonomy)
  
  number_of_entries_in_genomic_db <- number_of_entries %>% 
    filter(taxonomy == taxa) %>% 
    pull(n)
  
  size_gc_status <- df %>% 
    filter(taxonomy == taxa) %>% 
    select(size_mb,gc_percent,status, number_organism_name) %>%
    mutate(number_of_entries_in_NCBI = as.numeric(number_of_entries_in_genomic_db))
  
  Higher_status <- Assembly_level[min(match(size_gc_status$status,Assembly_level))]
  
  used_genomes <- size_gc_status %>% 
    filter(status == Higher_status)
  number_of_genomes_after_filter <- nrow(used_genomes)
  
  df_summ_entries <- used_genomes %>% 
    select(size_mb, gc_percent) %>% 
    summarize(across(everything(), ~summary(.) %>% as_tibble_row() %>% list())) %>% 
    unnest(cols = c(size_mb, gc_percent), names_sep = "_") %>%
    as_tibble(.name_repair = janitor::make_clean_names) %>%
    setNames(str_replace(names(.), "qu", "qtl"))
  
  genome_stats_table <- used_genomes %>% 
    select(-status) %>%
    mutate(df_summ_entries) %>%
    mutate(highest_assembly_level = Higher_status, .after = gc_percent ) %>%
    mutate(Taxonomy = taxa, .before = size_mb) %>%
    mutate(Taxonomy_from_dataset = number_organism_name, .after = Taxonomy, .keep = 'unused') %>%
    mutate(number_of_genomes_with_highest_assembly_level = number_of_genomes_after_filter, .before = size_mb_min)
  
  genome_stats_table %>% 
    slice_max(size_mb)
}


