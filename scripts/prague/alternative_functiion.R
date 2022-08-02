genomic_traits <- function(taxa){
  
  Assembly_level <- c("Complete Genome","Chromosome","Scaffold","Contig")
  first(Assembly_level)
  
  number_of_entries_in_genomic_db <- number_of_entries %>% 
                                          filter(taxonomy == taxa) %>% 
                                          pull(n)
  
  size_gc_status <- best_tax_prok %>% 
                        filter(taxonomy == taxa) %>% 
                        select(size_mb,gc_percent,status) %>%
                        mutate(number_of_entries_in_NCBI = as.numeric(number_of_entries_in_genomic_db))
  
  Higher_status <- Assembly_level[min(match(size_gc_status$status,Assembly_level))]
  
  used_genomes <- size_gc_status %>% 
    filter(status == Higher_status)

  number_of_genomes_after_filter <- nrow(used_genomes)
  
  df_summ_entries <- used_genomes %>% 
      select(size_mb, gc_percent) %>% 
      summarize(across(everything(), ~summary(.) %>% as_tibble_row() %>% list())) %>% 
      unnest(cols = c(size_mb, gc_percent), names_sep = "_") %>%
      setNames(str_replace(names(.), "\\.", "")) %>%
      setNames(str_replace(names(.), " Qu", "_qtl")) %>%
      setNames(str_replace(names(.), " Qu", "_qtl")) %>%
      as_tibble(.name_repair = janitor::make_clean_names)
 
  genome_stats_table <- used_genomes %>% 
    select(-status) %>%
    mutate(df_summ_entries) %>%
    mutate(highest_assembly_level = Higher_status, .after = gc_percent ) %>%
    mutate(Taxonomy = taxa, .before = size_mb) %>%
    mutate(number_of_genomes_with_highest_assembly_level = number_of_genomes_after_filter, .before = size_mb_min)
  
  genome_stats_table %>% 
    slice_max(size_mb)
}


