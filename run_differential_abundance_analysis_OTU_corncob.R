

# Analysing and plotting ASV differential abundances using corncob


# Usage:
# 
# Tax_order <- order_taxa(ps_obj, rank = "Phylum", rel_abund = TRUE)
# 
# comparison_str <- unlist(str_split(comparison, " vs. ")) # e.g. comparison <- "X vs. Y" 
#
# ps_obj_pairwise <- subset_physeq4pairwise(ps_obj, var2sub = "Treatment", comparison_str = comparison_str) # `var2sub` is a column in `sample_data(ps_obj)`, `comparison_str` contains two levels in `var2sub` 
# 
# # Remove species with prevalence < X%
# ps_obj_pairwise <- filter_prevalence(ps_obj = ps_obj_pairwise, prevalence = 0.1)
# 
# DA_obj <- run_corncob_DA(ps_obj = ps_obj_pairwise, var2test = "Treatment")
# 
# DA_df <- make_da_df(ps_obj = ps_obj_pairwise, da_obj = DA_obj, tax_rank = "Phylum")
# 
# DA_df %>% 
#   filter(Significance == "Pass") %>% 
#   write.csv(., file = paste0("corncob", "_", paste0(comparison, collapse = "_"), ".csv"))
# 
# p_corncob <- plot_corncob(da_df = DA_df, 
#                           p_title = comparison, 
#                           tax_rank = "Phylum",
#                           y_val = "Differential abundance", 
#                           sig_level = 0.05, 
#                           tax_order = Tax_order, 
#                           ASV_labels = F)



devtools::install_github("bryandmartin/corncob",
                         lib = lib.loc,
                         dependencies = TRUE)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(corncob)


# Order your taxa according to abundance (for plotting only)
order_taxa <- function(ps_obj, rank = "Phylum", rel_abund = TRUE){
  require(dplyr)
  require(phyloseq)
  require(speedyseq)
  
  ps_obj %>%
    tax_glom(taxrank = rank) %>%                     # agglomerate at 'Rank' level
    {if (rel_abund) transform_sample_counts(., function(x) x / sum(x)) else .} %>% # convert to rel abundance 
    psmelt() %>%                                        # Melt to long format
    arrange(rank) %>%                                  # arrange by 'Rank'
    group_by(across(all_of(rank))) %>% 
    summarise(Abundance = sum(Abundance)) %>%
    arrange(desc(Abundance)) %>% 
    mutate(across(all_of(rank), ~factor(., levels = fct_inorder(.)))) %>%  
    mutate(across(all_of(rank), ~fct_expand(., "Rare"))) %>% # just in case it doesn't exist
    mutate(across(all_of(rank), ~fct_relevel(., "Rare", after = Inf))) ->
    Tax_order 
  
  return(Tax_order)
}

# Filter OTUs by prevalence
filter_prevalence <- function(ps_obj, rank = "Phylum", prevalence = 0.1) {
  prevdf <- apply(
    X = otu_table(ps_obj),
    MARGIN = ifelse(taxa_are_rows(ps_obj), yes = 1, no = 2),
    FUN = function(x) {sum(x > 0)}
  )
  # Add taxonomy and total read counts to this data.frame
  prevdf <- data.frame(
    Prevalence = prevdf,
    TotalAbundance = taxa_sums(ps_obj),
    tax_table(ps_obj)
  )
  
  # Define prevalence threshold as 0.X of total samples
  prevalenceThreshold <- prevalence * nsamples(ps_obj)
  prevalenceThreshold
  
  # Execute prevalence filter, using `prune_taxa()` function
  
  prevdf %>% 
    subset(.,
           get(rank) %in% get_taxa_unique(ps_obj, rank)) ->
    prevdf_phylum_filt
  
  keepTaxa <- row.names(prevdf_phylum_filt)[(prevdf_phylum_filt$Prevalence >= prevalenceThreshold)]
  ps_obj_filt <- prune_taxa(keepTaxa, ps_obj)
  sample_data(ps_obj_filt)$Lib.size <-
    rowSums(otu_table(ps_obj_filt))
  print(ps_obj)
  print(ps_obj_filt)
  return(ps_obj_filt)
}

# Subset your phyloseq object according to the "levels" you want to test (e.g. labelled vs. control) var2sub is the variable that holds these levels.
# You can use glom_taxa and rank2glom to glom your data for testing at higher taxonomic levels (also faster) 
subset_physeq4pairwise <- function(ps_obj = Ps_obj_filt, var2sub = "Treatment", comparison_str = c("level1", "level2"), glom_taxa = FALSE, rank2glom = "Order") {
  ps_obj %>%
    subset_samples(get(var2sub) %in% comparison_str) %>% 
    {if (glom_taxa) tax_glom(., rank2glom) else .} %>%
    identity() -> # for debugging
    ps_obj_pairwise
  return(ps_obj_pairwise)
}

# This runs corncob. Please read the manual before you change anything here.
run_corncob_DA <- function(ps_obj = ps_obj_pairwise, var2test = "Treatment"){ 
  da_obj <- differentialTest(formula = ~ get(var2test),
                             phi.formula = ~ get(var2test),
                             formula_null = ~ 1,
                             phi.formula_null = ~ get(var2test), 
                             test = "Wald", 
                             boot = FALSE,
                             data = ps_obj,
                             fdr_cutoff = 0.05,
                             full_output = TRUE)
  which(is.na(da_obj$p)) %>% names
  return(da_obj)
}

# Convert the corncob object to a data.frame you can plot
make_da_df <- function(ps_obj = ps_obj_pairwise, da_obj = DA_obj, tax_rank = "Phylum") {
  DA_intervals <- plot(da_obj, data_only = TRUE)
  
  ps_obj %>%
    transform_sample_counts(., function(x) x / sum(x) * 100) %>% 
    taxa_sums(.) %>% 
    map_dbl(~(.x / nsamples(ps_obj))) %>% 
    enframe(name = "OTU", value = "Mean abundance (%)") -> 
    baseMean
  
  # Extract model coefficients
  # grab all mu.LocationSlope Estimates (differences in estimated population relative abundance)
  map(da_obj$all_models, 15) %>% # position 15 is the model coefficients
    map(.,2) %>%
    map(., ~ifelse(is.null(.x), NA, .x)) %>% # this ensures that all Nulls are converted to NA (otherwise unlist() will drop them)
    unlist(.) %>% 
    bind_cols(OTU = taxa_names(ps_obj), 
              tax_table(ps_obj), 
              `Differential abundance` = .,
              Significance = fct_recode(as_factor(taxa_names(ps_obj) %in% da_obj$significant_taxa), Pass = "TRUE", Fail = "FALSE"),
              ymin = as.numeric(NA),
              ymax = as.numeric(NA)
    ) %>%
    left_join(., baseMean, by = "OTU") ->
    da_df
  
  # Add confidence intervals (only works on the significant OTUs)
  da_df %<>% rows_update(., tibble(ymin = DA_intervals$xmin, OTU = da_obj$significant_taxa), by = "OTU")
  da_df %<>% rows_update(., tibble(ymax = DA_intervals$xmax, OTU = da_obj$significant_taxa), by = "OTU")
  return(da_df)
}

# Plot the results. 
plot_corncob <- function(da_df, p_title =  "X vs Y", tax_rank = "Phylum", y_val = "Differential abundance", sig_level = 0.05, tax_order = Tax_order, ASV_labels = FALSE) {
  # Plot differential abundance model results 
  
  pos <- position_jitter(width = 0.1, seed = 1)
  da_df %<>% 
    mutate_at(vars(matches(tax_rank)), as_factor) %>% 
    mutate(across(tax_rank, ~fct_expand(., levels(pull(tax_order, !!tax_rank))))) %>%  # just in case it doesn't match
    mutate(across(tax_rank, ~fct_relevel(., levels(pull(tax_order, !!tax_rank))))) %>% # reorder taxa according to what found in taxa_order
    mutate(across(tax_rank, ~fct_relevel(., "Rare", after = Inf))) 
  
  corncob_summary <- tibble(Label = c(paste0("⬆", 
                                             sum(da_df$`Differential abundance` > 0 &  da_df$Significance == "Pass"), 
                                             " ⬇", 
                                             sum(da_df$`Differential abundance` < 0 &  da_df$Significance == "Pass"), 
                                             " (", 
                                             nrow(da_df), 
                                             ")")))
  
  p <-
    ggplot(da_df) +
    geom_point(aes(
      x = !!sym(tax_rank),
      y = !!sym(y_val),
      colour = !!sym("Significance"),
      size = !!sym("Mean abundance (%)")),
      position = pos, 
      alpha = 2 / 3, 
      stroke = 0) +
    geom_linerange(aes(x = !!sym(tax_rank),
                       y = !!sym(y_val),
                       ymin = `ymin`,
                       ymax = `ymax`,
                       colour = !!sym("Significance")),
                   position = pos,
                   alpha = 1/5, 
                   show.legend = FALSE) +
    geom_text(
      data    = corncob_summary,
      mapping = aes(x = Inf, y = Inf, label = Label),
      hjust   = 1.1,
      vjust   = 1.6
    ) +
    xlab("") +
    ylab("Differential abundance") +
    labs(colour = paste("Significance at \n p <", sig_level), size = "Mean abundance (%)") +
    theme_grey(base_size = 18, ) +
    theme(axis.text.x = element_text(angle = 45.0, vjust = 1, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    scale_colour_manual(values = c(ggpomological:::pomological_base[[7]], ggpomological:::pomological_palette[[1]])) +
    scale_size_continuous(name = "Mean abundance (%)",
                          range = c(2, 8),
                          breaks = c(round(seq(min(da_df$`Mean abundance (%)`), max(da_df$`Mean abundance (%)`), length.out = 5), 1)))
  
  if (ASV_labels) {
    p <- p + geom_label_repel(
      aes(x = !!sym(tax_rank), y = !!sym(y_val)),
      size = 6,
      label = sub("Seq_([0-9]+)", "\\1", pull(da_df[da_df$Significance == "Pass", ], "OTU")),
      position = pos,
      data = da_df[da_df$Significance == "Pass", ],
      # nudge_x = 0.4,
      colour = "#4a4a4a",
      label.size = NA, 
      alpha = 0.75, 
      # fontface = 'bold',
      force = 1,
      force_pull = 1,
      box.padding = 0.4,
      point.padding = 0.1,
      max.overlaps = 100
    )
  }
  p <- p + labs(title = p_title)
  return(p)
}
