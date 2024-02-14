
# Helper functions --------------------------------------------------------



remove_suffix_from_gene_id <- function(dataframe, column = "gene_id") {
  #' Function that removes the suffix number from the ENSG gene code
  #' @param dataframe with column for gene_id that contains the ENSG code
  #' @param column the column that contains the variables
  #' @return dataframe with column with suffix removed

  col_name = rlang::sym(column)
  dataframe <- dataframe %>%
    dplyr::mutate(!!col_name := str_remove(!!col_name, "\\.\\d+$"))
  return(dataframe)
}



change_CHR_gencode_to_ensembl <- function(gwas_df) {
  #' Function that converts formating of gencode CHR to ensembl
  #' @param gwas_df dataframe  with the column CHR to transfor
  #' @return dataframe with reformated columns
  gwas_df <- gwas_df %>%
    dplyr::mutate(CHR = str_remove(CHR, "chr")) %>%
    dplyr::mutate(CHR = str_replace(CHR, "M", "MT")) # remove chr and recode M to MT
  return(gwas_df)
}


change_CHR_ensembl_to_gencode <- function(gwas_df) {
  #' Function that converts formating of gencode CHR to ensembl
  #' @param gwas_df dataframe  with the column CHR to transfor
  #' @return dataframe with reformated columns
  gwas_df %>%
    dplyr::mutate(CHR = str_replace(CHR, "MT", "M")) %>%  #change MT to M
    dplyr::mutate(CHR = str_c("chr", CHR)) # Add in chr prefix
  return(gwas_df)

}



load_edit_sites_file <- function(file_path, metadata_df) {
  #' Function that filters metadata based on location and contrast to explore,
  #' and loads relevant editing sites into a dataframe
  #' @param file_path the path to where files are located
  #' @param metadata_df the metadata dataframe
  sites <- read_delim(file_path)
  sites <- sites %>%
    dplyr::filter(z.score <= -1.96 | z.score >= 1.96) %>%
    left_join(metadata,
              by = c("Sample" = "sample_id"))
  return(sites)
}


filter_and_load_editing_results <- function(metadata_df, edit_files_list, contrast_to_explore, region_to_explore){
  #' Function that filters metadata based on location and contrast to explore,
  #' and loads relevant editing sites into a dataframe
  #' @param metadata dataframe with columns for sample_id, brain_region, Disease_Group
  #' @param contrast_to_explore a vector of contrasts to explore: eg control vs disease
  #' @param edit_sites_files the full list of files to explore (with outliers removed)
  #' @param region_to_explore a region to explore
  #' @return dataframe with editing sites and metadata

  if(is.null(region_to_explore)){

    print(str_c("Exploring",  ": ", contrast_to_explore[1], " vs. ", contrast_to_explore[[2]]))

    # Find samplename of files to explore
    files_to_explore <- metadata_df %>%
      dplyr::filter(Disease_Group %in% contrast_to_explore) %>%
      dplyr::select(sample_id) %>%
      unlist()

    # Find full filename
    edit_files <- edit_files_list[str_detect(edit_files_list, str_c(files_to_explore, collapse = "|"))]

    # Load edits and combine to df
    edit_sites_list <- lapply(edit_files, load_edit_sites_file, metadata_df = metadata)
    comparison_df <- do.call(rbind, edit_sites_list)

    return(comparison_df)
  } else if (!is.null(region_to_explore)){

    print(str_c("Exploring ", region_to_explore, ": ", contrast_to_explore[1], " vs. ", contrast_to_explore[[2]]))

    # Find samplename of files to explore
    files_to_explore <- metadata_df %>%
      dplyr::filter(brain_region == region_to_explore,
                    Disease_Group %in% contrast_to_explore) %>%
      dplyr::select(sample_id) %>%
      unlist()

    # Find full filename
    edit_files <- edit_files_list[str_detect(edit_files_list, str_c(files_to_explore, collapse = "|"))]

    # Load edits and combine to df
    edit_sites_list <- lapply(edit_files, load_edit_sites_file, metadata_df = metadata)
    comparison_df <- do.call(rbind, edit_sites_list)

    return(comparison_df)
  }
}

filter_sites_by_locus_vector <- function(sites_df, vector_of_loci, factor_list) {
  #' Efficient function for filtering a large dataframe for loci in a vector, using datatable package,
  #' and then filter for sites within both disease groups of interest
  #' @param sites_df dataframe with the column for locus
  #' @param vector_of_loci a vector of loci to use in filtering
  #' @return filtered dataframe
  # Use datatable package to filter the comparison table for rows that are within the complete_site_loci
  setDT(sites_df)

  filtered_sites <- sites_df[locus %in% common_loci]

  # Filter the dataframe for sites in both groups
  filtered_sites <- filtered_sites %>%
    group_by(locus) %>%
    filter(length(unique(Disease_Group)) == length(factor_list)) %>%
    ungroup()

  return(filtered_sites)
}



find_sites_loci_common_to_x_samples <- function(sites_df, sample_counts_df, no_common_samples_per_site) {
  #' Function that finds editing loci common to x number of samples with a treatment/control group
  #' @param sites_df dataframe  with the column for locus
  #' @param sample_counts_df a dataframe with a column called `Disease_Group`, and a column called `n`
  #' @param no_common_samples_per_site an integer, the number of common sites per disease/treatment group required
  #' @return vector of shared loci for each group

  filtered_sites = list()

  for (i in 1:nrow(sample_counts_df)) {
    group = sample_counts_df[i, 'Disease_Group'] %>% unlist()
    group_count = sample_counts_df[i, 'n'] %>%  unlist()

    if (no_common_samples_per_site > group_count) {
      print("The number of samples requested for the site to appear is greater than the number of Disease groups")
      break
    } else {

      print(str_c(group, " group; ", group_count, " samples"))
      group_sites <- sites_df %>%
        dplyr::filter(Disease_Group == group) %>%
        group_by(locus) %>%
        summarise(num_samples = n_distinct(Sample)) %>% # Count the number of distinct samples per locus
        dplyr::filter(num_samples >= no_common_samples_per_site) %>%
        dplyr::mutate(Disease_Group = group) %>%
        dplyr::select(locus) %>%
        unlist() %>%
        as.vector()
      filtered_sites[[group]] <- group_sites #create list of group sites

    }
  }
  # Find the intersect of the sites in each control/disease group
  # Reduce performs actions on successive components of a list, here we are using it with the intersect function.
  near_complete_site_loci <- Reduce(intersect, filtered_sites)
  return(near_complete_site_loci)
}


find_new_lost_sites <- function(sites_df, sample_counts_df, no_common_samples_per_site) {
  #' Function that finds sites common to 2 samples in each control/disease group, and which of these sites are new/lost
  #' @param sites_df dataframe  with the column for locus
  #' @param sample_counts_df a dataframe with a column called `Disease_Group`, and a column called `n`
  #' @param no_common_samples_per_site an integer, the number of common sites per disease/treatment group required
  #' @return list of sites unique to disease or control

  filtered_sites = list()

  for (i in 1:nrow(sample_counts_df)) {
    group = sample_counts_df[i, 'Disease_Group'] %>% unlist() %>% as.character
    group_count = sample_counts_df[i, 'n'] %>%  unlist() %>%  as.character

    if (no_common_samples_per_site > group_count) {
      print("The number of samples requested for the site to appear is greater than the number of samples")
      break
    } else {

      print(str_c(group, " group; ", group_count, " samples"))
      group_sites <- sites_df %>%
        dplyr::filter(Disease_Group == group) %>%
        group_by(locus) %>%
        summarise(num_samples = n_distinct(Sample)) %>% # Count the number of distinct samples per locus
        dplyr::filter(num_samples >= no_common_samples_per_site) %>%
        dplyr::mutate(Disease_Group = group) %>%
        dplyr::select(locus) %>%
        unlist() %>%
        as.vector()
      filtered_sites[[group]] <- group_sites #create list of group sites

    }
  }
  # Find the intersect of the sites in each control/disease group
  # Reduce performs actions on successive components of a list, here we are using it with the intersect function.
  control_unique <- tibble(locus = setdiff(filtered_sites$Control, filtered_sites$PD),
                           new_lost = "Lost")
  pd_unique <- tibble(locus = setdiff(filtered_sites$PD, filtered_sites$Control),
                      new_lost = "New")

  new_lost_df <- bind_rows(control_unique, pd_unique)
  return(new_lost_df)
}


# Create editing summary metric -------------------------------------------

get_editing_summary_metric <- function(summary_df, gene_lengths_df, outcome_measure, covariates, factor_list) {
  #' Takes a dataframe of editing sites and creates a new column of a summary editing metric per gene
  #' The metric is editing no_editing_sites_per_gene * mean_editing_rate / total_exon_length
  #' @param dataframe, a dataframe of editing sites including the columns:
  #'   * locus: defining the editing site in format chr_bp_strand
  #'   * gene_id: ENSG code
  #'   * gene_name: gene symbol
  #'   * Edits: editing proportion at that site
  #' @param  dataframe with exon and gene lengths per gene
  #' @return dataframe with additional columns

  # Prepare the grouping variables
  grouping_vars <- c("Sample", "gene_id", outcome_measure, unlist(covariates))
  print("Grouping gene summary by:")
  print(grouping_vars)

  # Create summary table with the number of edits and mean edits per gene
  edits_per_gene <- summary_df %>%
    drop_na(gene_id) %>%
    group_by(across(all_of(grouping_vars))) %>%
    dplyr::summarise(no_sites = n(),
                     mean_rate = mean(Edits),
                     median_rate = median(Edits),
                     mode_rate = custom_mode(Edits),
                     .groups = "drop")

  # Add in the gene length and exon length
  edits_per_gene <- edits_per_gene %>%
    left_join(.,
              gene_lengths_df,
              by = "gene_id")

  number_of_genes_without_length <- edits_per_gene %>%
    group_by(gene_id) %>%
    dplyr::slice(1) %>%
    dplyr::filter(is.na(gene_length)) %>%
    nrow()


  total_number_of_genes <- edits_per_gene %>%
    group_by(gene_id) %>%
    dplyr::slice(1) %>%
    nrow()
  print(str_c(number_of_genes_without_length, " out of ", total_number_of_genes, " genes do not have a gene length."))


  # Complete summary metric
  edits_per_gene <- edits_per_gene %>%
    drop_na(gene_length) %>%
    dplyr::mutate(complete_mean_metric = (no_sites * mean_rate / exon_length),
                  complete_median_metric = (no_sites * median_rate / exon_length),
                  complete_mode_metric = (no_sites * mode_rate / exon_length),
                  mean_over_exon_length = (mean_rate / exon_length),
                  mean_x_no_sites = (mean_rate * no_sites),
                  complete_mean_genelength = (no_sites * mean_rate / gene_length))

  max_mean_edit <- max(edits_per_gene$complete_mean_metric)
  max_median_edit <- max(edits_per_gene$complete_median_metric)
  max_mode_edit <- max(edits_per_gene$complete_mode_metric)
  max_no_sites <- max(edits_per_gene$no_sites)
  max_mean_over_exon_length = max(edits_per_gene$mean_over_exon_length)
  max_mean_x_no_sites = max(edits_per_gene$mean_x_no_sites)
  max_mean = max(edits_per_gene$mean_rate)
  max_complete_genelength = max(edits_per_gene$complete_mean_genelength)

  edits_per_gene <- edits_per_gene %>%
    dplyr::mutate(complete_mean_adjusted = (complete_mean_metric / max_mean_edit),
                  complete_median_adjusted = (complete_median_metric / max_median_edit),
                  complete_mode_adjusted = (complete_mode_metric / max_mode_edit),
                  no_sites_adjusted = (no_sites / max_no_sites),
                  mean_over_exon_length_adjusted = (mean_over_exon_length / max_mean_over_exon_length),
                  mean_x_no_sites_adjusted = (mean_x_no_sites / max_mean_x_no_sites),
                  mean_adjusted = (mean_rate / max_mean),
                  complete_mean_genelength_adjusted = (complete_mean_genelength / max_complete_genelength))

  # Filter the dataframe to ensure the gene is present in both groups
  edits_per_gene <- edits_per_gene %>%
    group_by(gene_id) %>%
    filter(length(unique(Disease_Group)) == length(factor_list)) %>%
    ungroup()

  return(edits_per_gene)

}



get_editing_summary_metric_scaled <- function(summary_df, gene_lengths_df, outcome_measure, covariates, factor_list) {
  #' Takes a dataframe of editing sites and creates a new column of a summary editing metric per gene
  #' The metric is editing no_editing_sites_per_gene * mean_editing_rate / total_exon_length
  #' This function takes the log10 of no_sites_per_gene, and scales other continuous variables
  #' @param dataframe, a dataframe of editing sites including the columns:
  #'   * locus: defining the editing site in format chr_bp_strand
  #'   * gene_id: ENSG code
  #'   * gene_name: gene symbol
  #'   * Edits: editing proportion at that site
  #' @param  dataframe with exon and gene lengths per gene
  #' @return dataframe with additional columns

  # Prepare the grouping variables
  grouping_vars <- c("Sample", "gene_id", outcome_measure, unlist(covariates))
  print("Grouping gene summary by:")
  print(grouping_vars)

  # Create summary table with the number of edits and mean edits per gene
  edits_per_gene <- summary_df %>%
    drop_na(gene_id) %>%
    group_by(across(all_of(grouping_vars))) %>%
    dplyr::summarise(no_sites = n(),
                     mean_rate = mean(Edits),
                     median_rate = median(Edits),
                     mode_rate = custom_mode(Edits),
                     .groups = "drop") %>%
    dplyr::mutate(no_sites = log10(no_sites))

  # Add in the gene length and exon length
  edits_per_gene <- edits_per_gene %>%
    left_join(.,
              gene_lengths_df,
              by = "gene_id")

  number_of_genes_without_length <- edits_per_gene %>%
    group_by(gene_id) %>%
    dplyr::slice(1) %>%
    dplyr::filter(is.na(gene_length)) %>%
    nrow()


  total_number_of_genes <- edits_per_gene %>%
    group_by(gene_id) %>%
    dplyr::slice(1) %>%
    nrow()
  print(str_c(number_of_genes_without_length, " out of ", total_number_of_genes, " genes do not have a gene length."))


  # Complete summary metric
  edits_per_gene <- edits_per_gene %>%
    drop_na(gene_length) %>%
    dplyr::mutate(complete_mean_metric = (no_sites * mean_rate / exon_length),
                  complete_median_metric = (no_sites * median_rate / exon_length),
                  complete_mode_metric = (no_sites * mode_rate / exon_length),
                  mean_over_exon_length = (mean_rate / exon_length),
                  mean_x_no_sites = (mean_rate * no_sites),
                  complete_mean_genelength = (no_sites * mean_rate / gene_length))

  max_mean_edit <- max(edits_per_gene$complete_mean_metric)
  max_median_edit <- max(edits_per_gene$complete_median_metric)
  max_mode_edit <- max(edits_per_gene$complete_mode_metric)
  max_no_sites <- max(edits_per_gene$no_sites)
  max_mean_over_exon_length = max(edits_per_gene$mean_over_exon_length)
  max_mean_x_no_sites = max(edits_per_gene$mean_x_no_sites)
  max_mean = max(edits_per_gene$mean_rate)
  max_complete_genelength = max(edits_per_gene$complete_mean_genelength)

  edits_per_gene <- edits_per_gene %>%
    dplyr::mutate(complete_mean_adjusted = (complete_mean_metric / max_mean_edit),
                  complete_median_adjusted = (complete_median_metric / max_median_edit),
                  complete_mode_adjusted = (complete_mode_metric / max_mode_edit),
                  no_sites_adjusted = (no_sites / max_no_sites),
                  mean_over_exon_length_adjusted = (mean_over_exon_length / max_mean_over_exon_length),
                  mean_x_no_sites_adjusted = (mean_x_no_sites / max_mean_x_no_sites),
                  mean_adjusted = (mean_rate / max_mean),
                  complete_mean_genelength_adjusted = (complete_mean_genelength / max_complete_genelength))

  # Filter the dataframe to ensure the gene is present in both groups
  edits_per_gene <- edits_per_gene %>%
    group_by(gene_id) %>%
    filter(length(unique(Disease_Group)) == length(factor_list)) %>%
    ungroup()

  return(edits_per_gene)

}



#'
#' get_editing_summary_metric_no_siteno <- function(summary_df, gene_lengths_df, outcome_measure, covariates) {
#'   #' Takes a dataframe of editing sites and creates a new column of a summary editing metric per gene
#'   #' The metric is editing no_editing_sites_per_gene * mean_editing_rate / total_exon_length
#'   #' @param dataframe, a dataframe of editing sites including the columns:
#'   #'   * locus: defining the editing site in format chr_bp_strand
#'   #'   * gene_id: ENSG code
#'   #'   * gene_name: gene symbol
#'   #'   * Edits: editing proportion at that site
#'   #' @param  dataframe with exon and gene lengths per gene
#'   #' @return dataframe with additional columns
#'
#'   # Prepare the grouping variables
#'   grouping_vars <- c("Sample", "gene_name", outcome_measure, unlist(covariates))
#'   print("Grouping gene summary by:")
#'   print(grouping_vars)
#'
#'   # Create summary table with the number of edits and mean edits per gene
#'   edits_per_gene <- summary_df %>%
#'     drop_na(gene_name) %>%
#'     group_by(across(all_of(grouping_vars))) %>%
#'     dplyr::summarise(no_sites = n(),
#'                      mean_rate = mean(Edits),
#'                      median_rate = median(Edits),
#'                      mode_rate = custom_mode(Edits),
#'                      .groups = "drop")
#'
#'   # Add in the gene length and exon length
#'   edits_per_gene <- edits_per_gene %>%
#'     left_join(.,
#'               gene_lengths_df,
#'               by = "gene_name")
#'
#'   number_of_genes_without_length <- edits_per_gene %>%
#'     group_by(gene_name) %>%
#'     dplyr::slice(n=1) %>%
#'     dplyr::filter(is.na(gene_length)) %>%
#'     nrow()
#'
#'
#'   total_number_of_genes <- edits_per_gene %>%
#'     group_by(gene_name) %>%
#'     dplyr::slice(n=1) %>%
#'     nrow()
#'   print(str_c(number_of_genes_without_length, " out of ", total_number_of_genes, " genes do not have a gene length."))
#'
#'   edits_per_gene <- edits_per_gene %>%
#'     drop_na(gene_length) %>%
#'     dplyr::mutate(edit_summary_metric = ( mean_rate / exon_length))
#'
#'   max_mean_edit <- max(edits_per_gene$edit_summary_metric)
#'
#'   edits_per_gene <- edits_per_gene %>%
#'     dplyr::mutate(edit_mean_adjusted = (edit_summary_metric / max_mean_edit)
#'     )
#'
#'   # Filter the dataframe to ensure the gene is present in both groups
#'   edits_per_gene <- edits_per_gene %>%
#'     group_by(gene_name) %>%
#'     filter(length(unique(Disease_Group)) == length(factor_list)) %>%
#'     ungroup()
#'
#'   return(edits_per_gene)
#'
#' }
#'
#'
#'
#'
#' get_editing_summary_metric_mean_edit_per_gene <- function(summary_df, gene_lengths_df, outcome_measure, covariates) {
#'   #' Takes a dataframe of editing sites and creates a new column of a summary editing metric per gene
#'   #' The metric is editing no_editing_sites_per_gene * mean_editing_rate / total_exon_length
#'   #' @param dataframe, a dataframe of editing sites including the columns:
#'   #'   * locus: defining the editing site in format chr_bp_strand
#'   #'   * gene_id: ENSG code
#'   #'   * gene_name: gene symbol
#'   #'   * Edits: editing proportion at that site
#'   #' @param  dataframe with exon and gene lengths per gene
#'   #' @return dataframe with additional columns
#'
#'   # Prepare the grouping variables
#'   grouping_vars <- c("Sample", "gene_name", outcome_measure, unlist(covariates))
#'   print("Grouping gene summary by:")
#'   print(grouping_vars)
#'
#'   # Create summary table with the number of edits and mean edits per gene
#'   edits_per_gene <- summary_df %>%
#'     drop_na(gene_name) %>%
#'     group_by(across(all_of(grouping_vars))) %>%
#'     dplyr::summarise(no_sites = n(),
#'                      mean_rate = mean(Edits),
#'                      median_rate = median(Edits),
#'                      mode_rate = custom_mode(Edits),
#'                      .groups = "drop")
#'
#'   # Add in the gene length and exon length
#'   edits_per_gene <- edits_per_gene %>%
#'     left_join(.,
#'               gene_lengths_df,
#'               by = "gene_name")
#'
#'   number_of_genes_without_length <- edits_per_gene %>%
#'     group_by(gene_name) %>%
#'     dplyr::slice(n=1) %>%
#'     dplyr::filter(is.na(gene_length)) %>%
#'     nrow()
#'
#'
#'   total_number_of_genes <- edits_per_gene %>%
#'     group_by(gene_name) %>%
#'     dplyr::slice(n=1) %>%
#'     nrow()
#'   print(str_c(number_of_genes_without_length, " out of ", total_number_of_genes, " genes do not have a gene length."))
#'
#'   edits_per_gene <- edits_per_gene %>%
#'     drop_na(gene_length) %>%
#'     dplyr::mutate(edit_summary_metric = ( mean_rate))
#'
#'   max_mean_edit <- max(edits_per_gene$edit_summary_metric)
#'
#'   edits_per_gene <- edits_per_gene %>%
#'     dplyr::mutate(edit_mean_adjusted = (edit_summary_metric / max_mean_edit)
#'     )
#'
#'   # Filter the dataframe to ensure the gene is present in both groups
#'   edits_per_gene <- edits_per_gene %>%
#'     group_by(gene_name) %>%
#'     filter(length(unique(Disease_Group)) == length(factor_list)) %>%
#'     ungroup()
#'
#'   return(edits_per_gene)
#'
#' }



run_linear_regression_edits <- function(sites_df, editing_metric, covariates, outcome_measure){
  #' Function that finds editing loci common to x number of samples with a treatment/control group
  #' @param sites_df dataframe  with the column for locus (either )
  #' @param editing_metric column of df to use as outcome variable for regression
  #' @param covariates a vector of covariates
  #' @return the linear regression results


  form <- reformulate(c(outcome_measure, covariates), response = editing_metric)
  group_covariates_lm <- lm(form, data = sites_df)

  return(group_covariates_lm)

}

# Define a custom mode function
custom_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



retrieve_lm_residuals <- function(dataframe, linear_model) {
  #' This function is meant to run on two types of data.
  #'   * The first is a dataframe of editing sites with a column for editing_locus and mean_editing_rate_per_locus
  #'   * The second is a dataframe of a summarised editing metric per gene, with a column for gene and editing_metric.
  #'It will takes results of a linear regression model and find the residuals for each datapoint
  #' @param lm_summary, results of a linear model of either mean editing rate "Edits" per "Locus", or the adjusted editing metric "edits_adjusted" per "gene_name"
  #' @param dataframe of the original data input into the regression model. In addition to the columns above, should also have a column for "Disease_Group" which contains the "disease" and control variables
  #' @param disease, which disease is being explored, relative to a the "Control" variable



  sites_subset_fitted <- broom::augment(x = linear_model,
                                        data = dataframe)
}



retrieve_baseline_mean_residuals_per_feature <- function(residuals_df, feature, factor_list) {
  #' This function is meant to run on two types of data. The first is a dataframe of residuals of with a column for editing_locus and mean_editing_rate_per_locus,
  #' the second is a dataframe of a residuals of the summarised editing metric per gene, with a column for gene and editing_metric. It will:
  #'    a) Finds the  mean of residuals across control factorised variable,
  #'    b) Adds a column per feature/locus outlining whether the mean of residuals is changed in the disease variable.
  #' @param residuals_df, dataframe with residuals mapped to each feature of interest. The output of results of a linear model of either mean editing rate "Edits" per "Locus", or the adjusted editing metric "edits_adjusted" per "gene_name"
  #' @param feature, reflecting the column of interest, either "locus" or "gene_name".
  #' @param disease, which disease is being explored, relative to a the "Control" variable

  feature_sym <- sym(feature)
  control_metric <- sym(factor_list[[1]])
  disease_metric <- sym(factor_list[[2]])

  mean_residuals_per_site <- residuals_df %>%
    dplyr::filter(Disease_Group == control_metric) %>%
    group_by(!!feature_sym) %>%
    summarise(mean_residuals_per_feature = mean(.resid), .groups = "drop") %>%
    ungroup() %>%
    dplyr::mutate(Direction = "Baseline")


  return(mean_residuals_per_site)

}

retrieve_mean_residuals_per_feature <- function(residuals_df, feature, factor_list) {
  #' This function is meant to run on two types of data. The first is a dataframe of residuals of with a column for editing_locus and mean_editing_rate_per_locus,
  #' the second is a dataframe of a residuals of the summarised editing metric per gene, with a column for gene and editing_metric. It will:
  #'    a) Finds the  mean of residuals across a disease/control factorised variable,
  #'    b) Adds a column per feature/locus outlining whether the mean of residuals is changed in the disease variable.
  #' @param residuals_df, dataframe with residuals mapped to each feature of interest. The output of results of a linear model of either mean editing rate "Edits" per "Locus", or the adjusted editing metric "edits_adjusted" per "gene_name"
  #' @param feature, reflecting the column of interest, either "locus" or "gene_name".
  #' @param disease, which disease is being explored, relative to a the "Control" variable

  feature_sym <- sym(feature)
  control_metric <- sym(factor_list[[1]])
  disease_metric <- sym(factor_list[[2]])

  mean_residuals_per_site <- residuals_df %>%
    group_by(Disease_Group, !!feature_sym) %>%
    summarise(mean_residuals_per_feature = mean(.resid), .groups = "drop") %>%
    ungroup()


  mean_residuals_per_site <- mean_residuals_per_site %>%
    pivot_wider(names_from = Disease_Group,
                values_from = mean_residuals_per_feature) %>%
    dplyr::mutate(edit_difference_per_feature = !!disease_metric - !!control_metric, # PD - Control. Factor levels need to be properly defined: untreated/control first
                  Direction = case_when(edit_difference_per_feature > 0 ~ "increased",
                                        edit_difference_per_feature < 0 ~ "decreased",
                                        TRUE ~ "equal_sites"))

  mean_residuals_per_site <- mean_residuals_per_site %>%
    pivot_longer(cols = c(disease_metric, control_metric),
                 names_to = "Disease_Group",
                 values_to = "mean_residuals_per_feature") %>%
    # dplyr::filter(Direction != "equal_sites") %>% # Remove sites where mean editing equals each other
    arrange(!!feature_sym) %>%
    dplyr::relocate(mean_residuals_per_feature, Direction, edit_difference_per_feature,  .after = last_col())


  mean_residuals_per_site$Disease_Group <- factor(mean_residuals_per_site$Disease_Group)
  levels(mean_residuals_per_site$Disease_Group) <- factor_list

  return(mean_residuals_per_site)

}



# Create plots ------------------------------------------------------------

# Function to create a title plot
title_plot <- function(title) {
  ggplot() +
    theme_void() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = title, fontface = "bold", hjust = 0.5)
}

plot_site_number_boxplot <- function(sites_df){
  number_plot <- sites_df %>%
    count(Sample, Disease_Group) %>%
    ggplot(aes(x = Disease_Group, y = n, fill = Disease_Group)) +
    geom_boxplot() +
 #   geom_jitter(width = 0.2) +
    theme_aw +
 #   scale_fill_manual(values = named_palette) +
    labs(y = "Number of sites", x = NULL)

  return(number_plot)
}

plot_site_number_violin <- function(sites_df, plot_label){
    number_plot <- sites_df %>%
        count(Sample, Disease_Group) %>%
        ggplot(aes(x = Disease_Group, y = n, fill = Disease_Group)) +
        geom_violin() + #draw_quantiles = c(0.25, 0.5, 0.75)) +
geom_jitter(width = 0.1) +
        theme_aw +
        scale_fill_manual(values = named_palette,
                          labels = {{plot_label}}) +
        labs(y = "Number of sites", x = NULL, fill = "Disease Group")

    return(number_plot)
}


create_box_plot_editing_metric <- function(edit_summary_df, plotting_variable) {
  #' This function takes a lm result and creates a forest plot from summarising the result of interest. It shows only the dependent variable of interest, not all covariates.
  #' @param editing_summary_df, result of get_editing_summary_metric
  #' @param plotting_variable, variable of interest to plot, input in ""
  #' @return box_plot
  #'
  sym_plotting_variable = sym(plotting_variable)
  edit_summary_box <- edit_summary_df %>%
    ggplot(aes(x = Sample, y = !!sym_plotting_variable, fill = Disease_Group)) +
    geom_boxplot(outlier.alpha = 0.2) +
    theme_aw +
    scale_y_continuous(trans = "log10",
                       labels = scales::comma,
                       limits = c(1e-6, 1.1),
                       breaks = c(1e-6, 1e-4, 1e-2, 1)) +
    scale_fill_manual(values = named_palette)

  return(edit_summary_box)

}


create_box_plot_editing_metric_per_disease_group <- function(edit_summary_df, plotting_variable) {
  #' This function takes a lm result and creates a forest plot from summarising the result of interest. It shows only the dependent variable of interest, not all covariates.
  #' @param editing_summary_df, result of get_editing_summary_metric
  #' @param plotting_variable, variable of interest to plot, input in ""
  #' @return box_plot
  #'
  sym_plotting_variable = sym(plotting_variable)
  edit_summary_box <- edit_summary_df %>%
    ggplot(aes(x = Disease_Group, y = !!sym_plotting_variable, fill = Disease_Group)) +
    geom_boxplot(outlier.alpha = 0.2) +
    theme_aw +
    scale_y_continuous(trans = "log10") +
    scale_fill_manual(values = named_palette) +
    labs(x = NULL)

  return(edit_summary_box)

}


create_editing_summary_plot_series <- function(sites_subset_df, edit_summary_df){
  #' This function creates a series of summary plots takes a lm result and creates a forest plot from summarising the result of interest. It shows only the dependent variable of interest, not all covariates.
  #' @param sites_subset_df, summary datraframe - the input to get_editing_summary_metric
  #' @param editing_summary_df, result of get_editing_summary_metric
  #' @return patchwork object of plots for editing summary

  edit_mean_rate_box <- sites_subset_df %>%
    create_box_plot_editing_metric(.,  "Edits") +
    labs(subtitle = "Editing rates per site",
         y = "Editing rate")

  genic_type_plot <- sites_subset_df %>%
    ggplot(aes(x = Sample, fill = gene_type)) +
    geom_bar(stat = "count", position = "fill" ) +
    scale_fill_manual(values = named_palette) +
    theme_aw +
    labs(y = "Propotion", subtitle = "Proportion of genic type of sites")

  edit_mean_metric_box <- edit_summary_df %>%
    create_box_plot_editing_metric(., "edit_summary_metric") +
    labs(y = "Editing metric", subtitle = "Editing metric: mean")

  adjusted_mean_metric_box <- create_box_plot_editing_metric(edit_summary_df, "edit_mean_adjusted") +
    labs(y = "Adjusted editing metric", subtitle = "Adjusted editing metric: mean")

  adjusted_median_metric_box <- create_box_plot_editing_metric(edit_summary_df, "edit_median_adjusted") +
    labs(y = "Adjusted editing metric", subtitle = "Adjusted editing metric: median")

  adjusted_mode_metric_box <- create_box_plot_editing_metric(edit_summary_df, "edit_mode_adjusted") +
    labs(y = "Adjusted editing metric", subtitle = "Adjusted editing metric: mode")

  site_number_histogram <- edit_summary_df %>%
    ggplot(aes(x = no_sites)) +
    geom_histogram(fill = "#4DBBD5FF", bins = 30) +
    scale_x_continuous(trans = "log10", labels = scales::comma) +
    theme_aw +
    labs(y = "Count", x = "Number of sites per gene", subtitle = "Number of sites per gene")

  summary_plot <- (edit_mean_rate_box | genic_type_plot | edit_mean_metric_box | adjusted_mean_metric_box | adjusted_median_metric_box | adjusted_mode_metric_box | site_number_histogram | plot_layout(guides = "collect"))

  return(summary_plot)

}

create_forest_plot_from_lm <- function(linear_model, disease) {
  #' This function takes a lm result and creates a forest plot from summarising the result of interest. It shows only the dependent variable of interest, not all covariates.
  #' @param linear_model, result of linear model
  #' @param disease, disease being explored (relative to control)
  #' @return forest plot

  # Retrieve estimates + confidence intervals
  coefs <- tidy(linear_model, conf.int = TRUE) %>%
    tidy_categorical(m = linear_model) %>%
    dplyr::filter(term == str_c('Disease_Group', disease)) # Select only the disease line

  # Add nicer names
  # coefs$names <- c( "PD (ref: Control)")

  # Gather point estimates and CIs
  coefs$display_ests <-sprintf("%.3f [%.3f, %.3f]" , coefs$estimate, coefs$conf.low, coefs$conf.high)

  # Create plot
  #  forest_plot <- ggplot( coefs, aes(estimate, names )) +
  forest_plot <- ggplot( coefs, aes(estimate, term )) +
    geom_point(aes(), size = 2, colour = "#0C2852")+
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, size = 1, colour = "#0C2852")+
    geom_text(aes(x = Inf, label = display_ests), vjust = -0.5, colour = "black", hjust = -0.1, size = 3.5) +
    coord_cartesian(clip = "off") +
    geom_vline(xintercept = 0, linetype = 2)+
    theme_aw +
    #theme(text = element_text(size = 14))
    theme(axis.text.y = element_blank()) +
    labs(subtitle = coefs$term, x = NULL, y = NULL)

  return(forest_plot)

}


create_residuals_plot <- function(residuals_df) {
  #' This function takes a residuals_df for case and control and plots them by sample
  #' @param residuals_df, dataframe with residuals, as output by retrieve_mean_residuals_per_feature()
  #' @return boxplot of residuals per sample of a disease group

  residuals_plot <- ggplot(residuals_df, aes(x = reorder(Sample, Disease_Group),
                                             y = .resid, fill = Disease_Group)) +
    geom_boxplot(outlier.alpha = 0.2) +
    theme_aw +
    scale_fill_manual(values = named_palette) +
    scale_y_continuous(trans = "log10",
                       labels = scales::comma,
                       limits = c(1e-6, 1.1),
                       breaks = c(1e-6, 1e-4, 1e-2, 1))  +
    labs(title = "Residual values",
         y = element_blank(),
         x = element_blank())

  return(residuals_plot)

}

create_mean_residuals_plot <- function(mean_residuals_df) {
  #' This function takes a mean_residuals_df for case and control and plots them by disease group
  #' @param mean_residuals_df, dataframe with residuals, as output by retrieve_mean_residuals_per_feature()
  #' @return boxplot of residuals per sample of a disease group

  mean_residuals_plot <- ggplot(mean_residuals_df, aes(x = Disease_Group, y = mean_residuals_per_feature, fill = Disease_Group)) +
    geom_boxplot(outlier.alpha = 0.2) +
    theme_aw +
    scale_fill_manual(values = named_palette) +
    scale_y_continuous(trans = "log10",
                       labels = scales::comma,
                       limits = c(1e-6, 1.1),
                       breaks = c(1e-6, 1e-4, 1e-2, 1))  +
    labs(title = "Mean residuals per group",
         y = element_blank(),
         x = element_blank())

  return(mean_residuals_plot)
}



create_parallel_plot <- function(mean_residuals_df, disease, feature, feature_label, scale) {
  #' This function takes a mean_residuals_df for case and control and plots a parallel line plot comparing them
  #' @param residuals_df, dataframe with residuals, as output by retrieve_mean_residuals_per_feature()
  #' @param disease, disease being explored (relative to control)
  #' @param feature, reflecting the column of interest, either "locus" or "gene_name".
  #' @return parallel line plot

  feature_sym <- sym(feature)

  parallel_plot <-  mean_residuals_df %>%
    ggplot(aes(x = Disease_Group, y = mean_residuals_per_feature, group = !!feature_sym, colour = Direction)) +
    geom_path(alpha = 0.04, linewidth = 0.1) +
    theme_aw +
    # theme(legend.position = "bottom") +
    scale_colour_manual(values = named_palette,
                        labels = plot_labels) +
    scale_x_discrete(labels = plot_labels) +
    labs(title = "Mean residuals per group",
         y = str_c("Mean residuals per ", feature_label),
         x = element_blank()) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1))) # Change legend so can see the colours

  # Conditionally add y-scale transformations
  if (scale == "log10") {
      parallel_plot <- parallel_plot +
          scale_y_continuous(trans = "log10",
                             labels = scales::comma,
                             limits = c(1e-6, 1.1),
                             breaks = c(1e-6, 1e-4, 1e-2, 1))
  }  else {
      parallel_plot <- parallel_plot +
          scale_y_continuous(limits = c(NA, NA)) # Adjust as needed or remove for automatic scaling
  }

   return(parallel_plot)
}


create_sites_count_plot <- function(mean_residuals_df, feature, disease, feature_label = feature_label) {

  feature_sym <- sym(feature)

  sites_count_plot <- mean_residuals_df %>%
    group_by(!!feature_sym) %>%
    slice(1) %>%  # Select out one of the two sites per locus
    ggplot(aes(x = Direction, fill = Direction)) +
    geom_bar(stat = "count") +
    scale_fill_manual(values = named_palette,
                      labels = plot_labels) +
    scale_x_discrete(labels = plot_labels) +
    theme_aw +
    theme(#legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(title = "Mean residuals per group",
         y = str_c("Number of ",  feature_label),
         x = element_blank())
  return(sites_count_plot)
}



create_regression_summary_plot <- function(lm_summary_with_outcome,
                                           sites_subset_residuals,
                                           mean_residuals_df,
                                           disease,
                                           feature,
                                           edit_label,
                                           feature_label) {

  forest_plot <- create_forest_plot_from_lm(linear_model = lm_summary_with_outcome, disease, edit_label = edit_label, feature_label = feature_label)
  residuals_plot <- create_residuals_plot(sites_subset_residuals)
  mean_residuals_plot <- create_mean_residuals_plot(mean_residuals_df)
  parallel_plot <- create_parallel_plot(mean_residuals_df = mean_residuals_df, disease = disease, feature = feature, feature_label = feature_label)
  sites_plot <- create_sites_count_plot(mean_residuals_df = mean_residuals_df,
                                        feature = feature,
                                        disease = disease,
                                        feature_label = feature_label)
  print("start combined plot")
  combined_plot = (forest_plot | residuals_plot | mean_residuals_plot | parallel_plot | sites_plot | plot_layout(guides = "collect"))
  return(combined_plot)
}


linear_regression_and_plot <- function(sites_df,
                                       editing_metric,
                                       covariates,
                                       factor_list,
                                       outcome_measure,
                                       disease) {
  #' This function is a wrapper of other functions required to take an input table, run linear regression and output a summary plot
  #' @param sites_df, input dataframe, whether with site information or per-gene-summary-metric information
  #' @param editing_metric, a string of the editing metric (column of dataframe) to use as the outcome measure
  #' @param covariates, a vector of covariates to use in linear regression
  #' @param disease, disease/treatment being explored (relative to control)
  #' @param feature, the corresponding feature of the column of interest: either "locus" or "gene_name".
  #' @return summary plot


  if (editing_metric == "Edits"){
    feature <- "locus"
    edit_label <- "Mean editing rate"
    feature_label <- "editing sites"
  } else if (editing_metric %in% c( "edit_mean_adjusted",
                                    "edit_median_adjusted",
                                    "edit_mode_adjusted")) {
    feature <- "gene_name"
    edit_label <- editing_metric
    feature_label <- "genes"
  } else {
    print("Editing metric not recognised")
    return(NULL)
  }


  lm_results_with_outcome <- run_linear_regression_edits(sites_df = sites_df,
                                                         editing_metric,
                                                         covariates = covariates,
                                                         outcome_measure = outcome_measure)

  # Make condition for when there are no covariates (ie, when it is IPSC model without covariates)
  if (length(covariates) == 0) {
    sites_subset_residuals <- sites_df %>%
      dplyr::rename(`.resid` = editing_metric)
  } else {
    lm_results_without_outcome <- run_linear_regression_edits(sites_df = sites_df,
                                                              editing_metric,
                                                              covariates = covariates,
                                                              outcome_measure = NULL)
    sites_subset_residuals <- retrieve_lm_residuals(dataframe = sites_df,
                                                    linear_model = lm_results_without_outcome)
  }

  mean_residuals_df <- retrieve_mean_residuals_per_feature(residuals_df = sites_subset_residuals,
                                                           feature = feature,
                                                           factor_list = factor_list)

  regression_summary_plot <- create_regression_summary_plot(lm_summary = lm_results_with_outcome,
                                                            sites_subset_residuals,
                                                            mean_residuals_df,
                                                            disease,
                                                            feature = feature,
                                                            edit_label = edit_label,
                                                            feature_label = feature_label) +
    plot_annotation()
  return(regression_summary_plot)
}




regress_editing <- function(sites_df,
                            editing_metric,
                            covariates,
                            factor_list,
                            outcome_measure,
                            disease) {
  #' This function takes an input dataframe (of sites or genes) and runs linear regression, then plots a forest plot of results)
  #' @param sites_df, input dataframe, whether with site information or per-gene-summary-metric information
  #' @param editing_metric, a string of the editing metric (column of dataframe) to use as the outcome measure
  #' @param covariates, a vector of covariates to use in linear regression
  #' @param disease, disease/treatment being explored (relative to control)
  #' @param feature, the corresponding feature of the column of interest: either "locus" or "gene_name".
  #' @return summary plot

  if (editing_metric == "Edits"){
    feature <- "locus"
    edit_label <- "Mean editing rate"
    feature_label <- "editing sites"
  } else if (editing_metric %in% c( "complete_mean_adjusted",
                                    "complete_mean_adjusted",
                                    "complete_median_adjusted",
                                    "complete_mode_adjusted",
                                    "edit_summary_metric",
                                    "no_sites_adjusted",
                                    "mean_over_exon_length_adjusted",
                                    "mean_x_no_sites_adjusted",
                                    "mean_adjusted",
                                    "complete_mean_genelength_adjusted")) {
    feature <- "gene_name"
    edit_label <- editing_metric
    feature_label <- "genes"
  } else {
    print("Editing metric not recognised")
    return(NULL)
  }

  lm_results_with_outcome <- run_linear_regression_edits(sites_df = sites_df,
                                                         editing_metric,
                                                         covariates = covariates,
                                                         outcome_measure = outcome_measure)
  sites_subset_residuals_with_outcome <- retrieve_lm_residuals(dataframe = sites_df,
                                                  linear_model = lm_results_with_outcome)

  # Make condition for when there are no covariates (ie, when it is IPSC model without covariates)
  if (length(covariates) == 0) {
    sites_subset_residuals <- sites_df %>%
      dplyr::rename(`.resid` = editing_metric)
  } else {
    lm_results_without_outcome <- run_linear_regression_edits(sites_df = sites_df,
                                                              editing_metric,
                                                              covariates = covariates,
                                                              outcome_measure = NULL)
    sites_subset_residuals <- retrieve_lm_residuals(dataframe = sites_df,
                                                    linear_model = lm_results_without_outcome)
  }

  mean_residuals_df <- retrieve_mean_residuals_per_feature(residuals_df = sites_subset_residuals,
                                                           feature = feature,
                                                           factor_list = factor_list)
  lm_results <- list(lm_results_with_outcome = lm_results_with_outcome,
                     sites_subset_residuals_with_outcome = sites_subset_residuals_with_outcome,
                     sites_subset_residuals = sites_subset_residuals,
                     mean_residuals_df = mean_residuals_df)

  return(lm_results)

}


plot_residual_histogram <- function(regression_df, resid_column, fill_variable){

  histogram <- regression_df %>%
    ggplot(aes_string(x = resid_column, colour = fill_variable)) +
    geom_freqpoly(linewidth = 1.1, binwidth = 0.05) +
    scale_colour_manual(values = named_palette) +
    labs(x = "Residual", y = "Frequency") +
    theme_aw
  return(histogram)
}

compare_proportions_baseline_altered <- function(lm_results_list,
                                                 loci_consequences,
                                                 factor_list) {

  baseline_mean_resids <- retrieve_baseline_mean_residuals_per_feature(residuals_df = lm_results_list$sites_subset_residuals,
                                                                       feature = "locus",
                                                                       factor_list = factor_list) %>%
    dplyr::left_join(loci_consequences,
                     by = "locus") %>%
    dplyr::select(locus, Direction, gene_name, biotype_plots, gene_type)

  altered_mean_resids <- lm_results_list$mean_residuals_df %>%
    distinct(locus, Direction) %>%
    left_join(loci_consequences,
              by = "locus")

  consequence_plot_data <- bind_rows(baseline_mean_resids,
                                     altered_mean_resids) %>%
    dplyr::mutate(Direction = as_factor(Direction))
}


plot_compare_proportions_baseline <- function(compare_proportions_df,
                                              metric_to_plot,
                                              fill_label) {
  enquo_metric <- sym(metric_to_plot)

  compare_proportions_plot <- compare_proportions_df %>%
    drop_na(!!enquo_metric) %>%
    ggplot(aes(x=Direction, fill = fct_rev(!!enquo_metric))) +
    geom_bar(position = "fill") +
    scale_fill_paletteer_d("ggthemes::Classic_20",
                           na.value = "grey90",
                           name = fill_label,
                           labels = plot_labels) +
    theme_aw +
    labs(y = "Proportion of editing sites", x = NULL, fill = fill_label) +
    scale_x_discrete(labels = plot_labels)  +
    guides(fill = guide_legend(reverse = TRUE),
    ) +
    labs(subtitle = "Proportion change in disease")

  return(compare_proportions_plot)

}





