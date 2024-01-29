# Functions for analysing A to I editing
# Aaron Wagen
# Sept 2022
# Some functions adapted from James Bayne

library(here)
library(GenomicRanges)
library(data.table)
library(tidyverse)
library(plyranges)
library(paletteer)
library(patchwork)
library(parallel)
library(gridExtra)
library(ggtranscript)
library(cowplot)
library(viridis)


# 1.  Handling bed files -------------------------------------------------------
# 1a. Reading bed files
list_beds <- function(dir_path, ag = F) {
  #' Lists all files that end in .bed in provided directory
  #' @param dir_path, path to dir containing bed files
  #' @param ag, if the file contains '.AG.bed' (old format)
  #' @return files, list of files

    if (ag) {
  files <- files[grep('.AG.bed$', files)]

  } else {
    files <- list.files(dir_path)
    files <- files[grep('.bed$', files)]
  }
  return(files)
}

read_bed <- function(file) {
  #' Reads a single bed file as columns and outputs a dataframe
  #' @param file, full path single .bed file
  #' @return data, data frame of bed file
  data <- read_delim(file = file, delim='\t',
                     col_names=c('chromosome','start','finish','edits','z.score','strand'),
                     col_types = 'cnncnc',
                     skip=1)  # Remove first line in read_delim phase, to stop read error
  data <- data[ , c('chromosome', 'start', 'edits', 'z.score', 'strand')]
  strand <- apply(data, 1, function(x) {ifelse(x[5] == "+",1,-1)})  # reformat strand from +,1 to 1,-1
  data$Strand_JB <- strand  #add in 2 strand columns of different formats
  return(data)
}


read_bed_QC <- function(file) {
  #' Reads a single bed file from a filepath for purpose of QC script
  #' @param file, full path of a single .bed file
  #' @return data, dataframe of bed file
  data <-  read_delim(file,
                      delim = '\t',
                      col_names=c('chromosome','start','finish','edits','z.score','strand'),
                      col_types = 'cnncnc',
                      skip = 1)
  data <- data %>%
    filter(chromosome %in% c(1:22, "X", "Y", "MT")) %>%
    mutate(Sample = annotation$Sample[sample_number], # Append sample number to each row of dataframe
           locus = str_c(chromosome, start, strand, sep = "_"),
           chromosome = factor(chromosome, levels = c(1:22, "X", "Y", "MT"))
    )
  return(data)
}


read_multi_beds <- function(dir_path, files) {
  #' Read multiple bed files as specified as a list
  #' @param dir_path, path to the dir containing the bed files
  #' @param files, list of files names (not full path)
  #' @return, list of dataframes with names as the file names
  full_paths <- lapply(files, function(x) {paste0(dir_path,x)})

  multi_data <- mclapply(full_paths, function(x) {
    read_bed(x)
  }, mc.cores = CORES)
  names(multi_data) <- files
  return(multi_data)
}


read_bed_original_strand <- function(file) {
  #' Reads a single bed file as columns and outputs a dataframe
  #' @param file, full path single .bed file
  #' @return data, data frame of bed file
  data <- read_delim(file = file, delim='\t',
                     col_names=c('chromosome','start','finish','edits','z.score','strand'),
                     col_types = 'cnncnc',
                     skip=1)  # Remove first line in read_delim phase, to stop read error
  data <- data[ , c('name', 'chromosome', 'start', 'edits', 'z.score', 'strand')]
  return(data)
}

read_bed_original_strand_names <- function(dir_path, file) {
  #' Reads a single bed file as columns and outputs a dataframe
  #' @param file, full path single .bed file
  #' @param dir_path, path to directory with .bed files
  #' @return data, data frame of bed file
  full_paths <-  paste0(dir_path, "/",file)
  data <- read_delim(file = full_paths, delim='\t',
                     col_names=c('chromosome','start','finish','edits','z.score','strand'),
                     col_types = 'cnncnc',
                     skip=1)  # Remove first line in read_delim phase, to stop read error
  Sample <- sub("\\..*", "", file) #Add in sample name derived from file name
  Strand_JB <- data$strand
  data <- cbind(Sample, data, Strand_JB)  #add in 2 strand columns of different formats
  data <- data[ , c('Sample', 'chromosome', 'start', 'edits', 'z.score', 'strand', 'Strand_JB')]
  return(data)
}

read_multi_beds_original_strand_names <- function(dir_path, files) {
  #' Read multiple bed files as specified as a list
  #' @param dir_path, path to the dir containing the bed files
  #' @param files, list of files names (not full path)
  #' @return, list of dataframes with names as the file names
  #full_paths <- lapply(files, function(x) {paste0(dir_path,x)})

  multi_data <- mclapply(files, function(x){
    read_bed_original_strand_names(dir_path, x)
  }, mc.cores = CORES)
  names(multi_data) <- files
  return(multi_data)
}

read_multi_beds_original_strand <- function(dir_path, files) {
  #' Read multiple bed files as specified as a list
  #' @param dir_path, path to the dir containing the bed files
  #' @param files, list of files names (not full path)
  #' @return, list of dataframes with names as the file names
  full_paths <- lapply(files, function(x) {paste0(dir_path,x)})

  multi_data <- mclapply(full_paths, function(x) {
    read_bed_original_strand(x)
  }, mc.cores = CORES)
  names(multi_data) <- files
  return(multi_data)
}

# 1b. Processing bed files
process_bed_info <- function(data, time_point, negatives =F) {
  #' Same as \code(process_bed_info) but doesn't calculate relative site
  #' @param data, dataframe whose contents are the read in from a bed file
  #' @param time_point, used as a label but must be supplied
  #' @return, a dataframe of the same num_rows but more columns

  # Get desired columns
  data <- data[,c('chromosome','start','edits','z.score','strand', 'Strand_JB')]

  # Extract edits, read and edits/read info from 'edits' column
  diff_info <- mclapply(data$edits, function(x) {substring(x,regexpr('[^_]+$',x))},
                        mc.cores = CORES) # get rid of 'AG'
  reads <- as.vector(unlist(lapply(diff_info, function(x) {as.integer(substring(x,regexpr('[0-9]+$',x)))})))
  raw_edits <- as.vector(unlist(lapply(diff_info, function(x) {as.integer(substring(x, 1, regexpr('/',x) - 1))})))
  edits <- raw_edits / reads

  # Cbind information and rename columns and remove -ve Z scores
  time_point <- rep(time_point, dim(data)[1])
  data <- cbind(time_point, data, reads, raw_edits, edits)
  names(data) <- c('Time Point', 'Chromosome','Locus','Diff Info', 'Z score', 'strand', 'Strand_JB',
                   'Reads', 'Raw Edits', 'Edits')
  data <- data %>% mutate_at('Time Point', as.character)
  if (!negatives) data <- data[data$`Z score` > 0,]  # remove negatives
  return(data)
}

process_bed_info_names <- function(data, annotations, negatives =F) {
  #' Same as \code(process_bed_info)
  #' @param data, dataframe whose contents are the read in from a bed file
  #' @param annotations, dataframe containing metadata (eg columns of timepoint,condition etc) with a column labelled with file_name
  #' @return, a dataframe of the same num_rows but more columns

  # Get desired columns
  data <- data[,c('Sample', 'chromosome','start','edits','z.score','strand')]

  # Extract edits, read and edits/read info from 'edits' column
  diff_info <- parallel::mclapply(data$edits, function(x) {substring(x,regexpr('[^_]+$',x))},
                        mc.cores = CORES) # get rid of 'AG'
  reads <- as.vector(unlist(lapply(diff_info, function(x) {as.integer(substring(x,regexpr('[0-9]+$',x)))})))
  raw_edits <- as.vector(unlist(lapply(diff_info, function(x) {as.integer(substring(x, 1, regexpr('/',x) - 1))})))
  Edits <- raw_edits / reads
  data <- cbind(data, reads, raw_edits, Edits)

  # Add metadata from annotations dataframe
  data <- left_join(data, annotations, by = "Sample")

  # Remove -ve Z scores
 data <- data %>%
   filter('Z-score' >= 0) # Remove zeroes
#  data <- data %>% mutate_at('Time Point', as.character)
#  if (!negatives) data <- data[data$`Z-score` > 0,]  # remove negatives
  return(data)
}


process_bed_QC <- function(data, annotations, negatives =F) {
  #' @param data, dataframe whose contents are the read in from a bed file for purpose of QC script
  #' @param annotations, dataframe with column for file_path, column called Sample with unique name for each sample, and other annotations as required
  #' @return, a dataframe of the same num_rows but more columns

  # Get desired columns
  data <- data[,c('Sample', 'chromosome','start','edits','z.score','strand')]

  # Extract edits, read and edits/read info from 'edits' column
  diff_info <- parallel::mclapply(data$edits, function(x) {substring(x,regexpr('[^_]+$',x))},
                                  mc.cores = CORES) # get rid of 'AG'
  reads <- as.vector(unlist(lapply(diff_info, function(x) {as.integer(substring(x,regexpr('[0-9]+$',x)))})))
  raw_edits <- as.vector(unlist(lapply(diff_info, function(x) {as.integer(substring(x, 1, regexpr('/',x) - 1))})))
  Edits <- raw_edits / reads
  data <- cbind(data, reads, raw_edits, Edits)

  # Add metadata from annotations dataframe
  data <- left_join(data, annotations, by = "Sample")

  # Remove -ve Z scores
  data <- data %>%
    filter('Z-score' >= 0) # Remove zeroes
  #  data <- data %>% mutate_at('Time Point', as.character)
  #  if (!negatives) data <- data[data$`Z-score` > 0,]  # remove negatives
  return(data)
}


process_multi_beds_names <- function(data, annotations) {
  #' Takes in list of pre_read bed files and processes all of them. outputs as a list
  #' @param data, output of \code(read_multi_beds) e.g. list of bed file data
  #' @param annotations, dataframe with metadata, and 1 column named 'Sample' with sample name
  #' @return list of data frames containing processed bed data

  multi_bed <- mapply(function(x) {process_bed_info_names(x, annotations)},
                      data, SIMPLIFY = F)

  return(multi_bed)
}

process_multi_beds <- function(data, time_points, negatives = F) {
  #' Takes in list of pre_read bed files and processes all of them. outputs as a list
  #' @param data, output of \code(read_multi_beds) e.g. list of bed file data
  #' @param time_points, timepoints for each dataframe in the list, input as a vector
  #' @return list of data frames containing processed bed data

  multi_bed <- mapply(function(x, y) {as.list(process_bed_info(x, y, negatives))},
                      data, time_points, SIMPLIFY = F)
  multi_bed <- lapply(multi_bed, function(x) {
    x <- data.frame(x, stringsAsFactors = F)
    names(x) <- c('Time Point', 'Chromosome','Locus','Diff Info', 'Z score', 'Strand_JB',
                  'Reads', 'Raw Edits', 'Edits')
    return(x)})
  return(multi_bed)
}

process_bed_info_conditions <- function(data, condition, time_point, negatives =F) {
  #' Same as \code(process_bed_info) but doesn't calculate relative site
  #' @param data, dataframe whose contents are the read in from a bed file
  #' @param time_point, used as a label but must be supplied
  #' @param condition, used as a label but must be supplied
  #' @return, a dataframe of the same num_rows but more columns

  # Get desired columns
  data <- data[,c('chromosome','start','edits','z.score','strand', 'Strand_JB')]

  # Extract edits, read and edits/read info from 'edits' column
  diff_info <- mclapply(data$edits, function(x) {substring(x,regexpr('[^_]+$',x))},
                        mc.cores = CORES) # get rid of 'AG'
  reads <- as.vector(unlist(lapply(diff_info, function(x) {as.integer(substring(x,regexpr('[0-9]+$',x)))})))
  raw_edits <- as.vector(unlist(lapply(diff_info, function(x) {as.integer(substring(x, 1, regexpr('/',x) - 1))})))
  edits <- raw_edits / reads

  # Cbind information and rename columns and remove -ve Z scores
  time_point <- rep(time_point, dim(data)[1])
  condition <- rep(condition, dim(data)[1])
  data <- cbind(condition, time_point, data, reads, raw_edits, edits)
  names(data) <- c('Condition', 'Time Point', 'Chromosome','Locus','Diff Info', 'Z score', 'strand', 'Strand_JB',
                   'Reads', 'Raw Edits', 'Edits')
  data <- data %>% mutate_at('Time Point', as.character)
  if (!negatives) data <- data[data$`Z score` > 0,]  # remove negatives
  return(data)
}

process_multi_beds_conditions <- function(data, conditions, time_points, negatives = F) {
  #' Takes in list of pre_read bed files and processes all of them. outputs as a list
  #' @param data, output of \code(read_multi_beds) e.g. list of bed file data
  #' @param conditions, conditions for each dataframe in the list, input as a vector
  #' @param time_points, timepoints for each dataframe in the list, input as a vector
  #' @return list of data frames containing processed bed data

  multi_bed <- mapply(function(x, y, z) {as.list(process_bed_info_conditions(x, y, z, negatives))},
                      data, conditions, time_points, SIMPLIFY = F)
  multi_bed <- lapply(multi_bed, function(x) {
    x <- data.frame(x, stringsAsFactors = F)
    names(x) <- c('Conditions', 'Time Point', 'Chromosome','Locus','Diff Info', 'Z score', 'strand', 'Strand_JB',
                  'Reads', 'Raw Edits', 'Edits')
    return(x)})
  return(multi_bed)
}



# 1c. Processing differential beds

# Read differential bed
read_dif_bed <- function(bed_dir, bed_file) {
  #' Reads a single bed file as columns and outputs a dataframe, adding .bed title as a column
  #' @param bed_dir, path to directory containing differential beds
  #' @param bed_file, the file to read in
  #' @return data, data frame of bed file
  data <- read_delim(file = str_c(bed_dir, "/", bed_file),
                     delim='\t',
                     col_names=c('chromosome','start','finish','edits','z.score','strand'),
                     col_types = 'cnncnc',
                     skip=1)  # Remove first line in read_delim phase, to stop read error
  data <- data[ , c('chromosome', 'start', 'edits', 'z.score', 'strand')]
  strand <- apply(data, 1, function(x) {ifelse(x[5] == "+",1,-1)})  # reformat strand from +,1 to 1,-1
  data$Strand_JB <- strand  #add in 2 strand columns of different formats
  data$SampleName <- str_replace(bed_file, "\\..*","")
  return(data)
}


# Read multiple differential beds
read_multi_dif_beds <- function(bed_dir, bed_files) {
  #' Read multiple bed files as specified as a list
  #' @param bed_dir, path to the dir containing the bed files
  #' @param bed_files, list of files names (not full path)
  #' @return, list of dataframes with names as the file names

  multi_data <- mclapply(bed_files, function(x) {
    read_dif_bed(bed_dir, x)
  }, mc.cores=CORES)
  names(multi_data) <- str_replace(bed_files, "\\..*","")
  return(multi_data)
}



# Process differential bed
process_dif_bed <- function(read_data, z_score_cuttoff=1.96) {
  #' Processes output of red_dif_bed, separating out columns of edits in treated and untreated, and filtering by submitted z score
  #' @param read_data, data output from red_dif_bed
  #' @param z_score_cuttoff, cuttoff to use for filtering significant differential editing
  #' @return data, data frame of bed file
  data <- read_data %>%
    dplyr::filter(z.score >=z_score_cuttoff | z.score <= -z_score_cuttoff) %>%
    tidyr::separate(edits, into = c("AG", "Untreated_edit", "Treated_edit"), sep="_") %>%
    dplyr::select(-AG) %>%
    tidyr::separate(Untreated_edit, into = c("untreated_raw_edits", "untreated_reads"), sep="/") %>%
    tidyr::separate(Treated_edit, into = c("treated_raw_edits", "treated_reads"), sep="/") %>%
    mutate(across(.cols=c(untreated_raw_edits, untreated_reads, treated_raw_edits, treated_reads), .fns=as.numeric)) %>%
    mutate(untreated_edits = untreated_raw_edits/untreated_reads,
           treated_edits = treated_raw_edits/treated_reads,
           treat_untreat_ratio = treated_edits/untreated_edits,
           locus = str_c(chromosome, "_", start, "_", strand)) %>%
    relocate(chromosome, start, strand, locus, untreated_edits, treated_edits, treat_untreat_ratio, z.score, SampleName)
  return(data)
}

#Process multiple differential beds
process_multi_dif_beds <- function(data, z_score_cuttoff=1.96) {
  #' Takes in list of pre_read bed files and processes all of them. outputs as a list
  #' @param data, output of \code(read_multi_beds) e.g. list of bed file data
  #' @param z_score_cuttoff, cuttoff to use for filtering significant differential editing
  #' @param annotations, dataframe with metadata, and 1 column named 'Sample' with sample name
  #' @return list of data frames containing processed bed data

  multi_bed <- mapply(function(x) {process_dif_bed(x, z_score_cuttoff)},
                      data, SIMPLIFY = F)

  return(multi_bed)
}



# 2. Annotating edits and filtering --------------------------------------------------------

annotate_edits <- function(data_range, reduced_ensembl_range, alu_annotations_range) {
  #' Takes data in genomic range form and annotates with ensembl gtf and alu_annotations supplied. Note this is without filtering.
  #' It does this via a left join with the data range, then selecting the first match in the case of multiple mapping
  #' @param data, data in genomic ranges format, in this case editing data
  #' @param reduced_ensembl_range, ensembl gtf previously processed, with exons collapsed, and intronic and intergenic regions
  #' @param alu_annotations_range, repeat annotations previously processed in genomic range form
  #' @return list of 3 dataframes: all annotated edits, edits with a z score filter>1.96, z score filter + minimum read depth >20
  annotated_edits <-
    plyranges::join_overlap_left_directed(x = data_range,
                                          y = reduced_ensembl_range) %>%
    plyranges::mutate(gene_type_embl = replace_na(gene_type_embl, "intergenic"))

  annotated_edits <- unique(annotated_edits) # Randomly select the first match with the edits (not the second)
  annotated_edits <- plyranges::join_overlap_left_directed(x = annotated_edits,
                                                           y = alu_annotations_range)
  annotated_edits <- unique(annotated_edits) %>%
    as_tibble() %>%
    mutate(across(c(gene_biotype_embl, gene_type_embl, repeat_type, repeat_class), factor))

  return(annotated_edits)

}


annotate_edits_gencode <- function(data_range, reduced_gencode_range, alu_annotations_range, all_genes) {
  #' Takes data in genomic range form and annotates with ensembl gtf and alu_annotations supplied. Note this is without filtering.
  #' It does this via a left join with the data range, then selecting the first match in the case of multiple mapping
  #' Note: this is different from annotate_edits in that the ensembl range has different column names
  #' Note: A filtering step is required here where multiple genes map to an editing site.
  #'    In this case any gene in the genes of interest is selected, then randomly selecting any other gene
  #' @param data, data in genomic ranges format, in this case editing data
  #' @param reduced_gencode_range, gencode gtf previously processed, with exons collapsed, and intronic and intergenic regions
  #' @param alu_annotations_range, repeat annotations previously processed in genomic range form
  #' @param all_genes, a vector of all genes of interest to preference in filtering
  #' @return list of 3 dataframes: all annotated edits, edits with a z score filter>1.96, z score filter + minimum read depth >20
  annotated_edits <-
    plyranges::join_overlap_left_directed(x = data_range,
                                          y = reduced_gencode_range) %>%
    plyranges::mutate(gene_type = replace_na(gene_type, "intergenic")) %>%
    plyranges::join_overlap_left_directed(.,
                                          y = alu_annotations_range)

  annotated_edits <- annotated_edits %>%
    as_tibble() %>%
    dplyr::mutate(gene_of_interest = if_else(gene_name %in% all_genes, "yes", "no")) %>%
    dplyr::arrange(seqnames, start, gene_of_interest) %>%
    dplyr::distinct(seqnames, start, gene_of_interest, .keep_all = TRUE) %>%
    dplyr::mutate(across(c(gene_biotype, gene_type, repeat_type, repeat_class), factor),
                  gene_id = str_replace(gene_id, '\\..*$', ""))  %>%  #Remove the sub-number from gene-id
    dplyr::relocate(seqnames, start, mean_editing, gene_name, gene_id, gene_type, gene_biotype)

  return(annotated_edits)

}



filter_edits_qc <- function(edits_df) {
  #' Takes a dataframe of edits and filters by z score (on column 'zscore') and then by total reads (on column 'reads')
  #' @param edits_df, dataframe containing editing information including a column 'zscore' and a column 'reads'
  #' @return list of 3 dataframes: all annotated edits, edits with a z score filter>1.96, z score filter + minimum read depth >20
  edits_df <- edits_df %>%
    mutate(filter = "all edits",
           locus = str_c(seqnames, start, strand, sep = "_"))
  z_score_filter_edits <- edits_df %>%
    dplyr::filter(`z.score` >= 1.96) %>%
    mutate(filter = "zscore")

  z_score_read_number_filter <- edits_df %>%
    dplyr::filter(`z.score` >1.96 & reads >10) %>%
    mutate(filter = "zscore_reads")

  return(list(edits_df = edits_df,
              z_score_edits = z_score_filter_edits,
              z_score_reads_edits = z_score_read_number_filter))
}





filter_duplicates <- function(dataframe, identifier_column, filtering_column, filtering_levels) {
  #' This function takes a dataframe, identifies duplicates within the dataframe and filters those duplicates by the filtering levels of a certain column
  #' @param dataframe, dataframe with a column by which duplicates are identified (identifier column), and one by which they are selected (filtering column)
  #' @param identifier_column, the column by which duplicates are identified. In the case of genomics it is often locus: <chr_start_strand>
  #' @param filtering_column, the column by which duplicates are selected. Should be a factor column
  #' @param filtering_levels, the levels by which the factors in filtering_column should be ordered, from most preferenced to least
  #' @return dataframe where duplicates have been removed as per factor levels above

  identifier <- dplyr::ensym(identifier_column)
  filterer <- dplyr::ensym(filtering_column)

  dataframe[[filterer]] <- factor(dataframe[[filterer]], levels = filtering_levels)
  filtered_from_duplicates <- dataframe %>%
    dplyr::arrange(!!identifier, !!filterer) %>%
    distinct(!!identifier, .keep_all = T)

  return(filtered_from_duplicates)
}


# 3. Munging data ------------------------------------------------------------


make_genomic_range <- function(data) {
  #' Takes dataframe and converts to genomic ranges format
  #' @param data, output dataframe of processed_multi_beds_conditions
  #' @return granges object with metadata retained
  processed_range <- GenomicRanges::makeGRangesFromDataFrame(
    df = data,
    seqnames.field = "Chromosome",
    start.field = "Start",
    end.field = "Start",
    strand.field = "Strand",
    keep.extra.columns = T,
    ignore.strand = F)
  return(processed_range)
}

make_multi_genomic_range <- function(processed_data_list) {
  #' @param processed_data_list, list of dataframes output by process_multi-beds_conditions
  #' @return list of genomic ranges
  processed_ranges <- mapply(function(x) {make_genomic_range(x)},
                             processed_data_list)
}



to_intron <- function(exons, group_var = NULL) {
  .check_coord_object(exons)
  .check_group_var(exons, group_var)

  # TODO - switch this to using GenomicRanges::gaps()?

  if (!is.null(group_var)) {
    exons <- exons %>% dplyr::group_by_at(.vars = group_var)
  }

  # make sure exons are arranged by coord, so that dplyr::lag works correctly
  exons <- exons %>%
    dplyr::arrange(start, end)

  # obtain intron start and ends
  introns <- exons %>%
    dplyr::mutate(
      intron_start := dplyr::lag(end),
      intron_end := start,
      type = "intron"
    ) %>%
    dplyr::select(-start, -end)

  # remove the introduced artifact NAs
  introns <- introns %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(intron_start) & !is.na(intron_end))

  # filter out introns with a width of 1, this should only happen when
  # utrs are included and are directly adjacent to end of cds
  introns <- introns %>% dplyr::filter(abs(intron_end - intron_start) != 1)

  introns <- introns %>% dplyr::rename(start = intron_start, end = intron_end)

  return(introns)
}


extract_nested_df <- function(summary_list, sublist, anotherlist) {
  #' Takes the summary list from samples, and extracts a subsublist of choice, convering it into a single dataframe
  #' @param summary_list, the output of the QC script, the first level will have list of samplenames, then there will be further nested lists of interest
  #' @param sublist, the sublist to extact initially, eg, 'All_edits_lists', or 'Alu_edits_lists
  #' @param subsublist, the subsublist to extract, eg, 'z_score_edits' or 'z_score_reads_edits'
  #' @return single dataframe or rowbinded dataframes
  category_data_frame <- summary_list %>%
    purrr::map(., sublist)  %>%
    purrr::map_dfr(., anotherlist)

  return(category_data_frame)
}


# 3.  Filtering  --------------------------------------------------------




# Cluster profile enrichment analysis -------------------------------------


GOenrichment <- function(filtered_df) {
  #' Takes a dataframe that has been filtered to differentially edited sites of interest and runs GOenrich referencing the gene_universe supplied. Must have a column for `gene_id` as per ensembl coding
  #' @param filtered_df, as dataframe filtered to genes of interest. Must contain a column called `gene_id` as per ensembl genes.
  #' @return, dot-plot of enriched go terms

  target_gene_list <- filtered_df$gene_id %>%
    unique()

  target_enrichment <- enrichGO(gene = target_gene_list ,
                                universe = gene_universe ,
                                OrgDb = "org.Hs.eg.db",
                                ont = "BP" ,
                                keyType = "ENSEMBL",
                                readable = TRUE,
                                pvalueCutoff = 1,
                                qvalueCutoff = 1)

  plot <- dotplot(target_enrichment)
  return(plot)

}


# Explore data/tables -----------------------------------------------------

summary_by_variable <- function(dataframe, sample_metric = Sample, variable = NULL){
  #' Makes a summary table by count and percent of a variable, relative to the Sample)
  #' @param dataframe, the dataframe to explore
  #' @param sample_metric, the samples comparison to look at
  #' @param variable, the variable to compare across samples
  #' @return 2 tables: one is raw numbers, one as percentage

  raw_table_withNA <- dataframe %>%
    dplyr::select({{sample_metric}}, {{variable}}) %>%
    base::table(useNA = "always")

  raw_table <- dataframe %>%
    dplyr::select({{sample_metric}}, {{variable}}) %>%
    base::table()

  prop_table_withNA <-  dataframe %>%
    dplyr::select({{sample_metric}}, {{variable}}) %>%
    base::table(useNA = "always") %>%
    prop.table(margin=1) %>% `*`(100) %>% round(2)

  prop_table <- dataframe %>%
    dplyr::select({{sample_metric}}, {{variable}}) %>%
    base::table() %>%
    prop.table(margin=1) %>% `*`(100) %>% round(2)

  sum_list <- list(Count_with_NA = raw_table_withNA, Count = raw_table, Percent_withNA = prop_table_withNA, Percent = prop_table)

  return(sum_list)

}

print_list_of_tables_for_viewing <- function(sum_list) {
  # Takes a list of tables (eg output by summary_by_variable) and creates kable style tables for viewing while working on the document

  for (i in seq_along(sum_list)){
    kable(sum_list[[i]], caption = names(sum_list)[i], longtable = FALSE) %>%
      kable_styling() %>%
      print()
  }
}

print_list_of_tables_for_knitting <- function(sum_list) {
  for (i in seq_along(sum_list)){
    DT::datatable(sum_list[[i]],
                  rownames = FALSE,
                  options = list(scrollX = TRUE),
                  class = 'white-space: nowrap') %>%
      print()
  }
}

chi_square_test_genic_ratio <- function(dataframe, culture_type) {
  #' Function that takes genic location of a variable, and performs the chisquare ratio rest of exon relative to introns and intergenic editing sites
  dataframe %>%
    dplyr::filter(culture == {{culture_type}}) %>%
    dplyr::select(Sample, gene_type) %>%
    dplyr::mutate(gene_type = fct_recode(gene_type, exon = "exon", other = "intron", other = "intergenic")) %>%
    droplevels() %>%
    table() %>%
    chisq.test()
}

chi_square_factor_test <- function(dataframe, culture_type, variable_of_interest, factor_of_interest) {
  #' Function that takes genic location of a variable, and performs the chisquare ratio comparing a factor of interest in a variable
  #' It will take a jacusa detect dataframe, allow you to pick a culture of interest, a variable of interst, and a factor of interest, where it will refactor all other variables to 'other'.
  #' It will then run the chisquare test
  #' @param dataframe
  #' @param culture_type write with ""
  #' @param variable_of_interest write **without** ""
  #' @param factor_of_interest write with ""
  dataframe %>%
    dplyr::filter(culture == {{culture_type}}) %>%
    dplyr::select(Sample, {{variable_of_interest}}) %>%
    dplyr::mutate("{{variable_of_interest}}" := fct_other({{variable_of_interest}}, keep = {{factor_of_interest}})) %>% # Use the := to define a variable
    droplevels() %>%
    table() %>%
    chisq.test()
}


chi_square_factor_test_differential <- function(dataframe, culture_type, variable_of_interest, factor_of_interest) {
  #' Function that takes genic location of a variable, and performs the chisquare ratio comparing a factor of interest in a variable
  #' It will take a jacusa differential dataframe, allow you to pick a culture of interest, a variable of interst, and a factor of interest, where it will refactor all other variables to 'other'.
  #' It will then run the chisquare test
  #' @param dataframe
  #' @param culture_type write with ""
  #' @param variable_of_interest write **without** ""
  #' @param factor_of_interest write with ""
  dataframe %>%
    dplyr::filter(culture == {{culture_type}}) %>%
    dplyr::select(Sample_change, {{variable_of_interest}}) %>%
    dplyr::mutate("{{variable_of_interest}}" := fct_other({{variable_of_interest}}, keep = {{factor_of_interest}})) %>% # Use the := to define a variable
    droplevels() %>%
    table() %>%
    chisq.test()
}




# 4. Plots -------------------------------------------------------------------
barplot_by_gene_annotation <- function(list_of_annotated_edits, proportional = c("stack", "fill"), comparison = gene_type) {
  #' Displays a barchart of the number of editing events in each filtered state from annotated edits, coloured by genic annotations
  #' @param list_of_annotated_edits, the list of 3 annotated edits produced by the function annotate_edits
  #' @param proportional, can be 'stack' to show the difference in number of edits, or 'fill' to see the proportional difference
  #' @return, barchart comparing each of the lists
  ggplot() +
    geom_bar(data = annotated_edits_list[["edits_df"]], aes(x="all edits", fill = {{comparison}}), position = proportional) +
    geom_bar(data = annotated_edits_list[["z_score_edits"]], aes(x="z filter", fill = {{comparison}}), position = proportional) +
    geom_bar(data = annotated_edits_list[["z_score_reads_edits"]], aes(x="z filter + reads", fill = {{comparison}}), position = proportional) +
    theme_bw() +
    scale_fill_paletteer_d("ggthemes::Classic_20",
                           na.value = "grey90",
                           direction = 1) +
    #  scale_fill_brewer(palette = "Dark2") +
    labs(x = NULL, y = "Number of editing events") +
    theme(legend.title = element_blank())
}

barplot_by_gene_annotation_dropNA <- function(list_of_annotated_edits, proportional = c("stack", "fill"), comparison = gene_type) {
  #' Displays a barchart of the number of editing events in each filtered state from annotated edits, coloured by genic annotations
  #' @param list_of_annotated_edits, the list of 3 annotated edits produced by the function annotate_edits
  #' @param proportional, can be 'stack' to show the difference in number of edits, or 'fill' to see the proportional difference
  #' @return, barchart comparing each of the lists
  ggplot() +
    geom_bar(data = annotated_edits_list[["edits_df"]] %>% drop_na({{comparison}}), aes(x="all edits", fill = {{comparison}}), position = proportional) +
    geom_bar(data = annotated_edits_list[["edits_df"]] %>% drop_na({{comparison}}), aes(x="z filter", fill = {{comparison}}), position = proportional) +
    geom_bar(data = annotated_edits_list[["edits_df"]] %>% drop_na({{comparison}}), aes(x="z filter + reads", fill = {{comparison}}), position = proportional) +
    theme_bw() +
    scale_fill_paletteer_d("ggthemes::Classic_20",
                           na.value = "grey90",
                           direction = 1) +
    #  scale_fill_brewer(palette = "Dark2") +
    labs(x = NULL, y = "Number of editing events") +
    theme(legend.title = element_blank())
}



## Fix this plot
violin_plot_editing_rate <- function(list_of_annotated_edits) {
  #' Displays a violin plot of the editing rate for each annotated edits list
  #' @param list_of_annotated_edits, the list of 3 annotated edits produced by the function annotate_edits
  #' @return, violin plot comparing each of the lists
  ggplot() +
    geom_violin(data = annotated_edits_list[["edits_df"]] %>%
                  dplyr::filter(Edits !=1), aes(x=filter, y=Edits, fill=Sample),
                draw_quantiles = c(0.25, 0.5, 0.75),
                trim=T) +
    geom_violin(data = annotated_edits_list[["z_score_edits"]], aes(x=filter, y=Edits, fill= Sample),
                draw_quantiles = c(0.25, 0.5, 0.75),
                trim=T) +
    geom_violin(data = annotated_edits_list[["z_score_reads_edits"]], aes(x=filter, y=Edits, fill=Sample),
                draw_quantiles = c(0.25, 0.5, 0.75),
                trim=T) +
    theme_bw() +
    scale_fill_paletteer_d("ggthemes::Classic_20",
                           na.value = "grey90",
                           direction = -1) +
    #  scale_fill_brewer(palette = "Dark2") +
    labs(x = NULL, y = "Editing rate per site") +
    theme(legend.title = element_blank(),
          legend.position = "none")
}


violin_plot_editing_rate_trimmed_log <- function(list_of_annotated_edits) {
  #' Displays a violin plot of the editing rate for each annotated edits list
  #' @param list_of_annotated_edits, the list of 3 annotated edits produced by the function annotate_edits
  #' @return, violin plot comparing each of the lists
  ggplot() +
    geom_violin(data = annotated_edits_list[["edits_df"]] %>%
                  dplyr::filter(Edits <= 0.99 & Edits >=0.01),
                aes(x=filter, y=Edits, fill=Sample),
                draw_quantiles = c(0.25, 0.5, 0.75),
                trim=T) +
    geom_violin(data = annotated_edits_list[["z_score_edits"]] %>%
                  dplyr::filter(Edits <= 0.99 & Edits >=0.01),
                aes(x=filter, y=Edits, fill= Sample),
                draw_quantiles = c(0.25, 0.5, 0.75),
                trim=T) +
    geom_violin(data = annotated_edits_list[["z_score_reads_edits"]] %>%
                  dplyr::filter(Edits <= 0.99 & Edits >=0.01),
                aes(x=filter, y=Edits, fill=Sample),
                draw_quantiles = c(0.25, 0.5, 0.75),
                trim=T) +
    theme_bw() +
    scale_fill_paletteer_d("ggthemes::Classic_20",
                           na.value = "grey90",
                           direction = -1) +
    scale_y_continuous(trans='log10') +
    #  scale_fill_brewer(palette = "Dark2") +
    labs(x = NULL, y = "Editing rate per site") +
    theme(legend.title = element_blank(),
          legend.position = "none")
}



compare_plot <- function(comparison_df, proportional = c("stack", "fill"), comparison, plot_label) {


  ggplot(comparison_df, aes(x = treatment, fill = fct_rev({{comparison}}))) +
    geom_bar(position = proportional) +
    scale_fill_paletteer_d("ggthemes::Classic_20",
                           na.value = "grey90",
                           labels = {{plot_label}}) +
    theme_aw +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y = "Number of editing sites", x = NULL) +
    scale_x_discrete(labels = {{plot_label}}) +
    guides(fill = guide_legend(reverse = TRUE)) +
    facet_wrap(~ culture,
               labeller = plot_labeller)
}

compare_plot_no_facet <- function(comparison_df, proportional = c("stack", "fill"), comparison, plot_label) {


   ggplot(comparison_df, aes(x = Sample,fill = fct_rev({{comparison}}))) +
    geom_bar(position = proportional) +
    scale_fill_paletteer_d("ggthemes::Classic_20",
                           na.value = "grey90",
                           labels = {{plot_label}}) +
    theme_aw +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y = "Number of editing sites", x = NULL) +
    scale_x_discrete(labels = {{plot_label}}) +
    guides(fill = guide_legend(reverse = TRUE))
}



compare_plot_dropNa <- function(comparison_df, proportional = c("stack", "fill"), comparison, plot_label) {
  ggplot(comparison_df %>% drop_na({{comparison}}), aes(x = Sample, fill = {{comparison}} )) + #fct_rev({{comparison}}))) +
    geom_bar(position = proportional) +
    scale_fill_paletteer_d("ggthemes::Classic_20",
                           na.value = "grey90",
                           labels = {{plot_label}}) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y = "Number of editing sites", x = NULL) +
    scale_x_discrete(labels = {{plot_label}}) +
    guides(fill = guide_legend(reverse = TRUE))
}

compare_violin_plots <- function(comparison_df, plot_label, palette_colours) {
  ggplot(comparison_df, aes(x=treatment, y=Edits, fill=Sample)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),
                trim=T) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y = "Editing rate per site", x = NULL) +
    scale_fill_manual(values = palette_colours,
                           labels = {{plot_label}}) +
    scale_x_discrete(labels = {{plot_label}})  +
    facet_wrap(~ culture,
               labeller = plot_labeller)
}

compare_violin_plots_trimmed <- function(comparison_df, plot_label, palette_colours) {
  ggplot(comparison_df %>%
           dplyr::filter(Edits <= 0.99 & Edits >=0.01),
         aes(x=treatment, y=Edits, fill=Sample)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),
                trim=T) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y = "Editing rate per site", x = NULL) +
    scale_y_continuous(trans='log10') +
    scale_fill_manual(values = palette_colours,
                      labels = {{plot_label}}) +
    scale_x_discrete(labels = {{plot_label}}) +
    facet_wrap(~ culture,
               labeller = plot_labeller)
}






# Helper Functions --------------------------------------------------------

getmode <- function(v) {
  #' Calculates the mode value of a vector
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

mean_length_of_files <- function(file_paths) {
  #' Takes a list of file paths and finds the mean length of them
  number_of_edits <- list() # Create empty list

  file_content <- lapply(file_paths, read.delim, header=F, quote="\"") # Create list of the filecontents

  for (i in 1:length(file_paths)){
    number_of_edits[[i]] <- length(file_content[[i]]$V1) # Find length of files

  }

  mean_length <- mean(unlist(number_of_edits))
  return(mean_length)
}


`%notin%` <- function(x,y) !(x %in% y)

# testing -----------------------------------------------------------------
# load bed data
# dir_path <- '../210825_Rodriques_Age-by-Editing_v210914/jacusa_differential/'
# files <- list_beds(dir_path)[1]
# multi <- read_multi_beds(dir_path, files)
# processed <- process_bed_info(data = multi[[1]], 'T2')
# processed_m <- process_multi_beds(data = multi, time_points = c('T2'))
# site <- site_to_id(bed_data[[1]][1,], gtf)
# data <- bed_data[[1]][1:3,]
# genes <- gtf$`Gene ID`[1:100]
# system.time(sites <- sites_for_all_timepoints(genes, bed_data, gtf))
# progressive_sites <- progressive_edits(all_sites)
# data <- all_sites[1:14, ]
# s <- relative_loci(data)
# data <- all_sites[all_sites$`Gene ID` == "ENSG00000187608",]
# data <- all_sites[1:30,]
# s <- relative_loci(data)




