#### James Bayne November 2021 ####
# Contains functions for handling the data for the time stamps project #
# Adapted by Aaron Wagen #

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
  full_paths <-  paste0(dir_path,file)
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
  #' Same as \code(process_bed_info) but doesn't calculate relative site
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



relative_loci <- function(data) {
  #' For an annotated bed data.frame, gives each of the sites a relative locus
  #' within its gene
  #' @param data, data.frame containing information from bed and gtf
  #' @return same data.frame as went in but with Site Locus fixed

  # Get unique genes
  genes <- unique(data$`Gene ID`)

  # Assign every site its index within the gene
  idxs <- lapply(genes, function(x) {
    data.frame('Site ID' = str_sort(unique(data[data$`Gene Name` == x,'Site ID'])),
               'Site Locus' = match(unique(data[data$`Gene Name` == x, 'Site ID']), str_sort(unique(data[data$`Gene Name` == x, 'Site ID']))),
               stringsAsFactors = FALSE)
  }) %>% rbindlist
  names(idxs) <- c('Site ID', 'Site Locus')

  # Add indexes to dataframe
  left_join(data, idxs)
}

progressive_edits <- function(data, labels = c()) {
  #' Subsets a dataframe of bed info so that only sites that have editing in
  #' all timepoints are included
  #' @param data, dataframe to be subset
  #' @param labels, optional, function will try figure out timepoints from
  #' the values in 'Time Points' column but if the data.frame is very small
  #' this may not have datapoints for all the desired time points, in which
  #' case they should be supplied as a list

  # Get or set time point labels
  time_points <- unique(data$`Time Point`)
  if (!is_empty(labels)) time_points <- labels

  # split into chromosomes
  chromosomes <- mclapply(unique(data$Chromosome), function(x) {
    data.frame(data[data$Chromosome == x,])
  }, mc.cores = CORES)

  # function to be applied to each chromosome in turn
  ind_chrom <- function(chrom) {
    #' Takes a sites data.frame for a single chromosome and returns
    #' progressivly edited sites
    #' @param chrom, single chromosome data.frame
    uniqs <- unique(chrom$Locus)
    progs <- lapply(uniqs, function(x) {
      if (sum(chrom[chrom$Locus == x, "Locus"]) == length(time_points)*x) x})
    not_nulls <- sapply(progs, function(x) !is.null(x))
    progs <- progs[not_nulls]  # get rid of NULLs
    chrom <- chrom[chrom$Locus %in% progs, ]
    return(chrom)
  }

  # get the progressive loci for every chromosome and rbind
  progressive_loci <- lapply(chromosomes, ind_chrom)
  progressive_loci <- rbindlist(progressive_loci)
  names(progressive_loci) <- names(data)  # clean column names
  return(progressive_loci)
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
           locus = str_c(chromosome, "_", start)) %>%
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



# 2.  Ensembl annotation of sites -----------------------------------------
# 2a. Extract ensembl information from gtf
loci_from_gtf <- function(annotation, feature = 'gene') {
  #' Function that takes a GTF file, pulls out every ensembl gene and
  #' @param annotation, preloaded gtf file
  #' @param feature, what to extract e.g. gene, three_prime_utr
  #' @returns a dataframe with columns c(chromosome, start, finish, strand, ensembl_id, gene_name) with
  #' each row as a unique gene.

  annotation <- annotation[annotation$type == feature, ]  # get unique genes rows
  genes <- annotation[, c("seqid", "start", "end", "strand", "gene_id", "gene_name")]

  genes$strand <- sapply(genes$strand, function(x) {ifelse(x == "+",1,-1)})  # reformat strand from +,1 to 1,-1

  # rename and convert from factor to character
  names(genes) <- c('Chromosome', 'Start', 'Finish', 'Strand', 'Gene ID','Gene Name')
  genes <- genes %>% mutate_at('Chromosome', as.character)
  genes <- genes %>% mutate_at('Gene ID', as.character)
  genes <- genes %>% mutate_at('Gene Name', as.character)
  return(genes)
}

ends_transcripts <- function(gtf=raw_gtf) {
  #' Function that takes a local GTF annotation file, pulls out every ensembl gene and
  #' @param gtf, gff/gtf input
  #' @returns a dataframe with columns c(chromosome, start, finish, strand, ensembl_id, gene_name) with
  #' each row as a unique gene.

  # subset to columns of interest
  gtf <- gtf[!is.na(gtf$transcript_id),c('seqid','start','end','strand','transcript_id','transcript_biotype','transcript_name',
                                         'gene_id','gene_name')]
  gtf$strand <- sapply(gtf$strand, function(x) {ifelse(x == "+",1,-1)})  # reformat strand from +,1 to 1,-1

  # get all unique transcripts
  transcript_ids <- unique(gtf$transcript_id)
  transcript_ids <- transcript_ids[!is.na(transcript_ids)] # remove na

  # get the ends of the transcripts depending on whether on -1 or 1 strand
  transcript_ends <- sapply(transcript_ids, function(x) {
    idx <- match(x, gtf$transcript_id)
    ifelse(gtf[idx, 'strand'] == 1, max(gtf[idx, 'end']), min(gtf[idx, 'start']))
  })
  return(transcript_ends)
}

get_sites_in_gene <- function(gene_id, data, gtf) {
  #' Extracts the sites for a given gene ID
  #' @param gene_id, string containing an Ensembl gene ID
  #' @param data, a data.frame containing the information from a single bed file
  #' @param gtf, a data.frame containing a gtf annotation file
  #' @return data.frame of sites that are found in that gene

  gtf <- gtf[gtf$`Gene ID` == gene_id,]  # get row with info for chosen gtf
  sites <- data[data$Chromosome == gtf$Chromosome & data$Strand == gtf$Strand &
                  data$Locus >= gtf$Start & data$Locus <= gtf$Finish, ]
  if (!is_empty(sites$Chromosome)) { # return nothing if no sites found

    sites <- cbind(sites, rep(gtf[,c('Gene ID', 'Gene Name')]))
    sites <- sites %>% mutate_at('Gene ID', as.character) %>% mutate_at('Gene Name', as.character)
    return(sites)
  }
}

sites_in_all_genes <- function(gene_list, data, gtf) {
  #' Slow function. Builds a data.frame of all the edited sites and their gene ids
  #' @param gene_list, list of ensembl IDs
  #' @param data, data.frame of bed file
  #' @param gtf, data.frame of gtf file
  #' @return dataframe

  # Apply to all gene ids
  sites <- mclapply(gene_list, function(x) get_sites_in_gene(x, data, gtf),
                    mc.cores = CORES)

  # clean up
  sites <- rbindlist(sites, idcol = F)
  sites <- drop_na(sites)
  return(sites)

}

sites_for_all_timepoints <- function(gene_list, data, gtf, file_name=NA) {
  #' For a list of bed dataframes, gets all the sites that are in genes'
  #' @param gene_list, list of genes IDs
  #' @param data, list of dataframes holding bed file data
  #' @param gtf, data.frame holding gtf annotation data
  #' @param file_name, name of file to write output to. Defaults to NA and no write
  #' @return dataframe containing all the sites at all time frames that are in genes

  # Apply to all bed dfs in the list of bed dfs
  sites <- lapply(data, function(x) sites_in_all_genes(gene_list, x, gtf))

  # clean up
  sites <- rbindlist(sites, idcol = F)
  sites <- drop_na(sites)

  # optional write to file
  if (!is.na(file_name)) {write_tsv(all_genes_all_info, file_name)}
  return(sites)
}

site_id <- function(data) {
  #' Creates unique ID for every site and adds to the dataframe
  #' @param data, dataframe to make unique identifier for
  #' @return dataframe

  id <- paste0(data$Chromosome,'_', data$Locus)
  data <- cbind(data[,1:3], id, data[,4:dim(data)[2]]) %>%
    mutate_at("id", as.character)
  names(data)[match('id', names(data))] <- 'Site ID'
  return(data)
}

.check_coord_object <- function(x,
                                check_seqnames = FALSE,
                                check_strand = FALSE) {
  if (!is.data.frame(x)) {
    stop(
      "object must be a data.frame. ",
      "GRanges objects are currently not supported and must be converted ",
      "using e.g. as.data.frame()"
    )
  }

  if (!all(c("start", "end") %in% colnames(x))) {
    stop("object must have the columns 'start' and 'end'")
  }

  if (check_seqnames) {
    if (!("seqnames" %in% colnames(x))) {
      stop("object must have the column 'seqnames'")
    }
  }

  if (check_strand) {
    if (!("strand" %in% colnames(x))) {
      stop("object must have the column 'strand'")
    }
  }
}

.check_group_var <- function(x, group_var) {
  if (!is.null(group_var)) {
    if (!all(group_var %in% colnames(x))) {
      stop(
        "group_var ('", group_var, "') ",
        "must be a column in object"
      )
    }
  }
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

# 3.  Filtering  --------------------------------------------------------
genes_x_sites <- function(data, min_sites) {
  #' Filters the dataframe to sites in genes that have at least a minimum
  #' number of sites. Genes with only 2 or 3 sites will not have a dynamic
  #' range that can give sufficient information on their editing rate.
  #' @param data, dataframe of bed data
  #' @param min_sites, minimum number of sites that all the genes have to hvae
  #' @return processed dataframe

  # Filter out genes with fewer than 5 sites
  unique_site_ids <- unique(data$`Site ID`) # list of Site IDs
  unique_sites <- data[match(unique_site_ids, data$`Site ID`),] # subseted df
  unique_genes <- unique(data$`Gene Name`)

  # Get all genes that have at least five sites
  genes_5_sites <- mclapply(unique_genes, function(x) {
    length(which(x == unique_sites$`Gene Name`)) >= min_sites
  }, mc.cores = CORES)
  genes_5_sites <- unlist(genes_5_sites)
  genes_5_sites <- unique_genes[genes_5_sites]

  # Filter all detected data to include only genes with > 5 sites
  data <- data[data$`Gene Name` %in% genes_5_sites,]
  return(data)
}

filter_median <- function(data, reads) {
  #' Filters dataframe to genes that have a median of at least `reads`
  #' @param data, dataframe to filter
  #' @param reads, int - minimum median of reads e.g. 80
  #' @return filtered dataframe

  # Get list of genes that meet filter
  genes <- lapply(unique(data$`Gene ID`), function(x) {
    median <-median(data[data$`Gene ID` == x,]$Reads)
    if (median > reads) return(x)
  })
  genes <- genes[sapply(genes, function(x) !is.null(x))]

  # Filter data
  filtered <- data[data$`Gene ID` %in% genes, ]
  return(filtered)
}

maxmin_utr <- function(data, raw_utr=raw_utr) {
  #' For every Gene ID, gets the biggest range of coordinates that describe the UTRs of the
  #' transcripts. Does not take into account information on transcripts and where their UTRs
  #' start. Just finds the max and min values across all the transcripts of that gene and
  #' returns.
  #' @param data, dataframe containing bed file data. Gene names will be drawn from here
  #' @param raw_utr, pre-processed utr data e.g. the output from loci_from_gtf('three_prime_utr')
  #' @return maxmin coordinates for UTRs (for genes that have them)

  utr <- lapply(unique(data$`Gene ID`), function(x) {
    y <- raw_utr[raw_utr$`Gene ID` == x,]
    if (!dim(y)[1] == 0) { # handle the fact that some genes don't have 3-UTRs
      # start <- ifelse(y$Strand == 1, min(y$Start), max(y$Start))  # min max is reversed for -1 strand
      # finish <- ifelse(y$Strand == 1, max(y$Finish), min(y$Finish))

      start <- min(y$Start)
      finish <- max(y$Finish)

      y <- y[1, ] # trim to first row to return
      y$Start[1] <- start
      y$Finish[1] <- finish
      y} else {
        NA
      }
  })
  utr <- utr[!is.na(utr)]
  utr <- rbindlist(utr)
  return(utr)
}


# 4.  Detect data -----------------------------------------------------------
get_detect_sites <- function(detect_data = detect_data, diff_data = all_sites, sites) {
  #' Given a list of filtered sites, pull all the data from detect files
  #' @param detect_data, list of dataframes containing detect_data
  #' @param diff_data, data frame containing annotated differential data
  #' @param sites, list of sites as Site IDs
  #' @return single df containing the detect data at all non-T0 timepoints
  #' for the shortlisted sites. Contains duplicates where genes overlap
  #' e.g. for Site ID == 1_1749346

  # For each time point, get the sites from the detect df that are in the shortlist
  detect_sites <- mclapply(detect_data, function(x) {
    x[x$`Site ID` %in% sites,]
  }, mc.cores = CORES)
  detect_sites <- rbindlist(detect_sites, idcol=F)
  detect_sites <- detect_sites %>% arrange(Chromosome, Locus)
  # detect_sites <- subset(detect_sites, select = -`Site Locus`) # don't want this

  # Get gene info from df with annotated sites for the shortlist of sites
  gene_info <- diff_data[diff_data$`Site ID` %in% detect_sites$`Site ID`,
                         c('Site ID', 'Gene ID', 'Gene Name','Site Locus')]
  gene_info <- gene_info %>% distinct()  # remove duplicate ROWS

  # merge (annotate detect data)
  detect_sites <- merge(detect_sites, gene_info, all.x = T, all.y = F)
  detect_sites <- detect_sites[, colnames(diff_data)]
  detect_sites <- detect_sites %>% arrange(Chromosome, Locus)
  return(detect_sites)
}



# 5. Fitting --------------------------------------------------------------

# Define models
exp_fit <- function(x, l, a) {1 - a * exp(1)^(- l * x)}
lin_fit <- function(x, l, a) {l * x + (1-a)}

# Fit models
site_to_xy <- function(site, ids = TRUE) {
  #' Subset dataframe to just Time Point and Edits and convert Time Point to int
  #' @param site, rows for a df for a single site
  if (ids) {
    site <- site[,c('Time Point', 'Edits','Site ID')]
    names(site) <- c('x','y','Site ID')
  } else {
    site <- site[,c('Time Point', 'Edits')]
    names(site) <- c('x','y')
  }

  site$x <- as.integer(substr(site$x, 2, 2))
  site <- as.data.frame(site)
  return(site)
}

fit_model <- function(site, FUN=exp_fit, r2 = TRUE) {
  #' Fits single site and returns lambda value and goodness of fit.
  #' Hard coded model to avoid slow down from conditional fitting
  #' @param site, a single site as a dataframe
  #' @param r2, if T then return R^2 value
  #' @return lambda and RSS as list

  # detect whether df has already been converted with `site_to_xy()`
  if (names(site)[1] == 'Time Point') {
  x <- as.integer(substr(site$`Time Point`, 2, 2)) # convert time point to int
  y <- site$Edits
  } else {
    x <- site$x
    y <- site$y
  }
  a <- (1- y[1]) # y intercept from data

  # Fit model
  nls_fit <- nls(y~FUN(x, l, a = a), start = c(l=0.5))
  print(rsquare(site, site))
  return(nls_fit)
}

fit_site_r2 <- function(site, FUN=exp_fit, scaled=FALSE) {
  #' For a single site, fit model and return lambda and r2 values
  #' @param site, site in xy format
  #' @return lambda and R2 values for sites

  x <- site$x
  y <- site$y
  a <- (1- y[1])
  if (scaled) {
    a <- 1
  }

  nls_fit <- nls(y~FUN(x, l, a = a), start = c(l=0.5))
  r2 <- rsquare(nls_fit, site)
  return(c(coefficients(nls_fit)[[1]], r2))
}

fit_gene_r2 <- function(gene, ..., FUN=exp_fit) {
  #' Fits model to all sites in a gene and returns lambda estimate and r2
  #' @param gene, df subset to a single gene
  #' @param FUN, function to fit
  #' @return dataframe of lambda, R2 with Site ID as the row names

  # get list of sites
  sites <- unique(gene$`Site ID`)

  # convert to xy format
  gene <- site_to_xy(gene, ids = TRUE)

  # fit
  fits <- mclapply(sites, function(x) {
    fit_site_r2(gene[gene$`Site ID` == x,])
    }, mc.cores = CORES)

  # turn into df
  fits <- t(data.frame(fits))
  fits <- data.frame(fits, row.names = sites)
  names(fits) <- c('Lambda','R2')
  return(fits)
}

fit_all_r2 <- function(data, FUN=exp_fit) {
  #' Fits the model to all sites in the datafarme, returning lambda estiamte and
  #' R2
  #' @param data, df to fit
  #' @param FUN, function to fit
  #' @return df with the fit for every site in the df
  data <- site_to_xy(data, ids=TRUE)

  sites <- unique(data$`Site ID`)

  fits <- mclapply(sites, function(x) {
    fit_site_r2(data[data$`Site ID` == x, ])
  }, mc.cores = CORES)

  # turn into df
  fits <- t(data.frame(fits))
  fits <- data.frame(fits, row.names = sites)
  names(fits) <- c('Lambda','R2')
  return(fits)
}

fit_site <- function(site, FUN=exp_fit) {
  #' Fits single site and returns lambda value and goodness of fit.
  #' Hard coded model to avoid slow down from conditional fitting
  #' @param site, a single site as a dataframe
  #' @return lambda and RSS as list

  x <- as.integer(substr(site$`Time Point`, 2, 2)) # convert time point to int
  y <- site$Edits
  a <- (1- y[1]) # y intercept from data

  # Fit model
  nls_fit <- nls(y~FUN(x, l, a = a), start = c(l=0.5))
  return(c(coefficients(nls_fit)[[1]], nls_fit$m$deviance()))

}

fit_gene <- function(gene, FUN=exp_fit) {
  #' Fits a specified model to all site in a gene and returns lambda and fit
  #' @param gene, a gene name
  #' @return dataframe named as the gene containing lambda, RSS and the site ID
  #' as the row name

  sites <- unique(gene$`Site ID`)

  fits <- mclapply(sites, function(x) {
    fit_site(gene[gene$`Site ID` == x,], FUN=FUN)
  }, mc.cores = CORES)

  # Turn into df
  fits <- t(data.frame(fits))
  fits <- data.frame(fits, row.names = sites)
  names(fits) <- c('Lambda','RSS')
  return(fits)
}

fit_all_sites <- function(data, FUN=exp_fit) {
  #' Fits a specified model to all site in a gene and returns lambda and fit
  #' @param data, data to fit to
  #' @param FUN, function to fit. Defaults to exp_cdf
  #' @return dataframe named as the gene containing lambda, RSS and the site ID
  #' as the row name

  sites <- unique(data$`Site ID`)

  fits <- mclapply(sites, function(x) {
    fit_site(data[data$`Site ID` == x,], FUN=FUN)
  }, mc.cores = CORES)

  # Turn into df
  fits <- t(data.frame(fits))
  fits <- data.frame(fits, row.names = sites)
  names(fits) <- c('Lambda','RSS')
  return(fits)
}

fit_site_interactive <- function(site, model='exp') {
  #' For a single site, fit a specified model and print goodness of fit metrics and show
  #' plot
  #' @param site, a single site as a dataframe
  #' @param model, if 'exp' then exponential otherwise does linear

  # Get values
  x <- as.integer(substr(site$`Time Point`, 2, 2)) # convert time point to int
  y <- site$Edits
  a <- (1- y[1]) # y intercept from data

  # Pick model
  if (model == 'exp') {
    fit_fun <- function(x, l, a) {1 - a * exp(1)^(- l * x)} # exponential model
  } else fit_fun <- function(x, l, a) l*x + (1-a) # linear model

  # Fit model
  nls_fit <- nls(y~fit_fun(x, l, a = a), start = c(l=0.5))
  lambda <- coefficients(nls_fit)[1] # extract lambda estimate

  # Print results
  print(paste0("lambda = ", lambda))
  print(paste0("RSS = ", nls_fit$m$deviance()))
  plot_fit(data = site, x = 'Time Point', y = 'Edits', lambda = lambda, intercept = a)
}

gene_violin <- function(gene, fits, data) {
  #' For a gene, plot summaries of the R2
  #' @param gene gene name
  #' @param fits fit data to pull from
  #' @param data data to get the Site IDs for a given gene
  #' @return boxplot

  # Unique sites
  sites <- unlist(unique(data[data$`Gene Name` == gene, 'Site ID']))

  # Subset R2 values to sites
  fits <- fits[fits$`Site ID` %in% sites,]
  fits$Gene <- 'SYT11'

  # Violin Plot
  fits %>% ggplot(aes(x=Gene, y=R2)) +
    geom_violin(fill=colours[3],trim = TRUE) +
    theme_classic() +
    ggtitle(paste0('R2 values for ', gene)) +
    ylim(0,1) +
    theme(text = element_text(size=12),
          plot.title = element_text(hjust=0.5))

  # For multiplot
  # fits %>% ggplot(aes(x=Gene, y=R2, fill=Gene)) +
  #   geom_violin(trim = TRUE) +
  #   scale_fill_viridis(option='plasma', end=0.9, discrete=T) +
  #   theme_classic() +
  #   ggtitle('R2 values for gene') +
  #   ylim(0,1) +
  #   theme(text = element_text(size=12),
  #         plot.title = element_text(hjust=0.5))
}


# 6. Extracting BAM data --------------------------------------------------

sites_to_granges <- function(g_data) {
  #' For a list of sites, converts them to a GRange object
  #' @param data, bed data df (pac bio data) subset to just sites of GOI
  #' @return GRanges object

  # Build GRange object of the sites to use in mapping from reference to query space
  gwhich <- GRanges(
    seqnames = g_data$Chromosome,
    ranges = IRanges(g_data$Locus + 1, width = rep(1L, length(g_data$Locus))), # need +1 to change to correct index
    strand = g_data$Strand
  )
  return(gwhich)
}

bam_for_gene <- function(g_data=g_data, bam_file=bam_file) {
  #' From a g_data frame (unique Site IDs), reads specified bam file and returns
  #' a Genomic Alignments object
  #' @param g_data, data with information to read
  #' @param bam_file, path to bam_file
  #' @return GAlignments object containing reads that align to region

  sites <- unique(g_data$`Locus`)

  # Get BAM data
  which <- IRangesList(IRanges(sites[1], sites[length(sites)]))
  names(which) <- g_data$Chromosome[1] # name range after chromosome
  what <- c('qname','strand','pos','qwidth','seq','cigar','qual')
  names(what) <- what
  param <- ScanBamParam(which=which, what=what)
  gbam <- readGAlignments(file=bam_file, use.names=T, param=param)
  gbam <- gbam[gbam@strand == g_data$Strand[1]] # remove if wrong strand
  gbam <- gbam[width(mcols(gbam)$seq) != 0]  # remove if no sequence

  return(gbam)
}

map_sites_reads <- function(site, gwhich=gwhich, gbam=gbam) {
  #' For a given site, gwhich and gbam. Retrieves the base for every read that
  #' maps to that site
  #' @param site, chr site ID
  #' @param gwhich, Granges object of the sites
  #' @param gbam, GAlignments object that contains the BAM information for each read, including sequence and cigar
  #' @return data.table of base calls and corresponding qnames

  # Build mapping
  idx <- match(as.integer(str_split(site, '_')[[1]][2]) + 1, gwhich@ranges@start) # find row of gwhich to use
  mapped <- mapToAlignments(gwhich[idx], gbam) # map first edited site onto the reads
  mapped <- mapped[mapped@ranges@width != 0] # remove reads that fail to map
  gbam <- gbam[match(mapped@seqnames, gbam@NAMES)]  # also remove from the gbam object

  # For a single site, get the base for every read that aligns
  all_reads <- lapply(1:length(mapped), function(x) {
    as.character(unlist(mcols(gbam)$seq[x])[mapped[x]@ranges@start])
  }) %>% unlist()
  all_reads <- data.table('calls'=all_reads, 'qname'=gbam@NAMES)
  names(all_reads) <- c(site,'qname')

  return(all_reads)
}

calls_to_counts <- function(data, simplify = FALSE) {
  #' Convert base calls to counts of Gs per read for histogram. Very dumb function
  #' @param data df containing the base calls for individual reads for a gene
  #' @param simplify whether to coerce anything not 'A' or 'G' to na
  #' @return counts of Gs, As and other e.g. NAs or SNPs for each read
  bases <- c('A','G','T','C')
  counts <- lapply(bases, function(x) {
    sapply(1:dim(data)[1], function(y) {sum(data[y,] == x, na.rm=TRUE)})
  })
  names(counts) <- bases
  strand <- ifelse(sum(counts$A) > sum(counts$T), '+','-')

  if (!simplify) {
    data.frame(
      'qname' = data$qname,
      'G' = counts$`G`,
      'A' = counts$`A`,
      'T' = counts$`T`,
      'C' = counts$`C`,
      'blank' = sapply(1:dim(data)[1], function(x) {sum(data[x,] == '',na.rm=TRUE)}),
      'na' = sapply(1:dim(data)[1], function(x) {sum(is.na(data[x,]))})
    )} else { if (strand == '+') {
      data.frame(
        'qname' = data$qname,
        'G' = counts$G,
        'A' = counts$A,
        'na' = rep(dim(data)[2],dim(data)[1]) - (counts$G + counts$A)
      )
    } else {
      data.frame(
        'qname' = data$qname,
        'C' = counts$C,
        'T' = counts$`T`,
        'na' = rep(dim(data)[2],dim(data)[1]) - (counts$C + counts$`T`)
      )
    }
    }

}

get_reads_for_bams <- function(bam_path=bam_path, bam_files=bam_files,
                               gene='SYT11', pb_data=pb_data) {
  #' For a list of bam files, pulls out the editing of each read that maps to
  #' sites.
  #' @param bam_path path to dir containing bam files
  #' @param bam_files list of bam_file names
  #' @param gene Gene Name to investigate
  #' @param pb_data processed pacbio data_df
  #' @return list of read_dfs, same length as @param bam_files

  # Construct paths
  # bam_files <- paste0(bam_path, bam_files)
  sites <- sort(unique(pb_data[pb_data$`Gene Name` == gene, 'Site ID'])) # edited sites in gene
  g_data <- pb_data[match(sites, pb_data$`Site ID`),] %>% arrange(Locus) # get data for the unique sites in gene
  g_data$Strand <- ifelse(g_data$Strand == 1,'+','-') # change strand to +/- notation
  gwhich <- sites_to_granges(g_data) # Convert sites to GRanges


  # Loop
  res <- lapply(bam_files, function(x) {
    print(x)
    # For each read get the base calls at the site loci
    gbam <- bam_for_gene(g_data, paste0(bam_path,x))
    read_df <- mclapply(1:length(sites), function(y) {
      map_sites_reads(sites[y], gwhich, gbam)
    }, mc.cores=CORES) %>% purrr::reduce(full_join, by = 'qname')

    read_df <- setcolorder(read_df,c("qname",colnames(read_df)[!(colnames(read_df) %in% c("qname"))])) # reorder (view)
  })

  # Name the list items after the bam files from which they were made
  names(res) <- bam_files
  return(res)
}

transpose_rdfs <- function(rdf) {
  #' Transposes the multi bam read_dfs
  #' @param rdf list of read dataframes
  #' @return same list but with the dataframes transposed

  # Apply to each read_df in turn
  res <- lapply(rdf, function(x) {
    reads <- x$qname
    x <- x[, 2:dim(x)[2]]
    x <- lapply(x, function(y) replace_na(y, '.')) # make all chr
    x <- data.frame(x, stringsAsFactors = FALSE)
    names(x) <- str_sub(names(x),2)

    x <- data.frame(t(x), stringsAsFactors = FALSE) # transpose
    names(x) <- reads

    # Get rid of reads with 0 edits
    x <- apply(x, 2, function(y) {
      if ('G' %in% y) {
        return(y) } else return(NA)
    })
    x <- x[!is.na(x)]

    # Return transposed df
    data.frame(x, stringsAsFactors = FALSE)
  })
  return(res)
}

# 7. Maximum Likehood Estimation ------------------------------------------

log_like <- function(t, l_e=lambda_e, l_n=lambda_n) {
  #' Builds log (ln) likelihood model for t given read r
  # print(t)

  # Edited sites
  G_list <- c()
  for (i in 1:length(l_e)) {
    G_list <- append(G_list, log(1 - exp(-l_e[i] * t)))
  }

  sum(G_list) - t * sum(l_n)
}

mle_read <- function(read, fits, FUN=log_like, strand = '+') {
  #' For a given read character vector, generates the MLE for t.
  #' Have to use super assignemnt operator in order to control scope
  #' @param read as a character vector
  #' @param fits fits for the sites in the gene
  #' @param FUN function for building the log_likelihood model
  #' @strand strand of the gene
  # Get lambda values
  if (strand == '+') {
    lambda_e <<- fits[which(read == 'G'), 'Lambda']
    lambda_n <<- fits[which(read == 'A'), 'Lambda']
  } else {
    lambda_e <<- fits[which(read == 'C'), 'Lambda']
    lambda_n <<- fits[which(read == 'T'), 'Lambda']
  }

  # MLE
  MLE_t <- maxLik::maxLik(FUN, start=c(t=25))
  return(MLE_t)
}

mle_gene <- function(read_df, fits, strand = '+') {
  #' @param read_df df containing the editing profile for all reads.
  #' rows = Site ID, cols = reads
  #' @param fits fits for the sites in the gene
  #' @param strand strand of the gene
  #' @return list of maxLik objects for each read

  res <- lapply(1:dim(read_df)[2], function(x) {
    mle_read(read_df[,x], fits, FUN=log_like)
  })

  names(res) <- g_fits$`Site ID`
  return(res)
}



# Plotting ----------------------------------------------------------------
# Colours
colours <- viridis(5, end=0.9, option='plasma')
names(colours) <- c('T0','T2', 'T4', 'T6', 'T8')

plot_sites <- function(data, x='Site Locus', y='Edits', fill='Time Point', xcoords = NA) {
  #' Customasiable bar plot for all time_points simultaneously
  #' @param data, should be a dataframe containing the sites for a single gene
  #' @param x, variable to plot on x axis
  #' @param y, variable to polot on y aix
  #' @param fill, variable to colour by
  #' @retrun plots the ggplot object

  p <- ggplot(data, aes_string(x = as.name(x), y = as.name(y), fill = as.name(fill))) +
    geom_bar(position=position_dodge2(width = 1, preserve = 'single', padding=0), stat='identity') +
    scale_fill_viridis(option='plasma', end = 0.9, discrete = T) +
    theme_classic() +
    theme(text = element_text(size=12))

  if (x == 'Site Locus') p <- p + xlim(min(data$`Site Locus`) - 1, max(data$`Site Locus`) + 1) # scale
  if (y == 'Edits') p <- p + ylim(0,1)  # scale
  if (!is.na(xcoords)) p <- p + xlim(xcoords)
  p
}

plot_point <- function(data, x='Locus', y='Edits', fill = 'Time Point', colours=colours, xcoords=NA) {
  #' Used for scatter plotting, expects a vector of hex codes for colours
  #' @param colours, should be a named vector of hex codes where names are e.g
  #' 'T8'

  # Select colours
  colours <- colours[names(colours) %in% unique(data$`Time Point`)]
  colours <- colours[c(order(colours))] # sort

  # Build plot
  p <- ggplot(data, aes_(x = sym(x), y = sym(y), colour = sym(fill))) +
    geom_point(size=3) +
    theme_classic() +
    scale_color_manual(values = colours) +
    theme(text = element_text(size=12))

  # if (x == 'Locus') p <- p + xlim(155883150, 155883278) # custom for seeing length of gene
  if (x == 'Site Locus') p <- p + xlim(0, max(data$`Site Locus`) + 1) # scale
  if (y == 'Edits') p <- p + ylim(0,1)  # scale
  if (!is.na(xcoords)) p <- p + xlim(xcoords)
  p
}

plot_fit <- function(data, x = 'Time Point', y, fill = 'Time Point', lambda = l, intercept = a) {
  #' Used to plot a single site with the x axis being time point (continuous)

  # Convert time point to int
  data$`Time Point` <- as.integer(substr(data$`Time Point`, 2, 2))

  p <- ggplot(data, aes_string(x = as.name(x), y = as.name(y), fill = as.name(fill))) +
    geom_bar(position=position_dodge2(width = 1, preserve = 'single', padding=0), stat='identity') +
    # geom_point(size=3) +
    scale_fill_viridis(option='plasma', end = 0.9, discrete = F) +
    theme_classic() +
    theme(text = element_text(size=12))

  # add model fit
  p <- p + stat_function(fun=function(x, l=lambda, a=intercept) exp_fit(x,l=lambda,a=intercept),
                         colour = '#3A3B3C', size=2)

  # adjust axes
  if (x == 'Site Locus') p <- p + xlim(min(data$`Site Locus`) - 1, max(data$`Site Locus`) + 1) # scale
  if (x == 'Time Point') p <- p + scale_x_continuous(breaks = c(0,2,4,6,8)) # scale
  if (y == 'Edits') p <- p + ylim(0,1)  # scale

  p
}

plot_rss <- function(data, ylims=NA) {
  #' Plots the RSS values for each Site Locus along a gene
  #' @param data, a dataframe of: Lambda, RSS e.g. output of fit_gene()
  #' @param ylims, list of ylims if desired

  p <- ggplot(data = data, aes(x = 1:dim(data)[1], y = RSS, colour=RSS)) +
    geom_point(size=2) +
    theme_classic() +
    ylim(0,0.08)

  if (!is.na(ylims)) {
    p <- p + ylim(ylims[1],ylims[2])
  }
  p
}

plot_density <- function(data, x=NA, filter=NA) {
  #' Makes a density plot of a variable
  #' @param data, a dataframe
  #' @param x, the variable to plot

  # Filter to values below threshold
  if (!is.na(filter)) {
    data <- data[data < filter]
  }

  # If no name supplied, try to take the first column's name
  if (is.na(x)) {
    x <- names(data)[1]
    p <- ggplot(data, aes(x = !!sym(x)))
  } else {
    p <- ggplot(data, aes_(x = sym(x)))
  }
  p <- p +
    # geom_histogram(fill=colours[5], binwidth = 1) +
    geom_density(colour=colours[1], lwd=1.5) +
    theme_classic() +
    theme(text = element_text(size=12))

  # p <- p + ylim(0,400) +
    # geom_vline(xintercept = 50, colour = 'grey')
  p
}

plot_histogram <- function(data, x=NA, t_point = "T0") {
  #' Makes a density plot of a variable
  #' @param data, a dataframe
  #' @param x, the variable to plot

  # If no name supplied, try to take the first column's name
  if (is.na(x)) {
    x <- names(data)[1]
    p <- ggplot(data, aes(x = !!sym(x)))
  } else {
    p <- ggplot(data, aes_(x = sym(x)))
  }
  p <- p +
    geom_histogram(fill=colours[match(t_point,DETECT_TIME_POINTS)], binwidth = 1) +
    theme_classic() +
    theme(text = element_text(size=12))

  p <- p + ylim(0,400) +
    geom_vline(xintercept = 50, colour = 'grey')
  p
}

plot_g_hist <- function(data, x='G', title='SYT11', no_zero = TRUE) {
  #' Makes a density plot of a variable
  #' @param data, a dataframe
  #' @param x, the variable to plot

  if (no_zero) {
    data <- data[data[match(x,names(data))] != 0,]
  }

  p <- ggplot(data, aes_(x = sym(x))) +
  geom_histogram(colour='white', fill=blues9[7], binwidth = 1) +
  theme_classic() +
  labs(x='Edits per read') +
  ggtitle(title) +
  # ylim(0,700) +
  # xlim(-1,30) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust=0.5))
  p
}


plot_cdf <- function(data, x=NA) {
  #' Makes a cumulative density plot of a variable
  #' @param data, a dataframe
  #' @param x, the variable to plot

  # If no name supplied, try to take the first column's name
  if (is.na(x)) {
    x <- names(data)[1]
    p <- ggplot(data, aes(x = !!sym(x)))
  } else {
    p <- ggplot(data, aes_(x = sym(x)))
  }
  p <- p +
    stat_ecdf(colour=colours[1], lwd=1.5) +
    # geom_density(colour=colours[1], lwd=1.5) +
    theme_classic() +
    theme(text = element_text(size=12)) +
    labs(y='Cumulative Density') +
  p
}

plot_multi_violins <- function(data, num, title) {
  #' Plot, for each gene, a violin for the R2 value
  #' @param data dataframe containing fits and a Gene Name for each fit
  #' @param num number of violins to plot
  #' @return plots multiple violins on same plot

  # Order by median R^2 value
  medians <- lapply(unique(data$`Gene Name`), function(x) {
    median(data[data$`Gene Name` == x, 'R2'])
  }) %>% unlist()
  medians <- data.frame('Gene Name'=unique(data$`Gene Name`), 'median'=medians,
                        stringsAsFactors = F)
  names(medians) <- c('Gene Name', 'median')

  # Count number of sites in each gene
  counts <- lapply(unique(data$`Gene Name`), function(x) {
    sum(data$`Gene Name` == x)
  }) %>% unlist()
  counts <- data.frame('Gene Name'=unique(data$`Gene Name`), 'counts'=counts,
                       stringsAsFactors = F)
  names(counts) <- c('Gene Name','counts')

  # Attach medians to fits df and sort desc
  data <- left_join(data, medians)
  data <- left_join(data, counts) %>% arrange(desc(median))
  data <- data[data$`Gene Name` %in% unique(data$`Gene Name`)[1:num],]

  # Rename column
  names(data)[4] <- 'Gene'

  # plot
  data %>%
    mutate(Gene=factor(Gene,levels=unique(data$Gene))) %>%
    ggplot(aes(x=Gene, y=R2, fill=Gene)) +
    geom_violin(trim = TRUE) +
    geom_point(aes(x=Gene, y=median), colour='black',size=4) +
    geom_point(aes(x=Gene, y=median), colour='white',size=3) +
    # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
    #              geom = "crossbar",
    #              width = 1,
    #              color = "black") +
    # stat_summary(fun.y  = 'median',
    #              geom = "point",
    #              size = 4,
    #              color = "black") +
    # stat_summary(fun.y  = 'median',
    #              geom = "point",
    #              size = 3,
    #              color = "white") +
    # geom_text(aes(label= sprintf("%.2f", median), y=0.25)) +
    geom_text(aes(label= counts), y=1.05, check_overlap = TRUE) +
    scale_fill_viridis(option='plasma', end=0.8, discrete=T) +
    theme_classic() +
    ggtitle(title) +
    ylim(0,1.05) +
    theme(text = element_text(size=12),
          plot.title = element_text(hjust=0.5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
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




