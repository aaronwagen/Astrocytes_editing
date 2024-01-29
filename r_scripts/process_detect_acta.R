####
#title: "Process acta jacusa detect files"
#author: "Aaron Wagen"
###

# The purpose of this script is to take a BED file from jacusa_to_bed.py, annotate it with genic and repeat information, and output QC metrics and plots.
# Inputs: 
#     - a dataframe called 'annotation' with column for file_path, column called Sample with unique name for each sample, and other annotations as required
#     - preprocessed genic, VEP, and alu annotations


# Can run this script from the terminal, AIediting folder with the code: `nohup Rscript r_scripts/process_detect_acta.R &> logs/process_detect_acta.log &`
# Or can run on slurm with: 
# ml R/4.2.0-foss-2021b
# sbatch --time=48:00:00 --mem=64G --output logs/process_detect_acta.log --wrap="Rscript r_scripts/process_detect_acta.R"

# Libraries and arguments -------------------------------------------------



# Load Library
options(datatable.showProgress = FALSE)
library(here) 
library(GenomicRanges)
library(data.table)
library(tidyverse)
library(plyranges)
library(paletteer)
library(patchwork)
library(gridExtra)
library(ggtranscript)
library(cowplot)

CORES=1

# Set relative path with here::here
here::i_am("r_scripts/process_detect_acta.R")

source(here::here("r_scripts",
                  "r_functions.R"))


# Args
args <- list(
  # vep_results_range_file = here::here("raw_data",
  #                                     "acta_neuropath_bulk",
  #                                     "editing",
  #                                     "vep",
  #                                     "vep_results_unique_annotations_range.RDS"),
  # gtf_file = "/camp/lab/gandhis/home/users/wagena/AIediting/raw_data/astro_neuron_bulk/references/Homo_sapiens.GRCh38.93.gtf.gz",
  # fasta_file = "/camp/lab/rodriquess/home/users/wagena/timestamps/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
  alu_annotations_range_file = here::here("references",
                                          "alu_processed_range.RDS"),
  reduced_ensembl_range_file = here::here("references",
                                          "ensembl93_reduced_gtf_range.RDS"),
  annotation_file = here::here("raw_data",
                               "acta_neuropath_bulk",
                               "SampleAnnot.txt"),
  # collapsed_annotation_file = here::here("raw_data",
  #                                        "acta_neuropath_bulk",
  #                                        "editing",
  #                                        "SampleAnnot.txt"),
  annotation_differential = here::here("raw_data",
                                       "acta_neuropath_bulk",
                                       "editing",
                                       "SampleAnnot_differential.txt"),
  detect_dir = here::here("raw_data",
                          "acta_neuropath_bulk",
                          "editing",
                          "jacusa_detect"),
  differential_dir = here::here("raw_data",
                                "acta_neuropath_bulk",
                                "editing",
                                "jacusa_differential"),
  output_dir = here::here("processed_data",
                          "acta_neuropath_bulk",
                          "jacusa_detect")
)




# Load data and define filtering levels -----------------------------------

source(here::here("r_scripts",                  
                  "filtering_args.R"))

print("Loading annotation files")

# 1) Load preprocessed annotation data
# These have been processed as per the script: "process_vep_results.R"

alu_annotations_range <- readRDS(args$alu_annotations_range_file)
reduced_ensembl_range <- readRDS(args$reduced_ensembl_range_file) 


#vep_results_range <- readRDS(args$vep_results_range_file)



## 2) Create annotation file
file_paths <- list.files(args$detect_dir, pattern=".AG.bed", full.names= TRUE)

paths_annotation <- list.files(args$detect_dir, pattern=".AG.bed", full.names= TRUE) %>% 
  as_tibble() %>% 
  mutate(Sample = basename(file_paths) %>% 
           str_replace("_detect.AG.bed","")) %>% 
  dplyr::rename(file_path = value)

# Match annotations with filepaths
annotation <- read_tsv(args$annotation_file) %>% 
  dplyr::rename(Sample = sample_name)  %>% 
  left_join(., paths_annotation, by = 'Sample') #%>% 
#mutate(sample_desc = sample_name) # to add in when multiplexing

# Select out columns of interest from annotation
annotation <- annotation %>% 
  dplyr::select(Sample, Disease_Group, Sex, AoO, AoD, DD, PMI, RIN = RINe_bulkRNA_Tapestation, aSN, TAU, AB = `thal AB`, aCG_aSN_score = `aCG aSN score`, sample_Id, file_path)


# Annotating detect data loop ---------------------------------------------


# Create summary list

summary_list <- list()

# Run script

for (sample_number in 1:nrow(annotation)) {
  sample_path <- annotation$file_path[sample_number]
  sample_name = annotation$Sample[sample_number]
  list_name = str_c(sample_name, "_list")
  
  print(sample_name)
  print(sample_path)
  
  data <- read_bed_QC(sample_path)
  data <- process_bed_QC(data, annotation)
  data <- data %>% 
    mutate(locus = str_c(chromosome, "_", start, "_", strand))
  data_range <- make_genomic_range(data)
  
  
  #   annotated_vep_edits <- plyranges::join_overlap_left_directed(x = data_range,
  #                                                                y = vep_results_range)
  #   
  #   
  #   # Annotate with vep 
  #   annotated_vep_gtf_edits <- plyranges::join_overlap_left_directed(x = annotated_vep_edits,
  #                                                                    y = reduced_ensembl_range)  %>%     
  #     plyranges::mutate(gene_type_embl = replace_na(gene_type_embl, "intergenic"))
  #   
  #   
  #   # For results with vep annotation, collapse the annotation into categories (exon, intron, intergenic)
  #   annotated_vep_gtf_edits_vepd <- annotated_vep_gtf_edits %>% 
  #     as_tibble() %>% 
  #     tidyr::drop_na(locus_vep) %>% 
  #     mutate(type = ifelse(Consequence %in% c("stop_lost",
  #                                             "start_lost",
  #                                             "stop_gained",
  #                                             "missense_variant",
  #                                             "synonymous_variant",
  #                                             "coding_sequence_variant",
  #                                             "5_prime_UTR_variant",
  #                                             "3_prime_UTR_variant",
  #                                             "stop_retained_variant",
  #                                             "non_coding_transcript_exon_variant"), "exon",
  #                          ifelse(Consequence %in% c("splice_acceptor_variant",
  #                                                    "splice_donor_variant",
  #                                                    "splice_region_variant",
  #                                                    "splice_donor_region_variant",
  #                                                    "incomplete_terminal_codon_variant",
  #                                                    "intron_variant",
  #                                                    "upstream_gene_variant",
  #                                                    "downstream_gene_variant"), "intron",
  #                                 ifelse(Consequence == "intergenic", "intergenic",
  #                                        "other")
  #                                 )
  #                          )
  #            ) 
  # 
  #   # For results with no vep annotation, take the gtf annotation for Gene, Gene name (symbol), and biotype
  # annotated_vep_gtf_edits_gtfd <- annotated_vep_gtf_edits %>% 
  #   as_tibble() %>% 
  #   dplyr::filter(is.na(locus_vep)) %>% 
  #   mutate(type = gene_type_embl,
  #          Gene = gene_id_embl,
  #          SYMBOL = gene_name_embl,
  #          BIOTYPE = gene_biotype_embl)
  # 
  # 
  # 
  # # Combine and tidy columns
  #  annotated_edits_combined <- bind_rows(annotated_vep_gtf_edits_vepd,
  #                                        annotated_vep_gtf_edits_gtfd) %>% 
  #    dplyr::arrange(seqnames, start) %>% 
  #   dplyr::select(-gene_type_embl, gene_id_embl, gene_biotype_embl, locus_embl, locus_vep) %>% 
  #   dplyr::rename(gene_type = type) %>% 
  #    dplyr::relocate(gene_type, .after = SYMBOL) %>% 
  #    as_granges()
  
  
  annotated_edits_combined <- plyranges::join_overlap_left_directed(x = data_range,
                                                                    y = reduced_ensembl_range)
  
  # Annotate with alus
  annotated_edits_complete <-
    plyranges::join_overlap_left_directed(x = annotated_edits_combined,
                                          y = alu_annotations_range)
  
  

  
  # Relevel factors
  annotated_edits_complete <- annotated_edits_complete %>% 
    as_tibble() %>% 
    dplyr::mutate(gene_type = factor(gene_type, levels = filtering_args$gtf_genic_type),
                  gene_biotype = factor(gene_biotype, levels = filtering_args$gtf_genic_biotype),
                 # IMPACT = factor(IMPACT, levels = filtering_args$vep_impact),
                #  Consequence = factor(Consequence, levels = filtering_args$vep_consequence),
                #  BIOTYPE = factor(BIOTYPE, levels = filtering_args$vep_biotype),
                #  TSL = factor(TSL, levels = filtering_args$vep_TSL),
                #  CANONICAL = factor(CANONICAL, levels = filtering_args$vep_canonical),
                  repeat_class = factor(repeat_class, levels = rev(filtering_args$alu_repeat_levels))
    )
  
  
  # Remove duplicated rows
  annotated_edits_complete <- annotated_edits_complete %>% 
    arrange(locus, gene_type, gene_biotype) %>% 
    group_by(locus) %>% 
    dplyr::slice(1)
  
  
  # Create summary dataframes with annotations
  annotated_edits_list <- filter_edits_qc(annotated_edits_complete)
  
  
  # number_of_edits_barplot <- barplot_by_gene_annotation(annotated_edits_list,
  #                                                       proportional="stack",
  #                                                       comparison = gene_type) +
  #   labs(title = "Number of edits")
  # 
  # proportion_of_edits_barplot <- barplot_by_gene_annotation(annotated_edits_list,
  #                                                           proportional="fill",
  #                                                           comparison = gene_type) +
  #   labs(title = "Proportion of edits", y = "Proportion of editing events" )
  # 
  # proportion_of_edits_gtf_biotype_barplot <- barplot_by_gene_annotation_dropNA(annotated_edits_list,
  #                                                                              proportional="fill",
  #                                                                              comparison = gene_biotype_embl) +
  #   scale_fill_paletteer_d("ggsci::default_igv", #"colorBlindness::SteppedSequential5Steps",
  #                          direction = 1,
  #                          na.value = "grey70") +
  #   theme(legend.text = element_text(size=8)) +
  #   labs(title = "Proportion of edits by biotype from GTF", subtitle = "NA's removed", y = "Proportion of editing events" )
  # 
  # 
  # proportion_of_edits_vep_biotype_barplot <- barplot_by_gene_annotation_dropNA(annotated_edits_list,
  #                                                                              proportional="fill",
  #                                                                              comparison = BIOTYPE) +
  #   scale_fill_paletteer_d("ggsci::default_igv", #"colorBlindness::SteppedSequential5Steps",
  #                          direction = 1,
  #                          na.value = "grey70") +
  #   theme(legend.text = element_text(size=8)) +
  #   labs(title = "Proportion of edits by biotype from VEP", subtitle = "NA's removed", y = "Proportion of editing events" )
  # 
  # 
  # 
  # proportion_of_edits_alus <- barplot_by_gene_annotation(annotated_edits_list,
  #                                                        proportional="fill",
  #                                                        comparison = repeat_class) +
  #   labs(title = "Proportion of repeat classes", y = "Proportion of editing events" )
  # 
  # proportion_of_edits_alus_dropna <- barplot_by_gene_annotation_dropNA(annotated_edits_list,
  #                                                                      proportional="fill",
  #                                                                      comparison = repeat_class) +
  #   theme(legend.position = "none") +
  #   labs(title = "Proportion of repeat classes, no NAs", y = "Proportion of editing events" )
  # 
  # alu_edits_annotated <- purrr::map(annotated_edits_list, ~dplyr::filter(., repeat_class == "SINE/Alu"))
  # 
  # proportion_of_alus_genic <- barplot_by_gene_annotation(alu_edits_annotated,
  #                                                        proportional="fill",
  #                                                        comparison = gene_type) +
  #   labs(title = "Genic annotations of alu edits", y = "Proportion of editing events" )
  # 
  # #Create summary plot
  # genic_legend <- cowplot::get_legend(number_of_edits_barplot)
  # alu_legend <- cowplot::get_legend(proportion_of_edits_alus)
  # sample_summary_plot <- grid.arrange(arrangeGrob(number_of_edits_barplot + theme(legend.position = "none"),
  #                                                 proportion_of_edits_barplot + theme(legend.position = "none"),
  #                                                 genic_legend,
  #                                                 ncol=3, widths = c(.9,.9,.25)),
  #                                     proportion_of_edits_gtf_biotype_barplot,
  #                                     proportion_of_edits_vep_biotype_barplot,
  #                                     arrangeGrob(proportion_of_edits_alus + theme(legend.position = "none"),
  #                                                 proportion_of_edits_alus_dropna + theme(legend.position = "none"),
  #                                                 alu_legend,
  #                                                 ncol=3, widths = c(.9,.9,.4)),
  #                                     arrangeGrob(proportion_of_alus_genic +theme(legend.position = "none"),
  #                                                 genic_legend,
  #                                                 ncol = 3, widths = c(0.9, 0.25, 0.9)), # make figure same width
  #                                     ncol = 1,
  #                                     top = str_c("Sample ", sample_name)
  # )
  
  
  # Save sample summary plot
  # ggsave(plot = sample_summary_plot, width = 10, height = 18, dpi = 300,
  #        filename = str_c(args$output_dir, "/", sample_name, ".pdf"))
  
  # Create sample summary list
  assign(list_name,
         list("All_edits_lists" = annotated_edits_list#,
              # "Number_of_edits_plot" = number_of_edits_barplot,
              # "Genic_proportions_plot" = proportion_of_edits_barplot,
              # "VEP_biotype_proportions_plot" = proportion_of_edits_vep_biotype_barplot,
              # "GTF_biotype_proportions_plot" =  proportion_of_edits_gtf_biotype_barplot,
              # "Repeat_proportions_plot" = proportion_of_edits_alus,
              # "Repeat_proportions_dropna+plot" = proportion_of_edits_alus_dropna,
              # "Alu_edits_lists" = alu_edits_annotated,
              # "Alu_genic_proportions_plot" = proportion_of_alus_genic,
              # "Sample_summary_plot" = sample_summary_plot
              )
  )
  
  
  # Combine lists to allow for comparisons
  summary_list[[list_name]] <- eval(as.name(list_name))
  
}


saveRDS(summary_list, str_c(args$output_dir, "/detect_summary_list.RDS"))




## Comparison across samples from QC data (Filtered for z score)

# A map function that will return all dataframes in the all_edits_lists list, including all edits, zscore_filter, and z_score_reads filter. Given the edits have all be labelled by individual filtering, I will likely take these and filter the dataframes
comparison_all_df <- summary_list %>% 
  purrr::map(., "All_edits_lists") %>% 
  purrr::map_dfr(., "edits_df")

comparison_sig_df <- summary_list %>% 
  purrr::map(., "All_edits_lists") %>% 
  purrr::map_dfr(., "z_score_edits")

alu_comparison_sig_df <-  summary_list %>% 
  purrr::map(., "Alu_edits_lists") %>% 
  purrr::map_dfr(., "z_score_edits")

saveRDS(comparison_all_df, str_c(args$output_dir, "/detect_summary_comparison_all_edits.RDS"))
saveRDS(comparison_sig_df, str_c(args$output_dir, "/detect_summary_comparison_sig_edits.RDS"))
saveRDS(alu_comparison_sig_df , str_c(args$output_dir, "/detect_summary_comparison_sig_alu_edits.RDS"))



# Useful numbers (To be completed_)

# 
#   
# # Number of loci with duplicate gtf annotations
# annotated_gtf_edits %>% 
#   as_tibble() %>% 
#   drop_na(gene_id) %>% 
#   dplyr::arrange(locus) %>% 
#   group_by(locus) %>% 
#   tidytable::filter(base::n()>1) %>% 
#   distinct(locus) %>% 
#   nrow()
# 
# # Number of unique editing sites annotated by gtf
# annotated_gtf_edits_df %>% 
#   tidyr::drop_na(gene_id) %>%  # Remove sites that gtf didn't annotate
#   distinct(locus) %>% # Find distinct loci
#   nrow() # Count rows
# 
# 
# 
# annotated_gtf_vep_edits %>% 
#   as_tibble() %>% 
#   dplyr::summarise(annotated_edits = length(unique(locus_vep)),
#             unmatched_edits = sum(is.na(annotated_gtf_vep_edits$locus_vep)),
#             duplicated_edis = sum(duplicated(locus_vep))
#   )
# 
# unmatched_vep_editing_sites <- annotated_gtf_vep_edits %>% 
#   dplyr::filter(is.na(locus_vep)) %>% 
#   as_tibble()
# 
# matched_vep_editing_sites <-  annotated_gtf_vep_edits %>% 
#   as_tibble() %>% 
#   arrange(locus_vep) %>% 
#   slice(1:100)
# 
# annotated_gtf_vep_edits %>% 
#   as_tibble() %>% 
#   dplyr::select(locus_vep) %>% 
#   is.na() %>% 
#   sum()
# 
# duplicated_rows %>% 
#   distinct(locus) %>% 
#   nrow()
# 
# 
# annotated_edits %>% 
#   as_tibble() 
# annotated_edits %>% 
#   dplyr::select(locus) %>% 
#   is.na() %>% 
#   sum()
# 
# nrow(dplyr::distinct(annotated_edits$locus))
# lewunique(annotated_edits$locus)
# 
# annotated_edits %>% 
#   as_tibble() %>% 
#   dplyr::n_distinct(seqnames, start, strand)
# 
# annotated_edits_df <- as_tibble(annotated_edits)
# n_distinct
# 
# annotated_edits_df %>% 
#   length(unique(seqnames))
# 
# annotated_edits$se
# 
# distinct(seqnames, start, strand) %>% 
#   nrow()



