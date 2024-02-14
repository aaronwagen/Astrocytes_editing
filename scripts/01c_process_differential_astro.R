# This script processes jacusa differential data from astrocyte neuron coculture experiments
# Aaron Wagen
# Nov 2022


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
library(parallel)
library(gridExtra)
library(ggtranscript)
library(cowplot)

CORES=1

here::i_am("r_scripts/process_detect_astro.R")

source(here::here("r_scripts",
                  "r_functions.R"))

# Args
args <- list(alu_annotations_range_file = here::here("raw_data",
                                                     "astro_neuron_bulk",
                                                     "16_alignment",
                                                     "vep",
                                                     "alu_annotations_distinct_range_230117.rds"),
             reduced_ensemble_range_file = here::here("raw_data",
                                                      "astro_neuron_bulk",
                                                      "16_alignment",
                                                      "vep",
                                                      "reduced_ensembl_38-93_range_230117_stranded_distinct.rds"),
             vep_results_range_file = here::here("raw_data",
                                                 "astro_neuron_bulk",
                                                 "16_alignment",
                                                 "vep",
                                                 "vep_results_unique_annotations_range.RDS"),
             are_results_range_file = here::here("references",
                                                 "ARE_locations_reduced_range_GRCh38.RDS"),
             gene_expression_file = here::here("raw_data",
                                               "astro_neuron_bulk",
                                               "16_alignment",
                                               "vep",
                                               "astr_cocult_expressed_genes.xlsx"),
             gtf_file = "/camp/lab/gandhis/home/users/wagena/AIediting/raw_data/astro_neuron_bulk/references/Homo_sapiens.GRCh38.93.gtf.gz",
             fasta_file = "/camp/lab/rodriquess/home/users/wagena/timestamps/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
             annotation_file = here::here("raw_data",
                                          "astro_neuron_bulk",
                                          "16_alignment",
                                          "SampleAnnot.txt"),
             annotation_differential = here::here("raw_data",
                                                  "astro_neuron_bulk",
                                                  "16_alignment",
                                                  "SampleAnnot_differential.txt"),
             detect_dir = here::here("raw_data",
                                     "astro_neuron_bulk",
                                     "16_alignment",
                                     "jacusa_detect"),
             differential_dir = here::here("raw_data",
                                           "astro_neuron_bulk",
                                           "16_alignment",
                                           "jacusa_differential"),
             output_dir = here::here("processed_data",
                                     "astro_neuro_coculture")
)

## 1B) Define levels for filtering

source(here::here("r_scripts",
                  "filtering_args.R"))


# Load data and filter  ---------------------------------------------------------------

print("Loading editing files")


file_paths <- list.files(args$differential_dir, pattern=".AG.bed", full.names= TRUE)
differential_files <- basename(file_paths)
CORES = 1

z_score_cuttoff = 1.96



read_multi_dif <- read_multi_dif_beds(args$differential_dir, differential_files)



processed_multi_dif_data <- process_multi_dif_beds(read_multi_dif, z_score_cuttoff= 0)
all_dif_data <- data.table::rbindlist(processed_multi_dif_data)

#processed_multi_dif_data_zscore <- process_multi_dif_beds(read_multi_dif, z_score_cuttoff= z_score_cuttoff)
#all_dif_data_zscore <- data.table::rbindlist(processed_multi_dif_data_zscore) 

# Add in column for significance
all_dif_data <- all_dif_data %>% 
  mutate(sig_edit = ifelse(z.score >=z_score_cuttoff | z.score <= -z_score_cuttoff, "significant", "not_significant"),
         across(chromosome, as_factor),
         SampleNumber = 1:nrow(all_dif_data)
         ) %>% 
  dplyr::filter(chromosome %in% c(1:22, "X", "Y", "MT")) %>% 
  droplevels() 

ggplot(all_dif_data, aes(x=z.score)) +
  geom_histogram(binwidth = 0.1) +
  xlim(-3,3)
  

sig_dif_data <- all_dif_data %>% 
  dplyr::filter(sig_edit == "significant")

all_dif_range <- make_genomic_range(all_dif_data)

head(all_dif_data)

# Annotated data ----------------------------------------------------------

print("Loading annotation files")

## 1) Load preprocessed annotation data
alu_annotations_range <- readRDS(args$alu_annotations_range_file)
reduced_ensembl_range <- readRDS(args$reduced_ensemble_range_file)
vep_results_range <- readRDS(args$vep_results_range_file)
are_results_range <- readRDS(args$are_results_range_file)

print("Annotating files")

annotated_vep_edits <- plyranges::join_overlap_left_directed(x = all_dif_range,
                                                             y = vep_results_range)


# Annotate with vep 
annotated_vep_gtf_edits <- plyranges::join_overlap_left_directed(x = annotated_vep_edits,
                                                                 y = reduced_ensembl_range) %>% 
  plyranges::mutate(gene_type_embl = replace_na(gene_type_embl, "intergenic"))



# For results with vep annotation, collapse the annotation into categories (exon, intron, intergenic)
annotated_vep_gtf_edits_vepd <- annotated_vep_gtf_edits %>% 
  as_tibble() %>% 
  tidyr::drop_na(locus_vep) %>% 
  mutate(type = ifelse(Consequence %in% c("stop_lost",
                                          "start_lost",
                                          "stop_gained",
                                          "frameshift_variant",
                                          "start_retained_variant",
                                          "stop_retained_variant",
                                          "transcript_amplification",
                                          "inframe_insertion",
                                          "inframe_deletion",
                                          "protein_altering_variant",
                                          "missense_variant",
                                          "synonymous_variant",
                                          "coding_sequence_variant",
                                          "5_prime_UTR_variant",
                                          "3_prime_UTR_variant",
                                          "non_coding_transcript_exon_variant"), "exon",
                       ifelse(Consequence %in% c("splice_acceptor_variant",
                                                 "splice_donor_variant",
                                                 "splice_region_variant",
                                                 "splice_donor_region_variant",
                                                 "incomplete_terminal_codon_variant",
                                                 "intron_variant",
                                                 "NMD_transcript_variant",
                                                 "non_coding_transcript_variant",
                                                 "TFBS_ablation",
                                                 "TFBS_amplification",
                                                 "TF_binding_site_variant",
                                                 "upstream_gene_variant",
                                                 "downstream_gene_variant",
                                                 "regulatory_region_ablation",
                                                 "regulatory_region_amplification",
                                                 "feature_elongation",
                                                 "regulatory_region_variant",
                                                 "feature_truncation"), "intron",
                              ifelse(Consequence == "intergenic", "intergenic",
                                     "other")
                       )
  )
  ) 

# For results with no vep annotation, take the gtf annotation for Gene, Gene name (symbol), and biotype
annotated_vep_gtf_edits_gtfd <- annotated_vep_gtf_edits %>% 
  as_tibble() %>% 
  dplyr::filter(is.na(locus_vep)) %>% 
  mutate(type = gene_type_embl,
         Gene = gene_id_embl,
         SYMBOL = gene_name_embl,
         BIOTYPE = gene_biotype_embl)

# Combine and tidy columns
annotated_edits_combined <- bind_rows(annotated_vep_gtf_edits_vepd,
                                      annotated_vep_gtf_edits_gtfd) %>% 
  dplyr::arrange(seqnames, start) %>% 
  dplyr::select(-variant.y) %>%  # Remove extra variant column
  dplyr::rename(gene_type = type) %>% # Take the  type column (combined from vep and gtf) as the final gene_ty
  dplyr::relocate(gene_type, .after = SYMBOL) %>% 
  as_granges()




# Annotate with alus
annotated_edits_alu <-
  plyranges::join_overlap_left_directed(x = annotated_edits_combined,
                                        y = alu_annotations_range)



# # Annotate with ares
annotated_edits_complete <-
  plyranges::join_overlap_left_directed(x = annotated_edits_alu,
                                        y = are_results_range %>%
                                          plyranges::select(-ARE_pattern)) # plyranges doesn't like the concatenated ARE_pattern column



# Relevel factors
annotated_edits_complete<- annotated_edits_complete %>% 
  as_tibble() %>% 
  dplyr::mutate(gene_type = factor(gene_type, levels = rev(filtering_args$gtf_genic_type)),
                gene_biotype_embl = factor(gene_biotype_embl, levels = filtering_args$gtf_genic_biotype),
                IMPACT = factor(IMPACT, levels = filtering_args$vep_impact),
                Consequence = factor(Consequence, levels = filtering_args$vep_consequence),
                BIOTYPE = factor(BIOTYPE, levels = filtering_args$vep_biotype),
                TSL = factor(TSL, levels = filtering_args$vep_TSL),
                CANONICAL = factor(CANONICAL, levels = filtering_args$vep_canonical),
                repeat_class = factor(repeat_class, levels = rev(filtering_args$alu_repeat_levels)),
                ARE_binary = case_when(!is.na(ARE_cluster) ~ "ARE",
                                              TRUE ~ NA)
  )


print("Saving files")

saveRDS(annotated_edits_complete, str_c(args$output_dir, "/differential_summary_comparison_all_edits.RDS"))


sig_annotated_edits <- annotated_edits_complete %>% 
  dplyr::filter(sig_edit == "significant")

saveRDS(sig_annotated_edits, str_c(args$output_dir, "/differential_summary_comparison_sig_edits.RDS"))










