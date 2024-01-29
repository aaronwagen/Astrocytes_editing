# Generate vep input
# Aaron Wagen
# Nov 2022

# This script parses results from all samples (detect and differential), and creates a single file for input to VEP containing all editing sites
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


# Args
args <- list(alu_annotations_range_file = here::here("references",
                                                     "sine_line_dfam2_hg38_range_220817.rds"),
             reduced_ensemble_range_file = here::here("references", 
                                                      "reduced_ensembl_38-93_range_221120_stranded.rds"),
             vep_results_range_file = here::here("raw_data",
                                                 "astro_neuron_bulk",
                                                 "16_alignment",
                                                 "vep",
                                                 "vep_results_all_samples_strand_filtered_range.RDS"),
             vep_dir = here::here("raw_data",
                                  "astro_neuron_bulk",
                                  "16_alignment",
                                  "vep"),
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



# Create annotation file

## Detect data
detect_file_paths <- list.files(args$detect_dir, pattern=".AG.bed", full.names= TRUE)
detect_files <- basename(detect_file_paths)

paths_annotation <- list.files(args$detect_dir, pattern=".AG.bed", full.names= TRUE) %>% 
  as_tibble() %>% 
  mutate(Sample = basename(detect_file_paths) %>% 
           str_replace(".AG.bed","")) %>% 
  dplyr::rename(file_path = value)

annotation <- read_tsv(args$annotation_file,
                       col_types = "ccff") %>% 
  dplyr::rename(Sample = sample_Id)  %>% 
  left_join(., paths_annotation, by = 'Sample') #%>% 
#mutate(sample_desc = sample_name) # to add in when multiplexing

annotation

read_multi_detect <- read_multi_beds_original_strand_names(args$detect_dir, detect_files)

processed_multi_detect <- process_multi_beds_names(read_multi_detect, annotation)

all_detect_data <- data.table::rbindlist(processed_multi_detect)

all_detect_data <- all_detect_data %>% 
  mutate(locus = str_c(chromosome, "_", start, "_", strand)) %>% 
  dplyr::select(Sample, chromosome, start, strand, locus)

all_detect_data

## Differential data

dif_file_paths <- list.files(args$differential_dir, pattern=".AG.bed", full.names= TRUE)
differential_files <- basename(dif_file_paths)
dif_file_paths
             
read_multi_dif <- read_multi_dif_beds(args$differential_dir, differential_files)

processed_multi_dif_data <- process_multi_dif_beds(read_multi_dif, z_score_cuttoff=0)

all_dif_data <- data.table::rbindlist(processed_multi_dif_data) 

all_dif_data <- all_dif_data %>% 
  dplyr::select(Sample = SampleName, chromosome, start, strand, locus)

## Combine and filter for unique locus

all_editing_sites <- bind_rows(all_detect_data, all_dif_data)


all_distinct_editing_sites <- all_editing_sites %>% 
  dplyr::distinct(locus, .keep_all = T) %>% 
  dplyr::arrange(chromosome, start)


## Put in vep format
vep_input <- all_distinct_editing_sites %>% 
  mutate(mutation = "A/G",
         end = start) %>% 
  dplyr::select(chromosome, start, end, mutation, strand)

# Save file
write_tsv(vep_input, 
          file = file.path(args$vep_dir, "all_samples_vep_input.tsv"),
          col_names = FALSE)            
             