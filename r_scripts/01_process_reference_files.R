# "Processing reference files for use
# author: "Aaron Wagen"
# date: "March 2023"





# Intro and prep ----------------------------------------------------------


# This script takes raw reference files and processes them for use subsequence workflows
# For the purpose of the parts of the project that use ensembl data (ie, all astrocyte analysis), will use their conventions:
# chromosome names are without 'chr'. 
# Will only use chromosomes 1-22, X, Y and MT (not M)

# GTEX data uses Gencode, which includes a chr in the chromosome name. For matching with the gtex editing matrices, this will be kept.


# Prep
library(tidyverse)
library(here)
library(data.table)
library(plyranges)


# Define arguments
# Need to check lines to skip each time and set as lines_to_skip

args <- list(reference_directory = here::here("references"),
             alu_annotations_raw_file = here::here("references",
                                                   "line_sine_dfam2.0_hg38.fa.out.gz"),
             alu_processed_file_path = here::here("references",
                                                  "alu_processed_range.RDS"),
             ensembl93_gtf = "/camp/lab/gandhis/home/users/wagena/references/GRCh38/annotation/Homo_sapiens.GRCh38.93.gtf",
             reduced_ensembl93_path = here::here("references",
                                                 "ensembl93_reduced_gtf.RDS"),
             reduced_ensembl93_range_path = here::here("references",
                                                 "ensembl93_reduced_gtf_range.RDS"),
             ensembl93_lengths_summary_path = here::here("references",
                                                         "ensembl93_gene_lengths_summary.RDS"),
             gencode26_gtf = "/camp/lab/gandhis/home/users/wagena/references/GRCh38/annotation/gencode.v26.annotation.gtf.gz",
             reduced_gencode26_path = here::here("references",
                                                       "gencode26_reduced_gtf.RDS"),
             reduced_gencode26_range_path = here::here("references",
                                                     "gencode26_reduced_gtf_range.RDS"),
             gencode26_lengths_summary_path = here::here("references",
                                                                 "gencode26_gene_lengths_summary.RDS")
)



# Process ensembl v93 gtf for use with astrocyte data -------------------------

# Note: ensembl has sequence names that are just numeric (no "chr"), and slightly different column names compared to gencode.

#Load ensembl file
ensembl_gtf <- plyranges::read_gff(args$ensembl93_gtf)

ensembl_annotation <- ensembl_gtf %>%
  dplyr::filter(seqnames %in% c(1:22, "X", "Y", "MT")) %>% 
  dplyr::select(type, gene_id, gene_name,  gene_biotype)  

ensembl_annot_df_redux <- ensembl_annotation %>% # Create dataframe with just gene_id and biotype information
  as_tibble() %>% 
  dplyr::select(-c(seqnames, start, end, width, strand, type)) %>%
  distinct(., gene_id,
           .keep_all = T)

# Create reduced exon to represent entire transcript
reduced_exons <- ensembl_annotation %>%
  plyranges::filter(type == 'exon') %>%
  plyranges::group_by(gene_id) %>%
  plyranges::reduce_ranges_directed()


# Add annotations to reduced exons
reduced_exons_df <- reduced_exons %>%
  as_tibble %>%
  left_join(.,
            ensembl_annot_df_redux,
            by = "gene_id") %>%
  mutate(gene_type = "exon")

reduced_introns <- ggtranscript::to_intron(reduced_exons_df,
                                           group_var = 'gene_id') %>%
  relocate(seqnames, start, end, width, strand, gene_id, gene_name, gene_biotype, gene_type) %>%
  mutate(gene_type = type) %>%
  dplyr::select(-type)


reduced_ensembl <- rbind(reduced_introns, reduced_exons_df) %>%
  arrange(seqnames, start)

# Find cumulative length of exon (for later use in editing summary metric)
exon_length <- reduced_ensembl %>% 
  group_by(gene_id) %>% 
  dplyr::filter(gene_type == "exon") %>% 
  summarise(exon_length = sum(width))

# Find cumulative length of gene (for later use in editing summary metric)
gene_length <- reduced_ensembl %>% 
  group_by(gene_id) %>% 
  summarise(gene_length = sum(width))

gene_exon_lengths_summary <- reduced_exons_df %>% 
  dplyr::distinct(gene_id, gene_name, gene_biotype) %>% 
  dplyr::left_join(.,
                   exon_length, by = "gene_id") %>% 
  dplyr::left_join(.,
                   gene_length, by = "gene_id")


saveRDS(gene_exon_lengths_summary,
        file = args$gencode26_lengths_summary_path)

reduced_ensembl_range <- reduced_ensembl %>%
  GenomicRanges::makeGRangesFromDataFrame(
    df = .,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = T,
    ignore.strand = F)

saveRDS(reduced_ensembl, file= args$reduced_ensembl93_path)

saveRDS(reduced_ensembl_range, args$reduced_ensembl93_range_path)


# Process gencode v26 gtf for use with GTEX8 data -----------------------------

#Load gencode file
gencode_gtf <- plyranges::read_gff(args$gencode26_gtf)


gencode_annotation <- gencode_gtf %>%
  plyranges::select(gene_type=type, gene_id, gene_name,  gene_biotype=gene_type)

gencode_annot_df_redux <- as_tibble(gencode_annotation) %>% # Create dataframe with just gene_id and biotype information
  dplyr::select(-c(seqnames, start, end, width, strand, gene_type)) %>%
  distinct(., gene_id,
           .keep_all = T)

# Create reduced exon to represent entire transcript (use directed to respect strand direction)
reduced_exons <- gencode_annotation %>%
  plyranges::filter(gene_type == 'exon') %>%
  plyranges::group_by(gene_id) %>%
  plyranges::reduce_ranges_directed()


# Add annotations to educed exons
reduced_exons_df <- reduced_exons %>%
  as_tibble %>%
  left_join(.,
            gencode_annot_df_redux,
            by = "gene_id") %>%
  mutate(gene_type = "exon")

reduced_introns <- ggtranscript::to_intron(reduced_exons_df,
                                           group_var = 'gene_id') %>%
  relocate(seqnames, start, end, width, strand, gene_id, gene_name, gene_biotype, gene_type) %>%
  mutate(gene_type = type) %>%
  dplyr::select(-type)


reduced_gencode <- rbind(reduced_introns, reduced_exons_df) %>%
  arrange(seqnames, start)

# Find cumulative length of exon (for later use in editing summary metric)
exon_length <- reduced_gencode %>% 
  group_by(gene_id) %>% 
  dplyr::filter(gene_type == "exon") %>% 
  summarise(exon_length = sum(width))

# Find cumulative length of gene (for later use in editing summary metric)
gene_length <- reduced_gencode %>% 
  group_by(gene_id) %>% 
  summarise(gene_length = sum(width))

gene_exon_lengths_summary <- reduced_exons_df %>% 
  dplyr::distinct(gene_id, gene_name, gene_biotype) %>% 
  dplyr::left_join(.,
                   exon_length, by = "gene_id") %>% 
  dplyr::left_join(.,
                   gene_length, by = "gene_id")
  
  
saveRDS(gene_exon_lengths_summary,
        file = args$ensembl93_lengths_summary_path)

reduced_gencode_range <- reduced_gencode %>%
  GenomicRanges::makeGRangesFromDataFrame(
    df = .,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = T,
    ignore.strand = F)

saveRDS(reduced_gencode, file= args$reduced_gencode26_path)

saveRDS(reduced_gencode_range, args$reduced_gencode26_range_path)





# Process raw alu results -------------------------------------------------

alu_annotations <- fread(args$alu_annotations_raw_file, 
                         fill=TRUE,
                         skip =2) 
column_names <- c("SW_score", "perc_div", "perc_del", "perc_ins", "chromosome", "start", "end", "chr_bp_remaining", "?strand", "repeat_type", "repeat_class", "repeat_position_start", "repeat_position_end", "repeat_bp_remaining", "ID") 
colnames(alu_annotations) <- column_names

alu_annotations_filtered <- alu_annotations %>% 
  dplyr::select(chromosome, start, end, repeat_type, repeat_class) %>% 
  mutate(chromosome = str_replace_all(chromosome, c("chr" = "", 
                                                    "M" = "MT")),
         across(c(chromosome, repeat_type, repeat_class), as_factor)) %>% 
  filter(chromosome %in% c(1:22, "X", "Y", "MT")) %>% 
  droplevels() 


# Combine classes of repeat - note these are the same classification as Zhongbo used in her work.
alu_annotations_grouped <- alu_annotations_filtered %>%
  mutate(repeat_class_extended = repeat_class,
         repeat_class = ifelse(repeat_class == "SINE/Alu", "SINE/Alu",
                               ifelse(str_detect(repeat_class, "LINE/L1"), "LINE/L1",
                                      ifelse(repeat_class == "Simple_repeat", "Simple_repeat", 
                                             ifelse(repeat_class == "SINE/MIR", "SINE/MIR",
                                                    ifelse(repeat_class == "LINE/L2", "LINE/L2",
                                                           ifelse(str_detect(repeat_class, "LTR"), "LTRs", 
                                                                  ifelse(str_detect(repeat_class, "RNA"), "RNAs",
                                                                         ifelse(str_detect(repeat_class, "DNA"), "DNAs",
                                                                                ifelse(repeat_class == "Retroposon/SVA", "Retroposon/SVA", 
                                                                                       ifelse(str_detect(repeat_class, "LINE") 
                                                                                              & repeat_class != "LINE/L1" 
                                                                                              & repeat_class !="LINE/L2", "Other LINE",
                                                                                              ifelse(str_detect(repeat_class, "SINE") 
                                                                                                     & repeat_class != "SINE/MIR" 
                                                                                                     & repeat_class !="SINE/Alu", 
                                                                                                     "Other SINE", 
                                                                                                     "Other")
                                                                                       )))))))))),
         repeat_class = as.factor(repeat_class)
  )



alu_annotations_range <- alu_annotations_grouped %>% 
  GenomicRanges::makeGRangesFromDataFrame(
    df = ., 
    seqnames.field = "chromosome", 
    start.field = "start",
    end.field = "end",
    keep.extra.columns = T,
    ignore.strand = T)


# Save complete alu range
saveRDS(alu_annotations_range, args$alu_processed_file_path)


