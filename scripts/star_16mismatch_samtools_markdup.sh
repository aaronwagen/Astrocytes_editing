#! /bin/bash

#####################################################

#SBATCH --job-name=star_samtools
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G

#SBATCH --time=0-24:00:00

######################################################

# README

# STAR SECOND PASS + SAMTOOLS
# Aaron Wagen, adapted from George Young
# Sept 2022

# This script runs the second pass of star alignment on a list of samples outlined in SampleAnnot, utilising the combined and filtered SJ.filtered.tab.
# It will then run samtools to further annotated reads in preparation for JACUSA editing detection.

# Example command to run script:
# tail -n+2 SampleAnnot.txt |
# parallel -kj1 --col-sep "\t" 'sbatch --parsable --export=sample=
# {2},R1=/camp/home/wagena/AIediting/raw_data/astro_neuron_bulk/preAlignmentQC/fastp/{1}*_R1_001_trimmed.fastq.gz,R2=/camp/home/wagena/AIediting/raw_data/astro_neuron_bulk/preAlignmentQC/fastp/{1}*_R3_001_trimmed.fastq.gz,index=/camp/home/wagena/gandhis/home/users/wagena/references/GRCh38/index/star/hg38.93_sjdboverhang149,bamdir=/camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/astro_neuron_bulk/16_alignment/star
# ~/scripts/star_16mismatch_samtools_markdup.sh' >> jobs



######################################################

# MODULES
ml purge
module load GCC/10.2.0 SAMtools/1.13-GCC-10.2.0 STAR/2.7.9a-GCC-10.3.0

######################################################


echo Start : `date`
echo Sample : ${sample}
echo Read1 : ${R1}

# read1=${R1}
# read2=${R2}
# prefix=${sample}

# echo read1 : ${read1}
# echo prefix : ${prefix}

cd ${bamdir}
tmp=`pwd -P | sed 's:/working/:/scratch/:'`
mkdir -p $tmp/samtools

## FLAG SUMMARY
# --outSAMattributesMD is crucial here for JACUSA to be run efficiently (if don't have this then JACUSA needs to run this itself with the fasta file. In hisat this is run with samtools calmd)
# intronMotif - adds some utility calling the strand of a transcript over a splice junction. George adds this as defauly.
# --outStd - 'standard out' sends the output to pipe. If want to use this option need a piped output to follow
# --outReadsUnmapped Fastx ", # output in separate fast/fastq files the unmapped/partially-mapped reads
# --outSAMtype BAM SortedByCoordinate ", # output as a sorted BAM,
# --outFilterType BySJout ", # removes spurious split reads
# --outFilterMultimapNmax 1 ", # only allows reads to be mapped to one position in the genome
# --outFilterMismatchNmax 999 ", # Maximum number of mismatches per pair. Large numbers switch off filter. Instead we filter by "--outFilterMismatchNoverReadLmax".
# --outFilterMismatchNoverReadLmax 0.04 is defauly. We are using 0.16 to not exlclude edits of interest
# --alignIntronMin 20 ", # min intron length. As per ENCODE options.
# --alignIntronMax 1000000 ", # max intron length. As per ENCODE options (currently from ensembl its 1,097,903 from KCNIP4).
# --alignMatesGapMax 1000000 ", # max gap between pair mates. As per ENCODE options.
# --alignSJoverhangMin 8 ", # minimum unannotated split read anchor. As per ENCODE options.
# --alignSJDBoverhangMin 3 " # minimum annotated split read anchor. Default is 3.
# MD is crucial here for JACUSA to be run efficiently (if don't have this then JACUSA needs to run this itself with the fasta file. In hisat this is run with samtools calmd)

STAR \
--runThreadN 32 \
--genomeDir ${index} \
--readFilesIn ${R1} ${R2} \
--outStd SAM \
--outSAMattributes NH HI AS nM NM MD \
--readFilesCommand zcat \
--outFileNamePrefix ${bamdir}/${sample}_ \
--outReadsUnmapped Fastx \
--sjdbFileChrStartEnd /camp/home/wagena/AIediting/raw_data/acta_neuropath_bulk/editing/star/SJ.filtered.tab \
--outSAMstrandField intronMotif \
--outFilterType BySJout \
--outFilterMultimapNmax 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.16 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 3 \
| samtools view -u - \
| samtools sort -n -u@4 -T "${tmp}/samtools/samtools_sort_${sample}" -o -  \
| samtools fixmate -mu@2 - - \
| samtools sort -u@4 -T "${tmp}/samtools/samtools_sort_${sample}_1" - \
| samtools markdup -@4 -f ${sample}"_samtools.md.log" -T "${tmp}/samtools/samtools_markdup_${sample}" - "${sample}_Aligned_fix_sort_markdup_calmd.bam" \
&& samtools index "${sample}_Aligned_fix_sort_markdup_calmd.bam"



# # set directory and temporary directory



# echo Running samtools

# echo Complete Samtools : `date`


# tail -n+2 SampleAnnot.txt
# tail -n+2 SampleAnnot.txt | parallel -kj1 --col-sep "\t" 'sbatch --parsable --export=sample={2},R1=/camp/home/wagena/AIediting/raw_data/astro_neuron_bulk/preAlignmentQC/fastp/{1}*_R1_001_trimmed.fastq.gz,R2=/camp/home/wagena/AIediting/raw_data/astro_neuron_bulk/preAlignmentQC/fastp/{1}*_R3_001_trimmed.fastq.gz,index=/camp/home/wagena/gandhis/home/users/wagena/references/GRCh38/index/star/hg38.93_sjdboverhang120,bamdir=/camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/astro_neuron_bulk/16_alignment/star ~/scripts/star_samtools_markdup.sh' >> jobs
