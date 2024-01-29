#! /bin/bash

#####################################################

#SBATCH --job-name=star_samtools
#SBATCH --output=/dev/null
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

#SBATCH --time=0-24:00:00

######################################################

# STAR FIRST PASS
# Aaron Wagen, adapted from George Young
# Sept 2022

# This script runs the first pass of star alignment on a list of samples outlined in SampleAnnot. 
# In the first pass it will save only the split junction .tab files, to collapse and filter into a combined SJ.filtered.tab file for second pass filtering.

# Example command to run on slurm:
# tail -n+2 SampleAnnot.txt |
# parallel -kj1 --col-sep "\t" 'sbatch --parsable --export=sample=
# {2},R1=preAlignmentQC/fastp/
# {1}*_R1_001_trimmed.fastq.gz,R2=preAlignmentQC/fastp/
# {1}*_R3_001_trimmed.fastq.gz,index=/camp/home/wagena/gandhis/home/users/wagena/references/GRCh38/index/star/hg38.93_sjdboverhang149,bamdir=/camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/astro_neuron_bulk/16_alignment/star
# ~/scripts/star_first_pass_16mismatch.sh' >> jobs

######################################################

# MODULES
ml purge
module load GCC/10.2.0 SAMtools/1.13-GCC-10.2.0 STAR/2.7.9a-GCC-10.3.0

######################################################


 # PARAMETERS

# #### Reference files
# genome="/camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/astro_neuron_bulk/JACUSA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# #bamfile= ${@1}
# #bamexample="/camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/astro_neuron_bulk/star_mismatch_trial/mismatch_04/NM6322_A1-BX0175-001_S37_Aligned.sortedByCoord.out.bam"
# sample=`basename $bamfile`
# prefix=`echo $sample | sed "s/_Aligned.sortedByCoord.out.bam//"` # replaces 'Aligned...out.bam' with empty space to get the prefix
# bamdir=`dirname $bamfile`
# mismatch=`basename $bamdir`

echo Start : `date`
echo Sample : $sample
echo: Read 1 : $R1


# tmp=`pwd -P | sed 's:/working/:/scratch/:'`
# mkdir -p hisat2 $tmp/{hisat2,trimmed}

# f_base="${R1##*/}"
# f_dir="${R1%$f_base}"
# f_base="${f_base%%_1.*}"


## FLAG SUMMARY
# --outSAMattributesMD is crucial here for JACUSA to be run efficiently (if don't have this then JACUSA needs to run this itself with the fasta file. In hisat this is run with samtools calmd)
# intronMotif - adds some utility calling the strand of a transcript over a splice junction. George adds this as defauly.
# --outStd - 'standard out' sends the output to pipe. If want to use this option need a piped output to follow
# --outReadsUnmapped Fastx ", # output in separate fast/fastq files the unmapped/partially-mapped reads
# --outSAMtype BAM SortedByCoordinate ", # output as a sorted BAM,
# --outFilterType BySJout ", # removes spurious split reads
# --outFilterMultimapNmax 1 ", # only allows reads to be mapped to one position in the genome
# --outFilterMismatchNmax 999 ", # Maximum number of mismatches per pair. Large numbers switch off filter. Instead we filter by "--outFilterMismatchNoverReadLmax".
# --outFilterMismatchNoverReadLmax 0.16 is defauly. We are using 0.16 to not exlclude edits of interest
# --alignIntronMin 20 ", # min intron length. As per ENCODE options.
# --alignIntronMax 1000000 ", # max intron length. As per ENCODE options (currently from ensembl its 1,097,903 from KCNIP4).
# --alignMatesGapMax 1000000 ", # max gap between pair mates. As per ENCODE options.
# --alignSJoverhangMin 8 ", # minimum unannotated split read anchor. As per ENCODE options.
# --alignSJDBoverhangMin 3 " # minimum annotated split read anchor. Default is 3. This is the whole purpose of using 2 passes, to allow previously unnanotated junctions to be considered annotated.
#--outSAMtype None : in the first pass don't need a sam output...

STAR \
--runThreadN 32 \
--genomeDir ${index} \
--readFilesIn ${R1} ${R2} \
--outSAMattributes NH HI AS nM NM MD \
--readFilesCommand zcat \
--outFileNamePrefix "${bamdir}/${sample}_" \
--outSAMtype None \
--outSAMstrandField intronMotif \
--outFilterMultimapNmax 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.16 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 3


# After this step can filter the splice junction SJ.tab files using the following:
# The SJ.filtered.tab was created by using the following command:
# From https://groups.google.com/g/rna-star/c/fdPiuJ3mR9w

# cat *.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > SJ.filtered.tab

# This filters out:
# - non canonical junctions (column 5 >0)
# - junctions supported by multimappers only (column 7 >0)
# - columns supported by too few reads (column 7>2)
# - remove annotated junctions as these are always included from the GTF (column 6==0)
# - sort the columns and remove duplicates (sort and uniq)




# tail -n+2 SampleAnnot.txt
# tail -n+2 SampleAnnot.txt | parallel -kj1 --col-sep "\t" 'sbatch --parsable --export=sample={2},R1=preAlignmentQC/fastp/{1}*_R1_001_trimmed.fastq.gz,R2=preAlignmentQC/fastp/{1}*_R3_001_trimmed.fastq.gz,index=/camp/home/wagena/gandhis/home/users/wagena/references/GRCh38/index/star/hg38.93_sjdboverhang120,bamdir=/camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/astro_neuron_bulk/16_alignment/star ~/scripts/star_first_pass.sh' >> jobs
