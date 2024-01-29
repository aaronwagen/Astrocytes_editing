#!/bin/bash
#SBATCH --job-name=jacusa_detect
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00

#SBATCH --mail-user=aaron.wagen@crick.ac.uk


#Script for running JACUSA detect on a single sample (submit sample name with script)

# tail -n 1 SampleAnnot.txt | cut -f2 | parallel -dryrun -kj1 'sbatch --parsable --export=EXAMPLE={} --output=/camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/acta_neuropath_bulk/editing_indiv/{}.log /camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/acta_neuropath_bulk/scripts/jacusa_detect_single_acta.sh' >> jobs
# tail -n 40 SampleAnnot.txt | cut -f2 | parallel -kj1 'sbatch --parsable --export=EXAMPLE={} --output=/camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/acta_neuropath_bulk/editing_indiv/{}.log /camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/acta_neuropath_bulk/scripts/jacusa_detect_single_acta.sh' >> jobs
# tail -n 24 SampleAnnot.txt | cut -f2 | parallel -kj1 'sbatch --parsable --export=EXAMPLE={} --output=/camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/acta_neuropath_bulk/editing/jacusa_detect/{}.log /camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/acta_neuropath_bulk/scripts/jacusa_detect_single_acta.sh' >> jobs

#load module
ml Java/13.0.2

BASEDIR="/camp/home/wagena/gandhis/home/users/wagena/AIediting/raw_data/acta_neuropath_bulk/editing/"
INDIR=${BASEDIR}"star/"
OUTDIR=${BASEDIR}"jacusa_detect/"


BAMFILE=`find $INDIR -maxdepth 1 -name "*${EXAMPLE}*Aligned_fix_sort_markdup_calmd.bam" -print`

echo Processing : $BAMFILE
echo Output in  : $OUTDIR

echo Start : `date`

# The following code regarding samples and prefixes is superfluous if the $EXAMPLE is the prefix, but useful if it is not
SAMPLE=`basename $BAMFILE`
PREFIX=`echo $SAMPLE | sed "s/_Aligned_fix_sort_markdup_calmd.bam//"`

echo Prefix : $PREFIX

java -jar ~/bin/JACUSA_v2.0.2-RC.jar call-1 -p 16 -P FR-SECONDSTRAND -a D,Y -F 3852 -r ${OUTDIR}${PREFIX}_detect.jacusa ${BAMFILE}
