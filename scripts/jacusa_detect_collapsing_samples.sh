#! /bin/bash
#SBATCH --job-name=jacusa_detect
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00


#Call the script with the arguments afterwards: outputfile_and_path, list_of_bam_files_to_collapse_separated_by_a_comma
# FR-SECONDSTAND is the correct way to call this script

ml Java/13.0.2

outfile=${1}
files_to_collapse=${2}

echo "Outfile :" $outfile
echo "Files to collapse : " $files_to_collapse

echo Running JACUSA

java -jar ~/bin/JACUSA_v2.0.2-RC.jar call-1 -p 8 -P FR-SECONDSTRAND -a D,Y -F 3852 -r ${outfile} ${files_to_collapse}

echo Complete JACUSA_v2