#!/bin/bash

#degraded and ancient DNA QC pipeline for Illumina NextSeq Data - mostly copied from Alisa's aDNA_qc.sh and adapted for redser4
#last edit KB 2021Apr8
		#modified TRBLQC_040521.sh added -S flag to seqprep
# for now, remove MIA analysis
# replaced $SAMTOOLS with samtools and ${MAPDAMAGE2} with mapDamage
#to make executable: chmod +x filename.sh

##to run ./SCRIPT.sh |& tee -a log.txt


############################################################
##General envelopes - EDIT THESE
PROJECT='TRBL_20xQC'

#cd /scratch/kim
DATE=`date +%Y-%m-%d`
mkdir TRBL/QC/${PROJECT}-${DATE}

HOME= /scratch/kim
RAW_READ_DIRECTORY=/redser4/raw/210317_NB551675_0092_AHMVMGBGXH/pgl_TRBL20x
PREFIX='KBSC201_' #raw read file prefix
PROCESSING_OUTPUT=/scratch/kim/TRBL/QC/${PROJECT}-${DATE}
TAXON=Agelaius
REFERENCE_SEQUENCE=/redser4/genomes/Agelaius_tricolor/TRBLgenome/TRBL-genome.fasta  # Remember to include the full path to reference genome -- need to unzip first?
MIA_REFERENCE_SEQUENCE=/redser4/genomes/mitogenomes/Agelaius_phoeniceus_NC_018801.1.fasta  # Remember to include the full path to mitogenome
REFERENCE_NAME=TRBL  # To allow naming of files, if multiple references are to be used
MIA_REFERENCE_NAME=AGPH  # To allow naming of files, if multiple references are to be used
GET_SAMPLES=/redser4/projects/kim/get_samplelist_reps.sh
CALC_STATS=/redser4/projects/chloe/Scripts/calculate_library_stats.R

### Initial setting up of directories
cd ${RAW_READ_DIRECTORY}
bash ${GET_SAMPLES} ${PREFIX} f > ${PROCESSING_OUTPUT}/${PROJECT}-${DATE}-samples.txt

cd ${PROCESSING_OUTPUT}

mkdir Raw_data_symlinks
mkdir BWA_analyses
mkdir MEGAN_analyses
mkdir MIA_analyses
mkdir SeqPrep_output
mkdir MapDamage_output
mkdir Scripts
mkdir Sample_lists_and_progress_files

echo "DATE=${DATE}"
echo "PREFIX=${PREFIX}"
echo "RAW_READ_DIRECTORY=${RAW_READ_DIRECTORY}"
echo "TAXON=${TAXON}"
echo "REFERENCE_SEQUENCE=${REFERENCE_SEQUENCE}"
echo "REFERENCE_NAME=${REFERENCE_NAME}"
echo "MIA_REFERENCE_SEQUENCE=${MIA_REFERENCE_SEQUENCE}"
echo "MIA_REFERENCE_NAME=${MIA_REFERENCE_NAME}"

##SeqPrep envelopes
SEQPREP_OUTPUT=${PROCESSING_OUTPUT}/SeqPrep_output # A directory to store SeqPrep output files
SEQPREP_MIN_LENGTH=28        # Removes unmappably short reads.
SEQPREP_OVERLAP=20           # Allows for confident merging of reads. Can be reduced to 10 if needed.
SEQPREP_LOCATION=/redser4/personal/jonas/bin      # To find SeqPrep, if using edser2

echo "Minimum read length is set to $SEQPREP_MIN_LENGTH."
echo "Merging reads overlapping at $SEQPREP_OVERLAP bases."

##BWA and SAMtools envelopes
BWA_OUTPUT=${PROCESSING_OUTPUT}/BWA_analyses     # A directory to store BWA intermediate files and output
# SAMTOOLS=/projects/redser2/usr/aesoares-bin/samtools019
INDEX_ALGORITHM=bwtsw           # If reference is <2Gb use 'is', if >2Gb use 'bwtsw'
SEED_DISABLE=1024            # Following ancient DNA data processing protocols
BWA_THREADS=10                # To speed up analysis edited (5 slow, 10 normal, 15 fast)
BAM_MIN_QUALITY=30           # Provides a fairly high cutoff

echo "Seed disabled for seedless mapping."
echo "Minimum mapping quality is set to $BAM_MIN_QUALITY phred."

##mapDamage envelopes - make sure REFERENCE_SEQUENCE, REFERENCE_NAME and BWA_OUTPUT are enabled
# MAPDAMAGE2=/soe/pheintzman/bin/mapDamage-master/bin/mapDamage		# To assess damage rates in the dataset
MAX_MISINCORP_FREQUENCY=0.3		# Use 0.3 if not too badly damaged, use 0.5 if badly damaged
READ_PLOT_LENGTH=25			# The number of nucleotides to plot at the 5' and 3' ends of the read
MAX_READ_LENGTH=150			# The maximum read length to consider

echo "Nucleotide misincorporation frequency is 0.5. Use 0.3 if not badly damaged."

##MEGAN envelopes
MEGAN_OUTPUT=${PROCESSING_OUTPUT}/MEGAN_analyses     # A directory to store MEGAN files and output
BLAST_DATABASE=/redser4/resources/blastdb/v5/nt   # Location of the BLAST database
MEGAN_THREADS=15                # To speed up analysis edited (10 normal)

# ##MIA envelopes
# MIA_OUTPUT=${PROCESSING_OUTPUT}/MIA_analyses     # A directory to store MIA output
# ANCIENT_DNA_MATRIX=/projects/redser3-notbackedup/projects/pheintzman/Scripts/ancient.submat.txt
# MIA_COVERAGE_FILTER=/projects/redser3-notbackedup/projects/pheintzman/Scripts/mia_consensus_coverage_filter.pl
# FASTX_TOOLKIT=/redser4/projects/chloe/Scripts/fastx_toolkit-0.0.13.2/src
# MIA_COVERAGE_FILTER_ANDRE=/projects/redser3-notbackedup/projects/common_jobs/coverage_filter_3.pl

##Prinseq envelopes
PRINSEQ_LITE=/redser4/projects/chloe/Scripts/prinseq-lite.pl
PRINSEQ_GRAPHS=/redser4/projects/chloe/Scripts/prinseq-graphs.pl
PRINSEQ_STATS=${PROCESSING_OUTPUT}/PRINSEQ_stats
COMPLEXITY_METHOD=dust			# dust is the standard used by BLAST. The entropy method is an alternative, but is not widely used.
COMPLEXITY_THRESHOLD=7			# Pretty low, but follows the PRINSEQ_LITE manual. Recommends 70 if using entropy.
COMBINE_PAIRED_END_READS=/redser4/projects/chloe/Scripts/combinePairedEndReads.pl
SPLIT_PAIRED_END_READS=/redser4/projects/chloe/Scripts/splitPairedEndReads.pl

##Other envelopes
GET_INSERT_SIZE=/redser4/projects/chloe/Scripts/getinsertsize.py
#############################################################################################

# #Initial setting up of files

for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
do
	gzip -d ${RAW_READ_DIRECTORY}/${SAMPLE}*.fastq.gz # unzips file in raw directory
	wait
	ln -s ${RAW_READ_DIRECTORY}/${SAMPLE}*_R1_001.fastq ${PROCESSING_OUTPUT}/Raw_data_symlinks/${SAMPLE}_R1_001.fastq
	wait
	ln -s ${RAW_READ_DIRECTORY}/${SAMPLE}*_R2_001.fastq ${PROCESSING_OUTPUT}/Raw_data_symlinks/${SAMPLE}_R2_001.fastq
	wait
done
wait
echo "...Initial data processing is complete" > ${PROJECT}_progress_file_${DATE}.txt

##################################################

# #SeqPrep -- removing adapters and merging reads. Overlap (-o) is set to 20, minimum length (-l) is set to 25, minimum quality (-q) is set to 15.

for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
do
	echo ${SAMPLE}
	${SEQPREP_LOCATION}/SeqPrep2 -f ${PROCESSING_OUTPUT}/Raw_data_symlinks/${SAMPLE}_R1_001.fastq -r ${PROCESSING_OUTPUT}/Raw_data_symlinks/${SAMPLE}_R2_001.fastq -1 ${SEQPREP_OUTPUT}/${SAMPLE}_R1_unmerged.fastq.gz -2 ${SEQPREP_OUTPUT}/${SAMPLE}_R2_unmerged.fastq.gz -S -q 15 -L ${SEQPREP_MIN_LENGTH} -A AGATCGGAAGAGCACACGTC -B AGATCGGAAGAGCGTCGTGT -s ${SEQPREP_OUTPUT}/${SAMPLE}_merged.fastq.gz -E ${SEQPREP_OUTPUT}/${SAMPLE}_readable_alignment.txt.gz -o ${SEQPREP_OVERLAP} -d 1 -C ATCTCGTATGCCGTCTTCTGCTTG -D GATCTCGGTGGTCGCCGTATCATT >& ${SEQPREP_OUTPUT}/${SAMPLE}_SeqPrep_output.txt
	wait
	gzip -d ${SEQPREP_OUTPUT}/${SAMPLE}*.gz #unzips files in seq prep output directory
	wait
	gzip ${RAW_READ_DIRECTORY}/${SAMPLE}*.fastq #zips raw data again
	wait
done
wait
echo "...Adapter trimming and merging of reads is complete" >> ${PROJECT}_progress_file_${DATE}.txt

##################################################

# #Filtering reads for low complexity sequences

for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
do

# #Remove low complexity reads from merged files
	perl ${PRINSEQ_LITE} -fastq ${SEQPREP_OUTPUT}/${SAMPLE}_merged.fastq -out_good ${SEQPREP_OUTPUT}/${SAMPLE}_merged.complexity_filtered -out_bad null -lc_method ${COMPLEXITY_METHOD} -lc_threshold ${COMPLEXITY_THRESHOLD} -line_width 0
	wait

# #Remove low complexity reads from unmerged files
	perl ${COMBINE_PAIRED_END_READS} ${SEQPREP_OUTPUT}/${SAMPLE}_R1_unmerged.fastq ${SEQPREP_OUTPUT}/${SAMPLE}_R2_unmerged.fastq ${SEQPREP_OUTPUT}/${SAMPLE}_unmerged_combined.fastq
	wait
	perl ${PRINSEQ_LITE} -fastq ${SEQPREP_OUTPUT}/${SAMPLE}_unmerged_combined.fastq -out_good ${SEQPREP_OUTPUT}/${SAMPLE}_unmerged_combined.complexity_filtered -out_bad null -lc_method ${COMPLEXITY_METHOD} -lc_threshold ${COMPLEXITY_THRESHOLD} -line_width 0
	wait
	perl ${SPLIT_PAIRED_END_READS} ${SEQPREP_OUTPUT}/${SAMPLE}_unmerged_combined.complexity_filtered.fastq
	wait
	mv ${SEQPREP_OUTPUT}/${SAMPLE}_unmerged_combined.complexity_filtered.fastq_1 ${SEQPREP_OUTPUT}/${SAMPLE}_R1_unmerged.complexity_filtered.fastq
	wait
	mv ${SEQPREP_OUTPUT}/${SAMPLE}_unmerged_combined.complexity_filtered.fastq_2 ${SEQPREP_OUTPUT}/${SAMPLE}_R2_unmerged.complexity_filtered.fastq

# #Compress raw merged and unmerged files and remove intermediates
	wait
	gzip ${SEQPREP_OUTPUT}/${SAMPLE}_*merged.fastq
	wait
	rm -f ${SEQPREP_OUTPUT}/${SAMPLE}_unmerged_combined*.fastq
	wait
done
wait
echo "...Low complexity reads have been filtered" >> ${PROJECT}_progress_file_${DATE}.txt

##################################################

# #BWA -- aligning reads to a reference sequence

# bwa index -p ${REFERENCE_SEQUENCE} -a ${INDEX_ALGORITHM} ${REFERENCE_SEQUENCE}

for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
do
	bwa aln -l ${SEED_DISABLE} -t ${BWA_THREADS} ${REFERENCE_SEQUENCE} ${SEQPREP_OUTPUT}/${SAMPLE}_merged.complexity_filtered.fastq > ${BWA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${REFERENCE_NAME}.sai
	bwa samse ${REFERENCE_SEQUENCE} ${BWA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${REFERENCE_NAME}.sai ${SEQPREP_OUTPUT}/${SAMPLE}_merged.complexity_filtered.fastq > ${BWA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${REFERENCE_NAME}.sam
done
wait
echo "...BWA alignments with merged reads are complete" >> ${PROJECT}_progress_file_${DATE}.txt

for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
do
	bwa aln -l ${SEED_DISABLE} -t ${BWA_THREADS} ${REFERENCE_SEQUENCE} ${SEQPREP_OUTPUT}/${SAMPLE}_R1_unmerged.complexity_filtered.fastq > ${BWA_OUTPUT}/${SAMPLE}_R1_unmerged.complexity_filtered.${REFERENCE_NAME}.sai
	bwa aln -l ${SEED_DISABLE} -t ${BWA_THREADS} ${REFERENCE_SEQUENCE} ${SEQPREP_OUTPUT}/${SAMPLE}_R2_unmerged.complexity_filtered.fastq > ${BWA_OUTPUT}/${SAMPLE}_R2_unmerged.complexity_filtered.${REFERENCE_NAME}.sai
	bwa sampe ${REFERENCE_SEQUENCE} ${BWA_OUTPUT}/${SAMPLE}_R1_unmerged.complexity_filtered.${REFERENCE_NAME}.sai ${BWA_OUTPUT}/${SAMPLE}_R2_unmerged.complexity_filtered.${REFERENCE_NAME}.sai ${SEQPREP_OUTPUT}/${SAMPLE}_R1_unmerged.complexity_filtered.fastq ${SEQPREP_OUTPUT}/${SAMPLE}_R2_unmerged.complexity_filtered.fastq > ${BWA_OUTPUT}/${SAMPLE}_unmerged.complexity_filtered.paired_end.${REFERENCE_NAME}.sam
done
wait
echo "...BWA alignments with unmerged reads are complete" >> ${PROJECT}_progress_file_${DATE}.txt

##################################################

# #SAMtools processing

for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
do

# #Convert SAM to BAM
	samtools view -q${BAM_MIN_QUALITY} -bSh ${BWA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${REFERENCE_NAME}.sam > ${BWA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${REFERENCE_NAME}.bam
	samtools view -q${BAM_MIN_QUALITY} -bSh ${BWA_OUTPUT}/${SAMPLE}_unmerged.complexity_filtered.paired_end.${REFERENCE_NAME}.sam > ${BWA_OUTPUT}/${SAMPLE}_unmerged.complexity_filtered.paired_end.${REFERENCE_NAME}.bam
	wait

# #Sort BAM files
	samtools sort ${BWA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${REFERENCE_NAME}.bam -o ${BWA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${REFERENCE_NAME}.sorted.bam
	samtools sort ${BWA_OUTPUT}/${SAMPLE}_unmerged.complexity_filtered.paired_end.${REFERENCE_NAME}.bam -o ${BWA_OUTPUT}/${SAMPLE}_unmerged.complexity_filtered.paired_end.${REFERENCE_NAME}.sorted.bam
	wait

# #Merge BAM files
	samtools merge ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.bam ${BWA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${REFERENCE_NAME}.sorted.bam ${BWA_OUTPUT}/${SAMPLE}_unmerged.complexity_filtered.paired_end.${REFERENCE_NAME}.sorted.bam
	wait

# #Sort and index all_reads BAM file
	samtools sort ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.bam -o ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.bam
	samtools index ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.bam
	wait

# #Remove duplicates and index
	samtools rmdup -S ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.bam ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.bam
	samtools index ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.bam
	wait

# #Generate statistics like number mapped, duplicates, and number aligned to each chromosome
	samtools flagstat ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.bam >& ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.flagstats.txt
	samtools flagstat ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.bam >& ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.flagstats.txt
	samtools idxstats ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.bam >& ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.idxstats.txt
	samtools idxstats ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.bam >& ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.idxstats.txt

done
wait
echo "...Alignment processing using SAMtools is complete" >> ${PROJECT}_progress_file_${DATE}.txt

##################################################

# #Fragment length distributions from raw merged files and from duplicate-removed BAM files

for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
do

# #Compute fragment length distributions of all reads in the library (mapped and unmapped)
	awk '{print length($10)}' ${BWA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${REFERENCE_NAME}.sam | sort -n | uniq -c | awk '{print $1"\t"$2}' | tail -n +2 > ${BWA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${REFERENCE_NAME}.sam.frag_lengths.txt
	wait

# #Compute fragment length distributions of all reads in the library (mapped)
	samtools view ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.bam | ${GET_INSERT_SIZE} -s ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.PE_insert_output.txt -r ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.SE_read_length_output.txt -
	wait
	cat ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.SE_read_length_output.txt ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.PE_insert_output.txt > ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.rmdup.read_and_insert_lengths.txt
	wait
	done
wait
echo "...Fragment length distribution data has been computed" >> ${PROJECT}_progress_file_${DATE}.txt

##################################################

# #mapDamage to assess damage rates from all aligned and duplicate-removed data and to draw fragment length distributions of aligned data

for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
do
	mapDamage -i ${BWA_OUTPUT}/${SAMPLE}_all_reads.complexity_filtered.${REFERENCE_NAME}.sorted.bam -r ${REFERENCE_SEQUENCE} --merge-reference-sequences -l 200 -d ${PROCESSING_OUTPUT}/MapDamage_output/mapDamage_${SAMPLE} -y 0.5 -m 25 -t ${SAMPLE} --no-stats
wait
done
wait
echo "...mapDamage plot created"  >> ${PROJECT}_progress_file_${DATE}.txt

##################################################

# #Remove intermediate BWA files to save space on edser

for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
do
	rm -f ${BWA_OUTPUT}/${SAMPLE}*${REFERENCE_NAME}.sam
	wait
	rm -f ${BWA_OUTPUT}/${SAMPLE}*${REFERENCE_NAME}.sai
	wait
	rm -f ${BWA_OUTPUT}/${SAMPLE}*merged*${REFERENCE_NAME}.bam
	wait
	rm -f ${BWA_OUTPUT}/${SAMPLE}*all_reads*${REFERENCE_NAME}.bam
	wait
done
wait
echo "...Intermediate BWA files have been removed" >> ${PROJECT}_progress_file_${DATE}.txt

##################################################

# #MEGAN and MIA -- Initial data processing

for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
do

# 3Concatenate merged and unmerged reads - MEGAN and MIA
	cat ${SEQPREP_OUTPUT}/${SAMPLE}_merged.complexity_filtered.fastq ${SEQPREP_OUTPUT}/${SAMPLE}_R1_unmerged.complexity_filtered.fastq ${SEQPREP_OUTPUT}/${SAMPLE}_R2_unmerged.complexity_filtered.fastq > ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fastq
	wait

# 3Convert FASTQ to FASTA - MEGAN and MIA (capture only)
  # ${FASTX_TOOLKIT}/fastq_to_fasta/fastq_to_fasta -n -Q33 -r -i ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fastq -o ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fasta
  fastq_to_fasta -n -Q33 -r -i ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fastq -o ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fasta
  # seqtk seq -Q33 -a ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fastq > ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fasta
  wait

# #Collapse identical reads - forward, reverse, and 5' duplicates - MEGAN and MIA (capture only)
	perl ${PRINSEQ_LITE} -fasta ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fasta -out_good ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed -out_bad null -derep 124 -line_width 0
	wait

# #Remove FASTA and compress FASTQ files
	rm -f ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fasta
	wait
# 	gzip ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fastq
done
wait
echo "...File setup for MEGAN and/or MIA is complete" >> ${PROJECT}_progress_file_${DATE}.txt

##################################################

#MEGAN -- metagenomic analysis. What is the taxonomic makeup of my sequence data?
export BLASTDB="/redser4/resources/blastdb/v5"
for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
do

#Compare reads to BLAST database
	blastn -num_threads ${MEGAN_THREADS} -query ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.fasta -db ${BLAST_DATABASE} -outfmt 6 -out ${MEGAN_OUTPUT}/${SAMPLE}_all_seqprep.duplicates_removed.BLAST.txt
done
wait
echo "...Comparing reads to BLAST database is done." >> ${PROJECT}_progress_file_${DATE}.txt

##################################################

# #MIA --- Option 1 -- creating a consensus using iterative mapping - use this with shotgun data
#
# for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
# do
#  	/soe/pheintzman/bin/mia-1.0/src/mia -r ${MIA_REFERENCE_SEQUENCE} -f ${SEQPREP_OUTPUT}/${SAMPLE}_merged.complexity_filtered.fastq -c -C -U -s ${ANCIENT_DNA_MATRIX} -i -F -k 14 -m ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln
#  	wait
#  	gzip ${SEQPREP_OUTPUT}/${SAMPLE}_merged.complexity_filtered.fastq
#   wait
# done
#
# for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
# do
#  	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.* -f 3 > ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.mia_stats.txt
#  	wait
#  	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.* -f 2 > ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.mia_coverage_per_site.txt
#  	wait
#  	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.* -f 5 > ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.mia_consensus.fasta
#  	wait
#  	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.* -f 41 > ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.inputfornext.txt
#  	wait
#  	perl ${MIA_COVERAGE_FILTER_ANDRE} -c 3 -p 0.67 -I ${SAMPLE}_3x_0.67 <${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.inputfornext.txt >${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.mia_consensus.3x_0.67_filtered.fasta
#   	wait
#  	perl ${MIA_COVERAGE_FILTER_ANDRE} -c 10 -p 0.9 -I ${SAMPLE}_10x_0.9 <${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.inputfornext.txt >${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.mia_consensus.10x_0.9_filtered.fasta
#
# done
# wait
# echo "...Comparing reads to BLAST database is done." >> ${PROJECT}_progress_file_${DATE}.txt

# #MIA --- Option 2 --- creating a consensus using iterative mapping - use this with capture data
#
# for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
# do
#	/soe/pheintzman/bin/mia-1.0/src/mia -r ${MIA_REFERENCE_SEQUENCE} -f ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.fasta -c -C -U -s ${ANCIENT_DNA_MATRIX} -i -F -k 14 -m ${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln
#	wait
# 	gzip ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.fasta
# 	wait
#done
#
# for SAMPLE in $(cat ${PROJECT}-${DATE}-samples.txt)
# do
#	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.* -f 3 > ${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.F.mia_stats.txt
#	wait
#	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.* -f 2 > ${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.F.mia_coverage_per_site.txt
#	wait
#	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.* -f 5 > ${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.F.mia_consensus.fasta
#	wait
#	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.* -f 41 > ${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.F.inputfornext.txt
#	wait
#	perl ${MIA_COVERAGE_FILTER_ANDRE} -c 3 -p 0.67 -I ${SAMPLE}_3x_0.67 <${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.F.inputfornext.txt >${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.F.mia_consensus.3x_0.67_filtered.fasta
# 	wait
#	perl ${MIA_COVERAGE_FILTER_ANDRE} -c 10 -p 0.9 -I ${SAMPLE}_10x_0.9 <${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.F.inputfornext.txt >${MIA_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed.${MIA_REFERENCE_NAME}.20M.maln.F.mia_consensus.10x_0.9_filtered.fasta
#done
#wait
# echo "...MIA analyses are complete" >> ${PROJECT}_progress_file_${DATE}.txt

cd ${PROCESSING_OUTPUT}

gzip ${SEQPREP_OUTPUT}/${SAMPLE}_*.fastq
gzip ${SEQPREP_OUTPUT}/${SAMPLE}_*.fasta
gzip ${MEGAN_OUTPUT}/${SAMPLE}_all_seqprep.duplicates_removed.BLAST.txt

echo "...Large files are gzipped." >> ${PROJECT}_progress_file_${DATE}.txt

Rscript ${CALC_STATS} BWA_analyses/ SeqPrep_output/ MapDamage_output/ # MIA_analyses/

echo "...Summary statistics are estimated!" >> ${PROJECT}_progress_file_${DATE}.txt

echo "... End." >> ${PROJECT}_progress_file_${DATE}.txt

mv ${PROJECT}-${DATE}-samples.txt ${PROCESSING_OUTPUT}/Sample_lists_and_progress_files
mv ${PROJECT}_progress_file_${DATE}.txt ${PROCESSING_OUTPUT}/Sample_lists_and_progress_files
