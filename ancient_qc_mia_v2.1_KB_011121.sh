#!/bin/bash

# 11/01/2021
#from Jonas
# Ancient QC script for nextseq data ONLY - naming scheme will not work for NovaSeq
# Script inspired by the Pete qc script, but a bit pared down
# Updated for redser4, with a few changes
# Hardcoded running mia with A. phoeniceus mt, for now
#arguments: prefix for naming files,
            #file with list of library names up to illumina sample S number (eg lib_S1), one on each line,
            #path to the raw directory
            #and the path to the reference genome to map to

#to make executable: chmod +x filename.sh

##to run ./SCRIPT.sh |& tee -a log.txt

# inputs

if [ "$#" -ne 4 ]; then
        >&2 echo "Args are PROJECT_NAME LIBRARY_LIST RAWDIR REFERENCE" && exit 1
fi

# file names and raw data directory
PROJECT=$1
shift
LIBS=$1 # file of library/sample names, before the lane information starts in the illumina sile naming scheme, but after the S info
shift
RAWDIR=$1
shift
REFERENCE=$1

mt_ref=/redser4/genomes/mitogenomes/Agelaius_phoeniceus_NC_018801.1.fasta

REF_BASE=${REFERENCE##*/}
REF_NAME=${REF_BASE%.*}

######
# Params and locations
######

#SeqPrep2 Variables
SEQPREP_MIN_LENGTH=30    #Min length to merge <-VARIABLE -L
SEQPREP_OVERLAP=15       #Overlap variable for merging <-VARIABLE -o
SEQPREP_MINQUAL=15       #Minumum quality score
MISMATCH_FRAC=0.05

#BWA and SAMtools Variables
BWA_THREADS=6
n=6

MAPDAMAGE2=/redser4/personal/jonas/bin/mapDamage
MAX_MISINCORP_FREQUENCY=0.5   # Use 0.3 if not too badly damaged, use 0.5 if badly damaged
READ_PLOT_LENGTH=25           # The number of nucleotides to plot at the 5' and 3' ends of the read

PICARD=/home/jooppenh/src/picard.jar

#####

mkdir -p $(date +%y%m%d)_${PROJECT}_${REF_NAME}

if [[ ! -f $(date +%y%m%d)_${PROJECT}_${REF_NAME}/$(basename "$LIBS") ]]; then
        mv "$LIBS" $(date +%y%m%d)_${PROJECT}_${REF_NAME}/$(date +%y%m%d)_${PROJECT}_${REF_NAME}.library_list.txt
fi

cd $(date +%y%m%d)_${PROJECT}_${REF_NAME} || exit

while read -r lib; do

    s=${lib%S*}
    samp=${s%_*}

    echo ${samp} >> sample_list.txt

    ###SeqPrep2 - Removing adapters and merging paired end reads
    mkdir -p ${samp} && cd ${samp} || exit

    /redser4/personal/jonas/bin/SeqPrep2 -f ${RAWDIR}/${lib}_R1_001.fastq.gz -r ${RAWDIR}/${lib}_R2_001.fastq.gz \
        -1 ${samp}_R1_unmerged.fastq.gz -2 ${samp}_R2_unmerged.fastq.gz \
        -q ${SEQPREP_MINQUAL} -L ${SEQPREP_MIN_LENGTH} \
        -A AGATCGGAAGAGCACACGTC -B AGATCGGAAGAGCGTCGTGT \
        -s ${samp}_merged.fastq.gz \
        -o ${SEQPREP_OVERLAP} -m ${MISMATCH_FRAC} \
        -C ATCTCGTATGCCGTCTTCTGCTTG -D GATCTCGGTGGTCGCCGTATCATT -S >& ${samp}.seqprep_output.txt

    # won't remove these yet
    #rm ${lib}_L00${lane}_R1_unmerged.fastq.gz ${lib}_L00${lane}_R2_unmerged.fastq.gz

    ###BWA - Aligning reads to a reference sequence
    merged_read_group='@RG\tID:'${samp}'_merged\tSM:'${samp}'\tPL:illumina\tLB:'${samp}'\tPU:'${samp}'_merged'
    unmerged_read_group='@RG\tID:'${samp}'_unmerged\tSM:'${samp}'\tPL:illumina\tLB:'${samp}'\tPU:'${samp}'_unmerged'

    # updated parameters to edit distance/gap penalty/disable seeding
    bwa samse ${REFERENCE} \
        <(bwa aln -l 16500 -n 0.01 -o 2 -t ${BWA_THREADS} ${REFERENCE} ${samp}_merged.fastq.gz) \
        ${samp}_merged.fastq.gz -r ${merged_read_group} | samtools sort -@ ${n} -o ${samp}_merged.bam


    bwa sampe ${REFERENCE} \
        <(bwa aln -l 16500 -n 0.01 -o 2 -t ${BWA_THREADS} ${REFERENCE} ${samp}_R1_unmerged.fastq.gz) \
        <(bwa aln -l 16500 -n 0.01 -o 2 -t ${BWA_THREADS} ${REFERENCE} ${samp}_R2_unmerged.fastq.gz) \
        ${samp}_R1_unmerged.fastq.gz ${samp}_R2_unmerged.fastq.gz -r ${unmerged_read_group} | samtools sort -@ ${n} -o ${samp}_unmerged.bam


    samtools flagstat -@ ${n} ${samp}_merged.bam > ${samp}_merged.flagstat.txt
    samtools flagstat -@ ${n} ${samp}_unmerged.bam > ${samp}_unmerged.flagstat.txt

    # just treating as implicit that the final bam is composed of both merged and unmerged reads- not indicating this with special filename handle like '_all' or something
    # essentially just viewing merged and unmerged fq/bams as temporary files to generated the combined bam
    # eventually should add option to merge/not merge and treat separately
    samtools merge -@ ${n} ${samp}.bam ${samp}_merged.bam ${samp}_unmerged.bam && rm ${samp}_merged.bam ${samp}_unmerged.bam

    samtools flagstat -@ ${n} ${samp}.bam > ${samp}.flagstat.txt
    samtools stats -@ ${n} ${samp}.bam > ${samp}.stats.txt

    java -jar -Xmx4g -Djava.io.tmpdir=$(pwd)/tmp -XX:ParallelGCThreads=8 $PICARD MarkDuplicates \
        INPUT="${samp}".bam OUTPUT="${samp}".rmdup.bam METRICS_FILE="${samp}".markdup_metrics.txt \
        REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 VALIDATION_STRINGENCY=LENIENT TMP_DIR=$(pwd)/tmp

    rmdir tmp

    ${MAPDAMAGE2} -i "${samp}".rmdup.bam -r ${REFERENCE} -l 150 -d mapdamage_${samp} -y ${MAX_MISINCORP_FREQUENCY} --merge-reference-sequences -m ${READ_PLOT_LENGTH} -t ${REF_NAME}_${samp}

    # mia
    /redser4/personal/jonas/scripts/run_mia.sh ${samp}_merged.fastq.gz ${mt_ref} Atricolor

    # plot read lengths for merged reads and insert size for the paired reads
    # this is only for mapped reads of course, maybe also want a read length histogram from the fastqs too
    samtools view -h -F 5 "${samp}".rmdup.bam | samtools stats -@ ${n} | /redser4/personal/jonas/scripts/readLengthPlotter.py /dev/stdin "${samp}"_merged_read_lengths.png

    samtools view -h -f 3 -F 4 "${samp}".rmdup.bam | samtools stats -@ ${n} | /redser4/personal/jonas/scripts/readLengthPlotter.py /dev/stdin "${samp}"_unmerged_insert_size.png --use_insert_size

    cd ..

done < $(date +%y%m%d)_${PROJECT}_${REF_NAME}.library_list.txt

# generate qc report
/redser4/personal/jonas/scripts/qcReportGenerator.py sample_list.txt $(date +%y%m%d)_${PROJECT}_${REF_NAME}.qcReport.txt && rm sample_list.txt
