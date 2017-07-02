#!/bin/bash
#$ -S /bin/bash
#$ -o rnaseq.outputlog.txt
#$ -e errorlog.txt
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=20G,scratch=10G
#$ -l h_rt=336:0:0
#$ -pe smp  4
#$ -t 2-50 #Skipping Job1 as it is the header line


#previously run for each:
# for REFERENCE in references/*.fna; do samtools faidx $REFERENCE; done
#previously run for each: samtools faidx $REFERENCE
# for REFERENCE in references/*.fna; do bowtie2-build $REFERENCE $REFERENCE.idx; done; mv references/*idx* index/      

FILE="VMR_ENB_samples.txt" #no spaces allowed in any category with tab delimited information about how to sequence

DIR=/scrapp2/VMR_ENB_RNASeq_June26 #location of a directory where all of your files will be
READDIR=$DIR/reads #location of your reads, I am assuming the directory will be called reads later, also assuming reads have been compressed with .gz

CURRJOB=$(sed "${SGE_TASK_ID}q;d" $FILE) # get the array job numbereth line of the file

SAMPLEID=$(echo $CURRJOB | cut -f1 -d' ')
REFERENCE=$(echo $CURRJOB | cut -f3 -d' ')
GTF=$(echo $CURRJOB | cut -f8 -d' ')
INDEX=$(echo $CURRJOB | cut -f9 -d' ')
FOR=$(echo $CURRJOB | cut -f4 -d' ')
REV=$(echo $CURRJOB | cut -f5 -d' ')
STRANDED=$(echo $CURRJOB | cut -f10 -d' ')
PAIRED=$(echo $CURRJOB | cut -f11 -d' ')
FORWARDTRIM=$(echo $CURRJOB | cut -f12 -d' ')
REVERSETRIM=$(echo $CURRJOB | cut -f13 -d' ')
MINQUAL=$(echo $CURRJOB | cut -f15 -d' ')


exec >logs/SGE${SGE_TASK_ID}.${SAMPLEID}.log 2>logs/SGE${SGE_TASK_ID}.${SAMPLEID}.errors

date
hostname

echo -e "Analysing: $SAMPLEID \nReference:$REFERENCE\nGTF:$GTF\nIndex:$INDEX\nForward Read:$FOR\nReverse Read:$REV\nStranded:$STRANDED\nPaired:$PAIRED\nForwardTrim:$FORWARDTRIM\nReverseTrim:$REVERSETRIM\nMinQual:$MINQUAL\n"

source /netapp/home/jbisanz/.bash_profile
source /netapp/home/jbisanz/qiime_for_cluster_1.9.1/settings/sourceme.sh

mkdir -p fastqc/$SAMPLEID

echo "$(date)  Doing QC and Aligning $SAMPLEID with bowtie2 2.2.6 and v0.11.4"
if [ $PAIRED == "yes" ]; then
	   /netapp/home/jbisanz/sftwr/FastQC/fastqc --outdir fastqc/$SAMPLEID $READDIR/$FOR
	   /netapp/home/jbisanz/sftwr/FastQC/fastqc --outdir fastqc/$SAMPLEID $READDIR/$REV
        bowtie2 -p $NSLOTS --trim5 $FORWARDTRIM --trim3 $REVERSETRIM --end-to-end --fr --sensitive -x index/$INDEX -1 $READDIR/$FOR -2 $READDIR/$REV -S alignments/$SAMPLEID.sam
elif [ $PAIRED == "no" ]; then
	   /netapp/home/jbisanz/sftwr/FastQC/fastqc --outdir fastqc/$SAMPLEID $READDIR/$FOR
        bowtie2 -p $NSLOTS --trim5 $FORWARDTRIM --trim3 $REVERSETRIM --end-to-end --sensitive -x index/$INDEX -U $READDIR/$FOR -S alignments/$SAMPLEID.sam
fi

echo "$(date)   Baming and sorting $SAMPLEID"
samtools view -bS --threads $NSLOTS alignments/$SAMPLEID.sam > alignments/$SAMPLEID.bam
samtools sort --threads $NSLOTS alignments/$SAMPLEID.bam > alignments/$SAMPLEID.sorted.bam
samtools index -b alignments/$SAMPLEID.sorted.bam alignments/$SAMPLEID.sorted.bai

echo "$(date)   Calling variants on $SAMPLEID"
samtools mpileup -uf references/$REFERENCE alignments/$SAMPLEID.sorted.bam | bcftools call -mv  > alignments/$SAMPLEID.vcf

echo "$(date)   Counting $SAMPLEID with HTSeq-0.8.0"
htseq-count -f bam --order pos --stranded $STRANDED --minaqual $MINQUAL --type CDS --order pos --idattr ID --nonunique all alignments/$SAMPLEID.sorted.bam references/$GTF > counts/$SAMPLEID.counts

echo "$(date)   Removing SAM file"
rm alignments/$SAMPLEID.sam 
