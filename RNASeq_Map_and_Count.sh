#!/bin/bash
#$ -S /bin/bash
#$ -o outputlog.txt
#$ -e errorlog.txt
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=20G,scratch=10G
#$ -l h_rt=12:0:0
#$ -pe smp  12

date
hostname
source /netapp/home/jbisanz/.bash_profile

########################################################
#Pipeline to: quality filter reads, build a mapping index using single ended reads
#	1. Build a mapping index (Bowtie2-build)
#	2. Map reads to new index (Bowtie2)
#	3. Pull feature counts (HTSeq)
#	4. Build a sorted BAM file with index to allow visualization of mapping with IGV
#
# Adjust settings below to your needs:
#
DIR=/scrapp2/YOURFOLDER #location of a directory where all of your files will be
READDIR=$DIR/comp_reads #location of your reads, I am assuming the directory will be called reads later, also assuming reads have been compressed with .gz
REFERENCE=166-Elenta-DSM2243.fna # a fasta file containing your contigs to map to
GTF=166-Elenta-DSM2243.gtf
######################################################


#Add in filtering at some point in the future when you have fastq files


######################################################
if [ ! -d "$DIR/decomp_reads" ]
then
    	echo "$(date)	Decompressing reads..."
        mkdir decomp_reads
        cp $READDIR/*.gz $DIR/decomp_reads/
        gunzip $DIR/decomp_reads/*
else
    	echo "Reads already decompressed"
fi
READDIR="$DIR/decomp_reads"

echo "$(date)	Making index from $REFERENCE"

mkdir index
bowtie2-build $REFERENCE index/$REFERENCE.idx

    	mkdir alignments
        mkdir counts


for curr_samp in decomp_reads/*.fastq
do

	SAMPLEID=$(echo $curr_samp | sed 's/\.fna//' | sed 's/decomp_reads\///')

  	echo "$(date)	Aligning $SAMPLEID"
                INDEX="index/$REFERENCE.idx"
                bowtie2 -q --met-file bowtie2log.txt -p $NSLOTS --trim5 5 --trim3 5 --end-to-end --sensitive -x $INDEX -U $curr_samp -S alignments/$SAMPLEID.sam
  	echo "$(date)	Counting $SAMPLEID"
                htseq-count --type=CDS --idattr=ID --stranded=yes --minaqual=10  alignments/$SAMPLEID.sam  $GTF > counts/$SAMPLEID.counts

  	echo "$(date)	Baming and sorting $SAMPLEID"
		samtools view -bS alignments/$SAMPLEID.sam > alignments/$SAMPLEID.bam
		samtools sort alignments/$SAMPLEID.bam  alignments/$SAMPLEID.sorted
        samtools index -b alignments/$SAMPLEID.sorted.bam alignments/$SAMPLEID.sorted.bai
done

rm -r $DIR/decomp_reads