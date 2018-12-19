#!/bin/bash
#$ -S /bin/bash
#$ -o assembly.log
#$ -e assembly.err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=16G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=20G
#$ -l h_rt=336:0:0
#$ -pe smp  8
#$ -t 3-34
#in the line above we do not use SGE TASK ID 1 because this is the header of the file, need 128Gb per job (8x16)

# Input is taken from a tsv file with the following format that should hopefully be intuitive. Up to 2 PE libraries are accepted in this version.
#Genome_Name	Directory	Genus	Species	Strain	GramStain	L1_F	L1_R	L2_F	L2_R
#Bacteroides_intestinalis_LPF21	GenomesForAssembly/Bacteroides_intestinalis_LPF21/	Bacteroides	intestinalis	LPF21	neg	Bacteroides_F_2_1_S285_R1_001.fastq.gz	Bacteroides_F_2_1_S285_R2_001.fastq.gz		
#Clostridioides_difficile_JBZPo1	GenomesForAssembly/Clostridioides_difficile_JBZPo1/	Clostridioides	difficile	JBZPo1	pos	CDGN_JBZPo1_S216_R1_001.fastq.gz	CDGN_JBZPo1_S216_R2_001.fastq.gz	JBZPo1_S1_L001_R1_001.fastq.gz	JBZPo1_S1_L001_R2_001.fastq.gz
#Ruminococcus_sp_LPE32	GenomesForAssembly/Ruminococcus_sp_LPE32/	Ruminococcus	sp	LPE32	pos	Ruminococcus_E_3_2_S292_R1_001.fastq.gz	Ruminococcus_E_3_2_S292_R2_001.fastq.gz		


#Import relevant variables
JOB=$( sed "${SGE_TASK_ID}q;d" genomes.txt )
Genome_Name=$( echo $JOB | cut -d' ' -f1 )
RDirectory=$( echo $JOB | cut -d' ' -f2 )
Genus=$( echo $JOB | cut -d' ' -f3 )
Species=$( echo $JOB | cut -d' ' -f4 )
Strain=$( echo $JOB | cut -d' ' -f5 )
GramStain=$( echo $JOB | cut -d' ' -f6 )
L1_F=$( echo $JOB | cut -d' ' -f7 )
L1_R=$( echo $JOB | cut -d' ' -f8 )
L2_F=$( echo $JOB | cut -d' ' -f9 )
L2_R=$( echo $JOB | cut -d' ' -f10 )
HomeDir=$( pwd )
TempDir=/scratch/${Genome_Name}_${SGE_TASK_ID}_${JOB_ID}


#Set up necessary directories
mkdir -p logs
mkdir -p assembled
mkdir $TempDir

echo $(date): Running $Genome_Name on $(hostname)

exec > ${HomeDir}/logs/${Genome_Name}.assembly.JID${JOB_ID}.SGE${SGE_TASK_ID}.log 2> ${HomeDir}/logs/${Genome_Name}.assembly.JID${JOB_ID}.SGE${SGE_TASK_ID}.err


#Import of relevant dependencies
source /turnbaugh/qb3share/shared_resources/export_python_2.7.14.sh
export PATH=/turnbaugh/qb3share/shared_resources/bin/:$PATH
export PATH=/turnbaugh/qb3share/shared_resources/sftwrshare/perl-5.28.0/bin:$PATH
export PATH=/turnbaugh/qb3share/shared_resources/sftwrshare/ncbi-blast-2.6.0+/bin/:$PATH
export PATH=/turnbaugh/qb3share/shared_resources/sftwrshare/hmmer-3.1b2-linux-intel-x86_64/binaries/$PATH
export PATH=/turnbaugh/qb3share/shared_resources/sftwrshare/aragorn1.2.36/:$PATH
export PATH=/turnbaugh/qb3share/shared_resources/sftwrshare/parallel-20170822/bin/:$PATH
export PATH=/turnbaugh/qb3share/shared_resources/sftwrshare/infernal-1.1.1-linux-intel-gcc/binaries/:$PATH
export PATH=/turnbaugh/qb3share/shared_resources/sftwrshare/barrnap-0.6/bin/:$PATH
export PATH=/turnbaugh/qb3share/shared_resources/sftwrshare/minced/:$PATH
export PATH=/turnbaugh/qb3share/shared_resources/sftwrshare/signalp-4.1/:$PATH
export PATH=/turnbaugh/qb3share/shared_resources/sftwrshare/bowtie2-2.3.4.1-linux-x86_64/:$PATH

echo $(date): Moving files to local scratch drive on $(hostname)
cp -r $RDirectory ${TempDir}/reads
cd $TempDir



####################################################################
#QC- using a new tool called fastp to do identification of adapters and quality control
echo $(date):Filtering reads
mkdir -p filtered

fastp \
 --in1 reads/$L1_F \
 --in2 reads/$L1_R \
 --out1 filtered/$L1_F \
 --out2 filtered/$L1_R \
 --trim_front1 15 \
 --trim_front2 15 \
 --cut_by_quality3 \
 --cut_window_size 4 \
 --cut_mean_quality 20 \
 --length_required 60 \
 --overrepresentation_analysis \
 --detect_adapter_for_pe \
 --json filtered/L1.json \
 --html filtered/L1.html \
 --thread $NSLOTS
 
if [[ $L2_F =~ .*fastq.gz ]]; then
 echo $(date): a second library is being used for assembly
 	fastp \
	 --in1 reads/$L2_F \
	 --in2 reads/$L2_R \
	 --out1 filtered/$L2_F \
	 --out2 filtered/$L2_R \
  --trim_front1 15 \
  --trim_front2 15 \
	 --cut_by_quality3 \
	 --cut_window_size 4 \
	 --cut_mean_quality 20 \
	 --length_required 60 \
	 --overrepresentation_analysis \
	 --detect_adapter_for_pe \
	 --json filtered/L2.json \
	 --html filtered/L2.html \
	 --thread $NSLOTS
fi

####################################################################
#QC- remove any PhiX
echo $(date): Removing PhiX

bowtie2 --local \
	--very-fast-local \
  --threads $NSLOTS \
  -x /turnbaugh/qb3share/shared_resources/databases/phix_bowtie2/phix174 \
  -1 filtered/$L1_F \
  -2 filtered/$L1_R \
  --un-conc filtered/L1.fastq \
  -S phix.aligned.sam
rm phix.aligned.sam
  
if [[ $L2_F =~ .*fastq.gz ]]; then
	bowtie2 --local \
		--very-fast-local \
		--threads $NSLOTS \
		-x /turnbaugh/qb3share/shared_resources/databases/phix_bowtie2/phix174 \
		-1 filtered/$L2_F \
		-2 filtered/$L2_R \
		--un-conc filtered/L2.fastq \
		-S phix.aligned.sam
	rm phix.aligned.sam
fi

####################################################################
#QC- Search for overlaps for use as separate library
echo $(date): Finding overlapped reads
/turnbaugh/qb3share/shared_resources/sftwrshare/vsearch-2.4.4-linux-x86_64/bin/vsearch \
	--threads $NSLOTS \
	--fastq_minovlen 30 \
	--fastq_mergepairs filtered/L1.1.fastq \
	--reverse filtered/L1.2.fastq \
	--fastqout_notmerged_fwd filtered/L1_1.fastq \
 	--fastqout_notmerged_rev filtered/L1_2.fastq \
	--fastqout filtered/L1_3.fastq
	
if [[ $L2_F =~ .*fastq.gz ]]; then
	/turnbaugh/qb3share/shared_resources/sftwrshare/vsearch-2.4.4-linux-x86_64/bin/vsearch \
	--threads $NSLOTS \
	--fastq_minovlen 30 \
	--fastq_mergepairs filtered/L2.1.fastq \
	--reverse filtered/L2.2.fastq \
	--fastqout_notmerged_fwd filtered/L2_1.fastq \
 	--fastqout_notmerged_rev filtered/L2_2.fastq \
	--fastqout filtered/L2_3.fastq
fi
	
####################################################################
#QC- cleaning up
echo $(date): Cleaning up intermediates
rm -r reads
rm filtered/$L1_F
rm filtered/$L1_R
rm filtered/$L2_F
rm filtered/$L2_R
rm filtered/L1.1.fastq
rm filtered/L1.2.fastq

####################################################################
#Assembly
echo $(date):Assembling

if [[ $L2_F =~ .*fastq.gz ]]; then
 /turnbaugh/qb3share/shared_resources/sftwrshare/SPAdes-3.13.0-Linux/bin/spades.py \
	--pe1-1 filtered/L1_1.fastq \
	--pe1-2 filtered/L1_2.fastq \
	--pe1-m filtered/L1_3.fastq \
	--pe2-1 filtered/L2_1.fastq \
	--pe2-2 filtered/L2_2.fastq \
	--pe2-m filtered/L2_3.fastq \
	--cov-cutoff auto \
	--threads $NSLOTS \
	--memory 128 \
	--tmp-dir tmpdir \
	-o assembly
else
/turnbaugh/qb3share/shared_resources/sftwrshare/SPAdes-3.13.0-Linux/bin/spades.py \
	--pe1-1 filtered/L1_1.fastq \
	--pe1-2 filtered/L1_2.fastq \
	--pe1-m filtered/L1_3.fastq \
	--cov-cutoff auto \
	--threads $NSLOTS \
	--memory 128 \
	--tmp-dir tmpdir \
	-o assembly
fi

echo $(date): Cleaning up assembly
pigz filtered/*fastq
rm -r assembly/corrected
rm -r assembly/K*

####################################################################
#Annotating
echo $(date): Annotating with PROKKA 1.12
/turnbaugh/qb3share/shared_resources/sftwrshare/prokka-1.12/bin/prokka \
	--outdir annotation \
	--prefix $Genome_Name \
	--locustag $Genome_Name \
	--genus $Genus \
	--species $Species \
	--strain $Strain \
	--kingdom Bacteria \
	--gram $GramStain \
	--cpus $NSLOTS \
	--mincontiglen 200 \
	assembly/scaffolds.fasta

####################################################################
#Running QC
echo $(date): Running QC

/turnbaugh/qb3share/shared_resources/sftwrshare/quast-4.5/quast.py \
 --output-dir qc \
 --threads $NSLOTS \
 annotation/${Genome_Name}.fna


####################################################################
#Finishing

echo $(date): Moving files back to biggut and cleaning up $TempDir
cp -r $TempDir ${HomeDir}/assembled/${Genome_Name}
rm -r $TempDir
