#!/bin/bash
#$ -S /bin/bash
#$ -o smrtpipe.outputlog.txt
#$ -e smrtpipe.errorlog.txt
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=20G,scratch=1G
#$ -l h_rt=96:0:0
#$ -pe smp 12

#currently running coupled 12 cores

SEQFOLDER=/netapp/home/jbisanz/SMRT_joel/UCB02 #Folder containing your results from SMRT (particularily the *.bax.h5 files
#PROTOCOL=/netapp/home/jbisanz/SMRT/current/common/protocols/RS_HGAP_Assembly.3.xml #problem with supplied, instead using one generated from portal

PROTOCOL=/netapp/home/jbisanz/SMRT/current/common/protocols/JB_RS_HGAP_Assembly.3.xml


mkdir ${SEQFOLDER}/HGAP3_Assembly
find $SEQFOLDER|grep bax.h5 >smrtcells.fofn 
/netapp/home/jbisanz/SMRT/install/smrtanalysis_2.3.0.140936/analysis/bin/fofnToSmrtpipeInput.py smrtcells.fofn > ${SEQFOLDER}/input.xml

/netapp/home/jbisanz/SMRT/smrtcmds/bin/smrtpipe  -D NPROC=$NSLOTS --output=${SEQFOLDER}/HGAP3_Assembly --params=${PROTOCOL} xml:${SEQFOLDER}/input.xml

