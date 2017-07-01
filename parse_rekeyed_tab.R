suppressMessages(require(ShortRead))
suppressMessages(require(data.table))

#setwd("~/Desktop/fibs_fastqs/")
ReKey<-"rekeyed_tab.txt"
ForRead<-"Burton2_S1_L001_R1_001.fastq.gz"
RevRead<-"Burton2_S1_L001_R2_001.fastq.gz"
outdir<-"demultiplexed"

message(paste("Reading", ReKey))
keys<-fread(ReKey, header=F, select=c(1,2,6))
colnames(keys)<-c("SeqID","SampleID","BarCodes")

counts<-data.table(table(keys$SampleID))
colnames(counts)<-c("SampleID","NinReKeyed")
counts$ForReads<-0
counts$RevReads<-0
counts<-rbind(counts, data.table(SampleID="NotAssigned",NinReKeyed=0,ForReads=0,RevReads=0))

keys$SeqID<-gsub("^@|:1$","", keys$SeqID)#appears that in the rekeyed, an extra :1 is added on to the read IDs

#Get paired reads in blocks
dir.create(outdir, showWarnings = F)
i=0
NSTREAM=1e6 #work on NSTREAM reads at a time rather than trying to load them all into RAM at the same time
streamF<-FastqStreamer(ForRead, n=NSTREAM) 
streamR<-FastqStreamer(RevRead, n=NSTREAM) 
repeat{
  i<-i+1
  message(paste("Processing  Read Chunk:", i))
  fqF<-yield(streamF)
  fqR<-yield(streamR)
  
  if(length(fqF)==0){break}
  
  fors<-data.table(SeqID=as.character(id(fqF)), ForSeq=as.character(sread(fqF)), ForQual=as.character(quality(quality(fqF))))
  fors$SeqID<-gsub(" ..+","",fors$SeqID)
  
  revs<-data.table(SeqID=as.character(id(fqR)), RevSeq=as.character(sread(fqR)), RevQual=as.character(quality(quality(fqR))))
  revs$SeqID<-gsub(" ..+","",revs$SeqID)
  
  if(sum(fors$SeqID!=revs$SeqID)>0){stop("ERROR READS NOT MATCHED")}
  merge<-cbind(fors,revs[,-1]) #going straight cbind for speed after checking they match
  merge<-merge(keys, merge, all.x=F, all.y=T, by="SeqID")
  merge$NewHeader<-paste(merge$SeqID, merge$SampleID, merge$BarCodes)
  
  merge[is.na(merge$SampleID)]$SampleID<-"NotAssigned"
  persamp<-split(merge, list(merge$SampleID))

  counts$ForReads<-counts$ForReads + sapply(counts$SampleID, function(x) NROW(persamp[[x]]))

  for(sample in names(persamp)){
    current<-persamp[[sample]]
    writeFastq(ShortReadQ(sread=DNAStringSet(current$ForSeq), quality=BStringSet(current$ForQual), id=BStringSet(current$NewHeader)), mode="a", paste0(outdir,"/",sample,"_1.fastq"), compress=F)
    writeFastq(ShortReadQ(sread=DNAStringSet(current$RevSeq), quality=BStringSet(current$RevQual), id=BStringSet(current$NewHeader)), mode="a", paste0(outdir,"/",sample,"_2.fastq"), compress=F)
  }
  
  message(paste(date(),"Processed", i*NSTREAM, "reads..."))
  
}
close(streamF)
close(streamR)

fwrite(counts, "parsing_summary.txt", sep='\t')
