#!/usr/bin/env Rscript
message("16STrimme.R v0.1 jbisanz")
message("Overlap sanger 16S reads and assign taxonomy to identify cultured isolates")
message("Requires the following packages from bioconductor: dada2, ShortRead, sangerseqR, CrispRVariants")


suppressMessages(library(optparse))
option_list = list(
  make_option(c("-s", "--Sample"), type="character", help="Isolate's name", metavar="character"),
  make_option(c("-f", "--Forward"), type="character", help="Forward sanger read (.ab1 format)", metavar="character"),
  make_option(c("-r", "--Reverse"), type="character", help="Reverse sanger read (.ab1 format)", metavar="character"),
  make_option(c("-o", "--OutFile"), type="character", default="Merged.16S.fasta", help="file with overlapped oriented 16S [default= %default]", metavar="character"),
  make_option(c("-S", "--SeqStrategy"), type="character", default="in", help="Sequencing Strategy- in refers to sequencing in from the ends (ex 8F 1543R) while out refers to sequencing out from a small overlap internally (ex V4F V4R) [default= %default]", metavar="character"),
  make_option(c("-5", "--Trim5"), type="numeric", default="30", help="How many bases to trim off start of reads [default= %default]", metavar="numeric"),
  make_option(c("-Q", "--truncQ"), type="numeric", default="10", help="Truncate after quality score less than this [default= %default]", metavar="numeric"),
  make_option(c("-N", "--maxN"), type="numeric", default="0", help="Number of N's allowed in final sequence [default= %default]", metavar="numeric"),
  make_option(c("-M", "--Append"), type="logical", default=FALSE, help="Should ouput be appended to an existing file? [default= %default]", metavar="logical"),
  make_option(c("-T", "--IDTaxa"), type="logical", default=TRUE, help="Should taxonomy be assigned? [default= %default]", metavar="logical"),
  make_option(c("-G", "--TaxDBGenus"), type="character", default="/Volumes/turnbaughlab/qb3share/shared_resources/databases/dada2_training_sets/silva_nr_v128_train_set.fa.gz", help="Dada2-formatted taxonomy database. The default is available when biggut is mounted on OSX: [default= %default]", metavar="character"),
  make_option(c("-D", "--TaxDBSpecies"), type="character", default="/Volumes/turnbaughlab/qb3share/shared_resources/databases/dada2_training_sets/silva_species_assignment_v128.fa.gz", help="Dada2-formatted species database. The default is available when biggut is mounted on OSX: [default= %default]", metavar="character"),
  make_option(c("-z", "--tmp"), type="character", default=tempdir(), help="Temporary Directory: [default= %default]", metavar="character")

)
opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)
if (is.null(opts$Sample) | is.null(opts$Forward) | is.null(opts$Reverse)){
  print_help(opt_parser)
  message("Example usage: ./16STrimme.R --Sample ElSeco --Forward Elseco_8F.ab1 --Reverse Elseco_1543R.abi --OutFile Elseco.fa")
  stop("At least the Sample, and Forward and Reverse must be specified.\n", call.=FALSE)
}

suppressMessages(library(CrispRVariants))
suppressMessages(library(sangerseqR))
suppressMessages(library(dada2))
suppressMessages(library(ShortRead))

unlink(paste0(opts$tmp,"/",opts$Sample, "_R1.fastq"))
unlink(paste0(opts$tmp,"/",opts$Sample, "_R2.fastq"))

sink<-abifToFastq(seqname=opts$Sample, fname=opts$Forward,  outfname=paste0(opts$tmp,"/",opts$Sample, "_R1.fastq"))
sink<-abifToFastq(seqname=opts$Sample, fname=opts$Reverse,  outfname=paste0(opts$tmp,"/",opts$Sample, "_R2.fastq"))
rm(sink)

fastqPairedFilter(fn=c(paste0(opts$tmp,"/",opts$Sample, "_R1.fastq"),paste0(opts$tmp,"/",opts$Sample, "_R2.fastq")), 
                  fout=c(paste0(opts$tmp,"/",opts$Sample, "_R1.filt.fastq"),paste0(opts$tmp,"/",opts$Sample, "_R2.filt.fastq")),
                  trimLeft=opts$Trim5,
                  truncQ=opts$truncQ,
                  maxN=opts$maxN)

R1=readFastq(paste0(opts$tmp,"/", opts$Sample, "_R1.filt.fastq"))
R2=readFastq(paste0(opts$tmp,"/", opts$Sample, "_R2.filt.fastq"))


align<-pairwiseAlignment(pattern=R1@sread, subject=reverseComplement(R2@sread), type="local")

if(pid(align)!=100){stop("Error, the overlap between reads is not perfect. Please manually inspect and correct .abi files!")}

if(opts$SeqStrategy=="in"){
  OLSeq<-DNAStringSet(
    paste0(
      as.character(subseq(R1@sread, start = 1, width=(align@pattern@range@start+align@pattern@range@width-1))),
      as.character(subseq(reverseComplement(R2@sread), start = align@subject@range@width+1, end=nchar(R2@sread)))
        )
  )
} else {
  OLSeq<-DNAStringSet(
    paste0(
      as.character(subseq(reverseComplement(R2@sread), start = 1, end=align@subject@range@start)),
      as.character(subseq(R1@sread, start=1, end=nchar(R1@sread)))
    )
  )
}

if(opts$IDTaxa){
  message("Assigning taxonomy using ", basename(opts$TaxDBGenus))
  message("This may take some time and can be bypassed with --IDTaxa==FALSE")
  tax<-assignTaxonomy(as.character(OLSeq), refFasta=opts$TaxDBGenus)
  tax<-addSpecies(tax, refFasta=opts$TaxDBSpecies, allowMultiple = TRUE, tryRC=TRUE)
  names(OLSeq)<-paste0(opts$Sample, " Taxonomy=",paste(tax[1,], collapse=","))
} else {
  names(OLSeq)<-opts$Sample
}

writeXStringSet(OLSeq, opts$OutFile, append = opts$Append)


