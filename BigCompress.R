#BigCompress.R
#Find files of a certain type within a target directory, and either create a list, or compress
#Intended use is to save on storage on lab server
#Requires Samtools and pigz be installed and in path

# Example usage: 
# Rscript BigCompress.R --indir testcomp/ --task compress --extension all --minsize=10 --cores=8 --report testreport.txt 
# above finds all fastq, fasta, fna, fa, txt, and sam files above 10mb and compresses them using 8 cores.

##########################
#Get arguments
suppressMessages(library(optparse))
option_list = list(
  make_option(c("-i", "--indir"), type="character", help="input directory to be scanned for uncompressed files", metavar="character"),
  make_option(c("-t", "--task"), type="character", default="list", help="task to perform, options are list or compress [default= %default]", metavar="character"),
  make_option(c("-e", "--extension"), type="character", default="fastq", help="file extension to work on, options are fastq,fasta,fna,fa,txt,sam, or all [default= %default]", metavar="character"),
  make_option(c("-s", "--minsize"), type="numeric", default="10", help="number of Mb file must be to be compressed [default= %default]", metavar="numeric"),
  make_option(c("-c", "--cores"), type="numeric", default="6", help="number of cpus to use /threads to launch [default= %default]", metavar="numeric"),
  make_option(c("-o", "--report"), type="character", default="comp_report.csv", help="text file to output a summary of activity to [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$indir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input directory --indir).n", call.=FALSE)
}

opt$minsize<-opt$minsize*1e6 #convert mb to bytes

if(!dir.exists(opt$indir)){stop("Input directory ", opt$indir," not found.")}

if(!opt$task %in% c("list","compress","dereplicate")){stop("Task ", opt$task," not known.")}

if(opt$task %in% c("list","compress")){
  
  if(opt$extension=="all"){
    files<-list.files(opt$indir, pattern="fastq$|fasta$|fna$|fa$|txt$|sam$", recursive=TRUE, full.names=TRUE)
  } else {
    files<-list.files(opt$indir, pattern=paste0(opt$extension,"$"), recursive=TRUE, full.names=TRUE)
  }
  info<-file.info(files)
  info<-cbind(data.frame(File=rownames(info), Size_Gb=sapply(info$size, function(x) utils:::format.object_size(x, "Gb")), Size_Bytes=info$size), info[,-1])
  rownames(info)<-NULL
  info$ToCompress<-sapply(info$Size_Bytes, function(x) if(x>=opt$minsize){TRUE}else{FALSE})
  write.csv(info, opt$report, row.names=F)
  print(info)
  message("----------------------------------------------------------------------")
  message("See ", opt$report, " for report on ", nrow(info)," uncompressed files, of which ", sum(info$ToCompress), " have met the --minsize criteria to compress.")
  message("A total of ",utils:::format.object_size(sum(info$Size_Bytes), units="Gb"), " is currently being used by ", opt$extension, " target files")
}


info<-subset(info, ToCompress)

if(opt$task=="compress" & nrow(info)>0){
  message("Compressing...")
  for(i in 1:nrow(info)){
    
    if(gsub("..+\\.","", info$File[i]) %in% c("fastq","fasta","fna", "faa", "txt")){
      message(date(), " ---> ", i,"/",nrow(info)," Compressing ", info$File[i])
      if(file.exists(paste0(info$File[i],".gz"))){  
        message("-----> A pre-existing compressed version named ", paste0(info$File[i],".gz") ," has already been found. Please manually investigate.")
      } else {
        system(paste0("pigz --processes ",opt$cores," ", info$File[i]), intern=T)  
        }
    }
    
    if(gsub("..+\\.","", info$File[i])=="sam"){
      message(date(), " ---> ", i,"/",nrow(info)," Compressing ", info$File[i])
      if(file.exists(gsub("sam$","bam", info$File[i]))){ 
        message("-----> A pre-existing compressed version named ", gsub("sam$","bam", info$File[i]) ," has already been found. Please manually investigate.")
      } else {
        system(paste0("samtools view -bS --threads ",opt$cores," ", info$File[i], " > ", gsub("sam$", "bam",  info$File[i])), intern=T)
        unlink(info$File[i])
      }
    }
  }
}
