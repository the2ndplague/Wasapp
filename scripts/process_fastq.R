library(here)
library(readxl)
library(stringr)
library(fs)


run_name<-"test2"
#read in run name and sample list
args<-commandArgs(trailingOnly = TRUE)
run_name<-args[1]

sample_list<-read_excel(here("..","output",run_name,"SAMPLE_LIST.xlsx"))


identifiers<-sample_list$identifier
identifiers<-unique(identifiers[!is.na(identifiers)])



#rename a.basecalled
#read in samples
path<-here("..","output",run_name,
           "basecalling","a.basecalled")
a_basecalled_paths<- list.files(
          path = path,
          pattern = "\\.fastq$",
          recursive = TRUE,
          full.names = TRUE
          )
a_basecalled_paths<-a_basecalled_paths[grep("fastq_pass",a_basecalled_paths)]
a_basecalled_samples <- basename(a_basecalled_paths)

a_basecalled_new_names<-sub(".*_(barcode[0-9]+|unclassified).*\\.fastq$", "\\1.fastq", a_basecalled_samples)

for(identifier in identifiers){
  sample_list_temp<-sample_list[sample_list$identifier==identifier,]
  sample_list_temp<-sample_list_temp[!is.na(sample_list_temp$barcode),]
  
  paths<-a_basecalled_paths[grep(identifier,a_basecalled_paths)]
  

  
  samples<-basename(paths) 
  new_names<-sub(".*_(barcode[0-9]+|unclassified).*\\.fastq$", "\\1.fastq", samples)
  directory<-dirname(paths[identifier])

  for(n in c(1:nrow(sample_list_temp))){
    concentration<-sample_list_temp$concentration[n]
    version<-sample_list_temp$version[n]
    barcode<-sample_list_temp$barcode[n]
    id<-grep(barcode,samples)
    directory<-dirname(paths[id])

  new_name<-paste0(sample_list_temp[sample_list_temp$barcode==barcode,"id"])
  directory<-gsub("/basecalling",paste0("/",version,"/basecalling"),directory)
  
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file_copy(paths[id],
            paste0(directory,"/",new_name,".fastq"))
  }
  

}


#rename a.basecalled_demux
#read in samples
path<-here("..","output",run_name,
           "basecalling","b.basecalled_demux")
b_basecalled_paths<- list.files(
  path = path,
  pattern = "\\.fastq$",
  recursive = TRUE,
  full.names = TRUE
)

b_basecalled_samples <- basename(b_basecalled_paths)

b_basecalled_new_names<-sub(".*_(barcode[0-9]+|unclassified).*\\.fastq$", "\\1.fastq", b_basecalled_samples)
identifier<-identifiers[1]
for(identifier in identifiers){
  sample_list_temp<-sample_list[sample_list$identifier==identifier,]
  sample_list_temp<-sample_list_temp[!is.na(sample_list_temp$barcode),]
  paths<-b_basecalled_paths[grep(identifier,b_basecalled_paths)]
  samples<-basename(paths) 
  new_names<-sub(".*_(barcode[0-9]+|unclassified).*\\.fastq$", "\\1.fastq", samples)

  for(n in c(1:nrow(sample_list_temp))){
    concentration<-sample_list_temp$concentration[n]
    version<-sample_list_temp$version[n]
    barcode<-sample_list_temp$barcode[n]
    id<-grep(barcode,samples)
    directory<-dirname(paths[id])

    new_name<-paste0(sample_list_temp[sample_list_temp$barcode==barcode,"id"])
    directory<-gsub("/basecalling",paste0("/",version,"/basecalling"),directory)

    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    file_copy(paths[id],
              paste0(directory,"/",new_name,".fastq"))
  }
  
  
}

samples<-unique(sample_list$version)
writeLines(samples,sep="\n",here("..","output",run_name,"versions.txt"
                                 ))







