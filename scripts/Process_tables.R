library(here)
library(dplyr)
library(R.utils)
library(stringr)
library(tidyr)
library(openxlsx)
run_name<-"exp16_WGS"
primer_version<-"COVID_WGS"
args<-commandArgs(trailingOnly = TRUE)
run_name<-args[1]
primer_version<-args[2]

sample_summary<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","table_2.csv"))
per_amplicon<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","table_3.csv"))
per_sample<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","table_4.csv"))
#clean data

sample_summary$Reads<-gsub(" k","k",sample_summary$Reads)
sample_summary$Bases<-gsub(" k","k",sample_summary$Bases)
sample_summary$Reads<-gsub(" M","M",sample_summary$Reads)
sample_summary$Bases<-gsub(" M","M",sample_summary$Bases)

per_amplicon$Reads<-gsub(" k","k",per_amplicon$Reads)
per_amplicon$Bases<-gsub(" k","k",per_amplicon$Bases)
per_amplicon$Reads<-gsub(" M","M",per_amplicon$Reads)
per_amplicon$Bases<-gsub(" M","M",per_amplicon$Bases)

per_sample$Reads<-gsub(" k","k",per_sample$Reads)
per_sample$Bases<-gsub(" k","k",per_sample$Bases)
per_sample$Reads<-gsub(" M","M",per_sample$Reads)
per_sample$Bases<-gsub(" M","M",per_sample$Bases)

per_amplicon$Reads<-sub("^\\S+\\s+(\\S+).*", "\\1", per_amplicon$Reads)
per_amplicon$Bases<-sub("^\\S+\\s+(\\S+).*", "\\1", per_amplicon$Bases)
per_sample$Reads<-sub("^\\S+\\s+(\\S+).*", "\\1", per_sample$Reads)
per_sample$Bases<-sub("^\\S+\\s+(\\S+).*", "\\1", per_sample$Bases)
sample_summary$Reads<-sub("^\\S+\\s+(\\S+).*", "\\1", sample_summary$Reads)
sample_summary$Bases<-sub("^\\S+\\s+(\\S+).*", "\\1", sample_summary$Bases)
#Clean up "k" to 1000
per_amplicon$Reads<-ifelse(grepl("k",per_amplicon$Reads,ignore.case=TRUE),
                           as.numeric(sub("k","",per_amplicon$Reads,ignore.case=TRUE))*1000,
                           per_amplicon$Reads)
per_amplicon$Bases<-ifelse(grepl("k",per_amplicon$Bases,ignore.case=TRUE),
                           as.numeric(sub("k","",per_amplicon$Bases,ignore.case=TRUE))*1000,
                           per_amplicon$Bases)
per_amplicon$Reads<-ifelse(grepl("M",per_amplicon$Reads),
                           as.numeric(sub("M","",per_amplicon$Reads))*1000000,
                           as.numeric(per_amplicon$Reads))
per_amplicon$Bases<-ifelse(grepl("M",per_amplicon$Bases),
                           as.numeric(sub("M","",per_amplicon$Bases))*1000000,
                           as.numeric(per_amplicon$Bases))
per_sample$Reads<-ifelse(grepl("k",per_sample$Reads),
                         as.numeric(sub("k","",per_sample$Reads))*1000,
                         per_sample$Reads)
per_sample$Bases<-ifelse(grepl("k",per_sample$Bases),
                         as.numeric(sub("k","",per_sample$Bases))*1000,
                         per_sample$Bases)
per_sample$Reads<-ifelse(grepl("M",per_sample$Reads),
                         as.numeric(sub("M","",per_sample$Reads))*1000000,
                         as.numeric(per_sample$Reads))
per_sample$Bases<-ifelse(grepl("M",per_sample$Bases),
                         as.numeric(sub("M","",per_sample$Bases))*1000000,
                         as.numeric(per_sample$Bases))
sample_summary$Reads<-ifelse(grepl("k",sample_summary$Reads,ignore.case=TRUE),
                             as.numeric(sub("k","",sample_summary$Reads,ignore.case=TRUE))*1000,
                             sample_summary$Reads)
sample_summary$Bases<-ifelse(grepl("k",sample_summary$Bases,ignore.case=TRUE),
                             as.numeric(sub("k","",sample_summary$Bases,ignore.case=TRUE))*1000,
                             sample_summary$Bases)
sample_summary$Reads<-ifelse(grepl("M",sample_summary$Reads),
                             as.numeric(sub("M","",sample_summary$Reads))*1000000,
                             as.numeric(sample_summary$Reads))
sample_summary$Bases<-ifelse(grepl("M",sample_summary$Bases),
                             as.numeric(sub("M","",sample_summary$Bases))*1000000,
                             as.numeric(sample_summary$Bases))

#split per_sample into different sheets 

df_list<-split(per_sample,per_sample$Sample)
wb<-createWorkbook()


#rename worksheets if necessary
names(df_list)<-gsub("barcode","",names(df_list))
names(df_list)<-gsub("PlasmidControl","Plasmid",names(df_list))
names(df_list)
for (name in names(df_list)) {
  addWorksheet(wb, name)
  writeData(wb, name, df_list[[name]])
}
saveWorkbook(wb, here("..","output",run_name,primer_version,"summary_stats.xlsx"), overwrite = TRUE)
write.csv(sample_summary,here("..","output",run_name,primer_version,"per_amplicon.csv"),row.names=F)
#write.csv(per_sample,here("..","output",run_name,primer_version,"per_sample.csv"),row.names=F)
write.csv(per_amplicon,here("..","output",run_name,primer_version,"per_amplicon.csv"),row.names=F)


##make summary table


#clean sample names
#sample_summary$Sample.alias<-gsub(".*_","",sample_summary$Sample.alias)
#per_sample$Sample<-gsub(".*_","",per_sample$Sample)
samples<-unique(sample_summary$Sample.alias)

sample_summary<-data.frame(t(sample_summary))
names(sample_summary)<-sample_summary[1,]
sample_summary<-sample_summary[-1,]

#make reads file
full<-data.frame(Amplicon=unique(per_sample$Amplicon))
for(sample in samples){
  sample_reads<-per_sample[per_sample$Sample==sample,c(2,3)]
  names(sample_reads)<-c("Amplicon",sample)
  full<-merge(full,sample_reads,by="Amplicon",all=T)
}
row.names(full)<-full$Amplicon
full<-full[,-1]

full<-rbind(sample_summary,full)
write.csv(full,here("..","output",run_name,primer_version,"summarised_data.csv"))






#get report folder
trace<-read.csv(here("..","output",run_name,primer_version,"wf-amplicon","execution","trace.txt"))
trace<-trace[nrow(trace),1]
code<-strsplit(trace,"\t")[[1]][2]
folder<-gsub("/.*","",code)
subfolder<-gsub(".*/","",code)
#################################################################

report_folder<-list.dirs(here("work"))
report_folder<-report_folder[grep(paste0(here("work",folder),"/",subfolder),report_folder)]
report_folder<-report_folder[grep("data",report_folder)]
#report_folder<-report_folder[-1]
report_folder
n<-report_folder[1]
depth_full<-data.frame()
for(n in report_folder){
  barcode<-sub(".*data/","",n)
  #get depth info
  depth_file_gz <- paste0(n, "/per-window-depth.tsv.gz")
  bamstats_file <- paste0(n, "/bamstats.tsv")
  barcode <- sub(".*data/", "", n)
  
  if (file.exists(depth_file_gz) && file.exists(bamstats_file)) {
  gunzip(paste0(n,"/per-window-depth.tsv.gz"),here("..","output",run_name,primer_version,"scrapped_tables","window_depth",paste0(barcode,"_per_window_depth.tsv")),remove=FALSE)
  depth<-read.table(here("..","output",run_name,primer_version,"scrapped_tables","window_depth",paste0(barcode,"_per_window_depth.tsv")),sep="\t")
  colnames(depth)<-depth[1,]
  depth<-depth[-1,]
  depth$barcode<-barcode
  depth_full<-rbind(depth_full,depth)
  write.csv(depth,here("..","output",run_name,primer_version,"scrapped_tables","depth",paste0(barcode,".csv")),row.names=F)
  
  
  #get bamstats info to working folder
  
  bamstats<-read.table(paste0(n,"/bamstats.tsv"),sep="\t")
  names(bamstats)<-bamstats[1,]
  bamstats<-bamstats[-1,]
  write.csv(bamstats,here("..","output",run_name,primer_version,"scrapped_tables","bamstats",paste0(barcode,"_bamstats.csv")),row.names=F)
  amplicon_list<-unique(per_sample$Amplicon)
  write.csv(amplicon_list,here("..","output",run_name,primer_version,"amplicon_list.csv"),row.names=F)
  
  
  # Check if both required files exist

    #get depth info
    depth <- read.table(here("..","output", run_name,primer_version, "scrapped_tables", "window_depth", paste0(barcode, "_per_window_depth.tsv")), sep = "\t")
    colnames(depth) <- depth[1, ]
    depth <- depth[-1, ]
    depth$barcode <- barcode
    depth_full <- rbind(depth_full, depth)
    write.csv(depth, here("..","output", run_name,primer_version, "scrapped_tables", "depth", paste0(barcode, ".csv")), row.names = F)

    #get bamstats info to working folder
    bamstats <- read.table(bamstats_file, sep = "\t")
    names(bamstats) <- bamstats[1, ]
    bamstats <- bamstats[-1, ]
    write.csv(bamstats, here("..","output", run_name,primer_version, "scrapped_tables", "bamstats", paste0(barcode, "_bamstats.csv")), row.names = F)
    amplicon_list <- unique(per_sample$Amplicon)
    write.csv(amplicon_list, here("..","output", run_name,primer_version, "amplicon_list.csv"), row.names = F)
  } else {
    # If files are missing, print a warning and skip this barcode
    warning(paste("Skipping barcode", barcode, "due to missing files in", n))
  }
}

write.csv(depth_full,here("..","output",run_name,primer_version,"scrapped_tables","depth",paste0("all_depths.csv")),row.names=F)

#get depth summary (note this is very unclean, to redo when free)
depth_full<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","depth",paste0("all_depths.csv")))
refs<-unique(depth_full$ref)
for(ref in refs){

  depth_matrix <- depth_full %>%
    dplyr::summarise(mean_value = mean(depth), .by = c(barcode, ref)) %>%
    tidyr::pivot_wider(names_from = barcode, values_from = mean_value)
  
  depth_cv_matrix <- depth_full %>%
    dplyr::summarise(cv_value = sd(depth)/mean(depth), .by = c(barcode, ref)) %>%
    tidyr::pivot_wider(names_from = barcode, values_from = cv_value)
  
  write.csv(depth_matrix,here("..","output",run_name,primer_version,"scrapped_tables","depth","depth_summary.csv"),row.names=F)
  write.csv(depth_cv_matrix,here("..","output",run_name,primer_version,"scrapped_tables","depth","depth_cv_summary.csv"),row.names=F)
}


###########################################################################################################
#create read summary
barcodes<-list.files(here("..","output",run_name,primer_version,"scrapped_tables","bamstats"))
barcodes<-gsub("_bamstats.csv","",barcodes)
amplicon_list<-read.csv(here("..","output",run_name,primer_version,"amplicon_list.csv"),col.names=1)
amplicon_list<-data.frame(amplicon=amplicon_list)
names(amplicon_list)<-"amplicon"
#make empty df for reads and coverage
read_full<-data.frame(amplicon=amplicon_list)
coverage_full<-data.frame(amplicon=amplicon_list)
meanqual_full<-data.frame(amplicon=amplicon_list)

summary_stats<-data.frame(amplicon=c("Total_reads",
                                    "Total_bases",
                                    "no. amplicons",
                                    "Mean length",
                                    "Mean coverage",
                                    "Mean quality"
))

for(barcode in barcodes){
  bamstats<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","bamstats",paste0(barcode,"_bamstats.csv")))
  #remove NaN lines
  bamstats<-bamstats%>% drop_na(coverage)
  if(nrow(bamstats)<1){
    next
  }
  #make summary portion
  total_bases<-round(sum(as.numeric(bamstats$length)),digits=3)
  total_reads<-round(length(bamstats$sample_name),digits=3)
  num_amplicons<-round(length(unique(bamstats$ref)),digits=3)
  mean_length<-round(mean(bamstats$length),digits=3)
  mean_cov<-round(mean(as.numeric(bamstats$ref_coverage)),digits=3)
  mean_qual<-round(mean(as.numeric(bamstats$mean_quality)),digits=3)
  stats_temp<-data.frame(amplicon=c("Total_reads",
                                      "Total_bases",
                                      "no. amplicons",
                                      "Mean length",
                                      "Mean coverage",
                                      "Mean quality"
  ),
  stats=c(total_reads,
          total_bases,
          num_amplicons,
          mean_length,
          mean_cov,
          mean_qual)
  )
  names(stats_temp)<-c("amplicon",barcode)
  summary_stats<-merge(summary_stats,stats_temp,by="amplicon")
  #make reads portion
  read_temp2<-data.frame()
  coverage_temp2<-data.frame()
  meanqual_temp2<-data.frame()
  amplicons<-unique(bamstats$ref)

  for(amplicon in amplicons){
    #counts
    read_count<-nrow(bamstats[bamstats$ref==amplicon,])
    read_temp<-data.frame(amplicon=amplicon,count=read_count)
    read_temp2<-rbind(read_temp2,read_temp)
    
    #coverage
    coverage_temp<-bamstats[bamstats$ref==amplicon,c("ref","length","ref_coverage")]
    mean_cov<-mean(as.numeric(coverage_temp$ref_coverage))
    coverage_temp<-data.frame(amplicon=amplicon,mean_coverage=mean_cov)
    coverage_temp2<-rbind(coverage_temp2,coverage_temp)
    
    #meanquality
    meanqual_temp<-bamstats[bamstats$ref==amplicon,c("ref","length","mean_quality")]
    mean_qual<-mean(as.numeric(meanqual_temp$mean_quality))
    meanqual_temp<-data.frame(amplicon=amplicon,mean_qual=mean_qual)
    meanqual_temp2<-rbind(meanqual_temp2,meanqual_temp)
  }
  
  names(read_temp2)<-c("amplicon",barcode)
  names(coverage_temp2)<-c("amplicon",barcode)
  names(meanqual_temp2)<-c("amplicon",barcode)
  
  read_full<-merge(read_full,read_temp2,by="amplicon",all=T)
  coverage_full<-merge(coverage_full,coverage_temp2,by="amplicon",all=T)
  meanqual_full<-merge(meanqual_full,meanqual_temp2,by="amplicon",all=T)
}

read_full<-rbind(summary_stats,read_full)
coverage_full<-rbind(summary_stats,coverage_full)
meanqual_full<-rbind(summary_stats,meanqual_full)

write.csv(read_full,here("..","output",run_name,primer_version,paste0("reads_summary.csv")),row.names=F)
write.csv(coverage_full,here("..","output",run_name,primer_version,paste0("coverage_summary.csv")),row.names=F)
write.csv(meanqual_full,here("..","output",run_name,primer_version,paste0("meanquality_summary.csv")),row.names=F)


#get reads 
library(ShortRead)
library(dplyr)

a_basecalled<-list.files(here("..","output",run_name,primer_version,"basecalling","a.basecalled"),pattern="\\.fastq$",recursive=T,full.names=T)

for(bc in samples){
  reads<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","bamstats",paste0(bc,"_bamstats.csv")))
  if(grepl("untrimmed",bc)==FALSE){
#trimmed reads from (c)
    here("..","output",run_name,primer_version,"basecalling","c.basecalled_demux_trim",paste0(bc),paste0(bc,".fastq"))
  fastq<-readFastq(here("..","output",run_name,primer_version,"basecalling","c.basecalled_demux_trim",paste0(bc),paste0(bc,".fastq")))
  }else{fastq<-readFastq(here("..","output",run_name,primer_version,"basecalling","c.basecalled_demux_trim",paste0(bc,""),paste0(bc,".fastq")))
}
  fastq<-data.frame(name=ShortRead::id(fastq),c=sread(fastq))
  fastq$name<-gsub("\t.*","",fastq$name)
  reads<-merge(reads,fastq,by="name")

#  full reads from (b)
  full_reads<-readFastq(here("..","output",run_name,primer_version,"basecalling","b.basecalled_demux",paste0(gsub("-untrimmedall","",bc),".fastq")))
  full_reads<-data.frame(name=ShortRead::id(full_reads),b=sread(full_reads))
  full_reads$name<-gsub("\t.*","",full_reads$name)
  reads<-merge(reads,full_reads,by="name")

# full reads + adapters + barcodes from (a)
  full_reads<-a_basecalled[grepl("fastq_pass",a_basecalled)&grepl(bc,a_basecalled)]
  full_reads<-readFastq(here(full_reads))
  full_reads<-data.frame(name=ShortRead::id(full_reads),a=sread(full_reads))
  full_reads$name<-gsub("\t.*","",full_reads$name)
  reads<-merge(reads,full_reads,by="name")
  write.csv(reads,here("..","output",run_name,primer_version,paste0(bc,"_reads.csv")))
}





# process variants for all barcodes
library(Rsamtools)
library(GenomicAlignments)
#create variants summary:
table_5<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","table_5.csv"))
table_6<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","table_6.csv"))
variants_summary<-rbind(table_5,table_6)

mean_q<-c()
strand_counts<-data.frame()
n_reads<-c()
#iterate through each variant, read the bam file and perform the following checks:
for(variant in c(1:nrow(variants_summary))){
  sample<-variants_summary$Sample[variant]
  amplicon<-variants_summary$Amplicon[variant]
  position<-variants_summary$Position[variant]
  ref_base<-variants_summary$Ref..allele[variant]
  which <- GRanges(amplicon, IRanges(start = position, end = position))
  param <- ScanBamParam(
    which = which,
    what = c("seq", "qual", "pos", "cigar"),
    tag=c("NM")
  )
  #checks to prevent index error
  bam_file <- here("..","output",run_name,primer_version,"wf-amplicon",sample,"alignments",
                   paste0(sample,".aligned.sorted.bam"))
  bai_file <- paste0(bam_file, ".bai")
  if (!file.exists(bai_file) || file.info(bai_file)$mtime < file.info(bam_file)$mtime) {
    message("Re-indexing: ", sample)
    system(paste("samtools index", bam_file))
  }
  aln <- readGAlignments(here("..","output",run_name,primer_version,"wf-amplicon",sample,"alignments",
                        paste0(sample,".aligned.sorted.bam")),
                         param = param)
 
  
  n_reads<-c(n_reads,length(aln))
  strand_counts <- rbind(strand_counts,table(strand(aln)))
  
  #mean quality at position of interest
  quals_list <- lapply(seq_along(aln), function(i) {
    aln_i <- aln[i]
    # aligned positions (reference)
    ref_pos <- start(aln_i):(start(aln_i) + width(aln_i) - 1)
    # qualities (Phred)
    q <- as.integer(mcols(aln_i)$qual[[1]])
    # match region
    keep <- ref_pos >= position-10 & ref_pos <= position+10
    q[keep]
  })
  # flatten
  all_quals <- unlist(quals_list)
  mean_q<-c(mean_q,mean(all_quals,na.rm=TRUE))
  
}

names(strand_counts)<-c("+","-","*")
variants_summary<-cbind(variants_summary,mean_q,n_reads,strand_counts)

