library(here)
library(dplyr)
library(R.utils)
library(stringr)
library(tidyr)
library(openxlsx)
library(patchwork)
library(ShortRead)
library(dplyr)
library(ggplot2)

run_name<-"test2"
primer_group<-"v1.0.2"
args<-commandArgs(trailingOnly = TRUE)
run_name<-args[1]
version<-args[2]

#read in length_standard
length_standard<-read.csv(here("..","input","length_standards","length_distribution_standard_all.csv"))
length_standard<-length_standard[order(length_standard$id),] #sort by id
names(length_standard)<-gsub("_mastersheet.csv","",names(length_standard))#remove .csv from names
length_standard<-length_standard[,c("range",paste0(primer_group))]

dir.create(here("..","output",run_name,primer_group,"length_distribution"), recursive = TRUE, showWarnings = FALSE)

samples<-list.files(here("..","output",run_name,primer_group,"basecalling","b.basecalled_demux"))
samples<-gsub(".fastq","",samples)


breaks <- seq(0, 1000, by = 10)
labels <- paste(
  breaks[-length(breaks)],
  breaks[-1],
  sep = "-"
)
length(labels)
full<-data.frame(id=c(1:length(labels)),range=labels)
bc<-samples[2]

for(bc in samples){
#  full reads from (b)
  # Read the FASTQ
  full_reads <- readFastq(here("..","output", run_name,primer_group, "basecalling", "b.basecalled_demux", paste0(gsub("-untrimmedall","", bc), ".fastq")))
  full_reads <- data.frame(
    name = ShortRead::id(full_reads),
    b = sread(full_reads)
  )
  full_reads$name <- gsub("\t.*", "", full_reads$name)
  full_reads$length <- nchar(full_reads$b)
  
  # Create bins
  full_reads$range <- cut(
    full_reads$length,
    breaks = breaks,
    labels = labels,
    include.lowest = TRUE,
    right = TRUE,
    ordered_result = TRUE
  )
  head(full_reads$range)
  summary_temp<-as.data.frame(table(full_reads$range))
  colnames(summary_temp)<-c("range",bc)
  full<-merge(full,summary_temp,by="range",all=TRUE)
  
  write.csv(full_reads,here("..","output",run_name,primer_group,"length_distribution",paste0(bc,"_read_length_distrubition.csv")),row.names=F)
  # Plot
  p <- ggplot(full_reads, aes(x = range)) +
    geom_bar() +
    scale_x_discrete(drop = FALSE) +  # keep all factor levels
    labs(
      x = "Length range",
      y = "Frequency",
      title = paste0(bc)
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6) # rotate + smaller font for many bins)
    ) + scale_y_log10(limits=c(1,15000))
  # Save
  ggsave(
    here("..","output", run_name,primer_group, "length_distribution", paste0(bc, "_read_length_distribution.png")),
    p,
    width = 12, height = 5
  )
  #write.csv(reads,here("..","output",run_name,paste0(bc,"_reads.csv")))

}

#append standard
full<-merge(length_standard,full,by="range",all=TRUE)
full<-full[order(full$id),]
full[is.na(full)]<-0
write.csv(full,here("..","output", run_name,primer_group, "length_distribution","read_length_distribution_summary.csv"),row.names=F)

