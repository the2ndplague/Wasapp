library(here)
#get sequali table summary
########################################################################
# Install and load the jsonlite package if you haven't already
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
library(jsonlite)

# Define the path to your JSON file
# Note: R works best with forward slashes in file paths, even on Windows.
#run_name<-readLines(here("..","output","run_name.txt"))
run_name<-"exp9_new"
args<-commandArgs(trailingOnly = TRUE)
run_name<-args[1]

sequali_reports<-list.files(here("..","output",run_name,"sequali_report"))
sequali_reports<-sequali_reports[grep(".json",sequali_reports)]

report_dfs <- list()

for(report in sequali_reports){
  file_path <- here("..","output",run_name,"sequali_report",report)
  # Read the JSON file and parse it into an R list object
  summary_list <- jsonlite::fromJSON(file_path)$summary
  
  # unlist() flattens any nested lists into a simple named vector,
  # which prevents the "unimplemented type 'list'" error.
  # data.frame() then converts this into a single-column data frame.
  summary_df <- data.frame(unlist(summary_list))
  report_name <- gsub("\\.fastq\\.json", "", report)
  names(summary_df) <- report_name
  
  report_dfs[[report_name]] <- summary_df
}

all_reports <- do.call(cbind, report_dfs)
all_reports<-t(all_reports)




checklist<-data.frame()
reports<-list.files(here("..","output",run_name,"sequali_report"),pattern="\\.json$")
for(report in reports){
  id<-gsub(".fastq.json","",report)
  json<-fromJSON(here("..","output",run_name,"sequali_report",report))
  adapter_section<-json$adapter_content
  # Extract adapter names
  adapter_names <- sapply(adapter_section$adapter_content, `[`, 1)
  # Extract the numeric data
  adapter_values <- lapply(adapter_section$adapter_content, `[`, 2)
  adapter_df<-data.frame(labels=adapter_section$x_labels)
  for(n in c(1:length(adapter_names))){
    name<-adapter_names[n]
    values<-adapter_values[[n]][[1]]
    adapter_temp<-data.frame(values)
    names(adapter_temp)<-name
    adapter_df<-cbind(adapter_df,adapter_temp)
  }
  check<-adapter_df[,-1]
  check<-c(id,
           any(check>1),
           round(as.numeric(max(check)),2),
           sum(check>1))
  checklist<-rbind(checklist,check)
  
}

names(checklist)<-c("id","adapter_content_present?","highest_AC%","bins_with_AC%>1")
row.names(checklist)<-checklist$id
checklist<-checklist[,-1]

all_reports<-merge(all_reports,checklist,by="row.names")

write.csv(all_reports,here("..","output",run_name,"sequali_summary.csv"),row.names=T)
