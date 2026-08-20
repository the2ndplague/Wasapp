library(here)
library(readxl)
library(openxlsx)

run_name="exp15"
args<-commandArgs(trailingOnly = TRUE)
run_name<-args[1]




amplicons<-c()

wb<-createWorkbook()

#sequencing overview
######################################################################################
sample_list<-read_excel(here("SAMPLE_LIST.xlsx"))
dirs<-unique(sample_list$version)

addWorksheet(wb, sheetName = "Sequencing overview")
writeData(wb, sheet = "Sequencing overview", x = sample_list)
######################################################################################




#sequali QC stats
######################################################################################
sequali_all<-data.frame()

  sequali<-read.csv(here("sequali_summary.csv"))
  sequali<-sequali[sequali$Row.names!="unclassified",]
  sequali_all<-rbind(sequali_all,sequali)


sequali_all<-sequali_all[,c(2,6:9,12:14),]
sequali_all$'Q20_reads%'<-round((sequali_all$q20_reads/sequali_all$total_reads)*100,2)
sequali_all$'Q20_bases%'<-round((sequali_all$q20_bases/sequali_all$total_bases)*100,2)

sequali_all<-sequali_all[,c(1,2,9,10,6)]
names(sequali_all)<-c("Sample","Total_reads","Q20_reads%","Q20_bases%","Adapter_content_present?")
sequali_all<-sequali_all[order(sequali_all$Sample),]
sequali_all$'Adapter_content_present?'<-gsub("FALSE","No",sequali_all$'Adapter_content_present?')
sequali_all$'Adapter_content_present?'<-gsub("TRUE","Yes",sequali_all$'Adapter_content_present?')

Notes<-c("Sequali reports",
         "Q20 reads (Ideally above >50%)",
         "Size distribution of sequence length should be within the range of your amplicons",
         "No adapter content should be present",rep("",nrow(sequali_all)-4))

sequali_all$Notes<-Notes

addWorksheet(wb, sheetName = "Sequali QC stats")
writeData(wb, sheet = "Sequali QC stats", x = sequali_all)

# 1. Create the style you want to apply (Red background, white text)
format <- createStyle(bgFill = "#FFC7CE", fontColour = "#9C0006")
# 2. Apply the conditional formatting rule
conditionalFormatting(
  wb, 
  sheet = "Sequali QC stats", 
  cols = 3,           
  rows = 2:(nrow(sequali_all)+1), 
  rule = "<50",       
  style = format
)

######################################################################################





#reads summary page
######################################################################################
stats_all<-data.frame(amplicon=NA)
summary_all<-data.frame(amplicon=NA)
for(dir in dirs){
  summary<-read.csv(here(dir,"reads_summary.csv"))[-c(1:6),]
  stats<-read.csv(here(dir,"reads_summary.csv"))[c(1:6),]
  amplicons<-unique(c(amplicons,summary$amplicon))  
  
  stats_all<-merge(stats_all,stats,by="amplicon",all=T)
  summary_all<-merge(summary_all,summary,by="amplicon",all=T)
}

summary_all<-summary_all[names(summary_all) != "unclassified"]
stats_all<-stats_all[names(stats_all) != "unclassified"]

samples<-names(summary_all)[-c(grep("unclassified",names(summary_all)),1)]
represent_stats<-data.frame(amplicon=c("Average","Stdev","CV",">1.92","<1.92"))
represent_all<-summary_all

for(sample in samples){
  totalreads<-stats_all[stats_all$amplicon=="Total_reads",sample]
  totalreads<-totalreads[!is.na(totalreads)]
  represent_all[,sample]<-round(as.numeric(represent_all[,sample])/totalreads*100,2)

  mean<-round(mean(represent_all[,sample],na.rm=TRUE),2)
  sd<-round(sd(represent_all[,sample],na.rm=TRUE),2)
  cv<-round(sd/mean,2)
  more<-sum(represent_all[,sample]>=1.92,na.rm=TRUE)
  less<-sum(represent_all[,sample]<1.92,na.rm=TRUE)
  
  represent_stats<-cbind(represent_stats,c(mean,sd,cv,round(more,0),round(less,0)))
  names(represent_stats)[ncol(represent_stats)]<-sample
}


summary_all<-rbind(stats_all,summary_all)

represent_all<-rbind(represent_all,represent_stats)

addWorksheet(wb, sheetName = "Results")
writeData(wb, sheet = "Results", x = summary_all)


format <- createStyle(bgFill = "#FFC7CE")
conditionalFormatting(
  wb, 
  sheet = "Results", 
  cols = 2:ncol(summary_all),           
  rows = 5, 
  rule = "<52",       
  style = format
)

format <- createStyle( fontColour = "red")
conditionalFormatting(
  wb, 
  sheet = "Results", 
  cols = c(2:ncol(summary_all)),           
  rows = 9:(nrow(summary_all)+1), 
  rule = "<100",       
  style = format
)


addWorksheet(wb, sheetName = "Coverage")
writeData(wb, sheet = "Coverage", x = represent_all) #not actually coverage, so the dataframe is called represent_all

format <- createStyle( fontColour = "red")
conditionalFormatting(
  wb, 
  sheet = "Coverage", 
  cols = 2:ncol(represent_all),           
  rows = 2:(nrow(represent_all)+1), 
  rule = "<1.92",       
  style = format
)

######################################################################################


#make variants list
######################################################################################
variants_all<-data.frame()
for(dir in dirs){
  variant_file<-read.csv(here(dir,"scrapped_tables","table_5.csv"))
  #add unique identifier for variants
  variant_file$id<-paste0(variant_file$Amplicon,
                          variant_file$Position,
                          variant_file$Ref..allele,
                          variant_file$Alt..allele,
                          variant_file$Type)
  variants_all<-rbind(variants_all,variant_file)
  
}
variants_all<-variants_all[names(variants_all) != "unclassified"]
variants_all<-variants_all[order(variants_all$Sample),]
variants_all$Sample<-gsub("-",".",variants_all$Sample)


depth_percents<-c()
amplicon_reads<-c()
for(variant in c(1:nrow(variants_all))){
  amplicon<-variants_all[variant,"Amplicon"]
  sample<-variants_all[variant,"Sample"]  
  amplicon_read<-summary_all[which(summary_all$amplicon==amplicon),sample]
  depth_percent<-round((variants_all[variant,"Depth"]/amplicon_read)*100,2)
  depth_percents<-c(depth_percents,depth_percent)
  
  amplicon_reads<-c(amplicon_reads,amplicon_read)
}

variants_all$Amplicon_reads<-amplicon_reads
variants_all$Depth_percent<-depth_percents

addWorksheet(wb, sheetName = "Variants")
writeData(wb, sheet = "Variants", x = variants_all)

######################################################################################
saveWorkbook(wb, file = here("experimental_summary.xlsx"), overwrite = TRUE)














primer_versions<-read_excel(here("SAMPLE_LIST.xlsx"))
primer_versions<-unique(primer_versions$version)

#PART 2
##########################################################################################################

# Create the final consolidated output directory
final_output_dir <- here::here("final_report")
if (!dir.exists(final_output_dir)) {
  dir.create(final_output_dir, recursive = TRUE)
}

# Get all unique viruses across all primer versions first
all_viruses <- c()
for (primer_version in primer_versions) {
  mastersheet_path <- here::here("..","..", "input", "primer_files", paste0(primer_version, "_mastersheet.csv"))
  if (file.exists(mastersheet_path)) {
    mastersheet <- read.csv(mastersheet_path)
    all_viruses <- c(all_viruses, unique(mastersheet$Virus))
  }
}
viruses <- unique(all_viruses)

for(virus in viruses){
  # Initialize collectors for all samples for this virus
  all_div_classes <- ""
  all_constants <- ""
  all_samples_for_css <- c()
  
  header<-paste0('<!DOCTYPE html>
<html lang="en">
<head> 
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Interactive Genomic Plot Study - ', virus, '</title>
    
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    
    <style>
        body {
            font-family: system-ui, -apple-system, sans-serif;
            margin: 40px;
            background-color: #f8fafc;
        }
        .sample-card {
            background: white;
            padding: 20px;
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.05);
            max-width: 1200px;
            margin: 0 auto 40px auto; /* Added bottom margin to separate them */
        }
        h1 {
            font-size: 1.5rem;
            color: #1e293b;
            margin-bottom: 5px;
        }
        p {
            color: #64748b;
            margin-bottom: 20px;
            font-size: 0.95rem;
        }
    </style>
</head>
<body>')
  
  script_setup<-paste0("<script>
    // Global configurations shared by all plots
    const config = {
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: ['select2d', 'lasso2d', 'autoScale2d']
    };")
  
  # Loop through each primer version to collect samples for the current virus
  for (primer_version in primer_versions) {
    graphing_data_path <- here::here("..", "output", run_name, primer_version, "final_report", "graphing_data")
    if (!dir.exists(graphing_data_path)) next
    
    graphing_files <- list.files(graphing_data_path, pattern = "_info.csv")
    samples <- gsub("_info.csv", "", graphing_files)
    mastersheet <- read.csv(here::here("..", "input", "primer_files", paste0(primer_version, "_mastersheet.csv")))
    
    for(sample in samples){
      sample_file<-read.csv(here::here("..","output",run_name,primer_version,"final_report","graphing_data",paste0(sample,"_info.csv")))
      
      virus_df<-sample_file[,grep(paste0("^", virus, "_"),names(sample_file))]
      if (is.vector(virus_df)) {
        virus_df <- as.data.frame(virus_df)
      }
      virus_df<-virus_df[!is.na(virus_df[,2]) & virus_df[,2] != "",]
      
      if (nrow(virus_df) == 0) {
        next
      }
      
      js_sample_suffix <- gsub("[^a-zA-Z0-9_]", "_", sample)
      all_samples_for_css <- c(all_samples_for_css, paste0(".", js_sample_suffix))
      
      # Generate and collect the HTML div for this sample
      div_class <- paste0('<div class="sample-card ', js_sample_suffix, '">
    <h1>Interactive Sequence Coverage Landscape - ', sample, ' (', primer_version, ')</h1>
    <p>Hover over data points to read values. Click and drag to zoom. Double-click to reset view.</p>
    <div id="genomic-plot-', js_sample_suffix, '" style="width:100%; height:600px;"></div>
</div>
                  
')
      all_div_classes <- paste0(all_div_classes, div_class)
      
      # Generate and collect the JavaScript for this sample
      x_axis<-paste0("[",paste(as.numeric(row.names(virus_df)),collapse=", "),"]")
      depth_values<-paste0("[",paste(as.numeric(virus_df[,2]),collapse=", "),"]")
      nucleotides<-paste0("['",paste(virus_df[,1],collapse="', '"),"']")
      
      variant_x<-variant_text<-virus_df[virus_df[,3]>0,]
      if (nrow(variant_x) > 0) {
        variant_x_val<-row.names(variant_x)
        variant_x_str<-paste0("[",paste(as.numeric(variant_x_val),collapse=", "),"]")
        
        variant_y<-as.numeric(virus_df[,3])
        variant_y<-variant_y[variant_y>0]
        variant_y_str<-paste0("[",paste(variant_y,collapse=", "),"]")
        
        variant_text_val<-virus_df[virus_df[,3]>0,4]
        variant_text_str<-paste0("['",paste(variant_text_val,collapse="', '"),"']")
      } else {
        variant_x_str<-"[]"
        variant_y_str<-"[]"
        variant_text_str<-"[]"
      }
      
      #amplicon_region
      mastersheet_temp<-mastersheet[mastersheet$Virus==virus,]
      mastersheet_temp<-mastersheet_temp[order(mastersheet_temp$amp_start),]
      
      max_y<-max(virus_df[,2], na.rm = TRUE)
      if (is.na(max_y) || max_y <= 0) {
        max_y <- 1
      }
      mastersheet_temp$pos<-max_y* rep(c(0.95, 1), length.out = nrow(mastersheet_temp))
      
      max_x<-max(mastersheet_temp$amp_end, na.rm = TRUE)
      
      amplicon_regions<-""
      annotations<-""
      
      if (nrow(mastersheet_temp) > 0) {
        for(amplicon in c(1:nrow(mastersheet_temp))){
          amplicon_name<-mastersheet_temp$amplicon_name[amplicon]
          fwdP_sequence<-mastersheet_temp$Fwd_primer[amplicon]
          revP_sequence<-mastersheet_temp$Rev_primer[amplicon]
          amplicon_name<-gsub("Screen-","",amplicon_name)
          amplicon_start<-as.numeric(mastersheet_temp$amp_start[amplicon])
          amplicon_end<-as.numeric(mastersheet_temp$amp_end[amplicon])
          
          fwdP_start<-as.numeric(mastersheet_temp$amp_start[amplicon]-nchar(mastersheet_temp$Fwd_primer[amplicon]))
          revP_end<-as.numeric(mastersheet_temp$amp_end[amplicon]+nchar(mastersheet_temp$Rev_primer[amplicon]))
          
          y_pos<-mastersheet_temp$pos[amplicon]  
          
          amplicon_region<-paste0("{
                      type: 'rect',
                      xref: 'x',
                      yref: 'paper',
                      x0: ",amplicon_start,",            // Amplicon Start (bp)
                      x1: ",amplicon_end,",            // Amplicon End (bp)
                      y0: ",(y_pos/max_y),",           // Height start
                      y1: ",(y_pos/max_y)+0.01,",           // Height stop
                      fillcolor: 'rgba(234, 179, 8, 0.4)', // Translucent Gold
                      line: { width: 0 }
                  },
             {
                      type: 'rect',
                      xref: 'x',
                      yref: 'paper',
                      x0: ",fwdP_start,",            // fwd primer Start (bp)
                      x1: ",amplicon_start,",            // fwd primer End (bp)
                      y0: ",(y_pos/max_y),",           // Height start
                      y1: ",(y_pos/max_y)+0.01,",           // Height stop
                      fillcolor: 'rgba(16, 185, 129, 0.4)', // Translucent Emerald Green
                      line: { width: 0 }
                  },
             {
                      type: 'rect',
                      xref: 'x',
                      yref: 'paper',
                      x0: ",amplicon_end,",            // Amplicon Start (bp)
                      x1: ",revP_end,",            // Amplicon End (bp)
                      y0: ",(y_pos/max_y),",           // Height start
                      y1: ",(y_pos/max_y)+0.01,",           // Height stop
                      fillcolor: 'rgba(239, 68, 68, 0.4)',  // Translucent Red
                      line: { width: 0 }
                  }")
          
          if (amplicon_regions == "") {
            amplicon_regions <- amplicon_region
          } else {
            amplicon_regions <- paste0(amplicon_regions, " , ", amplicon_region)
          }
          
          annotation<-paste0("{
                      xref: 'x',
                      yref: 'paper',
                      x: ",amplicon_start+((amplicon_end-amplicon_start)/2),",             // amplicon
                      y: ",(y_pos/max_y)+0.0005,",           // 
                      text: '",amplicon_name,"',     // The actual label text!
                      showarrow: false,   // 
                      font: {
                          size: 10,
                          color: '#71717a', // Dark gray text so it's easy to read
                          weight: 'bold'
                      }
                  },
                         {
                      xref: 'x',
                      yref: 'paper',
                      x: ",fwdP_start+((amplicon_start-fwdP_start)/2),",             // fwd primer
                      y: ",(y_pos/max_y)+0.0005,",           // 
                      text: '\",\"',     // 
                      showarrow: false,   // 
                      font: {
                          size: 10,
                          color: '#71717a', // Dark gray text so it's easy to read
                          weight: 'bold'
                      }
                  },
                         {
                      xref: 'x',
                      yref: 'paper',
                      x: ",amplicon_end+((revP_end-amplicon_end)/2),",             // rev primer
                      y: ",(y_pos/max_y)+0.0005,",           // 
                      text: '\",\"',     // 
                      showarrow: false,   // 
                      font: {
                          size: 10,
                          color: '#71717a', // Dark gray text so it's easy to read
                          weight: 'bold'
                      }
                  }")
          
          if (annotations == "") {
            annotations <- annotation
          } else {
            annotations <- paste0(annotations, " , ", annotation)
          }
        }
      }
      
      js_sample_suffix<-gsub("[^a-zA-Z0-9_]","_",sample)
      
      constant<-paste0('
    // ==========================================
    // ',sample,' PLOT SETUP
    // ==========================================
       
       const x_axis_',js_sample_suffix,' = ',x_axis,';
       
       const y_axis_',js_sample_suffix,' = ',depth_values,';
       
       const nucleotides_',js_sample_suffix,' = ',nucleotides,';
       
        const depthTrace_',js_sample_suffix,' = {
        x: x_axis_',js_sample_suffix,',
        y: y_axis_',js_sample_suffix,',
        type: \'scatter\',
        mode: \'lines\',
        name: \'Read Depth\',
        fill: \'tozeroy\',
        fillcolor: \'rgba(59, 130, 246, 0.15)\',
        line: { color: \'#3b82f6\', width: 2 },
        hovertemplate: \'Position: %{x}<br>Depth: %{y}<extra></extra>\'
    };
       
        const sequenceTrace_',js_sample_suffix,' = {
        x: x_axis_',js_sample_suffix,',
        y: x_axis_',js_sample_suffix,'.map(()=>4200),
        yref: \'paper\',
        type: \'scatter\',
        mode: \'text\',
        text: nucleotides_',js_sample_suffix,',
        textposition: \'top center\',
        name: \'Reference Base\',
        font: {
            family: \'Courier New, monospace\',
            size: 13,
            color: \'#0f172a\',
            weight: \'bold\'
        },
        hovertemplate: \'<b>Position:</b> %{x} bp<br><b>Base:</b> %{text}<extra></extra>\'
    };
       
        const variantTrace_',js_sample_suffix,' = {
        x: ',variant_x_str,',
        y: ',variant_y_str,',
        type: \'bar\',
        mode: \'markers+text\',
        name: \'Variants Found\',
        marker: {
          color: \'#ef4444\', 
          line: { color: \'#ffffff\', width: 1 } 
        },
        text: ',variant_text_str,',
        textposition: \'outside\',
        textfont: {
          size: 30,
          color: \'#1e293b\',
          weight: \'bold\'
        },
        hovertemplate: \'<b>Mutation Point</b><br>Position: %{x}<br>Variant Depth: %{y}<extra></extra>\'
    };
    
        const layout_',js_sample_suffix,' = {
        title: { text: \'Sequence Visualizer Asset Panel - ',sample,' - ',virus,'\', font: { color: \'#334155\', size: 16 } },
        hovermode: \'closest\',
        xaxis: { title: \'Genome Coordinate Position (bp)\', gridcolor: \'#f1f5f9\', zeroline: false,
        showspikes: true,
        spikemode: \'across\',
        spikedash: \'dash\',
        spikethickness: 1,
        spikecolor: \'#94a3b8\'},
        yaxis: { title: \'Sequence Amplification Scale\', gridcolor: \'#f1f5f9\', zerolinecolor: \'#cbd5e1\',
        range:[0 , 5000]},
        margin: { t: 50, b: 50, l: 60, r: 30 },
        plot_bgcolor: \'#ffffff\',
        paper_bgcolor: \'#ffffff\',
        legend: { orientation: \'h\', y: -0.15, x: 0.5, xanchor: \'center\' },
        
        shapes: [',amplicon_regions,'],
        annotations: [',annotations,']
    };
    
        Plotly.newPlot(\'genomic-plot-',js_sample_suffix,'\', [depthTrace_',js_sample_suffix,', sequenceTrace_',js_sample_suffix,', variantTrace_',js_sample_suffix,'], layout_',js_sample_suffix,', config);
')
      all_constants <- paste0(all_constants, constant)
    }
  }
  
  # Inject the sample-specific CSS classes after they have all been collected
  sample_names_css <- paste(unique(all_samples_for_css), collapse = ", ")
  header <- sub('<style>', paste0('<style>\n        ', sample_names_css, ' '), header)
  
  finish<-paste0("
</script>
</body>
</html>") 
  
  html<-c(header,all_div_classes,script_setup,all_constants,finish)
  writeLines(html,here::here(final_output_dir,paste0(virus,".html")))
}
