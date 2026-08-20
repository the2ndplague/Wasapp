
#this script consists of 2 parts:
#part 1: preparation of data for plotting and pdf files
#part 2: placing these data into .html files for viewing


library(here)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(grid)
library(ggrepel)
library(gridExtra)
library(patchwork)
run_name<-"exp16_WGS"
primer_version<-"COVID_WGS"
args<-commandArgs(trailingOnly = TRUE)
run_name<-args[1]
primer_version<-args[2]

dir.create(here("..","output",run_name,primer_version,"final_report"),showWarnings = FALSE)
dir.create(here("..","output",run_name,primer_version,"final_report","plots"),showWarnings = FALSE)
dir.create(here("..","output",run_name,primer_version,"final_report","graphing_data"),showWarnings = FALSE)

#PART 1
##########################################################################################################
#read in genomes and create dataframe
genomes<-readLines(here("..","input","panel_genomes","panel_genomes_FWD.fasta"))
genome_names <- genomes[seq(1, length(genomes), by = 2)]
genome_names <- sub("^>", "", genome_names)
sequence_strings <- genomes[seq(2, length(genomes), by = 2)]
genome_sequences_df <- data.frame(
  genome_name = genome_names,
  sequence    = sequence_strings,
  stringsAsFactors = FALSE
)
gsub("_FWD","",genome_sequences_df$genome_name)->genome_sequences_df$genome_name



#read in variants
variants<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","table_5.csv"))


mastersheet<-read.csv(here("..","input","primer_files",paste0(primer_version,"_mastersheet.csv")))

samples<-list.files(here("..","output",run_name,primer_version,"basecalling","b.basecalled_demux"))
samples<-gsub(".fastq","",samples)



for(sample in samples){
  if (!file.exists(here("..","output",run_name,primer_version,"scrapped_tables","depth",paste0(sample,".csv")))) {
    message(paste("Skipping, file not found:", sample))
    next # Skips the rest of this loop iteration
  }
  
  depth<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","depth",paste0(sample,".csv")))
  
  mastersheet_temp <- mastersheet[order(mastersheet$amplicon_name), ]
  mastersheet_temp$wfamplicon_name<-sort(unique(depth$ref))
  
  viruses<-unique(mastersheet_temp$Virus)
  
  virus_plots<-list()
  i<-0
  giant_df<-data.frame(start=c(1:3000))
  for(virus in viruses){
    i<-i+1
    amplicons<-mastersheet_temp$wfamplicon_name[mastersheet_temp$Virus==virus]
    
    variants[variants$Sample==sample&variants$Amplicon%in%amplicons,]
    
    variant_temp<-variants[variants$Sample==sample&variants$Amplicon%in%amplicons,]
    genome_positions<-c()
    for(variant in c(1:nrow(variant_temp))){
      amplicon<-variant_temp$Amplicon[variant]  
      #start pos of amplicon
      variant_pos_genome<-mastersheet_temp$amp_start[mastersheet_temp$wfamplicon_name==amplicon]+variant_temp$Position[variant]  
      genome_positions<-c(genome_positions,variant_pos_genome)
    }
    genome_positions<-genome_positions[!is.na(genome_positions)]
    variant_temp$genome_position<-genome_positions
   
    virus_sequence<-genome_sequences_df$sequence[genome_sequences_df$genome_name==virus]
    virus_df<-data.frame(start=c(1:nchar(virus_sequence)),
                         base=unlist(strsplit(virus_sequence, split = "")),
                         depth=0,
                         variant_depth=0,
                         variant_description="")
    
    idx <- match(variant_temp$genome_position, virus_df$start)
    virus_df$variant_depth[idx] <- virus_df$variant_depth[idx] + variant_temp$Depth
    virus_df$variant_description[idx] <- paste0(variant_temp$Type, ": ", variant_temp$Ref..allele, "->", variant_temp$Alt..allele)
    
    for(amplicon in amplicons){ #get depth
      depth_temp<-depth[depth$ref==amplicon,] #correct for position on genome
      depth_temp$start<-depth_temp$start+mastersheet_temp$amp_start[mastersheet_temp$wfamplicon_name==amplicon]
      depth_temp$end<-depth_temp$end+mastersheet_temp$amp_start[mastersheet_temp$wfamplicon_name==amplicon]
      
      incoming_depths <- depth_temp %>%
        group_by(end) %>%
        summarize(summed_depth = sum(depth, na.rm = TRUE))
      
      virus_df <- virus_df %>%
        left_join(incoming_depths, by = c("start" = "end")) %>%
        mutate(depth = depth + replace_na(summed_depth, 0)) %>%
        select(-summed_depth)
      
    }#amplicon
    
    #trim off unused ends
    virus_df <- virus_df[c(min(mastersheet_temp[mastersheet_temp$Virus==virus,"amp_start"]):
                             max(mastersheet_temp[mastersheet_temp$Virus==virus,"amp_end"])),]
    
    current_amplicons <- mastersheet[mastersheet$Virus == virus, ]
    current_amplicons<-current_amplicons[order(current_amplicons$amp_start), ]
    max_depth <- max(virus_df$depth, na.rm = TRUE)
    
    current_amplicons$y_pos <- -max_depth * rep(c(0.1, 0.2), length.out = nrow(current_amplicons))

    
    
    virus_df$primer   <- ""
    virus_df$amplicon <- ""
    
    for (row in 1:nrow(current_amplicons)) {
      
      name  <- current_amplicons$amplicon_name[row]
      a_st  <- current_amplicons$amp_start[row]
      a_en  <- current_amplicons$amp_end[row]
      f_len <- nchar(current_amplicons$Fwd_primer[row])
      r_len <- nchar(current_amplicons$Rev_primer[row])
      
      p_fwd_start <- a_st - f_len
      p_fwd_end   <- a_st
      p_rev_start <- a_en
      p_rev_end   <- a_en + r_len
      
      idx_fwd_start <- match(p_fwd_start, virus_df$start)
      idx_amp_start <- match(p_fwd_end,   virus_df$start) # amp_start
      idx_amp_end   <- match(p_rev_start, virus_df$start) # amp_end
      idx_rev_end   <- match(p_rev_end,   virus_df$start)
      
      if (!is.na(idx_amp_start)) {
        virus_df$amplicon[idx_amp_start] <- ifelse(virus_df$amplicon[idx_amp_start] == "", 
                                                   paste0("amplicon start: ", name), paste0(virus_df$amplicon[idx_amp_start], "; amplicon start: ", name))
      }
      if (!is.na(idx_amp_end)) {
        virus_df$amplicon[idx_amp_end] <- ifelse(virus_df$amplicon[idx_amp_end] == "", 
                                                 paste0("amplicon end: ", name), paste0(virus_df$amplicon[idx_amp_end], "; amplicon end: ", name))
      }
      
      if (!is.na(idx_fwd_start)) {
        virus_df$primer[idx_fwd_start] <- ifelse(virus_df$primer[idx_fwd_start] == "", 
                                                 paste0("primer start: ", name), paste0(virus_df$primer[idx_fwd_start], "; primer start: ", name))
      }
      if (!is.na(idx_amp_start)) {
        virus_df$primer[idx_amp_start] <- ifelse(virus_df$primer[idx_amp_start] == "", 
                                                 paste0("primer end: ", name), paste0(virus_df$primer[idx_amp_start], "; primer end: ", name))
      }
      if (!is.na(idx_amp_end)) {
        virus_df$primer[idx_amp_end] <- ifelse(virus_df$primer[idx_amp_end] == "", 
                                               paste0("primer start: ", name), paste0(virus_df$primer[idx_amp_end], "; primer start: ", name))
      }
      if (!is.na(idx_rev_end)) {
        virus_df$primer[idx_rev_end] <- ifelse(virus_df$primer[idx_rev_end] == "", 
                                               paste0("primer end: ", name), paste0(virus_df$primer[idx_rev_end], "; primer end: ", name))
      }
    }
    
    
    
    
    names(virus_df)[c(2:length(names(virus_df)))]<-paste0(virus,"_",names(virus_df))[c(2:length(names(virus_df)))]
    giant_df<-merge(giant_df,virus_df,by="start",all.x=TRUE)
    giant_df[is.na(giant_df)] <- ""
  }#virus
  write.csv(giant_df, here("..","output",run_name,primer_version,"final_report","graphing_data",paste0(sample,"_info.csv")), row.names = FALSE)
  
} #sample



#PART 2
##########################################################################################################
samples<-list.files(here("..","output",run_name,primer_version,"final_report","graphing_data"), pattern = "_info.csv") -> graphing_files
samples<-gsub("_info.csv","",samples)


for(sample in samples){
sample_file<-read.csv(here("..","output",run_name,primer_version,"final_report","graphing_data",paste0(sample,"_info.csv")))

viruses<-unique(gsub("_.*","",names(sample_file)))
viruses<-viruses[!viruses==c("start")]



virus_names<-paste0(".",paste(viruses,collapse=". "))
header<-paste0('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Interactive Genomic Plot Study</title>
    
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    
    <style>
        body {
            font-family: system-ui, -apple-system, sans-serif;
            margin: 40px;
            background-color: #f8fafc;
        }
        /* Fixed: Updated to match both card classes */
        ',virus_names,' {
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

#div class loop

div_classes<-paste0()

for(virus in viruses){
div_class<-paste0('<div class="',virus,'">
    <h1>Interactive Sequence Coverage Landscape -', virus,'</h1>
    <p>Hover over data points to read values. Click and drag to zoom. Double-click to reset view.</p>
    <div id="genomic-plot-',virus,'" style="width:100%; height:600px;"></div>
</div>
                  
')
div_classes<-paste0(div_classes,div_class)
}


script_setup<-paste0("<script>
    // Global configurations shared by both plots
    const config = {
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: ['select2d', 'lasso2d', 'autoScale2d']
    };")



mastersheet<-read.csv(here("..","input","primer_files",paste0(primer_version,"_mastersheet.csv")))
#constants_setup
constants<-paste0()
for(virus in viruses){
virus_df<-sample_file[,grep(virus,names(sample_file))]
virus_df<-virus_df[!is.na(virus_df[,2]),]

x_axis<-paste0("[",paste(as.numeric(row.names(virus_df)),collapse=", "),"]")
depth_values<-paste0("[",paste(as.numeric(virus_df[,2]),collapse=", "),"]")


nucleotides<-paste0("['",paste(virus_df[,1],collapse="', '"),"']")


variant_x<-variant_text<-virus_df[virus_df[,3]>0,]
variant_x<-row.names(variant_x)
variant_x<-paste0("[",paste(as.numeric(variant_x),collapse=", "),"]")

variant_y<-as.numeric(virus_df[,3])
variant_y<-variant_y[variant_y>0]
variant_y<-paste0("[",paste(variant_y,collapse=", "),"]")

variant_text<-virus_df[virus_df[,3]>0,4]
variant_text<-paste0("['",paste(variant_text,collapse="', '"),"']")


#amplicon_region
mastersheet_temp<-mastersheet[mastersheet$Virus==virus,]
mastersheet_temp<-mastersheet_temp[order(mastersheet_temp$amp_start),]
mastersheet_temp$pos<-max(virus_df[,2])* rep(c(0.95, 1), length.out = nrow(mastersheet_temp))

max_x<-max(mastersheet_temp$amp_end)
max_y<-max(virus_df[,2])

#needs to be in decimal place
amplicon_regions<-paste0()
annotations<-paste0()

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
                y0: ",(y_pos/max_y),",           // Height start (1% from bottom of grid)
                y1: ",(y_pos/max_y)+0.01,",           // Height stop (4% from bottom of grid)
                fillcolor: 'rgba(234, 179, 8, 0.4)', // Translucent Gold
                line: { width: 0 }
            },
       {
                type: 'rect',
                xref: 'x',
                yref: 'paper',
                x0: ",fwdP_start,",            // fwd primer Start (bp)
                x1: ",amplicon_start,",            // fwd primer End (bp)
                y0: ",(y_pos/max_y),",           // Height start (1% from bottom of grid)
                y1: ",(y_pos/max_y)+0.01,",           // Height stop (4% from bottom of grid)
                fillcolor: 'rgba(16, 185, 129, 0.4)', // Translucent Emerald Green
                line: { width: 0 }
            },
       {
                type: 'rect',
                xref: 'x',
                yref: 'paper',
                x0: ",amplicon_end,",            // Amplicon Start (bp)
                x1: ",revP_end,",            // Amplicon End (bp)
                y0: ",(y_pos/max_y),",           // Height start (1% from bottom of grid)
                y1: ",(y_pos/max_y)+0.01,",           // Height stop (4% from bottom of grid)
                fillcolor: 'rgba(239, 68, 68, 0.4)',  // Translucent Red
                line: { width: 0 }
            },

       ")

amplicon_regions<-paste0(amplicon_regions," , ",amplicon_region)



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
                text: '","","',     // 
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
                text: '","","',     // 
                showarrow: false,   // 
                font: {
                    size: 10,
                    color: '#71717a', // Dark gray text so it's easy to read
                    weight: 'bold'
                }
            }")

annotations<-paste0(annotations, " , ",annotation)

}

constant<-paste0('
    // ==========================================
    // ',virus,' PLOT SETUP
    // ==========================================','
       
       const x_axis',virus,' = ',x_axis,';
       
       const y_axis',virus,' = ',depth_values,";
       
       const nucleotides",virus,"= ",nucleotides,";
       
       
        const depthTrace",virus," = {
        x: x_axis",virus,",
        y: y_axis",virus,",
        type: 'scatter',
        mode: 'lines',
        name: 'Read Depth',
        fill: 'tozeroy',
        fillcolor: 'rgba(59, 130, 246, 0.15)',
        line: { color: '#3b82f6', width: 2 },
        hovertemplate: 'Position: %{x}<br>Depth: %{y}<extra></extra>'
    };
       
        const sequenceTrace",virus," = {
        x: x_axis",virus,",
        y: x_axis",virus,".map(()=>4200),
        yref: 'paper',
        type: 'scatter',
        mode: 'text',
        text: nucleotides",virus,",
        textposition: 'top center', // Pushes the letters slightly above the 0 line so they aren't sliced in half
        name: 'Reference Base',
        font: {
            family: 'Courier New, monospace', // 🌟 Crucial: Keeps widths perfectly uniform
            size: 13,
            color: '#0f172a',
            weight: 'bold'
        },
        hovertemplate: '<b>Position:</b> %{x} bp<br><b>Base:</b> %{text}<extra></extra>'
    };
       
        const variantTrace",virus," = {
        x: ",variant_x,",
        y: ",variant_y,",
        type: 'bar',
        mode: 'markers+text',
        width: 10,
        name: 'Variants Found',
        marker: {
         color: '#ef4444', 
          line: { color: '#ffffff', width: 1 } 
    },

        text: ",variant_text,",
        textposition: 'outside',
        textfont: {
        size: 30,           // Set your font size (default is 12)
        color: '#1e293b',   // Optional: Makes the text dark slate grey
        weight: 'bold'      // Optional: Makes the text bold
    },
        hovertemplate: '<b>Mutation Point</b><br>Position: %{x}<br>Variant Depth: %{y}<extra></extra>'
    };
    
        const layout",virus," = {
        title: { text: 'Sequence Visualizer Asset Panel - ",sample," - ",virus,"', font: { color: '#334155', size: 16 } },
        hovermode: 'closest',
        xaxis: { title: 'Genome Coordinate Position (bp)', gridcolor: '#f1f5f9', zeroline: false,
        showspikes: true,
        spikemode: 'across',
        spikedash: 'dash',
        spikethickness: 1,
        spikecolor: '#94a3b8'},
        yaxis: { title: 'Sequence Amplification Scale', gridcolor: '#f1f5f9', zerolinecolor: '#cbd5e1',
        range:[0 , 5000]},
        margin: { t: 50, b: 50, l: 60, r: 30 },
        plot_bgcolor: '#ffffff',
        paper_bgcolor: '#ffffff',
        legend: { orientation: 'h', y: -0.15, x: 0.5, xanchor: 'center' },
        
        // ==========================================
        // SHAPES TRACKS SETUP (Amplicons)
        // ==========================================
        shapes: [",amplicon_regions,"],
        annotations: [",annotations,"]
    };
    
        Plotly.newPlot('genomic-plot-",virus,"', [depthTrace",virus,", sequenceTrace",virus,",variantTrace",virus,"], layout",virus,", config);

")

constants<-paste0(constants,constant)
}

finish<-paste0("
</script>
</body>
</html>")



html<-c(header,div_classes,script_setup,constants,finish)

writeLines(html,here("..","output",run_name,primer_version,"final_report",paste0(sample,".html")))
}

