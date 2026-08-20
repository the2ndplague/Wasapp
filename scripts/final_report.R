library(here)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(grid)
library(ggrepel)
library(gridExtra)
run_name<-"exp16_WGS"
primer_version<-"Chikungunya_WGS"
args<-commandArgs(trailingOnly = TRUE)
run_name<-args[1]
primer_version<-args[2]

dir.create(here("..","output",run_name,primer_version,"final_report"),showWarnings = FALSE)
dir.create(here("..","output",run_name,primer_version,"final_report","plots"),showWarnings = FALSE)

#create plots
# ==================================================================================

# ==================================================================================
# get adapter_content
checklist<-read.csv(here("..","output",run_name,"sequali_summary.csv"),row.names=1)
checklist<-checklist[,c(1,11,12,13)]
names(checklist)<-c("id","adapter content present?","highest AC%","bins with AC% > 1")

#plot grid
explanation<-paste0("Explanation: Adapter content (AC%) from sequali
              Sequali flags any samples with AC% > 1%
              If none are present, AC graph shows no data\n
              highest AC%: highest AC% among all bins
              bins with AC%>1: no. of bins with AC%>1
")

#change the color for certain columns
value_matrix <- as.matrix(checklist)
text_colors <- matrix(
  "black",
  nrow = nrow(value_matrix),
  ncol = ncol(value_matrix)
)
col_idx <- 3
rows_to_red <- value_matrix[, col_idx] > 5


text_colors[rows_to_red,] <- "red"
#set theme for tablegrob
my_theme <- ttheme_default(
  core = list(
    fg_params = list(cex = 1, fontface = "plain", col = text_colors),  # table body font
    bg_params = list(fill = c("white", "#f9f9f9"), col = "gray90") # alternating row color + grid line color
  ),
  colhead = list(
    fg_params = list(cex = 1.1, fontface = "bold", col = "white"),
    bg_params = list(fill = "steelblue")
  ),
  rowhead = list(
    fg_params = list(cex = 1, fontface = "bold")
  )
)
#automatically set table size limits
n_rows <- nrow(checklist)
n_cols <- ncol(checklist)
col_width <- 50   # pixels per column
row_height <- 50   # pixels per row
header_height <- 70
margin <- 20       # extra margin
png_width  <- n_cols * col_width + margin
png_height <- n_rows * row_height + header_height + margin

#plot table
table_plot <- tableGrob(
  checklist,
  rows = NULL,
  theme = my_theme
)

ggsave(here("..","output",run_name,primer_version,"final_report","plots","summary_table.png"),
       plot=table_plot,
       width=png_width,
       height=png_height,
       units="mm",limitsize=FALSE)

#clean environment
rm(list = setdiff(ls(),c("run_name","primer_version")))
# ==================================================================================
# Q20 reads bar chart
explanation<-paste0("explanation: see title")
sequali_summary<-read.csv(here("..","output",run_name,"sequali_summary.csv"),row.names="Row.names")
q20_data<-sequali_summary[,c("q20_reads","total_reads")]

q20_data <- q20_data %>%
  # 1. Turn your row names (which are the actual barcodes) into a proper column
  tibble::rownames_to_column("barcode") %>%
  
  # 2. Since your columns are ALREADY the metrics, do the math directly!
  mutate( 
    `Q20 reads` = q20_reads,
    `Non-Q20 reads` = total_reads - q20_reads
  ) %>%
  
  # 3. Select only the columns needed for the chart
  select(barcode, `Q20 reads`, `Non-Q20 reads`) %>%
  
  # 4. Melt it down for ggplot2
  pivot_longer(
    cols = -barcode,
    names_to = "type",
    values_to = "reads"
  )


#add percentage
q20_data <- q20_data %>%
  group_by(barcode) %>%
  mutate(
    total_reads = sum(reads),
    percentage = reads / total_reads * 100
  ) %>%
  ungroup()
#add textcolumn
q20_data$geom_text<-paste0(q20_data$reads," (",sprintf("%.1f%%",q20_data$percentage),")")
q20_data <- q20_data %>%
  mutate(label_repel = total_reads <= 2000)

#create dataframe for total reads
q20_totals <- q20_data %>%
  group_by(barcode) %>%
  summarise(total_reads = unique(total_reads))

#plot
q20_barplot<-ggplot(q20_data, aes(
  x = barcode,
  y = reads,
  fill = type
)) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = c(
      "Q20 reads" = "steelblue",
      "Non-Q20 reads" = "grey80"
    )
  )+
  geom_text( #add total reads on top of bars
    data = q20_totals,
    inherit.aes=FALSE,
    aes(
      x = barcode,
      y = total_reads,
      label = paste0("Total: ", total_reads)
    ),
    vjust = -0.5,
    size = 2,
    fontface = "bold"
    )+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  ylab("Read count") +
  #labels
  ## 🔽 repel only for small stacked bars
  geom_text_repel(
    data = subset(q20_data, label_repel),
    aes(label = geom_text),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  geom_text(
    data = subset(q20_data, !label_repel),
    aes(label = geom_text),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  labs(title="Q20 vs non-Q20 reads per barcode (sequali)",
       subtitle=explanation)

q20_barplot
ggsave(
  here("..","output",run_name,primer_version,"final_report","plots","q20_barplot.png"),
  plot = q20_barplot,
  width = 210,
  height = 220,   # 👈 taller
  units = "mm",
  dpi = 300
)



#clean environment
rm(list = setdiff(ls(),c("run_name","primer_version")))




# ==================================================================================
# heatmap for sequencing depth
explanation<-paste0("depth (cv) for each amplicon. cv threshold = 0.1")
heatmap_data<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","depth","depth_summary.csv"))
heatmap_data <- heatmap_data %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))
cv_data<-read.csv(here("..","output",run_name,primer_version,"scrapped_tables","depth","depth_cv_summary.csv"))
cv_data <- cv_data %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))

#clean data for heatmap
heatmap_data$ref<-gsub("_et_Z_primer_sk_Z_SQK_NBD114_96","",heatmap_data$ref)
heatmap_data[is.na(heatmap_data)]<-0
cv_data$ref<-gsub("_et_Z_primer_sk_Z_SQK_NBD114_96","",cv_data$ref)
cv_data[is.na(heatmap_data)]<-0
#process heatmap data to long version
heatmap_data <- heatmap_data %>%
  pivot_longer(
    cols = -ref,        # all barcode columns
    names_to = "barcode",
    values_to = "depth"
  )

cv_data <- cv_data %>%
  pivot_longer(
    cols = -ref,        # all barcode columns
    names_to = "barcode",
    values_to = "cv"
  )

heatmap_data<-left_join(heatmap_data,cv_data,by=c("ref","barcode"))
cv_data$label<-sprintf("%.1f (%.1f)", heatmap_data$depth, heatmap_data$cv)

depth_heatmap<-ggplot(cv_data, aes(barcode, ref, fill = cv)) +
  labs(title="Heatmap: depth per amplicon and barcode",
       subtitle=explanation)+
  geom_tile() +
  geom_text(aes(label = label), size = 2) + #add text to heatmap
  scale_fill_gradientn(colours = c("green", "white", "lightcoral"), #
                       values  = rescale(c(0, 0.1, 1)),
                       limits  = c(0, 1),
                       oob     = squish
  ) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.subtitle = element_text(size=7,lineheight=0.9)
  )


ggsave(
  here("..","output",run_name,primer_version,"final_report","plots","depth_heatmap.png"),
  plot = depth_heatmap,
  width = 210,
  height = 300,   # 👈 taller
  units = "mm",
  dpi = 300
)



#clean environment
rm(list = setdiff(ls(),c("run_name","primer_version")))
# ==================================================================================
# heatmap for reads
explanation<-paste0("explanation: each cell represents reads count for a specific amplicon (row) and barcode (column)")
summary<-read.csv(here("..","output",run_name,primer_version,"reads_summary.csv"))
#clean data for heatmap
heatmap_data<-summary[-c(1:6),]
heatmap_data$amplicon<-gsub("_et_Z_primer_sk_Z_SQK_NBD114_96","",heatmap_data$amplicon)
heatmap_data[is.na(heatmap_data)]<-0
#process heatmap data to long version
heatmap_data <- heatmap_data %>%
  pivot_longer(
    cols = -amplicon,        # all barcode columns
    names_to = "barcode",
    values_to = "count"
  )

reads_heatmap<-ggplot(heatmap_data, aes(barcode, amplicon, fill = count)) +
  labs(title="Heatmap: reads per amplicon and barcode",
       subtitle=explanation)+
  geom_tile() +
  geom_text(aes(label = count), size = 2) + #add text to heatmap
  scale_fill_gradientn(colours = c("lightcoral", "white", "green"),
                       values  = rescale(c(0, 100, 4000)),
                       limits  = c(0, 4000),
                       oob     = squish
  ) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.subtitle = element_text(size=7,lineheight=0.9)
  )

reads_heatmap
ggsave(
  here("..","output",run_name,primer_version,"final_report","plots","reads_heatmap.png"),
  plot = reads_heatmap,
  width = 210,
  height = 300,   # 👈 taller
  units = "mm",
  dpi = 300
)
#clean environment
rm(list = setdiff(ls(),c("run_name","primer_version")))
# ==================================================================================
# heatmap for mean quality
explanation<-paste0("explanation: quality scores above 80% of all scores- green, below - red")
summary<-read.csv(here("..","output",run_name,primer_version,"meanquality_summary.csv"))
#clean data for heatmap
heatmap_data<-summary[-c(1:6),]
heatmap_data$amplicon<-gsub("_et_Z_primer_sk_Z_SQK_NBD114_96","",heatmap_data$amplicon)
heatmap_data[is.na(heatmap_data)]<-0
#process heatmap data to long version
heatmap_data <- heatmap_data %>%
  pivot_longer(
    cols = -amplicon,        # all barcode columns
    names_to = "barcode",
    values_to = "count"
  )
min<-min(heatmap_data$count)
max<-max(heatmap_data$count)
reads_heatmap<-ggplot(heatmap_data, aes(barcode, amplicon, fill = count)) +
  labs(title="Heatmap: quality per amplicon and barcode",
       subtitle=explanation)+
  geom_tile() +
  geom_text(aes(label = round(count,2)), size = 2) + #add text to heatmap
  scale_fill_gradientn(colours = c("red", "white", "green"),
                       values  = rescale(c(min, max*0.8, max)),
                       limits  = c(min, max),
                       oob     = squish
  ) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.subtitle = element_text(size=7,lineheight=0.9)
  )


ggsave(
  here("..","output",run_name,primer_version,"final_report","plots","meanquality_heatmap.png"),
  plot = reads_heatmap,
  width = 210,
  height = 300,   # 👈 taller
  units = "mm",
  dpi = 300
)
#clean environment
rm(list = setdiff(ls(),c("run_name","primer_version")))
# ==================================================================================
# heatmap for mean ref.coverage
explanation<-paste0("explanation: coverage values above 90% of all scores- green, below - red. Mean coverage is obtained by averaging all coverage values for each amplicon and barcode")
summary<-read.csv(here("..","output",run_name,primer_version,"coverage_summary.csv"))
#clean data for heatmap
heatmap_data<-summary[-c(1:6),]
heatmap_data$amplicon<-gsub("_et_Z_primer_sk_Z_SQK_NBD114_96","",heatmap_data$amplicon)
heatmap_data[is.na(heatmap_data)]<-0
#process heatmap data to long version
heatmap_data <- heatmap_data %>%
  pivot_longer(
    cols = -amplicon,        # all barcode columns
    names_to = "barcode",
    values_to = "count"
  )
min<-min(heatmap_data$count)
max<-max(heatmap_data$count)

reads_heatmap<-ggplot(heatmap_data, aes(barcode, amplicon, fill = count)) +
  labs(title="Heatmap: mean coverage per Amplicon and Barcode ",
       subtitle=explanation)+
  geom_tile() +
  geom_text(aes(label = round(count,2)), size = 2) + #add text to heatmap
  scale_fill_gradientn(colours = c("lightcoral", "white", "green"),
                       values  = rescale(c(0, 90, 100)),
                       limits  = c(0, 100),
                       oob     = squish
  ) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.subtitle = element_text(size=7,lineheight=0.9)
  )

ggsave(
  here("..","output",run_name,primer_version,"final_report","plots","coverage_heatmap.png"),
  plot = reads_heatmap,
  width = 210,
  height = 300,   # 👈 taller
  units = "mm",
  dpi = 300
)
#clean environment
rm(list = setdiff(ls(),c("run_name","primer_version")))
# ===========================================================================





# ===========================================================================
#get samples in length distribution
files<-list.files(here("..","output",run_name,primer_version,"length_distribution"))
files<-files[!files%in%c("read_length_distribution_summary.csv")]
files<-files[-c(grep(".csv",files))]
number<-length(files)
#get the charts
charts<-paste0()
for(file in files){
  chart<-paste0('<img src="../length_distribution/',file,'" style="width:90%;">')
  charts<-paste0(charts,
                 chart)
}
length_charts<-paste0('<h2 style="text-align:center;"><br><br>Length distribution charts<br><br></h2>
  <div style="
  display:grid;
  grid-template-columns: repeat(1,',number, 'fr);
  gap: 20px;
  justify-items: center;
">',charts
  ,'</div>')


#clean environment
rm(list = setdiff(ls(),c("run_name","length_charts","primer_version")))

# ===========================================================================
title<-paste0('<div style="width: 50%; text-align: center;">
                  <h3>Run performed for version: ',run_name,primer_version,'</h3>
                  <p></p>
                  </div>')
# read in and make starting page with notebooklm summary
bot_link<-"https://notebooklm.google.com/notebook/0982a2a4-eb35-4329-b96a-ddf38e769f0b?authuser=2"

prompt<-readChar(here("..","output",run_name,primer_version,"final_report","prompt","prompt.txt"),file.info(here("..","output",run_name,primer_version,"final_report","prompt","prompt.txt"))$size)
prompt_info<-paste0("chat with bot here:\n", bot_link)
prompt<-paste0('<p style="text-align:center;">',prompt, prompt_info,"</p>")
prompt <- gsub("\n", "<br>", prompt)
prompt<-gsub("\\*\\*(.*?)\\*\\*", "<strong>\\1</strong>", prompt, perl = TRUE)
# ===========================================================================
# low SNR
snr<-readChar(here("..","output",run_name,primer_version,"final_report","prompt","low_SNR.txt"),file.info(here("..","output",run_name,primer_version,"final_report","prompt","low_SNR.txt"))$size)
snr <- gsub("\n", "<br>", snr)
snr<-gsub("\\*\\*(.*?)\\*\\*", "<strong>\\1</strong>", snr, perl = TRUE)

#combine notebooklm outputs:
notebooklm<-paste0('<div style="display: flex; gap: 20px;">
  <div style="width: 50%; text-align: center;">
  <h3>Notebooklm summary</h3>
  <p>',prompt,'</p>
  </div>
  
  <div style="width: 50%; text-align: center;">
  <h3>Prompt: compare each sample and identify which performed the best across all amplicons</h3>
  <p>',snr,'</p>
  </div>
  </div>')

html<-readLines(here("..","output",run_name,primer_version,"wf-amplicon","wf-amplicon-report.html"))
new_block <- paste0(
  '
    <h2 style="text-align:center;"><br><br>Sequali metrics<br><br></h2>
    <div style="
  display:grid;
  grid-template-columns: repeat(2, 2fr);
  gap: 20px;
  justify-items: center;
">
    <img src="../final_report/plots/q20_barplot.png" style="width:90%;">
    <img src="../final_report/plots/cs_heatmap.png" style="width:90%;">
    <img src="../final_report/plots/summary_table.png" style="width:90%;">
    </div>',length_charts,'
  
  <h2 style="text-align:center;">Amplicon/sample specific heatmaps<br><br></h2>
    <div style="
  display:grid;
  grid-template-columns: repeat(2, 2fr);
  gap: 20px;
  justify-items: center;
">
    <img src="../final_report/plots/reads_heatmap.png" style="width:90%;">
    <img src="../final_report/plots/depth_heatmap.png" style="width:90%;">
    <img src="../final_report/plots/coverage_heatmap.png" style="width:90%;">
    <img src="../final_report/plots/meanquality_heatmap.png" style="width:90%;">
    </div>
  
'
)

# Collapse html lines into a single string
html_content <- paste(html, collapse = "\n")

# Combine summaries into a section
summary_section <- paste0('
        <section class="container p-4 mb-4 bg-white border rounded" id="Section_run_summary" tagname="section">
          <h2 class="h5 mb-0 pb-3">Run Summary (NotebookLM)</h2>
          ', title, '
          ', notebooklm, '
        </section>
')

# Wrap sequali and heatmaps block inside proper Bootstrap containers
visuals_section <- paste0('
        <section class="container p-4 mb-4 bg-white border rounded" id="Section_sequali_and_heatmaps" tagname="section">
          ', new_block, '
        </section>
')

# Insert them after the main-content section starts
main_content_start <- '<section id="main-content" role="region">'
if (grepl(main_content_start, html_content, fixed = TRUE)) {
  html_content <- sub(main_content_start, paste0(main_content_start, "\n", summary_section, "\n", visuals_section), html_content, fixed = TRUE)
}

# Add jump-to items in the dropdown menu
intro_pattern <- '<a class="dropdown-item" href="#Section_679359762d0a4728ba5fe7548092b6b3">Intro.</a>'
if (grepl(intro_pattern, html_content, fixed = TRUE)) {
  html_content <- sub(intro_pattern, paste0('<a class="dropdown-item" href="#Section_run_summary">Summary</a><a class="dropdown-item" href="#Section_sequali_and_heatmaps">QC & Heatmaps</a>', intro_pattern), html_content, fixed = TRUE)
}

# Write back to file
writeLines(html_content, here("..","output",run_name,primer_version,"wf-amplicon","wf-amplicon-report_modified.html"), useBytes = TRUE)





 