#!/usr/bin/env bash

# Configuration parameters (overridden by environment variables if set)
run_name="${RUN_NAME:-exp12}"
pod5_dir="${POD5_DIR:-../input/pod5}"
dorado_model_path="${DORADO_MODEL_PATH:-../input/packages/dorado-1.4.0-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0}"
kit_name="${KIT_NAME:-SQK-NBD114-96}"
min_qscore="${MIN_QSCORE:-10}"
primer_files_dir="${PRIMER_FILES_DIR:-../input/primer_files}"
threads="${THREADS:-8}"
notebooklm_exec="${NOTEBOOKLM_EXEC:-notebooklm}"

# Save the absolute path of the scripts directory to return to it reliably
scripts_dir="$(pwd)"

# Extract basename of model path for wf-amplicon override configuration
dorado_model_name=$(basename "${dorado_model_path}")

echo $run_name > ../output/run_name.txt
echo "exp15 cutadapt" > ../output/$run_name/README.txt

#basecalling and demultiplex using dorado
###############################################################################
dorado basecaller "${dorado_model_path}" "${pod5_dir}" -r --kit-name "${kit_name}" --min-qscore "${min_qscore}" --no-trim  --emit-fastq --barcode-both-ends -o ../output/$run_name/basecalling/a.basecalled
#
dorado demux ../output/$run_name/basecalling/a.basecalled -r --kit-name "${kit_name}" --emit-fastq  --barcode-both-ends -o ../output/$run_name/basecalling/b.basecalled_demux
# 


##############################################################################
#move all the fastq files out
base_dir="../output/${run_name}/basecalling/b.basecalled_demux"
# Find all .fastq files recursively
find "$base_dir" -type f -name "*.fastq" | while read -r file; do
    # Define destination path (keep original filename)
    dest="${base_dir}/$(basename "$file")"  
    # Move the file
    mv "$file" "$dest"
    echo "Moved: $file -> $dest"
done

rm -rf ../output/${run_name}/basecalling/b.basecalled_demux/no_sample
#rename fastq files in (b) and (a)
Rscript ./process_fastq.R $run_name


#wasPPipe######################################################################
mkdir -p ../output/$run_name/basecalling/a1.basecalled_wasPPipe/ 
cp -r ../output/$run_name/basecalling/a.basecalled/ ../output/$run_name/basecalling/a1.basecalled_wasPPipe/ 
find ../output/$run_name/basecalling/a1.basecalled_wasPPipe/a.basecalled -type f -name "*.fastq" -exec gzip {} +
nextflow run ./WaSPPipe-1.0.1/main.nf --fastq ../output/$run_name/basecalling/a1.basecalled_wasPPipe/a.basecalled/*/*/fastq_pass/

mv ./output ../output/$run_name/wasPPipe
rm -rf ../output/$run_name/basecalling/a1.basecalled_wasPPipe/ 






while IFS= read -r primer_version; do
    echo " processing data for version: $primer_version"
#=================================================
#=================================================
#=================================================
cp -r ../output/${run_name}/${primer_version}/basecalling/b.basecalled_demux/. ../output/${run_name}/${primer_version}/basecalling/c.basecalled_demux_trim

#cutadapt
for f in ../output/$run_name/${primer_version}/basecalling/b.basecalled_demux/*.fastq; do
    base=$(basename "$f")
    
    cutadapt -g file:"${primer_files_dir}/${primer_version}_linked_primers.fasta" \
        -q 20 \
        -o ../output/$run_name/${primer_version}/basecalling/c.basecalled_demux_trim/${base%.fastq}.fastq \
        "$f" \
        --discard-untrimmed \
    	--minimum-length 20 
done


cd ../output/${run_name}/${primer_version}/basecalling/c.basecalled_demux_trim
#put the fastq files in their own folder 
for file in *.fastq; do
    base="${file}"    
    # Remove the .fastq extension from base
    name="${base%.fastq}"    
    # Create folder
    mkdir -p "$name"    
    # Rename the file to name.fastq and move into the folder
    mv "$file" "$name/$name.fastq"
done
cd "${scripts_dir}"


#sequali qc
###############################################################################
find ../output/${run_name}/${primer_version}/basecalling/c.basecalled_demux_trim -name "*.fastq" | while read fastq_file; do
    # Get the file name and its directory
    filename=$(basename "$fastq_file")
    # Create an output	 directory for the report
    mkdir -p ../output/$run_name/sequali_report

    # Run sequali on the file
    sequali "$fastq_file" --outdir ../output/$run_name/sequali_report
done

Rscript ./Process_sequali.R $run_name
Rscript ./sequali_replacement.R $run_name $primer_version

#run amplicon-wf
###############################################################################
nextflow run epi2me-labs/wf-amplicon \
    --fastq ../output/${run_name}/${primer_version}/basecalling/c.basecalled_demux_trim \
    --analyse_unclassified True \
    --reference "${primer_files_dir}/${primer_version}.fasta" \
    -profile standard \
    --reads_downsampling_size 0 \
    --min_read_length 200 \
    --max_read_length 600 \
    --min_read_qual 10 \
    --drop_frac_longest_reads 0 \
    --min_n_reads 40 \
    --take_longest_remaining_reads False \
    --min_coverage 20 \
    --force_spoa_length_threshold 2000 \
    --spoa_minimum_relative_coverage 0.15 \
    --minimum_mean_depth 30 \
    --override_basecaller_cfg "${dorado_model_name}"\
    --igv \
    --number_depth_windows 100000 \
    --medaka_target_depth_per_strand 10000 \
    --out_dir ../output/$run_name/$primer_version/wf-amplicon \
    --threads "${threads}"
#downstream cleaning of data using R and python
###############################################################################
mkdir -p ../output/$run_name/$primer_version/scrapped_tables/bamstats
mkdir -p ../output/$run_name/$primer_version/scrapped_tables/depth
mkdir -p ../output/$run_name/$primer_version/scrapped_tables/window_depth

python tablemaker.py ../output/$run_name/$primer_version/wf-amplicon/wf-amplicon-report.html ../output/$run_name/$primer_version/scrapped_tables

Rscript ./Process_tables.R $run_name $primer_version

###############################################################################
#get variant coordinates



###############################################################################
#get summary from notebooklm
mkdir -p ../output/$run_name/$primer_version/final_report/prompt #make prompt directory

"${notebooklm_exec}" delete -y
#make new workbook
"${notebooklm_exec}" create "waspp read report"
#use new notebook
"${notebooklm_exec}" list > ../output/$run_name/$primer_version/final_report/prompt/workbook_list.txt
id=$(grep "waspp" ../output/$run_name/$primer_version/final_report/prompt/workbook_list.txt \
  | tr -d '│ ' \
  | grep -oE '[a-f0-9-]{6,}' \
  | head -n 1)
id="${id:0:6}"
echo "$id"
"${notebooklm_exec}" use $id

#upload sources
"${notebooklm_exec}" source add ../output/$run_name/$primer_version/reads_summary.csv
"${notebooklm_exec}" source add ../output/$run_name/$primer_version/meanquality_summary.csv
"${notebooklm_exec}" source add ../output/$run_name/$primer_version/coverage_summary.csv
"${notebooklm_exec}" source add ../output/$run_name/$primer_version/sequali_summary.csv
"${notebooklm_exec}" source add ../output/$run_name/$primer_version/length_distribution/read_length_distribution_summary.csv

#prompt 1:
"${notebooklm_exec}" ask "you are a bioinformatician. read all the .csv files and describe any differences between amplicons and samples. ignore the unclassified group and deliver any insights you might find on the following topics, in 100 words less for each:

1. percentage of q20 reads over total reads should be above 50% (look in sequali summary)
2. the bulk of the read sizes of sequence length should be between 200 to 500 bp long. identify any outlier peaks for any samples (look in read_length_distribution_summary)
3. coverage of amplicons should ideally be >99%. Identify samples and amplicons that do not meet this criteria. (look in Coverage_summary) 
4. the targeted reads for each amplicon should be at least 100. Identify samples and amplicons that do not meet this criteria (look in Reads_summary) 
5. identify any amplicon noise in the low concentration/NTC samples, where reads should be low (look in Reads_summary)

provide 1 point for each topic only and nothing else." > ../output/$run_name/$primer_version/final_report/prompt/prompt.txt


#prompt 2:
"${notebooklm_exec}" ask "search for any patterns within and between the samples" > ../output/$run_name/$primer_version/final_report/prompt/low_SNR.txt

Rscript ./final_report.R $run_name $primer_version
#Rscript ./final_report_v2.R $run_name $primer_version
Rscript ./final_report_WGS.R $run_name $primer_version

done < ../output/${run_name}/versions.txt


#make big summary
cp -r ../output/make_big_summary_v3.R ../output/${run_name}/
Rscript ../output/${run_name}/make_big_summary_v3.R $run_name
