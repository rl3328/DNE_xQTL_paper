# qc result
sos_template=./sos_template_rss_gwas.sh
QC_method="rss_qc"
FINEMAP_method="susie_rss"
sos_submit=./sos_submit_QC_${QC_method}_finemap_${FINEMAP_method}.sh

# Input the region list, note: this may subject to change of the format 
# Here it serves convert the first three columns to format of chr:start-end
region_list=$(cat ./ld_meta_file.tsv \
    | sed '1d' \
    | awk -F'\t' '{gsub(/^chr/, "", $1); print $1 ":" $2 "-" $3}')

echo "$region_list" | while IFS= read -r line; do 
    processed_line=$(sed "s/fill_in_region/${line}/g" "${sos_template}" | sed "s/QC_METHOD/${QC_method}/g" | sed "s/FINEMAP/${FINEMAP_method}/g")
    echo "$processed_line" >> "${sos_submit}"
done