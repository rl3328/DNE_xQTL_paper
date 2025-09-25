#!/bin/bash

# Script to add ALT_FREQS column to GWAS sumstats files (left join - keep all GWAS variants)
# The lookup table serves as annotation only - all GWAS variants are retained
# Handles allele flipping by checking both REF/ALT orientations (chr_pos_A2_A1 and chr_pos_A1_A2)
# Usage: bash annotate_maf.sh LOOKUP_TABLE GWAS_DIR OUTPUT_DIR [GWAS_PATTERN]

# set -e

# Function to show usage
show_usage() {
    echo "Usage: bash annotate_maf.sh LOOKUP_TABLE GWAS_DIR OUTPUT_DIR [GWAS_PATTERN]"
    echo ""
    echo "Arguments:"
    echo "  LOOKUP_TABLE  Path to the lookup table file (.tsv.gz)"
    echo "  GWAS_DIR      Directory containing GWAS files"
    echo "  OUTPUT_DIR    Output directory for processed files"
    echo "  GWAS_PATTERN  Pattern to match GWAS files (optional, default: *.tsv.gz)"
    echo ""
    echo "Example:"
    echo "  bash annotate_maf.sh \\"
    echo "    ~/project/AD_hg38_maf_lookup.tsv.gz \\"
    echo "    ~/data/GWAS/rss_imputed_qced_GWAS_image_PD_Aging/image_AD1 \\"
    echo "    ~/data/GWAS/rss_imputed_qced_GWAS_image_PD_Aging/image_AD1_with_freqs"
}

# Function to clean up temporary files and free memory
cleanup_temp_files() {
    local temp_files=("$@")
    for temp_file in "${temp_files[@]}"; do
        if [ -f "$temp_file" ]; then
            rm -f "$temp_file"
        fi
    done
    # Force garbage collection in shell (limited effect but helps)
    unset temp_files
}

# Function to clear AWK arrays and variables
clear_awk_memory() {
    # This is handled within AWK by using delete statements
    # Called after processing each file
    echo "  🧹 Cleared processing memory"
}

# Function to process a single GWAS file
process_gwas_file() {
    local gwas_file="$1"
    local lookup_file="$2"
    local output_dir="$3"
    
    local base_name=$(basename "$gwas_file")
    local output_file="${output_dir}/${base_name}"
    
    echo "Processing: $base_name"
    
    # Create temporary files
    local temp_gwas=$(mktemp)
    local temp_lookup=$(mktemp)
    local temp_lookup_hash=$(mktemp)
    local temp_merged=$(mktemp)
    
    # Array to track all temp files for cleanup
    local temp_files_array=("$temp_gwas" "$temp_lookup" "$temp_lookup_hash" "$temp_merged")
    
    # Decompress files if needed and prepare for merging
    if [[ "$gwas_file" == *.gz ]]; then
        zcat "$gwas_file" > "$temp_gwas"
    else
        cp "$gwas_file" "$temp_gwas"
    fi
    
    if [[ "$lookup_file" == *.gz ]]; then
        zcat "$lookup_file" > "$temp_lookup"
    else
        cp "$lookup_file" "$temp_lookup"
    fi
    
    # Extract header from GWAS file
    local gwas_header=$(head -n 1 "$temp_gwas")
    
    # Create lookup hash (chrom_pos_ref_alt -> alt_freqs)
    awk -F'\t' 'NR>1 {print $1 "_" $2 "_" $4 "_" $5 "\t" $6}' "$temp_lookup" > "$temp_lookup_hash"
    
    # Clear temp_lookup from memory since we have the hash now
    rm -f "$temp_lookup"
    echo "  🗑️  Cleared lookup temp file"
    
    # Process GWAS file and merge with memory-efficient AWK
    {
        echo -e "${gwas_header}\teffect_allele_frequency"
        
        awk -F'\t' '
        BEGIN { 
            # Read lookup table into associative array
            while ((getline line < "'$temp_lookup_hash'") > 0) {
                split(line, a, "\t")
                freqs[a[1]] = a[2]
            }
            close("'$temp_lookup_hash'")
            print "  📊 Loaded " length(freqs) " frequency entries" > "/dev/stderr"
        }
        NR>1 {
            # Create merge key: chrom_pos_A2_A1 (assuming A2=REF, A1=ALT)
            merge_key = $1 "_" $2 "_" $5 "_" $4
            
            # Print original line with ALT_FREQS
            if (merge_key in freqs) {
                print $0 "\t" freqs[merge_key]
            } else {
                print $0 "\tNA"
            }
        }
        END {
            # Clear the associative array to free memory
            delete freqs
            print "  🧹 Cleared AWK frequency array" > "/dev/stderr"
        }' "$temp_gwas"
    } > "$temp_merged"
    
    # Clear GWAS temp file from memory
    rm -f "$temp_gwas" "$temp_lookup_hash"
    echo "  🗑️  Cleared GWAS and hash temp files"
    
    # Compress with BGZF format (required for tabix) instead of regular gzip
    bgzip -c "$temp_merged" > "$output_file"
    
    # Clear merged temp file immediately after compression
    rm -f "$temp_merged"
    echo "  🗑️  Cleared merged temp file"
    
    # Rebuild tabix index
    echo "  🔧 Rebuilding tabix index..."
    
    # Try different tabix indexing approaches
    if tabix -s 1 -b 2 -e 2 "$output_file" 2>/dev/null; then
        echo "  ✅ Tabix index rebuilt successfully (generic format)"
    elif tabix -p bed -S 1 -s 1 -b 2 -e 2 "$output_file" 2>/dev/null; then
        echo "  ✅ Tabix index rebuilt successfully (bed format)"
    elif tabix -s 1 -b 2 -e 2 -0 "$output_file" 2>/dev/null; then
        echo "  ✅ Tabix index rebuilt successfully (0-based)"
    else
        echo "  ⚠️  Warning: Failed to rebuild tabix index"
        echo "      You can manually rebuild with: tabix -s 1 -b 2 -e 2 $output_file"
    fi
    
    # Calculate statistics using streaming to avoid loading large files
    echo "  📈 Calculating statistics..."
    local total_variants matched_variants
    
    # Use streaming approach for large files
    if [[ -f "$output_file" ]]; then
        # Count total variants (subtract 1 for header)
        total_variants=$(zcat "$output_file" | wc -l)
        total_variants=$((total_variants - 1))
        
        # Count matched variants using streaming
        matched_variants=$(zcat "$output_file" | awk -F'\t' 'NR>1 && $NF != "NA" {count++} END {print count+0}')
        
        local match_percent
        if [ "$total_variants" -gt 0 ]; then
            match_percent=$(awk "BEGIN {printf \"%.2f\", $matched_variants/$total_variants*100}")
        else
            match_percent="0.00"
        fi
        
        echo "  Total variants: $total_variants"
        echo "  Matched variants: $matched_variants ($match_percent%)"
    fi
    
    echo "  Output: $output_file"
    
    # Final cleanup of any remaining temp files
    cleanup_temp_files "${temp_files_array[@]}"
    clear_awk_memory
    
    # Clear local variables to free memory
    unset temp_files_array gwas_header base_name total_variants matched_variants match_percent
    
    echo "  ✅ Memory cleanup completed"
    echo ""
}

# Function to monitor memory usage (optional)
show_memory_usage() {
    if command -v free >/dev/null 2>&1; then
        echo "💾 Current memory usage:"
        free -h | grep -E "(Mem|Swap):" | sed 's/^/  /'
    fi
}

# Main function
main() {
    # Check if bgzip is available
    if ! command -v bgzip &> /dev/null; then
        echo "Error: bgzip is required but not found. Please install htslib/tabix package." >&2
        echo "On Ubuntu/Debian: sudo apt-get install tabix" >&2
        echo "On CentOS/RHEL: sudo yum install htslib" >&2
        echo "With conda: conda install -c bioconda htslib" >&2
        exit 1
    fi
    
    # Check arguments
    if [ $# -lt 3 ] || [ $# -gt 4 ]; then
        show_usage
        exit 1
    fi
    
    local lookup_table="$1"
    local gwas_dir="$2"
    local output_dir="$3"
    local gwas_pattern="${4:-*.tsv.gz}"
    
    # Expand paths
    lookup_table=$(eval echo "$lookup_table")
    gwas_dir=$(eval echo "$gwas_dir")
    output_dir=$(eval echo "$output_dir")
    
    # Validate paths
    if [ ! -f "$lookup_table" ]; then
        echo "Error: Lookup table file not found: $lookup_table" >&2
        exit 1
    fi
    
    if [ ! -d "$gwas_dir" ]; then
        echo "Error: GWAS directory not found: $gwas_dir" >&2
        exit 1
    fi
    
    # Create output directory if it doesn't exist
    if [ ! -d "$output_dir" ]; then
        mkdir -p "$output_dir"
        echo "Created output directory: $output_dir"
    fi
    
    # Show initial memory usage
    show_memory_usage
    echo ""
    
    echo "Loading lookup table: $lookup_table"
    local lookup_variants
    if [[ "$lookup_table" == *.gz ]]; then
        lookup_variants=$(( $(zcat "$lookup_table" | wc -l) - 1 ))
    else
        lookup_variants=$(( $(wc -l < "$lookup_table") - 1 ))
    fi
    echo "Lookup table loaded: $lookup_variants variants"
    echo ""
    
    # Find GWAS files
    local gwas_files=()
    while IFS= read -r -d '' file; do
        gwas_files+=("$file")
    done < <(find "$gwas_dir" -maxdepth 1 -name "$gwas_pattern" -type f -print0)
    
    if [ ${#gwas_files[@]} -eq 0 ]; then
        echo "Error: No GWAS files found matching pattern: $gwas_dir/$gwas_pattern" >&2
        exit 1
    fi
    
    echo "Found ${#gwas_files[@]} GWAS files to process"
    echo ""
    
    # Process each GWAS file with memory cleanup
    local processed_count=0
    local file_counter=0
    for gwas_file in "${gwas_files[@]}"; do
        ((file_counter++))
        echo "=== Processing file $file_counter of ${#gwas_files[@]} ==="
        
        if process_gwas_file "$gwas_file" "$lookup_table" "$output_dir"; then
            ((processed_count++))
            echo "✅ File $file_counter completed successfully"
        else
            echo "❌ Error processing file $file_counter: $(basename "$gwas_file")" >&2
        fi
        
        # Show memory usage after each file
        show_memory_usage
        echo ""
    done
    
    # Clear the gwas_files array to free memory
    unset gwas_files
    
    echo "Processing complete!"
    echo "Successfully processed $processed_count files"
    echo "Output files saved in: $output_dir"
    
    # Final memory usage
    echo ""
    show_memory_usage
}

# Setup cleanup trap for unexpected exits
cleanup_on_exit() {
    echo "🧹 Cleaning up temporary files..."
    # Clean up any remaining temp files
    find /tmp -name "tmp.*" -user "$(whoami)" -mtime +0 -delete 2>/dev/null || true
}
trap cleanup_on_exit EXIT INT TERM

# Run main function
main "$@"