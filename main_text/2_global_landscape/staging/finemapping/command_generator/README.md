# Commands for genome-wide scan using RSS analysis pipeline for the AD image and MS GWAS

## Aim: This document recorded the procedure for performing genome-wide scans using the RSS analysis pipeline, specifically for Alzheimer's Disease (AD) imaging and Multiple Sclerosis (MS) GWAS data(2.8.2015). For comprehensive instructions, please refer to the [pipeline documentation](https://github.com/cumc/xqtl-analysis/blob/main/analysis/Wang_Columbia/GWAS/AD_GWAS_rss.ipynb).

## Input:
1. generator.sh: the command generator, which can add --region-name for each sos script and customize your finemapping method and qc method.
2. sos_template_rss_gwas.sh: the template whose path is used in the generator.sh.
3. submission.sh: the submission commands for mmjobman.
## Output:
sos_submit_QC_rss_qc_finemap_susie_rss.sh: a file containing multiple sos scipts that was produced by `bash generator.sh` in your terminal.

### Steps
1. customize the parameter and path in the generator.sh, I used "susie_rss" finemapping method, "rss_qc" qc method and RAISS imputed. I modified the `sos_template` and `sos_submit`.

**Note:**
**1. This RSS analysis is conducted by LD_block, so the `--region-name` for each sos scrip is from `s3://statfungen/ftp_fgc_xqtl/resource/ADSP_R4_EUR_LD/ld_meta_file.tsv`.**
**2. Use `--oem-packages` for batch jobs if you want to use packages in the shared folder. This is also the default option.**
