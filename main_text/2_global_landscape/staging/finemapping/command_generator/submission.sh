username=rl3328
jobname=rss_image_test
jobman="$HOME/Downloads/mmcloud/src/mm_batch.sh"
bash $jobman \
 --job-script  ./sos_submit_QC_rss_qc_finemap_susie_rss.sh \
 -c 2:32 -m 16:128 -jn $jobname \
 --job-size 87 \
 --parallel-commands 0 \
 --min-cores-per-command 1 \
 --min-mem-per-command 7 \
 --mount statfungen/ftp_fgc_xqtl:/home/$username/data statfungen/ftp_fgc_xqtl/analysis_result/image_MS_finemapping/:/home/$username/output \
 --mountOpt "mode=r" "mode=rw" \
 --cwd "/home/$username/data" \
 --ebs-mount "/home/$username/input=3" \
 --no-fail-fast --opcenter 3.82.198.55 