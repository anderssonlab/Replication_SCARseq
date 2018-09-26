#### Process SCAR-seq data ####
#### Commands (Maria Dalby 0310-2018) ####

## QC and mapping ##
# Install Kundaje pipeline and genome data from https://github.com/kundajelab/chipseq_pipeline#pipeline
# python chipseq.py [CHIPSEQ_TYPE] [FINAL_StegE] [CONF_JSON_FILE] -species [SPECIES] -nth [NUM_THREADS] -fastq1 ... -fastq2 ... -ctl_fastq1 ...

## Run pipeline
nice nohup chipseq.py histone files.json  -species mm10 -final_stage filt_bam > chipseq_files.out 2>&1 &
disown # important to add after command

# Split into forward and reverse from SE data
for file in K36_NS_r1_chr_filt.nodup.bam; do
    base=$( basename $file | sed -e "s/.nodup.bam//" )
    samtools view -F 20 -h ${file} | samtools view -Sb -h > ${base}_F.nodup.bam
    samtools view -f 16 -h ${file} | samtools view -Sb -h > ${base}_R.nodup.bam
done

# Split into forward and reverse from PE data
#for file in *.PE2SE.nodup.bam; do
#    base=$( basename $file | sed -e "s/.nodup.bam//" )
#    samtools view -f 99 -h ${file} | samtools view -Sb -h > ${base}_99_F.bam
#    samtools view -f 147 -h ${file} | samtools view -Sb -h > ${base}_147_F.bam
#    samtools merge -f ${base}_F.nodup.bam ${base}_99_F.bam ${base}_147_F.bam
#    samtools index ${base}_F.nodup.bam
#
#    samtools view -f 83 -h ${file} | samtools view -Sb -h > ${base}_83_R.bam
#    samtools view -f 163 -h ${file} | samtools view -Sb -h > ${base}_163_R.bam
#    samtools merge -f ${base}_R.nodup.bam ${base}_83_R.bam ${base}_163_R.bam
#    samtools index ${base}_R.nodup.bam
#done

# Convert bam to bg + bw files
ls *nodup.bam | parallel --gnu -j 10 "genomeCoverageBed -ibam {} -bg -g mm10.chrom.sizes" > bw_files/$(basename {} .bam).bdg
ls *.bdg | parallel --gnu -j 10 "bdg2bw {} /isdata/alab/people/maria/genome/mm10/mm10.chrom.sizes {}.bw"

# Calculate histone partition from SCAR-seq data
# nice nohup ./partition_smooth_find_breakpoints.sh ${DIR}/${FILE}_F.nodup.bw ${DIR}/${FILE}_R.nodup.bw ${WIN_size} ${RADIUS_BIN} ${DRADIUS_BIN} ${ZRADIUS} ${OUT_DIR} ${CHROM} &
nice nohup ./partition_smooth_find_breakpoints.sh ${DIR}/${FILE}_F.nodup.bw ${DIR}/${FILE}_R.nodup.bw 1000 30 30 1 res_${FILE} mm10.chrom.sizes &

# Calculate RFD from OK-seq data
# nice nohup ./RFD_smooth_find_breakpoints.sh ${FILE}_F.nodup.bw ${FILE}_R.nodup.bw ${WIN_size} ${RADIUS_BIN} ${DRADIUS_BIN} ${ZRADIUS} ${OUT_DIR} ${CHROM} &
nice nohup ./RFD_smooth_find_breakpoints.sh ${FILE}_F.nodup.bam.bw ${FILE}_R.nodup.bam.bw 1000 30 30 1 res_OKseq mm10.chrom.sizes &


# raw bam files correlations 
nice nohup multiBamSummary bins -b K36_313_r1.nodup.bam K36_313_r2.nodup.bam K36me3_316_r1.nodup.bam K36me3_316_r2.nodup.bam K36_NS_r1.nodup.bam K36_r3_merge.nodup.bam K36_r4_merge.nodup.bam -out bam_cor_K36_bs1kb.out -bs 1000 &
nice nohup multiBamSummary BED-file --BED gencode.vM17.annotation.bed -b K36_313_r1.nodup.bam K36_313_r2.nodup.bam K36me3_316_r1.nodup.bam K36me3_316_r2.nodup.bam K36_NS_r1.nodup.bam K36_r3_merge.nodup.bam K36_r4_merge.nodup.bam -out bam_cor_K36_genes.out &

plotCorrelation -in bam_cor_K36_bs1kb.out -o  bam_cor_K36_bs1kb_pearson.pdf -c pearson -p heatmap --outFileCorMatrix bam_pearson_cor_K36_bs1kb.matrix.txt
plotCorrelation -in bam_cor_K36_genes.out -o  bam_cor_K36_genes_pearson.pdf -c pearson -p heatmap --outFileCorMatrix bam_pearson_cor_K36_genes.matrix.txt
