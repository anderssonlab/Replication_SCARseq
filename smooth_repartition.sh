## ./find_breakpoints_smooth.sh ../Okazaki_data/Okazaki_mm10_r1_F.nodup.bw ../Okazaki_data/Okazaki_mm10_r1_R.nodup.bw 1000 20 20 results mm10.chrom.sizes
##Usage: find_breakpoints_smooth.sh F.bw R.bw window_size smooth_radius derivative_radius zero_crossing_radius outdir chrom_sizes

FWD=$1
REV=$2
WIN=$3
RADIUS=$4
DRADIUS=$5
ZRADIUS=$6
OUT=$7
CHROM=$8


mkdir -p results
mkdir -p ${OUT}

## Make windows (step the same distance as window sizes)
bedtools makewindows -i srcwinnum -g $CHROM -w $WIN -s $WIN > ${OUT}/windows.bed

## Split windows on chromosome
for CHR in $( cut -f 1 $CHROM | sort | uniq ); do
	awk -v c=$CHR '$1==c' ${OUT}/windows.bed > ${OUT}/${CHR}_windows.bed &
done
wait

## Aggregate counts in windows
for CHR in $( cut -f 1 $CHROM | sort | uniq ); do
	bigWigAverageOverBed ${FWD} ${OUT}/${CHR}_windows.bed ${OUT}/${CHR}_windows_F.tab &
	bigWigAverageOverBed ${REV} ${OUT}/${CHR}_windows.bed ${OUT}/${CHR}_windows_R.tab &
done
wait

## CPM normalise strands individually (with pseudocount 1)
CPM=`cat ${OUT}/*_windows_F.tab cat ${OUT}/*_windows_R.tab | awk '{SUM += ($4+1)} END {print SUM/10^6}'`

echo "CPM factors: $CPM"
for CHR in $( cut -f 1 $CHROM | sort | uniq ); do
	awk -v c=$CPM 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${OUT}/${CHR}_windows_F.tab > ${OUT}/${CHR}_windows_F_CPM.tab &
	awk -v c=$CPM 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${OUT}/${CHR}_windows_R.tab > ${OUT}/${CHR}_windows_R_CPM.tab &
done
wait

## Calculate RFD, derivative and boundary scores
for CHR in $( cut -f 1 $CHROM | sort | uniq ); do
	./RFD_smooth.repart.pl ${OUT}/${CHR}_windows_F_CPM.tab ${OUT}/${CHR}_windows_R_CPM.tab $RADIUS $DRADIUS $ZRADIUS > ${OUT}/${CHR}_windows_RFD.txt &
done
wait

## Gather output
if [ -f ${OUT}/${OUT}_smooth_results_w${WIN}_s${RADIUS}_d${DRADIUS}_z${ZRADIUS}.txt ]; then rm ${OUT}/${OUT}_smooth_results_w${WIN}_s${RADIUS}_d${DRADIUS}_z${ZRADIUS}.txt; fi
touch ${OUT}/${OUT}_smooth_results_w${WIN}_s${RADIUS}_d${DRADIUS}_z${ZRADIUS}.txt
for CHR in $( cut -f 1 $CHROM | sort | uniq ); do
	paste ${OUT}/${CHR}_windows.bed ${OUT}/${CHR}_windows_F.tab ${OUT}/${CHR}_windows_R.tab ${OUT}/${CHR}_windows_F_CPM.tab ${OUT}/${CHR}_windows_R_CPM.tab ${OUT}/${CHR}_windows_RFD.txt | awk '$8>0 || $14>0' | cut -f -3,8,14,20,26,30- >> ${OUT}/${OUT}_smooth_results_w${WIN}_s${RADIUS}_d${DRADIUS}_z${ZRADIUS}.txt
done

rm ${OUT}/*_windows_RFD.txt
rm ${OUT}/*_windows_R.tab
rm ${OUT}/*_windows_F.tab
rm ${OUT}/*_windows_R_CPM.tab
rm ${OUT}/*_windows_F_CPM.tab
rm ${OUT}/*_windows.bed
rm ${OUT}/windows.bed
