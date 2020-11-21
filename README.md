# Replication_SCARseq

Bash, perl and R commands for:
"MCM2 promotes symmetric inheritance of modified histones during DNA replication" Petryk N. et. al. 2018. Science.

A few data processing steps are listed here: SCARseq_commands.sh

For OK-seq data,run RFD_smooth_find_breakpoints.sh with RFD_smooth.pl stored in the same folder.
For all SCAR-seq data, run partition_smooth_find_breakpoints.sh	with partition_smooth.pl stored in the same folder.

Input signal files to both RFD_smooth_find_breakpoints.sh and partition_smooth_find_breakpoints.sh should be strand seperated bigWig files.

Finally, the replication_SCARseq.R file stores the commands used to process the binned signal, including most analysis and plots for the paper. Worth noticing is that the smoothed SCAR-seq partition or smoothed OK-seq RFD values are all referred to as a common "RFD" collumn. Raw (non-smoothed) SCAR-seq partition or OK-seq RFD values are referred to as "RFD.raw"

Contact: Maria Dalby (mdalbydk@gmail.com) or Robin Andersson (Robin@bio.ku.dk)
