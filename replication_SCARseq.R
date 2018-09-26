## Replication and Chromatin RFD

require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(rtracklayer)
require(ggplot2)
require(reshape2)
require(ggthemes)
library(GenomicAlignments)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(stringr)
library(ChIPseeker)
library(stringi)
library(clusterProfiler)
library(ReactomePA)

server = ""
data_dir <- "/binf-isilon/alab/projects/rep_chromatin"
setwd(data_dir)

my.pal <- c(brewer.pal(6,"RdBu"),brewer.pal(6,"BrBG"),brewer.pal(6,"PRGn"))
myseqnames <- c(paste("chr",1:19,sep=""))

#load("Rfiles/data.gr_current_130418.RDat") # Initial submission
#load("dataGR_fult_K5incl_wInput.RDat") # Final submission

#### load data ####
Ok <- read.table(file.path(server,data_dir,"/RFD/results/Okazaki_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE) # Output from RFD_smooth_find_breakpoints.sh

## Histone/input SCAR-seq binned partition output from partition_smooth_find_breakpoints.sh
## input data
input_r1 <- read.table(file.path(server,data_dir, "RFD/results/res_input_r1_teg_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
input_r2 <- read.table(file.path(server,data_dir, "RFD/results/res_input_NS_r2_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
input_r3 <- read.table(file.path(server,data_dir, "RFD/results/res_input_NS_r3_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)

input_316_r1 <- read.table(file.path(server,data_dir, "RFD/results/res_input_316_r2_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
input_313_r1 <- read.table(file.path(server,data_dir, "RFD/results/res_input_313_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
input_m_r1 <- read.table(file.path(server,data_dir, "RFD/results/res_input_MS_r1_teg_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)

input_K5_316_r1 <- read.table(file.path(server,data_dir, "RFD/results/res_input_K5_316_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
input_K5_313_r1 <- read.table(file.path(server,data_dir, "RFD/results/res_input_K5_313_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)

## K5 data
K5_3 <- read.table( file.path(server,data_dir,"/RFD/results/res_K5_r3_teg_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K5_4 <- read.table( file.path(server,data_dir,"/RFD/results/res_K5_r4_teg_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)

# K5 mutant data - 316: mutant AA#1 (double KD), 313: mutant AA#2 (single KD)
K5m316_1 <- read.table(file.path(server,data_dir,"/RFD/results/res_K5_316_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K5m316_2 <- read.table(file.path(server,data_dir,"/RFD/results/res_K5_316_r2_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K5m313_1 <- read.table(file.path(server,data_dir,"/RFD/results/res_K5_313_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K5m313_2 <- read.table(file.path(server,data_dir,"/RFD/results/res_K5_313_r2_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)


## K20 data
K20_1 <- read.table(file.path(server,data_dir,"/RFD/results/res_K20_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K20_2 <- read.table(file.path(server,data_dir,"/RFD/results/res_K20_r2_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K20_3 <- read.table(file.path(server,data_dir,"/RFD/results/res_K20_r3_teg_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)

# K20 mutant data - 316: mutant AA#1, 313: mutant AA#2
K20m316_1 <- read.table(file.path(server,data_dir,"/RFD/results/res_K20_316_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K20m316_2 <- read.table(file.path(server,data_dir,"/RFD/results/res_K20_316_r2_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K20m313_1 <- read.table(file.path(server,data_dir,"/RFD/results/res_K20_313_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K20m313_2 <- read.table(file.path(server,data_dir,"/RFD/results/res_K20_313_r2_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)


## K36 data
K36_1 <- read.table(file.path(server,data_dir,"/RFD/results/res_K36_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K36_3 <- read.table(file.path(server,data_dir,"/RFD/results/res_K36_r3_merge_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K36_4 <- read.table(file.path(server,data_dir,"/RFD/results/res_K36_r4_merge_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)

# K36 mutant data - 316: mutant AA#1, 313: mutant AA#2
K36m316_1 <- read.table(file.path(server,data_dir,"/RFD/results/res_K36_316_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K36m316_2 <- read.table(file.path(server,data_dir,"/RFD/results/res_K36_316_r2_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K36m313_1 <- read.table(file.path(server,data_dir,"/RFD/results/res_K36_313_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K36m313_2 <- read.table(file.path(server,data_dir,"/RFD/results/res_K36_313_r2_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)

## parental strand
K5_par <- read.table( file.path(server,data_dir,"/RFD/results/res_K5_NS_par_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
K20_par <- read.table( file.path(server,data_dir,"/RFD/results/res_K20_NS_par_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)

## expression data
# CAGE mm10 ES cells (pooled signal)
#cage <- read.table(file.path(server,data_dir,"/RFD/results/res_mESCAGE_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
pro1 <- read.table(file.path(server,data_dir,"/RFD/results/res_proseq_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
pro2 <- read.table(file.path(server,data_dir,"/RFD/results/res_proseq_r2_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
#rna1 <- read.table(file.path(server,data_dir,"/RFD/results/res_rnaseq_r1_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)
#rna2 <- read.table(file.path(server,data_dir,"/RFD/results/res_rnaseq_r2_smooth_results_w1000_s30_d30_z1.txt"),sep="\t",header=FALSE,as.is=TRUE)

# store all data into data_df
data_df <- rbind(data.frame(Ok,type=rep("Ok",times=nrow(Ok))),
              data.frame(K5_3,type=rep("K5_3",times=nrow(K5_3))),
              data.frame(K5_4,type=rep("K5_4",times=nrow(K5_4))),
              data.frame(K5_par,type=rep("K5par",times=nrow(K5_par))),
              data.frame(K20m316_1,type=rep("K20m16_1",times=nrow(K20m316_1))),
              data.frame(K20m316_2,type=rep("K20m16_2",times=nrow(K20m316_2))),
              data.frame(K36m316_1,type=rep("K36m16_1",times=nrow(K36m316_1))),
              data.frame(K36m316_2,type=rep("K36m16_2",times=nrow(K36m316_2))),
              data.frame(K20m313_1,type=rep("K20m13_1",times=nrow(K20m313_1))),
              data.frame(K20m313_2,type=rep("K20m13_2",times=nrow(K20m313_2))),
              data.frame(K36m313_1,type=rep("K36m13_1",times=nrow(K36m313_1))),
              data.frame(K36m313_2,type=rep("K36m13_2",times=nrow(K36m313_2))),
              data.frame(K5m316_1,type=rep("K5m16_1",times=nrow(K5m316_1))),
              data.frame(K5m316_2,type=rep("K5m16_2",times=nrow(K5m316_2))),
              data.frame(K5m313_1,type=rep("K5m13_1",times=nrow(K5m313_1))),
              data.frame(K5m313_2,type=rep("K5m13_2",times=nrow(K5m313_2))),
              data.frame(K20_1,type=rep("K20_1",times=nrow(K20_1))),
              data.frame(K20_2,type=rep("K20_2",times=nrow(K20_2))),
              data.frame(K20_3,type=rep("K20_3",times=nrow(K20_3))),
              data.frame(K20_par,type=rep("K20par",times=nrow(K20_par))),
              data.frame(K36_1,type=rep("K36_1",times=nrow(K36_1))),
              data.frame(K36_3,type=rep("K36_3",times=nrow(K36_3))),
              data.frame(K36_4,type=rep("K36_4",times=nrow(K36_4))),
              data.frame(pro1,type=rep("proseq_r1",times=nrow(pro1))),
              data.frame(pro2,type=rep("proseq_r2",times=nrow(pro2))))
colnames(data_df) <- c("chr","start","end","F","R","F.cpm","R.cpm","RFD.raw","RFD","RFD.deriv","score","zero.deriv","type")

# keep only standart chromosomes
data_df <- subset(data_df, chr %in% myseqnames)

data_df$mid <- rowMeans(data_df[,c("start","end")]) # genomic midpoint
data_df$ID <- paste(data_df$chr,data_df$start,sep=":") %>% paste(.,data_df$end,sep="-") # genomic identifier
data_df$mark <- sapply(str_split(as.character(data_df$type),"_"),"[",1) # histone mod.
data_df$exprs <- rowSums(data_df[,c("F","R")]) # total raw counts
data_df$CPM <- data_df$F.cpm + data_df$R.cpm # total counts pr. million

# create GRange object
data.gr <- makeGRangesFromDataFrame(data,
                                    seqnames.field = "chr",
                                    start.field = "start",
                                    end.field = "end",
                                    keep.extra.columns = T)

# filter out ENCODE blacklist regions (mm10 https://sites.google.com/site/anshulkundaje/projects/blacklists)
neg_list <- read.table(file.path(server,"/binf-isilon/alab/people/maria/genome/mm10/mm10.blacklist.bed"))
neg_gr <-  makeGRangesFromDataFrame(df=neg_list,
                                    seqnames.field="V1",
                                    start.field="V2",
                                    end.field="V3")
neg.overlap <- findOverlaps(data.gr,neg_gr)
data.gr <- data.gr[-queryHits(neg.overlap)]

# filter based on coverage
data.gr <- subset(data.gr,CPM>=0.3) # & exprs>=500
## loaded file

#### define Ok breakpoints ####
data.Ok <- read.table("RFD/results/Okazaki_smooth_results_w1000_s30_d30_z1.txt",header=FALSE,as.is=TRUE,sep="\t",fill=TRUE) # output from RFD_smooth_find_breakpoints.sh
colnames(data.Ok) <- c("chr","start","end","F","R","F.cpm","R.cpm","RFD.raw","RFD","RFD.deriv","score","zero.deriv")
data.Ok <- subset(data.Ok, chr %in% myseqnames)
data.Ok$exprs <- rowSums(data.Ok[,c(4,5)])
data.Ok$CPM <- rowSums(data.Ok[,c(6,7)])
data.Ok$ID <- paste(data.Ok$chr,data.Ok$start,sep=":") %>% paste(.,data.Ok$end,sep="-") # remove again
data.Ok <- subset(data.Ok, CPM >= 0.3 ) # & exprs>500 # filter similar to SCAR-seq data

# identify regions with RFD derivative (only zero crossing) higher than 90% quantile of RFD derivatives
ids <- which(data.Ok[,12]>quantile(data.Ok[,10],probs=0.9,na.rm=TRUE))

# print regions and merge proximal regions using bedtools. IDs will be collapsed
write.table(cbind(data.Ok[ids,1:3],ids),"boundaries_cpm.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("bedtools merge -c 4 -d 3000 -o collapse -i boundaries_cpm.bed > boundaries_3000merged_cpm.bed")
b.data <- read.table("boundaries_3000merged_cpm.bed",sep="\t",as.is=TRUE,header=FALSE)

# For each collapsed region, select the bin with the highest derivative
filtered.data <- data.Ok[sapply(b.data[,4], function(x) {
  y <- as.numeric(strsplit(x,",")[[1]])
  y[which.max(data.Ok[y,10])]
}),]

# write OK-seq initiation zone file
#write.table(filtered.data,file="RFD/filtered_s30_3000boundaries0.9_cpm0.3.txt",quote = F,sep="\t",col.names = F,row.names = F)


#### signal around Ok-seq initiation zones ####
ok_break <- read.table(file.path(server,data_dir,"RFD/filtered_s30_3000boundaries0.9_cpm0.3.txt"),sep="\t",header=FALSE,as.is=TRUE) 
colnames(ok_break) <- c("chr","start","end","F","R","F.cpm","R.cpm","RFD.raw","RFD","RFD.deriv","score","zero.deriv","CPM")
ok_break$mid <- rowMeans(ok_break[,c("start","end")])
ok_break$ID <- paste(ok_break$chr,ok_break$start,sep=":") %>% paste(.,ok_break$end,sep="-")

ok.gr <- makeGRangesFromDataFrame(ok_break,
                                  seqnames.field="chr",
                                  start.field = "start",
                                  end.field = "end",
                                  keep.extra.columns = T)

# Add initiation zone label to full data table
data.gr$IZ <- ifelse(data.gr$ID %in% ok.gr$ID, TRUE,FALSE)

# Rank initiation zones for figure S3F, S7C
#ok.gr$rank <- findInterval(ok.gr$score,median(ok.gr$score,na.rm=T))

# add random initiation zone
#ranID <- sample(data.gr$ID[data.gr$mark=="Ok" & !data.gr$ID %in% ok.gr$ID],length(ok.gr))
#data.gr$RanIZ <- ifelse(data.gr$ID %in% ranID, TRUE,FALSE)

# signal in proximity to ok-breakpoints
#ran.gr <- sample(data.gr[data.gr$mark=="Ok" & data.gr$IZ == FALSE],length(ok.gr)) # DELETE?
ok.ext.gr <- resize(ok.gr,200000,fix="center") # figure 1C, 2D, 3A, S3A, S3C-F, S5E, S7A-C, S8G
ok.ext.gr$break_start <- start(ok.gr)

## signal in proximity to enhancer midpoints (Specific to figure 2C, S8G)
#ok.ext.gr <- resize(Enh_SE.gr,200000,fix="center") # center on enhancer midpoints
#ok.ext.gr$break_start <- mid(ranges(Enh_SE.gr))

# exclude overlapping initiation zone (within 200000 bp).
ok.dist <- distanceToNearest(ok.ext.gr)
ok.ext.gr <- ok.ext.gr[-queryHits(subset(ok.dist,ok.dist@elementMetadata$distance==0))]

# subset to data within search-space (ok.ext.gr)
overlap_pairs <- findOverlaps(ok.ext.gr, data.gr)
data.break.gr <- data.gr[subjectHits(overlap_pairs)]
data.break.gr$break_ID <- ok.ext.gr$ID[queryHits(overlap_pairs)] # initiation zone ID
#data.break.gr$rank <- ok.ext.gr$rank[queryHits(overlap_pairs)] # rank ID

# for enhancer midpoint plots
#data.break.gr$break_ID <- ok.ext.gr$No[queryHits(overlap_pairs)] # enhancer ID
#data.break.gr$enh_active <- ok.ext.gr$active[queryHits(overlap_pairs)]

data.break.gr$dist <- start(data.break.gr) - ok.ext.gr$break_start[queryHits(overlap_pairs)]
#data.break.gr$dist <- round_any(data.break.gr$dist,1000) # for enhancer midpoints

data.mean.df <- data.frame(data.break.gr %>% data.frame %>% 
                             dplyr::group_by(dist,mark,type) %>%  # rank, enh_active
                             dplyr::summarise(#RFD_sd = sd(RFD,na.rm = T),
                                  RFD.raw = mean(RFD.raw, na.rm = T),
                                  RFD = mean(RFD,na.rm = T)))

data.mean.df$dist <- data.mean.df$dist / 1000 # convert to kb

ggplot(data=subset(data.mean.df,!is.na(mut))) + 
  geom_line(aes(colour=factor(rank),dist,RFD)) + #alpha=rep
  #geom_smooth(aes(colour=mark),size=0.4) + #alpha=rep
  facet_wrap(~ mark + PTM,scale="free_y",ncol=3) +
  theme_minimal() +
  scale_alpha_manual(values = c(1,0.7,0.5)) +
  scale_colour_manual(values=my.pal[c(5,6,1,2)]) +
  #scale_colour_brewer(palette = my.cols) + #
  scale_linetype_manual(values=c(3,2,1)) +
  #coord_cartesian(ylim =  c(-0.003,0.003))  +
  xlab("Distance (kb) from IZ")
pdf("SCAR_RDF_100kb_enhancer_pr.rep_wOverlap.pdf") # ,height = 7,width = 5
print(plot)
dev.off()


#### identify initiation zone edges ####
# based on min/max OK-seq data
data.break.ok.df <- data.frame(data.break.gr[data.break.df$type=="Ok",]) # subset to OKseq
data.ok_F <- data.break.ok.df[data.break.ok.df$dist>0,] # downstream
data.ok_R <- data.break.ok.df[data.break.ok.df$dist<0,] # upstream

# highest RFD downstream of initiation zone
ok_max <- data.frame(data.ok_F %>% 
                       dplyr::group_by(break_ID,enh_active) %>%
                       dplyr::summarise(ID = ID[which.max(RFD)],
                                        RFDedge = max(RFD),
                                        dist = dist[which.max(RFD)]))
# lowest RFD upstream of initiation zone
ok_min <- data.frame(data.ok_R %>% 
                       dplyr::group_by(break_ID,enh_active) %>%
                       dplyr::summarise(ID = ID[which.min(RFD)],
                                        RFDedge = min(RFD),
                                        dist = dist[which.min(RFD)]))

ok_extreme <- rbind(ok_max,ok_min)
data.extreme <- data.break.df[which(data.break.df$ID %in% ok_extreme$ID),] # all data at initiation zone edges (two bins pr. initiation zone)

#### statistical test ####
data.break.df <- data.frame(data.break.gr)
data.extreme.min <- data.break.df[which(data.break.df$ID %in% ok_min$ID),]
data.extreme.max <- data.break.df[which(data.break.df$ID %in% ok_max$ID),]

RFD.mat.min <- dcast(data.extreme.min, break_ID ~ type, value.var = "RFD", fun.aggregate = mean)
RFD.mat.max <- dcast(data.extreme.max, break_ID ~ type, value.var = "RFD", fun.aggregate = mean)

# define samples to test
RFD.mat.min <- RFD.mat.min[,colnames(RFD.mat.min) %in% c("break_ID","K5_3","K5_4","K20_1","K20_2","K20_3","K36_1","K36_3","K36_4","K36m13_1","K36m16_1","K36m13_2","K36m16_2","K20m13_1","K20m16_1","K20m13_2","K20m16_2","K5m13_1","K5m16_1","K5m13_2","K5m16_2","Ok")]
RFD.mat.max <- RFD.mat.max[,colnames(RFD.mat.max) %in% c("break_ID","K5_3","K5_4","K20_1","K20_2","K20_3","K36_1","K36_3","K36_4","K36m13_1","K36m16_1","K36m13_2","K36m16_2","K20m13_1","K20m16_1","K20m13_2","K20m16_2","K5m13_1","K5m16_1","K5m13_2","K5m16_2","Ok")]

# perform a Wilcoxon Rank Sum test in samples partitions
w.pval <- list() 
for (nr_mark in 3:ncol(RFD.mat.min)) { # 
  w.pval[colnames(RFD.mat.min)[nr_mark]] <- wilcox.test(RFD.mat.min[,nr_mark],RFD.mat.max[,nr_mark],paired = T)$p.value #,alternative = c("g"))
}
w.test.df <- data.frame(pvalue=unlist(w.pval))
w.test.df$type <- rownames(w.test.df)
w.test.df$mark <- sapply(strsplit(as.character(w.test.df$type),"_"),"[",1)

# bar plot on p-values
ggplot(subset(w.test.df,!mark %in% c("Ok","proseq","input","K5par","K20par")), aes(type,-log10(pvalue))) +
  geom_bar(stat="identity", position="dodge",aes(fill=mark),alpha=0.8) +
  theme_minimal() + 
  geom_hline(yintercept = c(-log10(0.01),-log10(0.05)),colour="red",size=0.2) +
  scale_fill_manual(values = my.pal) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) 
ggsave("P.values.w.test_min_vs_max_smoothRFD.pdf")

# boxplot on smoothed partition at edges (figure 1D, S3C, S3E, S5E, S7B)
RFD.mat.min$stand_rep <- "lagging"
RFD.mat.max$stand_rep <- "leading"
RFD.mat <- rbind(RFD.mat.min, RFD.mat.max)

mat_m <- melt(RFD.mat)
mat_m$stand_rep <- factor(mat_m$stand_rep,levels=c("leading","lagging"))
ggplot(data=subset(mat_m,!variable %in% c("K5_3","K5_4","K20_1","K20_2","K20_3","K36_1","K36_3","K36_4")), aes(variable,value)) + # ,"K36_1","K36_3","K36_4"
  #geom_violin(aes(fill=stand_rep),alpha=0.8,outlier.shape = NA) + #outlier.shape = NA
  geom_boxplot(aes(fill=stand_rep),alpha=0.8,outlier.shape = NA,notch = 1) + #outlier.shape = NA
  scale_fill_manual(values=my.pal[c(1,5)]) +
  coord_cartesian(ylim=c(-0.55,0.55)) +
  theme_minimal() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("RFD.smooth_min.max.RFDsmooth.dist100kb.pdf")

#### IZ distance and size ####
# figure S2B, C
ok_extreme <- merge(ok_max,ok_min,by="break_ID",suffixes = c(".max",".min"))
ok_extreme$size <- ok_extreme$dist.max - ok_extreme$dist.min
ok_extreme$score <- ok.gr$score[match(ok_extreme$break_ID,ok.gr$ID)]
ok_extreme$DeltaRFD <- ok_extreme$RFDedge.max - ok_extreme$RFDedge.min


#### example regions ####
HOX_list <-list(c("HoxA","chr6",52149505,52274087), 
                c("HoxB","chr11",96226566,96390579),
                c("HoxC","chr15",102849312,103128870),
                c("HoxD","chr2",74651075,74773405))

focus.data <- data.frame(subset(data.gr,mid>32720000 & mid<36720000 & seqnames=="chr3")) # example region figure 1B, 3B, S5D, S6D 
focus.data$focus <- ok_break$score[match(focus.data$ID,ok_break$ID)] # add RFD initiation zones
focus.data$include <- ifelse(focus.data$CPM > 0.3, TRUE,FALSE) # & focus.data$exprs > 500

focus.data <- subset(focus.data,!mark %in% c("input","proseq","K5par","K20par")) # plot on subset of data
focus.data_m <- melt(focus.data[,c("ID","mid","type","focus","F.cpm","R.cpm")],id.vars = c("ID","mid","type","focus"))
ggplot(subset(focus.data_m),aes(mid,value)) + 
  geom_vline(data=subset(focus.data_m,!is.na(focus)),aes(xintercept=mid),color="black") + 
  geom_bar(stat="identity",aes(fill=variable)) +
  #geom_line(color="red") +
  #geom_vline(aes(xintercept=mid,color=abs(score))) +
  theme_minimal() +
  facet_wrap(~ type + variable,ncol=3,scale="free_y") +
  scale_alpha_manual(values = c(0.3,1)) +
  scale_colour_manual(values = c(my.pal[c(2,2,2,3,1,1,1,)],my.pal)) + 
  scale_fill_manual(values = c(my.pal[c(1,6,2,3,5,6,7,8,11,12,13)],my.pal))
#ggsave("K5_K20_type_mESC_RFD_fixY.pdf",height = 19,width = 5.98)

# Ok IZ demonstation
focus.data <- data.frame(subset(data.gr,mid>75060000 & mid<78560000 & seqnames=="chr10")) # RFD region (figure S2A)
focus.data$focus <- ok_break$score[match(focus.data$ID,ok_break$ID)] # add RFD initiation zones

# add min/max and IZ info, for RFD breakpoint detection
focus.data$IZ <- FALSE
# ok_extreme initiated in "Identify initiation zone edges"
focus.data$IZ[match(ok_extreme$break_ID,focus.data$ID)] <- "IZ"
focus.data$IZ[match(ok_extreme$ID.min,focus.data$ID)] <- "min"
focus.data$IZ[match(ok_extreme$ID.max,focus.data$ID)] <- "max"

focus.data$include <- ifelse(focus.data$exprs > 500, TRUE,FALSE) # focus.data$CPM > 0.3 & 
focus.data <- subset(focus.data,mark %in% c("Ok"))

plot <- ggplot(focus.data,aes(mid,score)) + 
  geom_vline(data=subset(focus.data,!is.na(focus)),aes(xintercept=mid),color="black") +
  geom_vline(data=subset(focus.data,IZ== "min" | IZ == "max"),aes(xintercept=mid),color="grey") +
  #geom_vline(aes(xintercept=mid,color=abs(score))) + 
  #geom_bar(aes(fill=mark,colour=mark),stat="identity",size=0.2) + # 
  #scale_colour_distiller(palette="Blues",direction=1) + xlab("position")
  geom_line(aes(fill=mark,colour=mark),size=0.5) +
  #coord_cartesian(ylim=c(0,0.7)) + 
  theme_minimal() +
  scale_fill_manual(values = my.pal[c(6)]) + scale_colour_manual(values = my.pal[c(6)]) 
pdf("mESC_score_chr1075Mb_Ok_IZ_fixY_min_max_line.pdf",width = 5.98,height = 1)
print(plot)
dev.off()

## include gene track for example regions
library(Gviz)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
gtr <- GenomeAxisTrack()
TAD_annot <- AnnotationTrack(TAD.gr,chromosome = "chr3", # [which(Enh.gr$active=="active")]
                             feature = TAD.gr$Compartment, name="TAD")
IZ_annot <- AnnotationTrack(ok.gr,chromosome = "chr3",name="IZ") # [which(Enh.gr$active=="active")]
shift_annot <- AnnotationTrack(shift.gr,chromosome = "chr3",name="shift") # [which(Enh.gr$active=="active")]

enh_annot <- AnnotationTrack(Enh.gr,chromosome = "chr3", # [which(Enh.gr$active=="active")]
                             feature = Enh.gr$active, name="enhancer")

genes.gr <- genes(txdb)
gene_annot <- AnnotationTrack(genes.gr,chromosome = "chr3",
                              feature = genes.gr$RNA_active, name="genes") # strand(genes_ex.gr)
plotTracks(list(gtr,gene_annot,enh_annot),chromosome = "chr3", from = 32720000, to = 36720000,  #  from = 32720000, to = 36720000, ,TAD_annot,IZ_annot,shift_annot
           showTitle = TRUE,showId=FALSE,shape = "fixedArrow",
           #"+"="blue","-"="red",  # fixedArrow
           "TRUE"="red","FALSE"="black",
           "A"="red","B"="blue",
           "active"="red","inactive"="black")

#genes.gr <- transcripts(txdb,filter=list(tx_chrom = myseqnames))
genes.gr$ENSEMBL <- ensem2$ENSEMBL[match(genes.gr$gene_id,ensem2$ENTREZID)]
genes.gr$ALIAS <- ensem2$ALIAS[match(genes.gr$gene_id,ensem2$ENTREZID)]
genes.gr$RNA_active <- RNA_df$active[match(genes.gr$ENSEMBL,RNA_df$gene_id)]
genes_ex.gr <- subset(genes.gr,seqnames=="chr3" & start > 32720000 & end < 36720000)
gene_annot <- AnnotationTrack(genes_ex.gr,chromosome = "chr3",
                              feature = genes_ex.gr$RNA_active, name="genes", id = genes_ex.gr$ALIAS) # strand(genes_ex.gr)
plotTracks(list(gtr,gene_annot),chromosome = "chr3", from = 34000000, to = 36000000,  #  from = 32720000, to = 36720000, ,TAD_annot,IZ_annot,shift_annot
           showTitle = TRUE,showId=T,shape = "fixedArrow",
           #"+"="blue","-"="red",  # fixedArrow
           "TRUE"="red","FALSE"="black",
           "A"="red","B"="blue",
           "active"="red","inactive"="black")
gene_annot@range[is.na(gene_annot@range$feature)]



#### input bias based on logLR track ####
# load files similar to "load data" and store into logLR_data
logLR_data <- rbind(data.frame(K5_3_LR,type=rep("K5_3",times=nrow(K5_3_LR))),
                    data.frame(K5_4_LR,type=rep("K5_4",times=nrow(K5_4_LR))),
                    data.frame(K20_1_LR,type=rep("K20_1",times=nrow(K20_1_LR))),
                    data.frame(K20_2_LR,type=rep("K20_2",times=nrow(K20_2_LR))),
                    data.frame(K20_3_LR,type=rep("K20_3",times=nrow(K20_3_LR))),
                    data.frame(K36_1_LR,type=rep("K36_1",times=nrow(K36_1_LR))),
                    data.frame(K36_3_LR,type=rep("K36_3",times=nrow(K36_3_LR))),
                    data.frame(K36_4_LR,type=rep("K36_4",times=nrow(K36_4_LR))),
                    data.frame(K20m3_1_LR,type=rep("K20m13_1",times=nrow(K20m3_1_LR))),
                    data.frame(K20m3_2_LR,type=rep("K20m13_2",times=nrow(K20m3_2_LR))),
                    data.frame(K36m3_1_LR,type=rep("K36m13_1",times=nrow(K36m3_1_LR))),
                    data.frame(K36m3_2_LR,type=rep("K36m13_2",times=nrow(K36m3_2_LR))),
                    data.frame(K20m6_1_LR,type=rep("K20m16_1",times=nrow(K20m6_1_LR))),
                    data.frame(K20m6_2_LR,type=rep("K20m16_2",times=nrow(K20m6_2_LR))),
                    data.frame(K36m6_1_LR,type=rep("K36m16_1",times=nrow(K36m6_1_LR))),
                    data.frame(K36m6_2_LR,type=rep("K36m16_2",times=nrow(K36m6_2_LR))),
                    data.frame(K5m3_1_LR,type=rep("K5m13_1",times=nrow(K5m3_1_LR))),
                    data.frame(K5m3_2_LR,type=rep("K5m13_2",times=nrow(K5m3_2_LR))),
                    data.frame(K5m6_1_LR,type=rep("K5m16_1",times=nrow(K5m6_1_LR))),
                    data.frame(K5m6_2_LR,type=rep("K5m16_2",times=nrow(K5m6_2_LR))))
colnames(logLR_data) <- c("chr","start","end","sumF","meanF","sumR","meanR","RFD.raw","RFD","RFD.deriv","score","zero.deriv","type")

logLR_data$mid <- rowMeans(logLR_data[,c("start","end")])
logLR_data$ID <- paste(logLR_data$chr,logLR_data$start,sep=":") %>% paste(.,logLR_data$end,sep="-")
logLR_data$mark <- sapply(strsplit(as.character(logLR_data$type),"_"),"[",1)
logLR_data$exprs <- rowSums(logLR_data[,c("meanF","meanR")])
logLR_data$ID_type <- paste(logLR_data$ID,logLR_data$type,sep="_")

data.gr$ID_type <- paste(data.gr$ID,data.gr$type,sep="_")
data.gr$F.LR <- logLR_data$meanF[match(data.gr$ID_type, logLR_data$ID_type)]
data.gr$R.LR <- logLR_data$meanR[match(data.gr$ID_type, logLR_data$ID_type)]
data.gr$RFD.LR <- logLR_data$RFD[match(data.gr$ID_type, logLR_data$ID_type)]
data.gr$RFD.raw.LR <- logLR_data$RFD.raw[match(data.gr$ID_type, logLR_data$ID_type)]

data.gr$LR_pass <- ifelse(data.gr$F.LR>0 | data.gr$R.LR>0, "pass","fail")


#### Enhancer active/inactive test ####
# mark specific min/max amplitude
data_F <- subset(data.break.df,(dist>0 & !mark %in% c("K5","K5m13","K5m16")) | dist<0 & mark %in% c("K5","K5m13","K5m16"))
data_R <- subset(data.break.df,(dist<0 & !mark %in% c("K5","K5m13","K5m16")) | dist>0 & mark %in% c("K5","K5m13","K5m16")) 

data_max <- data.frame(data_F %>% 
                       dplyr::group_by(break_ID,enh_active,type,mark) %>%
                       dplyr::summarise(ID = ID[which.max(RFD)],
                                        #RFDedge = RFD[which.max(abs(RFD))],
                                        RFDedge = max(RFD),
                                        dist = dist[which.max(RFD)]))
data_min <- data.frame(data_R %>% 
                       dplyr::group_by(break_ID,enh_active,type,mark) %>%
                       dplyr::summarise(ID = ID[which.min(RFD)],
                                        #RFDedge = RFD[which.max(RFD)],
                                        RFDedge = min(RFD),
                                        dist = dist[which.min(RFD)]))

data_max$break_type <- paste(data_max$break_ID,data_max$type,sep="-")
data_min$break_type <- paste(data_min$break_ID,data_min$type,sep="-")
data_edge <- merge(data_max,data_min,by="break_type",suffixes = c(".max",".min"))

data_edge$size <- data_edge$dist.max - data_edge$dist.min
data_edge$dRFD <- data_edge$RFDedge.max - data_edge$RFDedge.min

ggplot(subset(data_edge,mark.min %in% c("K5","K20","K36","K5m13","K20m13","K36m13")),aes(enh_active.min,dRFD)) +
  geom_boxplot(aes(colour=enh_active.min)) +
  facet_wrap(~type.min,scale="free_y")

w.pval_SE <- w.pval_active <- t.pval_SE <- t.pval_active <- no_SE <- no_active <- list()
for(type in type_list) { 
  no_SE[type] <- length(data_edge$dRFD[which(data_edge$enh_active.min=="SE" & data_edge$type.max == type)])
  no_active[type] <- length(data_edge$dRFD[which(data_edge$enh_active.min=="active" & data_edge$type.max == type)])
  w.pval_SE[type] <- wilcox.test(data_edge$dRFD[which((data_edge$enh_active.min=="SE" | data_edge$enh_active.min=="active") & data_edge$type.max == type)],data_edge$dRFD[which(data_edge$enh_active.min=="inactive" & data_edge$type.max == type)])$p.value
  #w.pval_active[type] <- wilcox.test(data_edge$dRFD[which(data_edge$enh_active.min=="active" & data_edge$type.max == type)],data_edge$dRFD[which(data_edge$enh_active.min=="inactive" & data_edge$type.max == type)])$p.value
  #t.pval_SE[type] <- t.test(data_edge$dRFD[which(data_edge$enh_active.min=="SE" & data_edge$type.max == type)],data_edge$dRFD[which(data_edge$enh_active.min=="inactive" & data_edge$type.max == type)])$p.value
  #t.pval_active[type] <- t.test(data_edge$dRFD[which(data_edge$enh_active.min=="active" & data_edge$type.max == type)],data_edge$dRFD[which(data_edge$enh_active.min=="inactive" & data_edge$type.max == type)])$p.value
}

pval_df <- data.frame(type = names(w.pval_SE),
                      w.SE = unlist(w.pval_SE),
                      #t.SE = unlist(t.pval_SE),
                      noSE = unlist(no_SE),
                      #w.active = unlist(w.pval_active),
                      #t.active = unlist(t.pval_active),
                      noActive = unlist(no_active))

ggplot(pval_df,aes(type,-log10(active))) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(0.05)) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  
# defined by OK edges
ok_max$break_edge <- paste(ok_max$break_ID,ok_max$ID,sep="-")
ok_min$break_edge <- paste(ok_min$break_ID,ok_min$ID,sep="-")

data.break.df$break_edge <- paste(data.break.df$break_ID,data.break.df$ID,sep="-")

data_max <- data.break.df[which(data.break.df$break_edge %in% ok_max$break_edge),]
data_min <- data.break.df[which(data.break.df$break_edge %in% ok_min$break_edge),]

data_max$break_type <- paste(data_max$break_ID,data_max$type,sep="-")
data_min$break_type <- paste(data_min$break_ID,data_min$type,sep="-")
data_edge <- merge(data_max,data_min,by="break_type",suffixes = c(".max",".min"))

data_edge$active <- Enh_SE.gr$active[match(data_edge$break_ID.max,Enh_SE.gr$name)]
data_edge$size <- data_edge$dist.max - data_edge$dist.min
data_edge$dRFD <- data_edge$RFD.max - data_edge$RFD.min

ggplot(subset(data_edge,type.max %in% c("K5_3","K5_4","K5m13_1","K5m13_2","K5m16_1","K5m16_2")),aes(active,dRFD)) +
  geom_boxplot(aes(colour=active)) +
  facet_wrap(~type.min)

w.pval_SE <- w.pval_active <- w.no_SE <- w.no_active <- list()
for(type in type_list) {
  w.no_SE[type] <- length(data_edge$dRFD[which(data_edge$active=="SE" & data_edge$type.max == type)])
  w.no_active[type] <- length(data_edge$dRFD[which(data_edge$active=="active" & data_edge$type.max == type)])
  w.pval_SE[type] <- wilcox.test(data_edge$dRFD[which(data_edge$active=="SE" & data_edge$type.max == type)],data_edge$dRFD[which(data_edge$active=="inactive" & data_edge$type.max == type)])$p.value
  w.pval_active[type] <- wilcox.test(data_edge$dRFD[which(data_edge$active=="active" & data_edge$type.max == type)],data_edge$dRFD[which(data_edge$active=="inactive" & data_edge$type.max == type)])$p.value
}

pval_df <- data.frame(type = names(w.pval_active),
                        SE = unlist(w.pval_SE),
                        noSE = unlist(no_SE),
                        active = unlist(w.pval_active),
                        noActive = unlist(no_active))

ggplot(w.pval_df,aes(type,-log10(active))) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(0.05)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))

# data.extreme <- data.trim[which(data.trim$ID %in% ok_top.gr$ID),]
exprs.mat <- dcast(data.extreme, ID ~ type, value.var = "exprs", fun.aggregate = mean)
RFD.mat <- dcast(data.extreme, ID ~ type, value.var = "RFD", fun.aggregate = mean)
mat <- data.frame(RFD.mat,exprs=exprs.mat[,-1])

ggplot(mat, aes(Ok,K20m13_1)) +
  geom_hex(bins=100) +
  #geom_point(aes(colour=exprs.proseq),alpha=0.5,size=1) +
  scale_fill_gradientn(colours=c("orange","red")) +
  geom_smooth(se=F,method="lm",size=0.4,colour="red") +
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE,size=3) +
  #scale_colour_continuous(low="blue",high="yellow",trans="log10") +
  #coord_cartesian(ylim=c (-1,1)) +
  #xlab("K5 rep1 Repartition") + ylab("K5 rep2 Repartition") +
  theme_minimal()

#source("/isdata/alab/people/maria/scripts/ggplot_hex.R")
df = data.frame(x=mat$Ok, y=mat$K5,var = log10(mat$exprs.Ok))
ggplot(df, aes(x=x, y=y, col=var)) + geom_point()
plot(hex_bin(x=df$x, y=df$y, var4=df$var))
plot(hex_bin(x=df$x, y=df$y, var4=df$var, frequency.to.area=TRUE))


#### heatmaps ####
# heatmap of SCAR-seq sample or mark (S3B)
heat_mark = "K5"
#data.break.df: signal around initiatin zones
#heat.mat <- dcast(data.break.df[data.break.df$type==heat_mark,],break_ID ~ dist,value.var = "RFD") # pr. samples
heat.mat <-  dcast(data.break.df[data.break.df$mark==heat_mark,],break_ID ~ dist,value.var = "RFD", fun.aggregate = mean) # average on histone mark
rownames(heat.mat) <- heat.mat$break_ID
heat.mat <- heat.mat[,-1]

# order according to RFD defined initiation zone width (distance between edges)
break_df <- data.frame(min=ok_min,max=ok_max)
break_order <- break_df$min.break_ID[order(break_df$min.dist-break_df$max.dist,decreasing = T)]
heat.mat <- heat.mat[match(break_order,rownames(heat.mat)),]

# filter out top outliers to improve signal (optional)
#high_qt <- quantile(heat.mat,0.9999,na.rm=T)
#low_qt <- quantile(heat.mat,0.0001,na.rm=T)
#
#heat.mat <- heat.mat[apply(heat.mat,1,function(x) max(x,na.rm=T))<high_qt,]
#heat.mat <- heat.mat[apply(heat.mat,1,function(x) min(x,na.rm=T))>low_qt,]

# filter out regions with too many NAs (optional)
#qt <- apply(heat.mat,1,function(x) length(which(is.na(x)))) %>% quantile(0.9)
#heat.plot.mat <- heat.mat[apply(heat.mat,1,function(x) length(which(is.na(x)))<=qt),]

# costumize heatmap colour and range (optional)
breaksList = c(seq(-max(abs(heat.plot.mat),na.rm=T),0, by = 0.01),
               rev(seq(-max(abs(heat.plot.mat),na.rm=T),0, by = 0.01)*-1))
white_space <- round((length(breaksList)/2)/4)
colList <- (length(breaksList)/2) + white_space
col_pal <- c(colorRampPalette(rev(brewer.pal(n=9,name="Blues")))(colList)[-c((colList-(white_space-1)):colList)],
             rev(colorRampPalette(rev(brewer.pal(n=9,name="Reds")))(colList)[-c((colList-(white_space-1)):colList)]))

heat.res <- pheatmap(heat.plot.mat,
                     color=col_pal,
                     #colorRampPalette(rev(brewer.pal(n=8,name="RdBu")))(length(breaksList)),
                     cluster_rows = F, cluster_cols = F,
                     breaks = breaksList,
                     filename = paste0(heat_mark,"_heatmap_dist300kb_IZ_Ok_order.pdf"),
                     show_rownames = F, show_colnames = F)



#### correlate regions ####
RFD.mat <- subset(data.gr, mark %in% c("Ok","K5","K20","K20m13","K20m16","K36","K36m13","K5m13","K5m16","K36m16","proseq")) %>% data.frame %>% 
  dcast(., ID ~ mark, value.var = "RFD", fun.aggregate = mean)

## correlate with Hi-C directionality index
#RFD.mat$DI <- data.gr$DI[match(RFD.mat$ID,data.gr$ID)] 

RFD.mat_m <- melt(RFD.mat) # long format data.frame
data.joint <- left_join(RFD.mat_m, RFD.mat_m, by=c("ID"))

# example: correlate K5 with K20 WT and mutants
data.joint <- data.joint[data.joint$variable.x %in% c("K5","K5m13","K5m16"),]
data.joint <- data.joint[data.joint$variable.y %in% c("K20","K20m13","K20m16"),]

# hex plot all data points (figure S7D, 3C-D)
breaks = seq(0,1300, by = 300) # adjust to fit marks
col_pal <- colorRampPalette(brewer.pal(n=7,name="Greens"))(7)
ggplot(subset(data.joint,variable.y == "K20" & variable.x == "K5"), aes(value.x,value.y)) +
    geom_hex(bins=120) +
    geom_smooth(se=F,method="lm",size=0.3,colour="red") +
     scale_fill_gradientn(colours = (brewer.pal(n=9,name="Blues")[2:8])) + # ,limits=c(0,1250)) + # ,trans="sqrt"
    #ggtitle(mark) +
    coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) +
    xlab("RFD") + ylab("Partition") +
    theme_minimal()
ggsave(paste0(mark1,mark2,"_RFDsmooth_blues_limits.pdf"),width = 3.39 ,height = 2.67)

# correlate all marks with RFD (OK-seq) 
data.joint <- left_join(RFD.mat_m, RFD.mat_m, by=c("ID"))
# example: correlate K5 with K20 WT and mutants
data.joint <- data.joint[data.joint$variable.x %in% "Ok",]

# loop over all marks and calculate correlations to OK RFD
cor_list <- pval_list <- list() 
for(mark in unique(data.joint$variable.y)) {
  cor_list[mark] <- cor.test(data.joint$value.x[data.joint$variable.y==mark],data.joint$value.y[data.joint$variable.y==mark],method = "spearman",exact=FALSE)$estimate
  pval_list[mark] <- cor.test(data.joint$value.x[data.joint$variable.y==mark],data.joint$value.y[data.joint$variable.y==mark],method = "spearman",exact=FALSE)$p.value
}

cor_df <- data.frame(estimate=unlist(cor_list),
                     pval=unlist(pval_list),
                     mark = names(pval_list))

# bar plot correlations (such as in figure S8D, E)
ggplot(cor_df,aes(mark,estimate)) +
  geom_bar(aes(fill=mark),stat="identity",position="dodge") +
  scale_fill_manual(values=c(my.pal[c(7,2,2,2,1,1,1,12,12,12,12,6)],my.pal)) +
  theme_minimal() +
  facet_wrap(~mark) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) 
ggsave("Correlation_spearman_mark_vs_OK.pdf")


## K36 bam correlation (figure S5C, S8A)
cor_mat <- read.table(file.path(data_dir,"footprint/bam_pearson_cor_K36_genes.matrix.txt"))
rownames(cor_mat) <- gsub(".nodup.bam","",rownames(cor_mat))
colnames(cor_mat) <- gsub(".nodup.bam","",colnames(cor_mat))

breaksList = c(seq((-1),1,by=0.005))
col_pal <- colorRampPalette(brewer.pal(n=7,name="RdBu"))(length(breaksList)+3)

pheatmap(cor_mat,
         color=col_pal,
         breaks = breaksList,
         cluster_rows = T, cluster_cols = F,
         filename = "K36wt_FE_K36encode_bigwig_pearson_cor_bs1kb_col2.pdf",
         show_rownames = T, show_colnames = T)

#### define mutant shifts ####
data.mut <- read.table("RFD/results/res_K36_r1_smooth_results_w1000_s30_d30_z1.txt",header=FALSE,as.is=TRUE,sep="\t",fill=TRUE)
colnames(data.mut) <- c("chr","start","end","F","R","F.cpm","R.cpm","RFD.raw","RFD","RFD.deriv","score","zero.deriv")
data.mut <- subset(data.mut, chr %in% myseqnames)
data.mut$exprs <- rowSums(data.mut[,c(4,5)])
data.mut$CPM <- rowSums(data.mut[,c(6,7)])

data.mut <- subset(data.mut, CPM >= 0.3 & exprs > 500)
#data.mut <- subset(data.mut, exprs > quantile(exprs,probs=0.1))

ids <- which(data.mut[,12]>quantile(data.mut[,10],probs=0.9,na.rm=TRUE))

write.table(cbind(data.mut[ids,1:3],ids),"shifts_cpm.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("bedtools merge -c 4 -d 3000 -o collapse -i shifts_cpm.bed > shifts_3000merged_cpm.bed")
b.data <- read.table("shifts_3000merged_cpm.bed",sep="\t",as.is=TRUE,header=FALSE)

shift_df <- data.mut[sapply(b.data[,4], function(x) {
  y <- as.numeric(strsplit(x,",")[[1]])
  y[which.max(data.mut[y,10])]
}),]

shift_df$mid <- rowMeans(shift_df[,c("start","end")])
shift_df$ID <- paste(shift_df$chr,shift_df$start,sep=":") %>% paste(.,shift_df$end,sep="-")
shift.gr <- makeGRangesFromDataFrame(shift_df,
                                     seqnames.field="chr",
                                     start.field = "start",
                                     end.field = "end",
                                     keep.extra.columns = T)
save(shift.gr,file="shiftGR_K36_r1.RDat")

#### RT and U domains ####
# load replication timing data
RTdata <- read.table(file.path(server,data_dir,"extern_data/RT/RT_binned_results_w1000.txt"),sep="\t",header=FALSE,as.is=TRUE)
#RTdata <- read.table(file.path(server,data_dir,"extern_data/RT/RT-D3-1.5.bed"),sep="\t",header=FALSE,as.is=TRUE)
RTdata <- RTdata[,c(1:3,6)]
colnames(RTdata) <- c("chr","start","end","log2FC")
RTdata$ID <- paste(RTdata$chr,RTdata$start,sep=":") %>% paste(.,RTdata$end,sep="-")
RT.gr <- makeGRangesFromDataFrame(RTdata,
                                  seqnames.field="chr",
                                  start.field = "start",
                                  end.field = "end",
                                  keep.extra.columns = T)
#RT.gr <- reduce(RT.gr,min.gapwidth=20000)
#RT.gr$log2FC[RT.gr$log2FC==0] <- NA
RT.gr <- subset(RT.gr,seqnames %in% myseqnames)

# add replication timing to data range
data.gr$RT <- RT.gr$log2FC[match(data.gr$ID,RT.gr$ID)]

# load U-domain coordinates (figure S4)
Udomains <- read.table(file.path(server,data_dir,"extern_data/U-domainBG01mapped_mm10.bed"),sep="\t",header=FALSE,as.is=TRUE)
colnames(Udomains) <- c("chr","start","end")
Udomains$mid <- rowMeans(Udomains[,c("start","end")])
Udomains$ID <- paste(Udomains$chr,Udomains$start,sep=":") %>% paste(.,Udomains$end,sep="-")
U.gr <- makeGRangesFromDataFrame(Udomains,
                                 seqnames.field="chr",
                                 start.field = "start",
                                 end.field = "end",
                                 keep.extra.columns = T)
border.func <- function(gr) {
  start.gr <-gr
  start.gr$pos <- "start"
  end(start.gr) <- start(start.gr)
  end.gr <-gr
  end.gr$pos <- "end"
  start(end.gr) <- end(end.gr)
  c(start.gr,end.gr)
  #reduce(c(start.gr,end.gr))
}

U.border.gr <- border.func(U.gr)

#### distantce enrichment ####

# nearest enrichments
x.gr <- ok.gr
#y.gr <- data.gr[which(data.gr$mark=="Ok" & !is.na(data.gr$RFD.extreme))][,"ID"]

#load(file.path(data_dir,"RFD/shiftGR_mcm2_mut313_r2.RDat"))
y.gr <- unique(TF.gr[,0]) # shift.gr

# define background set
#background.gr <- cage.gr[-(subjectHits(findOverlaps(y.gr,cage.gr)))]
background.gr <- unique(data.gr[data.gr$annotation=="Promoter" & !data.gr$ID %in% y.gr$ID][,0])
#background.gr <- data.gr[data.gr$type=="K20m13_2" & !data.gr$ID %in% y.gr$ID]

# actual distance
x_y_dist <- distanceToNearest(x.gr,y.gr)
dist_y <- x_y_dist@elementMetadata$distance+1
dist_y[start(x.gr) > start(y.gr)[subjectHits(x_y_dist)]] <- dist_y[start(x.gr) > start(y.gr)[subjectHits(x_y_dist)]] * (-1)

# random distance
dist_ran_list <- list()
for(i in 1:10) {
  ran_y.gr <- sample(background.gr,length(y.gr))
  #ran_y.gr <- sample(data.gr[data.gr$mark=="Ok" & !data.gr$ID %in% y.gr$ID],length(y.gr))
  x_ran_dist <- distanceToNearest(x.gr,ran_y.gr)
  
  dist_ran_list[i] <- list(x_ran_dist@elementMetadata$distance)
  dist_ran_list[[i]][start(x.gr) > start(ran_y.gr)[subjectHits(x_ran_dist)]] <- dist_ran_list[[i]][start(x.gr) > start(ran_y.gr)[subjectHits(x_ran_dist)]] * (-1)
}
ran_dist_df <- data.frame(do.call(cbind,dist_ran_list))

# quick plot
distnace <- data.frame(dist_y=unlist(dist_y),
                       dist_ran = unlist(dist_ran_list[2]))
distnace_m <- melt(distnace)
ggplot(distnace_m,aes(value)) + 
  #geom_density(alpha=0.4,aes(fill=variable,colour=variable)) +
  geom_histogram(alpha=0.2,aes(fill=variable,colour=variable),bins = 2000) +
  scale_fill_manual(values=c(my.pal[c(6,8)])) + scale_colour_manual(values=my.pal[c(6,8)]) + 
  #geom_vline(xintercept = c(150000,400000),colour="red3",size=0.2) +
  #facet_wrap(~variable,scale="free_y") +
  # coord_cartesian(xlim=c(-1000000,1000000)) +
  theme_minimal() 
#ggsave("distance_quick_IZ_random_histogram_coord1000kb_TAD_A_borders.pdf")

# add all TSSs to pool 
x_tss_dist <- distanceToNearest(x.gr,ES_annot.gr)
dist_tss <- x_tss_dist@elementMetadata$distance+1
dist_tss[start(x.gr) > start(ES_annot.gr)[subjectHits(x_tss_dist)]] <- dist_tss[start(x.gr) > start(ES_annot.gr)[subjectHits(x_tss_dist)]] * (-1)

# faction of overlap pr. tested regions pull.
frac_y_list <- frac_ran_list <- frac_tss_list <- list()
i_0 = -500000
for(i in seq(-490000,500000,by=10000)) {
  #print(length(which(dist==i)) )
  frac_y_list[as.character(i)] <- length(which(dist_y>=i_0 & dist_y<i)) / length(x.gr)
  #frac_y_list[as.character(i)] <- length(which(dist_y==i)) / length(y.gr)
  
  # all tss's 
  frac_tss_list[as.character(i)] <- length(which(dist_tss>=i_0 & dist_tss<i)) / length(x.gr)
  
  boot_frac <- list()
  for(j in 1:ncol(ran_dist_df)) {
    #frac_ran_list[as.character(i)] <- length(which(dist_ran>=i_0 & dist_ran<i)) / length(ran_y.gr)
    boot_frac[j] <- length(which(ran_dist_df[,j]>=i_0 & ran_dist_df[,j]<i)) / length(x.gr) #nrow(ran_dist_df)
  }
  frac_ran_list[as.character(i)] <- mean(unlist(boot_frac))
  #frac_ran_list[as.character(i)] <- length(which(dist_ran==i)) / length(ran_y.gr)
  i_0 <- i
}

dist_frac_df <- data.frame(dist = seq(-490000,500000,by=10000),
                           frac_y = as.numeric(unlist(frac_y_list)),
                           frac_tss = as.numeric(unlist(frac_tss_list)),
                           frac_ran = as.numeric(unlist(frac_ran_list)))

#t.test(dist_frac_df$frac_y,dist_frac_df$frac_ran)
dist_frac_m <- melt(dist_frac_df,measure.vars=c("frac_y","frac_ran","frac_tss")) #frac_tss
dist_frac_m$dist_num <- as.numeric(as.character(dist_frac_m$dist))
dist_frac_m$dist_num <- dist_frac_m$dist_num/1000
dist_frac_m$variable <- factor(dist_frac_m$variable,levels = c("frac_ran","frac_y","frac_tss")) # frac_tss
ggplot(dist_frac_m,aes(dist_num,value)) +
  #geom_bar(stat="identity",position = "dodge",aes(fill=variable,colour=variable),size=0.3,alpha=0.3) +
  geom_line(aes(colour=variable),size=0.5) +
  #geom_smooth(fill=NA,aes(colour=variable))  +
  geom_hline(yintercept = mean(dist_frac_df$frac_ran),linetype=3) + 
  #geom_vline(xintercept = 0,linetype=3) +
  #coord_cartesian(ylim=c(0,0.065)) +
  theme_minimal() +
  scale_color_manual(values = c("grey",my.pal[c(6,12)])) + scale_fill_manual(values = c("grey",my.pal[c(6,12)]))
ggsave("Nearest_distance_K20m13_r2_shift_vs_IZ_fractionOverlap_10perm_1000Kb.pdf")


# distance from TF TSss
TF.gr <- data.gr[data.gr$TF==TRUE]
ranID <- sample(data.gr$ID[data.gr$annotation=="Promoter" & !data.gr$ID %in% TF.gr$ID],length(TF.gr))
data.gr$RanTF <- ifelse(data.gr$ID %in% ranID, TRUE,FALSE)

y.gr <- ok.gr
x.gr <- unique(TF.gr[,0]) # shift.gr
background.gr <- unique(data.gr[data.gr$annotation=="Promoter" & !data.gr$ID %in% x.gr$ID][,0])

# actual distance
x_y_dist <- distanceToNearest(x.gr,y.gr)
dist_y <- x_y_dist@elementMetadata$distance+1
#dist_y[start(x.gr) > start(y.gr)[subjectHits(x_y_dist)]] <- dist_y[start(x.gr) > start(y.gr)[subjectHits(x_y_dist)]] * (-1)

# random distance
dist_ran_list <- list()
for(i in 1:10) {
  ran_y.gr <- sample(background.gr,length(x.gr))
  #ran_y.gr <- sample(data.gr[data.gr$mark=="Ok" & !data.gr$ID %in% y.gr$ID],length(y.gr))
  x_ran_dist <- distanceToNearest(ran_y.gr,y.gr)
  
  dist_ran_list[i] <- list(x_ran_dist@elementMetadata$distance)
  #dist_ran_list[[i]][start(y.gr) > start(ran_y.gr)[subjectHits(x_ran_dist)]] <- dist_ran_list[[i]][start(y.gr) > start(ran_y.gr)[subjectHits(x_ran_dist)]] * (-1)
}
ran_dist_df <- data.frame(do.call(cbind,dist_ran_list))

t.test(dist_ran_list[[1]],dist_y)

#### annotate ####
# annotate the full data range based on UCSCmm10
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
annot <- annotatePeak(data.gr, tssRegion=c(-1,1),
                      TxDb=txdb, annoDb="org.Mm.eg.db", level="transcript",
                      sameStrand = FALSE, overlap="all")
# currate UCSC annotations
annot.df <- data.frame(annot)
annot.df$annotation_cur <- annot.df$annotation
annot.df$annotation_cur[grep("exon",annot.df$annotation)] <- "exon"
annot.df$annotation_cur[grep("intron",annot.df$annotation)] <- "intron"
annot.df$annotation_cur[grep("Downstream",annot.df$annotation)] <- "Downstream <3kb"

# add annotations to data range
data.gr$annotation <- annot.df$annotation[match(data.gr$ID,annot.df$ID)]
data.gr$annotation_cur <- annot.df$annotation_cur[match(data.gr$ID,annot.df$ID)]
data.gr$SYMBOL <- annot.df$SYMBOL[match(data.gr$ID,annot.df$ID)]
data.gr$ENSEMBL <- annot.df$ENSEMBL[match(data.gr$ID,annot.df$ID)]
data.gr$geneStrand <- annot.df$geneStrand[match(data.gr$ID,annot.df$ID)]
data.gr$geneStrand <- ifelse(data.gr$geneStrand=="1","+","-")
data.gr$distToTSS <- annot.df$distanceToTSS[match(data.gr$ID,annot.df$ID)]
data.gr$geneLength <- annot.df$geneLength[match(data.gr$ID,annot.df$ID)]
data.gr$geneStart <- annot.df$geneStart[match(data.gr$ID,annot.df$ID)]
data.gr$geneEnd <- annot.df$geneEnd[match(data.gr$ID,annot.df$ID)]

# simplify annotation
data.gr$annotation_simple <- ifelse(data.gr$annotation_cur %in% c("intron_first","intron_other","exon_first","exon_other","3' UTR","5' UTR","Promoter"),"gene","intergenic")
data.gr$annotation_simple[data.gr$annotation_simple == "gene"] <- paste(data.gr$annotation_simple[data.gr$annotation_simple == "gene"],data.gr$RNA_active[data.gr$annotation_simple == "gene"],sep=" ")

# RNA-seq activity (M. Snyder RNA-seq 2 mES samples)
RNA_r1 <- read.table(file.path(data_dir,"/extern_data/M_Snyder/rep1_gene_ENCFF243VGH.tsv"),header = T,as.is = T)
RNA_r2 <- read.table(file.path(data_dir,"/extern_data/M_Snyder/rep2_gene_ENCFF476NAH.tsv"),header = T,as.is = T)
RNA_df <- merge(RNA_r1,RNA_r2,by="gene_id",suffixes = c("_r1","_r2"))
RNA_df$FPKM <- rowMeans(RNA_df[,c("FPKM_r1","FPKM_r2")])
RNA_df$active <- ifelse(RNA_df$FPKM>0.5,TRUE,FALSE)

RNA_df$gene_id <- sapply(str_split(RNA_df$gene_id,"\\."),"[[",1)
data.gr$FPKM <- RNA_df$FPKM[match(data.gr$ENSEMBL, RNA_df$gene_id)]

ensem2 <- bitr(RNA_df$gene_id, 
               fromType = "ENTREZID", 
               toType = c("ENSEMBL","ALIAS"),
               OrgDb = "org.Mm.eg.db")
RNA_df$ALIAS <- ensem2$ALIAS[match(RNA_df$gene_id,ensem2$ENTREZID)]
data.gr$RNA_active <- ifelse(data.gr$ENSEMBL %in% RNA_df$gene_id[RNA_df$active], TRUE,FALSE)


# CAGE enhancers
load("extern_data/FANTOM5/F5_enh_robust_nonTSSoverlap_GRange.RDat")
start(Enh.gr) <- Enh.gr$peakStart
end(Enh.gr) <- Enh.gr$peakEnd
Enh.gr <- Enh.gr[!is.na(Enh.gr$active)]
enh_ovarlap <- findOverlaps(data.gr,Enh.gr) #[which(Enh.gr$active=="active")]
data.gr$enhancer <- data.gr$enhancerNo <- FALSE
data.gr$enhancer[queryHits(enh_ovarlap)] <-  Enh.gr$active[subjectHits(enh_ovarlap)]
data.gr$enhancerNo[queryHits(enh_ovarlap)] <-  Enh.gr$name[subjectHits(enh_ovarlap)]

# Super enhancers
SEnh_meta <- read.csv(file.path(file.path(server,data_dir,"extern_data/enhancers/mmc4.csv")),skip = 3,sep=";",header = T)

SEnh <- read.table(file.path(file.path(server,data_dir,"extern_data/enhancers/cell_6822_mmc1_ESc.bed")),skip=1,header = F)
SEnh <- read.table(file.path(file.path(server,data_dir,"extern_data/enhancers/mmc2.csv")),skip=2,sep=";",header=T)

SEnh.gr <- makeGRangesFromDataFrame(SEnh[SEnh$isSuper=="YES",],seqnames.field = "chrom",
                                    start.field = "start",
                                    end.field = "end",
                                    keep.extra.columns = T)
chain <- import.chain("/isdata/alab/people/maria/genome/mm9ToMm10.over.chain")
SEnh.gr <- unlist(liftOver(SEnh.gr,chain))

ensem2 <- bitr(SEnh.gr$proximal_gene, 
               fromType = "ACCNUM",
               toType = c("ENTREZID"),
               OrgDb = "org.Mm.eg.db")

SEnh.gr$ENTREZID <- ensem2$ENTREZID[match(SEnh.gr$proximal_gene,ensem2$ACCNUM)]
SEnh.gr$name <- paste(seqnames(SEnh.gr),start(SEnh.gr),sep="-") %>% paste(.,end(SEnh.gr),sep=":")
SEnh.gr <- subset(SEnh.gr,seqnames %in% myseqnames)

SEnh.gr$EPU <- SEnh_meta$Enhancer_and_gene_overlap_EPU[match(SEnh.gr$ID,SEnh_meta$Enhancer_ID)]
SEnh.gr$EPtopo <- SEnh_meta$Enhancer_gene_overlap_topoDomain[match(SEnh.gr$ID,SEnh_meta$Enhancer_ID)]

ok_dist <- distanceToNearest(SEnh.gr,ok.gr)
SEnh.gr[ok_dist@elementMetadata$distance<10000]

overlap <- findOverlaps(ok.gr,SEnh.gr,maxgap = 20000)
SEnh.gr$IZ <- FALSE
SEnh.gr$IZ[subjectHits(overlap)] <- TRUE

#Enh_noCTCF.gr <- subsetByOverlaps(Enh.gr,TAD.border.gr,invert=T,maxgap = 100000)
#SEnh.gr$active <- "SE"
#Enh_SE.gr <- c(Enh.gr[,c(1,9)],SEnh.gr[,c(14,18)])

Senh_ovarlap <- findOverlaps(data.gr,SEnh.gr)
data.gr$Senhancer <- FALSE
data.gr$Senhancer[queryHits(Senh_ovarlap)] <-  TRUE

# TAD borders
TAD_ovarlap <- findOverlaps(data.gr,TAD.border.gr)
data.gr$TADborder <- FALSE
data.gr$TADborder[queryHits(TAD_ovarlap)] <-  TAD.border.gr$Compartment[subjectHits(TAD_ovarlap)]

#data.gr$TADborder <- FALSE
#data.gr$TADborder[queryHits(TAD_ovarlap)] <- TRUE

#### Compare genome background and data annotation ####
# figure S1C, S2D, S8B

# load 1kb binned genome coordinates
genome.df <- read.table(file.path(server,data_dir,"/macs2_files/mm10_w1000.bed"),sep="\t",header=FALSE,as.is=TRUE)
genome.gr <- makeGRangesFromDataFrame(genome.df,
                                      seqnames.field = "V1",
                                      start.field = "V2",
                                      end.field = "V3",
                                      keep.extra.columns = T)
# use standard chromosomes
genome.gr <- subset(genome.gr, seqnames %in% myseqnames)

# enclude non-mappable regions
genome.gr <- subset(genome.gr, start > 3000000) 
genome.gr <- subsetByOverlaps(genome.gr,neg_gr,invert=T)

# annotate genomic background
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genome_annot.gr <- annotatePeak(genome.gr, tssRegion=c(-1,1),
                         TxDb=txdb, annoDb="org.Mm.eg.db", level="transcript",
                         sameStrand = FALSE, overlap="all")@anno

genome_annot.gr$annotation_cur <- genome_annot.gr$annotation
genome_annot.gr$annotation_cur[grep("exon",genome_annot.gr$annotation)] <- "exon"
genome_annot.gr$annotation_cur[grep("intron",genome_annot.gr$annotation)] <- "intron"
genome_annot.gr$annotation_simple <- ifelse(genome_annot.gr$annotation_cur %in% c("intron","exon","exon_other","3' UTR","5' UTR","Promoter"),"gene","intergenic")

# initiation zone vs. data annotation (figure S2D)
# initiation zone (data.Ok) as gRange
IZ.gr$annotation_cur <- IZ.gr$annotation
IZ.gr$annotation_cur[grep("exon",IZ.gr$annotation)] <- "exon"
IZ.gr$annotation_cur[grep("intron",IZ.gr$annotation)] <- "intron"
IZ.gr$annotation_cur[grep("Downstream",IZ.gr$annotation)] <- "Downstream <3kb"
IZ.gr$annotation_simple <- ifelse(IZ.gr$annotation_cur %in% c("intron","exon","exon_other","3' UTR","5' UTR","Promoter"),"gene","intergenic")

annot_data <- data.frame(table(data.gr$annotation_cur,data.gr$mark))
annot.df$mark = "genome"
annot_genome <- data.frame(table(annot.df$annotation_cur,annot.df$mark))
annot_genome <- rbind(annot_genome,annot_data)

annot_IZ <- data.frame(table(IZ.gr$annotation_cur,IZ.gr$mark))
mark_count <- data.frame(table(IZ.gr$mark))
annot_IZ$Fraq <- annot_IZ$Freq / mark_count$Freq[match(annot_IZ$Var2,mark_count$Var1)]
annot_IZ$Var2 <- "IZ"
annot_genome <- rbind(annot_genome,annot_IZ)

mark_count <- rbind(data.frame(table(data.gr$mark)),
                    data.frame(Var1 = "genome",Freq = 2405535))
annot_genome$Fraq <- annot_genome$Freq / mark_count$Freq[match(annot_genome$Var2,mark_count$Var1)]
ggplot(subset(annot_genome,Var2 %in% c("IZ","Ok","genome")),aes(Var1,Fraq)) +
  geom_bar(aes(fill=Var2),stat="identity",position="dodge") +
  theme_minimal() +
  coord_flip() +
  scale_fill_manual(values = my.pal[c(6,5,7)])
ggsave("IZ_Ok_genome_annot_data_genome_frac.pdf")

# gene orientation
gene_st.gr <- resize(genes.gr, 1)
gene_en.gr <- resize(genes.gr, 1,fix="end")

data.gr$upstreamGenes <- follow(data.gr, unstrand(gene_st.gr))
data.gr$upstreamOrient <- NA
data.gr$upstreamOrient[!is.na(data.gr$upstreamGenes)] <- strand(gene_st.gr)[data.gr$upstreamGenes[!is.na(data.gr$upstreamGenes)]]
data.gr$upstreamDist <- NA
data.gr$upstreamDist[which(!is.na(data.gr$upstreamGenes) & data.gr$upstreamOrient == "+")] <- distance(data.gr[which(!is.na(data.gr$upstreamGenes) & data.gr$upstreamOrient == "+")],
                                                                                                       gene_en.gr[data.gr$upstreamGenes[which(!is.na(data.gr$upstreamGenes) & data.gr$upstreamOrient == "+")]])
data.gr$upstreamDist[which(!is.na(data.gr$upstreamGenes) & data.gr$upstreamOrient == "-")] <- distance(data.gr[which(!is.na(data.gr$upstreamGenes) & data.gr$upstreamOrient == "-")],
                                                                                                       gene_st.gr[data.gr$upstreamGenes[which(!is.na(data.gr$upstreamGenes) & data.gr$upstreamOrient == "-")]])
data.gr$downstreamGenes <- precede(data.gr, unstrand(gene_st.gr))
data.gr$downstreamOrient <- NA
data.gr$downstreamOrient[!is.na(data.gr$downstreamGenes)] <- strand(gene_st.gr)[data.gr$downstreamGenes[!is.na(data.gr$downstreamGenes)]]
data.gr$downstreamDist <- NA
data.gr$downstreamDist[which(!is.na(data.gr$downstreamGenes) & data.gr$downstreamOrient == "+")] <- distance(data.gr[which(!is.na(data.gr$downstreamGenes) & data.gr$downstreamOrient == "+")],
                                                                                                             gene_st.gr[data.gr$downstreamGenes[which(!is.na(data.gr$downstreamGenes) & data.gr$downstreamOrient == "+")]])
data.gr$downstreamDist[which(!is.na(data.gr$downstreamGenes) & data.gr$downstreamOrient == "-")] <- distance(data.gr[which(!is.na(data.gr$downstreamGenes) & data.gr$downstreamOrient == "-")],
                                                                                                             gene_en.gr[data.gr$downstreamGenes[which(!is.na(data.gr$downstreamGenes) & data.gr$downstreamOrient == "-")]])
data.gr$nbOrient <- NA
data.gr$nbOrient[which(data.gr$upstreamOrient=="-" & data.gr$downstreamOrient=="+")] <- "divergent" # & data.gr$downstreamDist < 50000
data.gr$nbOrient[which(data.gr$upstreamOrient=="+" & data.gr$downstreamOrient=="-")] <- "convergent" # & data.gr$downstreamDist < 50000
data.gr$nbOrient[which(data.gr$upstreamOrient=="+" & data.gr$downstreamOrient=="+")] <- "tandem" # & data.gr$downstreamDist < 50000
data.gr$nbOrient[which(data.gr$upstreamOrient=="-" & data.gr$downstreamOrient=="-")] <- "tandem" # & data.gr$downstreamDist < 50000

data.gr$max_nb_dist <- data.frame(data.gr)[,c(51,54)] %>% apply(.,1,max) %>% round_any(.,5000)
data.gr$round_nb <- NA
data.gr$round_nb[data.gr$max_nb_dist < 5000] <-  "<5"
data.gr$round_nb[data.gr$max_nb_dist >= 5000 & data.gr$max_nb_dist < 20000] <-  "5-20"
data.gr$round_nb[data.gr$max_nb_dist >= 20000 & data.gr$max_nb_dist < 50000] <-  "20-50"
data.gr$round_nb[data.gr$max_nb_dist >= 50000] <-  ">50"

# plot fractions
mark_df<- data.frame(table(data.gr$mark[data.gr$mark=="Ok" & data.gr$max_nb_dist < 50000 & data.gr$annotation_simple == "intergenic"]))
IZ_df<- data.frame(table(data.gr$TADborder[data.gr$mark=="Ok" & data.gr$max_nb_dist < 50000 & data.gr$annotation_simple == "intergenic"]))

orient_df<- data.frame(table(data.gr$nbOrient[data.gr$mark=="Ok" & data.gr$max_nb_dist < 50000 & data.gr$annotation_simple == "intergenic"]),
                       mark="Ok")
orient_df$Frac <- orient_df$Freq / mark_df$Freq

orientIZ_df <- data.frame(table(data.gr$nbOrient[data.gr$mark=="Ok" & data.gr$TADborder & data.gr$max_nb_dist < 50000 & data.gr$annotation_simple == "intergenic"]),
                          mark = "IZ")
orientIZ_df$Frac <- orientIZ_df$Freq / IZ_df$Freq[2]

enrich_df <- data.frame(orientIZ_df,
                        enrichment = log2(orientIZ_df$Frac / orient_df$Frac))
orient_df <- rbind(orient_df,orientIZ_df)
ggplot(orient_df,aes(Var1,Frac)) +
  geom_bar(aes(fill=mark),stat="identity",position="dodge") +
  theme_minimal() +
  scale_fill_manual(values = my.pal[c(7,12,13)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
ggsave("Ok_TADborder_annotations_Div.Con.Tan_50kb.max.dist_intergenic.pdf")

# gene annotation
mark_df<- data.frame(table(data.gr$mark[data.gr$mark=="K5par"]))
IZ_df<- data.frame(table(data.gr$RFD.extreme[data.gr$mark=="K5par"]))

orient_df<- data.frame(table(data.gr$annotation_simple[data.gr$mark=="K5par"]), mark="K5par")
orient_df$Frac <- orient_df$Freq / mark_df$Freq

orientIZ_df <- data.frame(table(data.gr$annotation_simple[data.gr$mark=="K5par" & !is.na(data.gr$RFD.extreme)]),mark = "extreme")
orientIZ_df$Frac <- orientIZ_df$Freq / sum(IZ_df$Freq)

orient_df <- rbind(orient_df,orientIZ_df)
ggplot(orient_df,aes(Var1,Frac)) +
  geom_bar(aes(fill=mark),stat="identity",position="dodge") +
  theme_minimal() +
  scale_fill_manual(values = my.pal[c(7,12,13)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
ggsave("K5par_extreme_gene_annotations.pdf")

ggplot(data.frame(subset(data.gr,mark %in% c("K36m16","K36"))),aes(RFD)) +
  geom_density(aes(colour=annotation_simple,linetype=mark)) +
  #geom_freqpoly(aes(colour=annotation_simple,linetype=mark),bins=100) +
  #geom_histogram(aes(colour=annotation_simple,fill=annotation_simple),position="dodge") +
  theme_minimal() +
  #facet_wrap(~mark,ncol=1) +
  #coord_cartesian(xlim=c(-0.5,0.5)) +
  scale_colour_manual(values = my.pal[c(3,4,1,6,8)]) +  scale_fill_manual(values = my.pal[c(3,4,1,6,8)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
ggsave("K36_WT_316mut_gene_annotations_density.pdf")

t.test(data.gr$RFD[data.gr$mark=="K36m13" & data.gr$annotation_simple == "gene + TRUE"],data.gr$RFD[data.gr$mark=="K36m13" & data.gr$annotation_simple == "gene + TRUE"])


##
nbTotal <- data.frame(table(data.gr$nbOrient,data.gr$RNA_active))
nbDist <- data.frame(table(data.gr$nbOrient[data.gr$IZ],data.gr$RNA_active[data.gr$IZ]))
nbDist$Fraq <- nbDist$Freq / nbTotal$Freq[match(nbDist$Var1,nbTotal$Var1)]
ggplot(nbDist,aes(Var1,Fraq)) +
  geom_bar(aes(fill= Var2),stat="identity",position="dodge") +
  theme_minimal() +
  scale_fill_manual(values = my.pal[c(7,12)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
ggsave("IZ_annotations_Div.Con.Tan_FracOfTotal_RNAactiveCol.pdf")

#chr.size <- 195471971/1000
#annot.table <- data.frame(table(annot.df$mark,annot.df$annotation_cur))
#annot.table$freq <- (annot.table$Freq / chr.size ) 
#
#ggplot(annot.df[annot.df$geneStrand=="1",],aes(annotation_cur,RFD.raw)) +
#  geom_hline(yintercept = 0,colour="black",size=0.2) + 
#  geom_boxplot(aes(colour=annotation_cur),outlier.shape = NA) +
#  #geom_bar(aes(fill=Var2),stat="identity",position = "dodge") +
#  facet_wrap(~mark,scale="free",ncol=3) +
#  coord_flip() +
#  theme_minimal() +
#  #theme(axis.text.x = element_text(angle = 0, hjust = 1,vjust=0.5)) +
#  scale_colour_manual(values = my.pal)
#ggsave("Annotation_box_outlier_genePlus_RFD.raw_mark_CPM0.3_flip_freeY_notch.pdf",height = 11,width = 12)
#
## test for significant RFD differencies btw. plus and minus strand
#pval_df <- data.frame(matrix(data=NA,
#                             nrow=length(unique(annot.df$mark)),
#                             ncol=length(unique(annot.df$annotation_cur)),
#                             dimnames=list(unique(annot.df$mark),
#                                           unique(annot.df$annotation_cur))))
#colnames(pval_df) <- gsub("\\.","", colnames(pval_df))
#colnames(pval_df) <- gsub(" ","", colnames(pval_df))
#
#estimate_df <- pval_df
#for(mark_i in unique(annot.df$mark)) {

for(annot_i in unique(annot.df$annotation_cur)) {
  data_df <- subset(annot.df,mark %in% mark_i & annotation_cur %in% annot_i)
  res <-  t.test(data_df$RFD.raw[data_df$geneStrand==1],
                 data_df$RFD.raw[data_df$geneStrand==2])
  pval_df[mark_i,annot_i] <- res$p.value
  estimate_df[mark_i,annot_i] <- res$statistic
}

#estimate_df$mark <- rownames(estimate_df)
#estimate_m <- melt(estimate_df)
#
#pval_df$mark <- rownames(pval_df)
#pval_m <- melt(pval_df)
#pval_m$statistic <- estimate_m$value
#ggplot(data=subset(pval_m,variable!="DistalIntergenic"),aes(variable,-log10(value))) +
#  geom_bar(stat="identity",position="dodge",aes(fill=statistic>0)) +
#  scale_fill_manual(values = my.pal) +
#  coord_flip() +
#  geom_hline(yintercept = -log10(0.01)) +
#  facet_wrap(~mark,scale="free_x",ncol=3) +
#  theme_minimal() +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
#ggsave("Repartition.raw_strand.diff_t.test_p.value_signStatistic_flip_freeX.pdf",height = 10)

# chromHMM
chromHMM <- read.table(file.path(server,data_dir,"extern_data/ChromHMM/mm9_ENCFF376DPH.bed.gz"))
#chromHMM <- read.table(file.path(server,data_dir,"extern_data/ChromHMM/mESC_E14_12_dense.annotated.bed.gz"),skip=1)
colnames(chromHMM) <- c("chr","start","end","seqment","score","strand","thickStart","thickEnd","edge")

chromHMM$segment_annot <- gsub("^1$","K4m3",chromHMM$seqment) %>%
  gsub("^2$","K4m1/3",.,perl = TRUE) %>%
  gsub("^3$","K4m1",.,perl = TRUE) %>%
  gsub("^4$","K4m1+K36m3",.,perl = TRUE) %>%
  gsub("^5$","K36m3",.,perl = TRUE) %>%
  gsub("^6$","Unmarked",.,perl = TRUE) %>%
  gsub("^7$","K27m3",.,perl = TRUE)

chromHMM.gr <- makeGRangesFromDataFrame(chromHMM,
                                        seqnames.field="chr",
                                        start.field = "start",
                                        end.field = "end",
                                        keep.extra.columns = T)
chain <- import.chain("/isdata/alab/people/maria/genome/mm9ToMm10.over.chain")
chromHMM.gr <- unlist(liftOver(chromHMM.gr,chain))

overlap <- findOverlaps(data.gr,chromHMM.gr)
data.gr$seqment <- data.gr$segment_annot <- NA
data.gr$seqment[queryHits(overlap)] <- chromHMM.gr$seqment[subjectHits(overlap)]
data.gr$segment_annot[queryHits(overlap)] <- chromHMM.gr$segment_annot[subjectHits(overlap)]

# DHS peaks
DHS <- read.table(file.path(server,data_dir,"extern_data/Dnase_UW/rep1_rep2_hotspot_broadPeak.bed"))
DHS.gr <- makeGRangesFromDataFrame(DHS,
                                   seqnames.field="V1",
                                   start.field = "V2",
                                   end.field = "V3",
                                   keep.extra.columns = T)
overlap <- findOverlaps(data.gr,DHS.gr)
data.gr$DHS <- FALSE
data.gr$DHS[queryHits(overlap)] <- TRUE

CTCF_bed <- read.table(file.path(server,data_dir,"extern_data/B_Bonev/CHIP_seq/GSE96107_ES_CTCF.IDR0.05.filt_ID.narrowPeak"),sep="\t")
CTCF_bed <- CTCF_bed[,c(1:4)]
colnames(CTCF_bed) <- c("chr","start","end","ID")
CTCF.gr <- makeGRangesFromDataFrame(CTCF_bed,
                                    seqnames.field = "chr",,
                                    start.field = "start",
                                    end.field = "end",
                                    keep.extra.columns = T)
CTCF.gr <- subset(CTCF.gr,seqnames %in% myseqnames)

overlap <- findOverlaps(data.gr,CTCF.gr)
data.gr$CTCF <- FALSE
data.gr$CTCF[queryHits(overlap)] <- TRUE


#### Hi-C DI ####
DI <- read.csv(file.path(server,data_dir,"extern_data/B_Bonev/HiC_10kb/B.Bonev_10Kb_res_DI_1Mb.txt"),sep="\t",as.is = T,header = F)
colnames(DI) <- c("chr","end","A","B","DI")
DI.gr <- makeGRangesFromDataFrame(DI,
                                  seqnames.field = "chr",
                                  start.field = "end",
                                  end.field = "end",
                                  keep.extra.columns = T)
start(DI.gr) <- start(DI.gr) - 10000

overlap <- findOverlaps(data.gr,DI.gr,minoverlap = 2)
data.gr$DI <- NA
data.gr$DI[queryHits(overlap)] <- DI.gr$DI[subjectHits(overlap)]

#### TADs ####
# Signal within TADs (figure 2A, S8)
# Coordinate can be changed to U-domain coordingate (figure S4A)
TAD <- read.csv(file.path(server,data_dir,"extern_data/B_Bonev/HiC_meta/mmc2.csv"),sep=";",as.is = T)
TAD <- TAD[TAD$end!=0,]
TAD <- TAD[!is.na(TAD$chrom),]

TAD.gr <- makeGRangesFromDataFrame(TAD,
                                   seqnames.field = "chrom",
                                   start.field = "start",
                                   end.field = "end",
                                   keep.extra.columns = T)
TAD.gr <- subset(TAD.gr, seqnames %in% myseqnames)

# Signal within TADs
overlap_pairs <- findOverlaps(TAD.gr,data.gr)
TAD.data.gr <- data.gr[subjectHits(overlap_pairs)] # take all data within TADs
TAD.data.gr$tad_width <- width(TAD.gr)[queryHits(overlap_pairs)] # TAD width
TAD.data.gr$compartment <- TAD.gr$Compartment[queryHits(overlap_pairs)] # TAD compartment

# distance to TAD border
TAD.data.gr$dist <- TAD.data.gr$mid - start(TAD.gr)[queryHits(overlap_pairs)]
TAD.data.gr$dist <- round_any(TAD.data.gr$dist,1000)

# relative distance to TAD border
TAD.data.gr <- TAD.data.gr[-which(TAD.data.gr$dist>TAD.data.gr$tad_width),] # restrict to bins within TADs
TAD.data.gr$rel.dist <- TAD.data.gr$dist / TAD.data.gr$tad_width # relative dist. = dist. / TAD width
TAD.data.gr$rel.dist <- round_any(TAD.data.gr$rel.dist,0.01) # round to clear signal

# summarize signal pr. binned relative distance pr. mark
rel.dist <- data.frame(TAD.data.gr %>% data.frame %>%
                         dplyr::group_by(rel.dist,mark) %>% 
                         dplyr::summarise(#sdRFD = sd(RFD),
                           RFD = mean(RFD,rm.na=T), # mean signal across TAD
                           RFD.raw = mean(RFD.raw,rm.na=T), # mean raw signal across TAD
                           #DI = mean(DI,na.rm=T), # trim = 0.5 # mean directionalty index
                           CPM = mean(CPM,rm.na=T))) # mean CPM counts


ggplot(data=subset(rel.dist,mark %in% c("Ok","K20m13","K20m16","K36m13","K36m16","K5m13","K5m16","K5","K20","K36")),aes(rel.dist,RFD)) + 
  #geom_vline(aes(xintercept=rel.dist,colour=IZ_frac)) +
  #geom_bar(aes(fill=mark),stat="identity",position="dodge",size=0.8) + # colour=mark
  geom_line(aes(colour=mark,linetype=mut),size=0.8) + # colour=mark,alpha=mark ,
  scale_alpha_manual(values=c(1,0.4)) + scale_linetype_manual(values=c(2,1,3)) +
  facet_wrap( ~ PTM,scale="free_y",ncol=2) +
  #coord_cartesian(ylim=c(-0.015,0.015)) +
  ##scale_color_brewer(palette = "Paired") +
  #scale_colour_distiller(palette="Blues",direction=1) + xlab("position") +
  scale_colour_manual(values=c(my.pal[c(2,2,2,1,1,1,12,12,12,6)])) +
  theme_minimal() 
# ylab("Smoothed repartition") + xlab("Relative distance within TADs")
ggsave("rel.TAD.dist0.01_line_RFD.smooth_wK5mut_WT_freeY.pdf")


#### heatmap TADs ####
TAD.gr$ID <- paste(seqnames(TAD.gr),start(TAD.gr),sep=":") %>% paste(.,end(TAD.gr),sep="-")
#U.gr$ID <- paste(seqnames(U.gr),start(U.gr),sep=":") %>% paste(.,end(U.gr),sep="-")

#TAD_A.gr <- subset(TAD.gr,Compartment=="A")
TAD.ext.gr <- resize(TAD.gr,2200000,fix="center") 
TAD.ext.gr$TAD_mid <- mid(ranges(TAD.ext.gr))
TAD.ext.gr$TAD_mid <- round_any(TAD.ext.gr$TAD_mid,1000)

overlap_pairs <- findOverlaps(TAD.ext.gr, data.gr)
data.TAD.gr <- data.gr[subjectHits(overlap_pairs)]
data.TAD.gr$TAD_ID <- TAD.ext.gr$ID[queryHits(overlap_pairs)]

data.TAD.gr$dist <- start(data.TAD.gr) - TAD.ext.gr$TAD_mid[queryHits(overlap_pairs)]
data.TAD <- data.frame(data.TAD.gr)

for(TAD_mark in unique(data.TAD$mark)) {
  TAD_mat <- dcast(subset(data.TAD,mark==TAD_mark),TAD_ID ~ dist,value.var = "RFD",fun.aggregate = mean) 
  rownames(TAD_mat) <- TAD_mat$TAD_ID
  TAD_mat <- TAD_mat[,-1]
  
  #TAD_order <- TAD.gr$ID[order(width(TAD.gr),decreasing = F)]
  TAD_order <- U.gr$ID[order(width(U.gr),decreasing = F)]
  TAD_mat_order <- TAD_mat[match(TAD_order,rownames(TAD_mat)),]
  
  # filter out top outliers
  high_qt <- quantile(TAD_mat_order,0.9999,na.rm=T)
  low_qt <- quantile(TAD_mat_order,0.0001,na.rm=T)
  
  # remove entire row
  #TAD_mat_order <- TAD_mat_order[apply(TAD_mat_order,1,function(x) max(x,na.rm=T))<high_qt,]
  #TAD_mat_order <- TAD_mat_order[apply(TAD_mat_order,1,function(x) min(x,na.rm=T))>low_qt,]
  
  # set NA for outliers
  TAD_mat_order[TAD_mat_order>high_qt] <- NA
  TAD_mat_order[TAD_mat_order<low_qt] <- NA
  
  #apply(TAD_mat_order,1,function(x) length(which(is.na(x)))) %>% summary
  qt <- apply(TAD_mat_order,1,function(x) length(which(is.na(x)))) %>% quantile(0.9)
  TAD_plot_mat <- TAD_mat_order[apply(TAD_mat_order,1,function(x) length(which(is.na(x)))<=qt),]
  
  breaksList = c(seq(-max(abs(TAD_plot_mat),na.rm=T),0, by = 0.005),rev(seq(-max(abs(TAD_plot_mat),na.rm=T),0, by = 0.005)*-1))
  #breaksList = seq(-max(abs(TAD_plot_mat),na.rm=T), max(abs(TAD_plot_mat),na.rm=T), by = 0.01)
  colList <- (length(breaksList)/2)+8
  col_pal <- c(colorRampPalette(rev(brewer.pal(n=9,name="Blues")))(colList)[-c((colList-7):colList)],
               rev(colorRampPalette(rev(brewer.pal(n=9,name="Reds")))(colList)[-c((colList-7):colList)]))
  
  heat.res <- pheatmap(TAD_plot_mat, #scale="row",
                       color=col_pal,
                       cluster_rows = F, cluster_cols = F,
                       breaks = breaksList,
                       filename = paste0(TAD_mark,"_heatmap_9qt_0.9999qt_setNA_filt_dist2.2Mb_Udomains_CPM0.3Filt_v2.pdf"),
                       show_rownames = F, show_colnames = F)
  
}

#### FANTOM5 TSSs ####
load("/seqdata/sandelin/projects/FANTOM5_mouse_enhancerome/0_dpi_tss_clusters/mm10.TC.exp_matrix_tpm.Rdata")
rownames(mm10TcExpMatrixTpm) <- mm10TcExpMatrixTpm$id
mm10TcExpMatrixTpm <- mm10TcExpMatrixTpm[,-1]

annot <- read.table("/seqdata/sandelin/projects/FANTOM5_mouse_enhancerome/0_dpi_tss_clusters/annotated_mm10_TC.cage_peak",as.is=T)
colnames(annot) <- c("chr","start","end","name","score","strand","peakStart","peakEnd","tssAnnoation","transcript_id","gene_id","gene_name")
annot_tss <- annot[which(annot$tssAnnoation %in% c("primary_tss","alternative_tss")),]

ES_data <- mm10TcExpMatrixTpm[,colnames(mm10TcExpMatrixTpm) %in% c("CNhs14098","CNhs14099","CNhs14100","CNhs14091","CNhs14092","CNhs14093")]

# expressed above thres 
above_one <- rowSums(ES_data > 1)
ES_data$activity <- ifelse(above_one>2,"active","inter")
ES_data$activity[above_one==0] <- "inactive"

#ES_data <- subset(ES_data, above_one > 3)
ES_tss_data <- ES_data[which(rownames(ES_data) %in% annot_tss$name),]
ES_annot <- annot_tss[which(annot_tss$name %in% rownames(ES_tss_data)),]
ES_annot$gene_id <- sapply(strsplit(ES_annot$gene_id,"\\."),"[[",1)
ES_annot$activity <- ES_data$activity[match(ES_annot$name,rownames(ES_data))]
#save(ES_annot,file="extern_data/FANTOM5/ES_annot.RDat")
load("extern_data/FANTOM5/ES_annot.RDat")

#ES_tss_data$strand <- ES_annot$strand[match(rownames(ES_tss_data),ES_annot$name)]
#ES_tss_data$name <- rownames(ES_tss_data)
#ES_tss_m <- melt(ES_tss_data)
#data.frame(ES_tss_m %>% 
#             dplyr::group_by(strand,variable) %>%
#             dplyr::summarise(mean = mean(value)))
#ggplot(ES_tss_m,aes(variable,log10(value))) + 
#  geom_boxplot(aes(colour=strand)) + 
#  scale_colour_manual(values = my.pal[c(1,6)]) +
#  theme_minimal()
#ggsave("mEScells_F5.strand.diff.pdf",width = 4,height = 3.5)

ES_annot_act <- ES_annot[ES_annot$activity=="active",c("chr","peakStart","peakEnd","name","score","strand")]
ES_annot_act$peakStart <- ES_annot_act$peakStart - 300
ES_annot_act$peakEnd <- ES_annot_act$peakEnd + 299

ES_annot_act <- ES_annot_act[ES_annot_act$peakStart>0,]
write.table(ES_annot_act,file=file.path("ECS_1tpm_3lib_cage.TSS_600bp.bed"),quote = F,sep="\t",col.names = F,row.names = F)
ES_annot.gr <- makeGRangesFromDataFrame(ES_annot,seqnames.field = "chr",
                                        start.field = "start",
                                        end.field = "end", 
                                        keep.extra.columns = T)
start(ES_annot.gr) <- ES_annot.gr$peakStart
end(ES_annot.gr) <- ES_annot.gr$peakEnd
ES_annot.gr <- promoters(ES_annot.gr,upstream = 500, downstream = 200)
#save(ES_annot.gr,file="extern_data/FANTOM5/ES_annot_GRange_up500_down200.RDat")
#load("extern_data/FANTOM5/ES_annot_GRange_up500_down200.RDat")

# subset according to early transcription
gene_set <- read_delim(file.path(server,data_dir,"extern_data/GO_term_summary_20180226_072440.txt"),"\t", escape_double = FALSE, trim_ws = TRUE) %>% as.data.frame
ES_early <- ES_annot[ES_annot$gene_name %in% gene_set$Symbol,]

#ES_matrix$cell_cycle <- FALSE
#ES_matrix$cell_cycle[which(ES_matrix$ID %in% ES_early$name)] <- "TRUE"

#### FANTOM5 enhancers ####
load("/seqdata/sandelin/projects/FANTOM5_mouse_enhancerome/0_data_freeze/enhancers/mm10.enhancers.wRep.freeze.robust.exp_matrix.Rdata")
#Enh_matrix <- read.table(file.path(server,data_dir,"/extern_data/FANTOM5/F5.mm10.enhancers.expression.matrix.gz"))
Enh_matrix <- read.table("/seqdata/sandelin/projects/FANTOM5_mouse_enhancerome/0_permissive_enhancers/mm10/permissive_enhancers_v2/F5.mm10.enhancers.expression.tpm.matrix")
Enh_annot <- read.table(file.path(server,data_dir,"/extern_data/FANTOM5/mm10.enhancers_full_freeze.sort.bed"))
colnames(Enh_annot) <- c("chr","start","end","name","score","strand","peakStart","peakEnd","Col9","Col10","Col11","Col12")

# use enhancers above noise thorshold
Enh_matrix <- Enh_matrix[which(rownames(Enh_matrix) %in% rownames(enh_robust_count_wRep)),]

# enhancer expression activity
Enh_matrix <- Enh_matrix[,colnames(Enh_matrix) %in% c("CNhs14098","CNhs14099","CNhs14100","CNhs14091","CNhs14092","CNhs14093")] # use only mES enhancers
above_one <- rowSums(Enh_matrix > 0.5) # TPM threshold
Enh_matrix$activity <- ifelse(above_one>2,"active","inter") # active: more than 2 samples above 0.5 TPM
Enh_matrix$activity[above_one==0] <- "inactive" # inactive: no samples above 0.5 TPM

Enh_annot$active <- Enh_matrix$activity[match(Enh_annot$name, rownames(Enh_matrix))]
Enh.gr <- makeGRangesFromDataFrame(Enh_annot,seqnames.field = "chr",
                                   start.field = "start",
                                   end.field = "end", keep.extra.columns = T)

Enh.gr <- subsetByOverlaps(Enh.gr,ES_annot.gr,invert=T) # exclude annotated promoters
Enh.gr <- subset(Enh.gr,!is.na(Enh.gr$active))
colnames(elementMetadata(Enh.gr))[1] <- "ID"
#Enh_annot <- subset(Enh_annot,name %in% Enh.gr$name)
#save(Enh.gr,file="extern_data/FANTOM5/F5_enh_robust_nonTSSoverlap_GRange.RDat")
load("extern_data/FANTOM5/F5_enh_robust_nonTSSoverlap_GRange.RDat")


#### footprints ####
CTCF_bed <- read.table(file.path(server,data_dir,"extern_data/B_Bonev/CHIP_seq/GSE96107_ES_CTCF.IDR0.05.filt_ID.narrowPeak"),sep="\t")
CTCF_bed <- CTCF_bed[,c(1:4)]
colnames(CTCF_bed) <- c("chr","start","end","ID")
CTCF.gr <- makeGRangesFromDataFrame(CTCF_bed,
                                    seqnames.field = "chr",
                                    start.field = "start",
                                    end.field = "end",
                                    keep.extra.columns = T)
CTCF.gr <- subset(CTCF.gr,seqnames %in% myseqnames)
#CTCF.gr <- subsetByOverlaps(CTSF.gr, ES_annot.gr,invert=T,maxgap = 200)
#CTCF.gr <- reduce(CTSF.gr)
#CTCF.gr$ID <- paste(seqnames(CTSF.gr),start(CTSF.gr),sep=":") %>% paste(.,end(CTSF.gr),sep="-")

# only CTCFs within TADs
CTSF.gr <- subsetByOverlaps(CTSF.gr,TAD.gr,minoverlap = 2)
overlap <- findOverlaps(CTSF.gr,TAD.gr,minoverlap = 2)

#TADdist <- distanceToNearest(CTSF.gr,TAD.start.gr)
CTSF.gr$TADstart <- start(TAD.gr)[subjectHits(overlap)]
CTSF.gr$TADend <- end(TAD.gr)[subjectHits(overlap)]
CTSF.gr$TADwidth <- width(TAD.gr)[subjectHits(overlap)]
CTSF.gr$TADno <- TAD.gr$No[subjectHits(overlap)]
CTSF.gr$Compartment <- TAD.gr$Compartment[subjectHits(overlap)]
CTSF.gr$Compartment <- TAD.gr$Compartment[subjectHits(overlap)]
CTSF.gr$start.dist <- start(CTSF.gr) - CTSF.gr$TADstart
CTSF.gr$end.dist <- abs(end(CTSF.gr) - CTSF.gr$TADend)
CTSF.gr$rel.dist <- CTSF.gr$start.dist / CTSF.gr$TADwidth

TSSoverlap <- findOverlaps(CTSF.gr, ES_annot.gr)
CTSF.gr$gene_active <- "NA"
CTSF.gr$gene_active[queryHits(TSSoverlap)] <- ES_annot.gr$activity[subjectHits(TSSoverlap)]

ggplot(data.frame(CTSF.gr),aes(rel.dist)) +
  geom_histogram(bins=50) +
  theme_minimal()
#ggsave("CTCF_rel.dist.TADs.pdf")

# nearest quantiles
startID <- unique(CTSF.gr$ID[CTSF.gr$start.dist < quantile(CTSF.gr$start.dist[CTSF.gr$rel.dist<0.5],0.1)] )
midID <- unique(CTSF.gr$ID[which(CTSF.gr$rel.dist>0.475 & CTSF.gr$rel.dist<0.525)] )
endID <- unique(CTSF.gr$ID[CTSF.gr$end.dist < quantile(CTSF.gr$end.dist[CTSF.gr$rel.dist>0.5],0.1)] )
CTSF_start.gr <- unique(CTSF.gr[which(CTSF.gr$ID %in% startID)])
CTSF_mid.gr <- unique(CTSF.gr[which(CTSF.gr$ID %in% midID)])
CTSF_end.gr <- unique(CTSF.gr[which(CTSF.gr$ID %in% endID)])

# nearest
CTSF.nearest <- data.frame(data.frame(CTSF.gr) %>% 
                             dplyr::group_by(TADno) %>%
                             dplyr::summarise(start = ID[which.min(start.dist)],
                                              end = ID[which.min(end.dist)],
                                              mid = ID[which.min(abs(rel.dist-0.5))],
                                              start.dist = min(start.dist),
                                              end.dist = min(end.dist),
                                              mid.dist = min(abs(rel.dist-0.5))))

CTSF_start.gr <- unique(CTSF.gr[which(CTSF.gr$ID %in% CTSF.nearest$start)])
CTSF_mid.gr <- unique(CTSF.gr[which(CTSF.gr$ID %in% CTSF.nearest$mid)])
CTSF_end.gr <- unique(CTSF.gr[which(CTSF.gr$ID %in% CTSF.nearest$end)])

write.table(data.frame(CTSF_start.gr)[,c(1:3,6)],file=file.path(server,data_dir,"extern_data/B_Bonev/CHIP_seq/GSE96107_ES_CTCF.IDR0.05.filt_ID_startTAD_nearest10qt_new.narrowPeak"),sep="\t",quote = F,row.names = F,col.names = F)

# load bwtools matrix
ES_matrix <- fread(file.path(data_dir,"footprint/CPM_ext_mm10ES_TSS.av10.matrix"),stringsAsFactors = F,sep="\t")
colnames(ES_matrix) <- c('Sample','ID','Position','Signal')
#ES_matrix$Signal[is.na(ES_matrix$Signal)] <- 0
ES_matrix <- separate(ES_matrix, Sample, c("mark","rep","strand"), sep = "_", remove = F)
#ES_matrix$mark <- sapply(str_split(ES_matrix$Sample,"_"),"[[",1)
#ES_matrix$strand <- stri_sub(ES_matrix$Sample,-1,-1)
#ES_matrix$rep <- stri_sub(ES_matrix$sample,-4,-3)
ES_matrix$tc.strand <- stri_sub(ES_matrix$ID,-1,-1)

#ES_matrix$Compartment <- CTSF.gr$Compartment[match(ES_matrix$ID,CTSF.gr$ID)]
#ES_matrix$gene_active <- CTSF.gr$gene_active[match(ES_matrix$ID,CTSF.gr$ID)]
ES_matrix$gene_active <- ES_annot.gr$activity[match(ES_matrix$ID,ES_annot.gr$name)]

# revert Ok-seq strandness
#ES_matrix$strand[c(which(ES_matrix$mark=="Ok" & ES_matrix$strand == "F"),which(ES_matrix$mark=="Ok" & ES_matrix$strand == "R"))] = ES_matrix$strand[c(which(ES_matrix$mark=="Ok" & ES_matrix$strand == "R"),which(ES_matrix$mark=="Ok" & ES_matrix$strand == "F"))] 

signal_mean <- data.frame(subset(ES_matrix) %>% #
                            dplyr::group_by(Position,mark,strand,tc.strand) %>% # tc.strand Compartment
                            dplyr::summarise(Signal = mean(Signal,na.rm=T)))

signal_mean$mark <- factor(signal_mean$mark,levels=c("cage","Pro","Ok","input","K20","K20par","K20m13","K20m16","K36","K36m13","K36m16","K5","K5par"))
signal_mean$same_strand <- paste(signal_mean$tc.strand,signal_mean$strand,sep=" ")
signal_zoom <- signal_mean[which(abs(signal_mean$Position)<=1000 & signal_mean$mark %in% c("Ok","K5","K20","K36","K36m13","K36m16","K20m13","K20m16","input","proseq")),] # 
ggplot(signal_zoom,aes(Position,Signal)) + 
  geom_line(aes(colour=same_strand),size=0.5) +  # alpha=Compartment
  theme_minimal() + 
  scale_alpha_manual(values=c(1,0.8,0.6,1,0.8,0.6,0.4)) +
  theme(text=element_text(family="Helvetica"),legend.position="top") + 
  guides(col=guide_legend(nrow=1)) + 
  ylab("Signal (CPM)") +
  scale_color_manual(values = my.pal[c(1,2,6,5,7,8,3,13,14)]) +
  facet_wrap(~mark ,scale="free_y",ncol=3) +
  xlab("Position relative to peak")
ggsave("ECs_F5_1tpm_3lib_TSSs_freeY_2kb_CPM_ext_10av_v2.pdf")

repart_matrix <- data.frame(subset(ES_matrix, mark %in% c("Ok","K5","K20","K36","K36m13","K36m16","K20m13","K20m16","input","proseq")) %>%
                              dplyr::group_by(ID,Position,rep,mark,tc.strand) %>%
                              dplyr::summarise(RFD = (Signal[strand=="F"] - Signal[strand=="R"]) / (Signal[strand=="F"] + Signal[strand=="R"] )))

repart_mean <- data.frame(repart_matrix %>% 
                            dplyr::group_by(Position,mark,tc.strand) %>%
                            dplyr::summarise(RFD = mean(RFD,na.rm=T)))

signal_zoom <- repart_mean[abs(repart_mean$Position)<=1000,]
ggplot(signal_zoom,aes(Position,RFD)) + 
  geom_line(aes(colour=tc.strand),size=0.5) + 
  theme_minimal() + 
  scale_alpha_manual(values=c(1,0.5)) +
  theme(text=element_text(family="Helvetica"),legend.position="top") + 
  guides(col=guide_legend(nrow=2)) + 
  ylab("signal around cage TSSs") +  xlab("Position relative to TSS") +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~mark ,scale="fixed") +
  coord_cartesian(ylim=c(-0.2,0.2)) 
ggsave("ECs_F5_TSS_repartition_fixY_1kbwindow_CPM_ext_10av.pdf")

ggplot(signal_zoom,aes(tc.strand,RFD)) + 
  geom_boxplot(aes(fill=rep),size=0.5,notch=T) + 
  theme_minimal() + 
  theme(text=element_text(family="Helvetica"),legend.position="top") + 
  guides(col=guide_legend(nrow=2)) + 
  scale_fill_brewer(type = 'div', palette = 4, direction = -1) +
  ylab("signal around cage TSSs") +
  #scale_fill_brewer(palette = "Accent") +
  facet_wrap(~mark,scale="free_y") +
  xlab("Position relative to TSS")
ggsave("CPM_footprint_TSS_1000kb_window.pdf")


# enhancers
enh <- read.table(file.path(data_dir,"extern_data/F5/F5.mm10.enhancers.bed.gz"),as.is=T)
colnames(enh) <- c("chr","start","end","id","score","strand","summit_start","summit_end","col9","col10","col11","col12")

enh_mat <- read.table(file.path(data_dir,"extern_data/F5/F5.mm10.enhancers.expression.matrix.gz"),as.is = T, header = T,row.names = 1)
enhES_mat <- enh_mat[,which(colnames(enh_mat) %in% c("CNhs14098","CNhs14099","CNhs14100"))]
enhES_mat <- enhES_mat[rowSums(enhES_mat) > 2,]

enhES <- enh[which(enh$id %in% rownames(enhES_mat)),]
enhES.gr <- makeGRangesFromDataFrame(enhES,seqnames.field = "chr",start.field = "start",end.field = "end",keep.extra.columns = T)

data.pro1$ID <- paste(data.pro1$V1,data.pro1$V2,sep=":") %>% paste(.,data.pro1$V3,sep="-")
data.pro2$ID <- paste(data.pro2$V1,data.pro2$V2,sep=":") %>% paste(.,data.pro1$V3,sep="-")
data.pro <- rbind(data.pro1,data.pro2)
data.pro.gr <- makeGRangesFromDataFrame(data.pro,seqnames.field = "V1",
                                        start.field = "V2",end.field = "V3")
overlaps <- findOverlaps(data.pro.gr,enhES.gr)
enhES.gr <- enhES.gr[unique(subjectHits(overlaps))]

overlaps <- findOverlaps(data.gr,Enh.gr)
enh_data.gr <- data.gr[queryHits(overlaps)]
enh_data.gr$enh_ID <- Enh.gr$ID[subjectHits(overlaps)]
enh_data.df <- data.frame(enh_data.gr) 
ok.enh.df <- merge(enh_data.df[enh_data.df$mark=="Ok",],enh_data.df[enh_data.df$mark=="proseq",],by="ID")
ok.enh.df$pro.bin <- round_any(ok.enh.df$F.cpm.y / ok.enh.df$CPM.y,0.05)
ok.enh.df$pro.bin <- round_any(log10(ok.enh.df$CPM.y),0.1)
ggplot(ok.enh.df,aes(factor(pro.bin),RFD.x)) +
  geom_boxplot() +
  geom_smooth(aes()) + # 
  #facet_wrap(~ gene_strand.x) +
  theme_minimal()
ggsave("Ok.raw_RFD_proseq_bin_WholeGenome.pdf",width = 5,height = 4)

ggplot(ok.enh.df,aes(RFD.y,RFD.x)) +
  geom_point() +
  geom_smooth(aes()) + # 
  #facet_wrap(~ gene_strand.x) +
  theme_minimal()

overlaps <- findOverlaps(ok.gr,enhES.gr,maxgap = 2000)
length(unique(subjectHits(overlaps)))

#### simulate full RFD ####
# scale mutant partition to max OKseq RFD range to test the effect of MCM2 mutant on RFD
dataOK.gr <- data.gr[data.gr$mark=="Ok"]
maxRFD <- max(max(dataOK.gr$RFD, abs(min(dataOK.gr$RFD))))
data.gr$RFDscale <- data.gr$RFD 
data.gr$RFDscale[-which(data.gr$mark %in% c("K20m16","K36m16","K20m13","K36m13"))] <- data.gr$RFD[-which(data.gr$mark %in% c("K20m16","K36m16","K20m13","K36m13"))] / maxRFD
