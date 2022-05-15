library(dplyr)
library(scales)
library(tidyverse)
library(reshape)
library(GGally)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

limitRange <- function(data, mapping, ...) {
    ggplot(data=data, mapping=mapping, ...) +
        geom_point(...) +
        geom_smooth(method = "lm", se = FALSE)
}

library(DESeq2)
library(gplots)
library(RColorBrewer)
library(IHW)
library(GGally)
library(ggplot2)
library(vidger)
library(DEGreport)
library(genefilter)
library(ggrepel)
library(ggbeeswarm)
library(gghighlight)
library(ggpubr)
plotPCA.hk <- function (object, intgroup = "condition", ntop = 500, pc_1 = 1, pc_2 = 2, returnData = FALSE, scree = FALSE)
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, pc_1], PC2 = pca$x[, pc_2], group = group,
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc_1:pc_2]
    return(d)
  }
  if (scree) {
    xx<- barplot(round(percentVar, digits = 2)*100, names.arg=c(1:length(percentVar)),xlab="PC",ylab="% Variance",ylim=c(0,100), main="Scree Plot")
    text(x = xx, y = round(percentVar, digits = 4)*100, label = round(percentVar, digits = 4)*100, pos = 3, cex = 0.8, col = "black")
  }
  else {
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color="condition", shape="type", label = "sample")) + scale_shape_manual(values=seq(0,127)) + geom_point(size = 3) + xlab(paste0("PC",pc_1,": ", round(percentVar[pc_1] * 100, digits = 2), "% variance")) + ylab(paste0("PC",pc_2,": ", round(percentVar[pc_2] * 100, digits=2), "% variance")) + coord_fixed() + geom_text_repel(size=3)
  }
}

degPlotCluster.hk <- function (table, time, color = NULL, min_genes = 10, process = FALSE,
    points = TRUE, boxes = TRUE, smooth = TRUE, lines = TRUE,
    facet = TRUE, cluster_column = "cluster", prefix_title = "Group:")
{
    stopifnot(class(table)[1] == "data.frame")
    if (cluster_column %in% colnames(table)) {
        table[["cluster"]] = table[[cluster_column]]
    }
    if (process) {
        table <- .process(table, time, color)
    }
    if ("cluster" %in% colnames(table)) {
        counts <- table(distinct(table, genes, cluster)[["cluster"]])
        counts <- counts[counts >= min_genes]
        if (length(counts) == 0)
            stop("No clusters with min_genes > ", min_genes)
        table <- inner_join(table, data.frame(cluster = as.integer(names(counts)),
            title = paste(prefix_title, names(counts), "-variants:",
                counts), stringsAsFactors = FALSE), by = "cluster")
        table %>% arrange(table$cluster)
        tempe <- table$title[!duplicated(table$title)]
        table$title = factor(table$title, levels=tempe)
        remove(tempe)
    }
    if (is.null(color)) {
        color = "dummy"
        table[[color]] = ""
        lines = FALSE
    }
    table[["line_group"]] = paste(table[["genes"]], table[[color]])
    splan <- length(unique(table[[time]])) - 1L
    p <- ggplot(table, aes_string(x = time, y = "value", fill = color,
        color = color))
    if (boxes)
        p <- p + geom_boxplot(alpha = 0, outlier.size = 0, outlier.shape = NA,
            )
    if (points)
        p <- p + geom_point(alpha = 0.4, size = 1, position = position_jitterdodge(dodge.width = 0.9))
    if (smooth)
        p <- p + stat_smooth(aes_string(x = time, y = "value",
            group = color, color = color), se = FALSE, method = "lm",
            formula = y ~ poly(x, splan))
    if (lines)
        p <- p + geom_line(aes_string(group = "line_group"),
            alpha = 0.1)
    if (facet)
        p <- p + facet_wrap(~title)
    p <- p + theme(strip.text = element_text(size=8), axis.text.x = element_text(angle = 90, hjust = 1,
        vjust = 0.5)) + ylab("Z-score of variant abundance") + xlab("")
    p
}

ggpairs_ext <- function(data, mapping, pts=list(), smt=list(), ...){
              ggplot(data = data, mapping = mapping, ...) +
                         do.call(geom_point, pts) +
                         do.call(geom_smooth, smt)
}

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), source = x)
  }))
}

appendDataFrameColumns <-function(df, prefix="", suffix="", sep="",drop="") {
  if (drop!="") {
  
  if (prefix=="") {
  newname <- paste(drop, suffix, sep=sep)
  colnames(df) <- paste(colnames(df), suffix, sep=sep)
  names(df)[names(df) == newname ] <- drop
  return(df)
  }
  
  if (suffix=="") {
    newname <- paste(prefix,drop, sep=sep)
    colnames(df) <- paste(prefix, colnames(df), sep=sep)
    names(df)[names(df) == newname ] <- drop
  return(df)
  }
  }
  else {
    if (prefix=="") {
  colnames(df) <- paste(colnames(df), suffix, sep=sep)
  return(df)
  }
    if (suffix=="") {
    colnames(df) <- paste(prefix, colnames(df), sep=sep)
  return(df)
  }
    
  }
}

# The palette with grey:
cbgPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
#scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
# scale_colour_manual(values=cbPalette)


read <- function(sg,rep,e,day) {
  file_all_count=paste0("./SGE_sg",sg,"_Rep",rep,"_Exon",e,"_Day",day,"_E",e,"_all_count.txt")
  all_count=paste0("all_sg",sg,"_E",e,"_D",day,"_R",rep)
  f_all <- read.csv(file_all_count ,sep="\t",header =FALSE)
  colnames(f_all) <- c("Count","Seq")
  f_all <- f_all %>% relocate(Seq)
  return(f_all)
  assign(all_count, f_all, envir = .GlobalEnv)
}

read_annotation <- function(sg,e) {

  file_VEP1=paste0("./Exon",e,"_sg",sg,"_combine.txt")
  vep1=paste0("Exon",e,"sg",sg,"_annotation")
  f_vep1 <- read.csv(file_VEP1 ,sep="\t",header =TRUE)
  assign(vep1, f_vep1, envir = .GlobalEnv)
  
}

selectmut <- function(counts, min.count=10, N=0.50){
 
  lib.size <- colSums(counts)
  MedianLibSize <- median(lib.size)
  CPM.Cutoff <- min.count / MedianLibSize*1e6
  CPM <- apply(counts,2, function(x) (x/sum(x))*1000000)
 
  min.samples <- round(N * ncol(counts))
  
    keep <- apply(counts, 1, function(x, n = min.samples){
    t = sum(x >= min.count) >= n
    t
  
    #keep <- apply(CPM, 1, function(x, n = min.samples){
    #t = sum(x >= CPM.Cutoff) >= n
    #t
    })
    
  ## the same as:
  #f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
  #flist <- genefilter::filterfun(f1)
  #keep <- genefilter::genefilter(CPM, flist)

  return(keep)
}

load_count <- function(sg,exon) {

out_cts = paste0("Exon",exon,"sg",sg,"_bind_filtered")
out_norm_mat = paste0("Exon",exon,"sg",sg,"_normFactors")

sg <- c(sg)
rep <-c(1,2,3)
e<- c(exon)
day <- c(4,7,11,15,21)

read_table <- expand.grid(exon=e,rep=rep,day=day, sg=sg)

E1_bind <- mapply(function(a,b,c,d) read (a,b,c,d), a=read_table$sg, b=read_table$rep, read_table$exon, read_table$day , SIMPLIFY = FALSE) %>% purrr::reduce(full_join,by=c("Seq")) %>% `colnames<-` (c("Seq", "D4R1","D4R2","D4R3","D7R1","D7R2","D7R3","D11R1","D11R2","D11R3","D15R1","D15R2","D15R3","D21R1","D21R2","D21R3")) %>% replace(is.na(.), 0)

cmt_E1 <- data.matrix(E1_bind[2:ncol(E1_bind)])
#Keep the sequence if 25% of the sample has 5 counts and above
keep.mut <- selectmut(cmt_E1, min.count=5, N=0.25)
cmt_E1_filtered <- cmt_E1[keep.mut,]

E1_bind_filtered <- E1_bind[keep.mut,]

E1_bind_sizeFactor <- apply(cmt_E1_filtered,2, function(x) (sum(x)/1000000))

#Consider starting 5% diploid and grow exponentially to 20% at day17 (==day21), 0.05 and 0.1 becoz the haploid and diploid
E1_bind_sizeFactor_d4 <- apply(cmt_E1_filtered[,1:3],2, function(x) (sum(x)/1000000))
E1_bind_sizeFactor_d7 <- apply(cmt_E1_filtered[,4:6],2, function(x) (sum(x)/(1000000/((1-0.05*exp((log(4)/17)*3)+0.1*exp((log(4)/17)*3))/1.05))))
E1_bind_sizeFactor_d11 <- apply(cmt_E1_filtered[,7:9],2, function(x) (sum(x)/(1000000/((1-0.05*exp((log(4)/17)*7)+0.1*exp((log(4)/17)*7))/1.05))))
E1_bind_sizeFactor_d15 <- apply(cmt_E1_filtered[,10:12],2, function(x) (sum(x)/(1000000/((1-0.05*exp((log(4)/17)*11)+0.1*exp((log(4)/17)*11))/1.05))))
E1_bind_sizeFactor_d21 <- apply(cmt_E1_filtered[,13:15],2, function(x) (sum(x)/(1000000/((1-0.05*exp((log(4)/17)*17)+0.1*exp((log(4)/17)*17))/1.05))))

E1_bind_sizeFactor_diploid_adj <- c(E1_bind_sizeFactor_d4,E1_bind_sizeFactor_d7,E1_bind_sizeFactor_d11,E1_bind_sizeFactor_d15,E1_bind_sizeFactor_d21)

#E1_bind_cpm <- apply(cmt_E1_filtered,2, function(x) (x/sum(x))*1000000) %>% cbind(E1_bind_filtered[1],.)

row.names(E1_bind_filtered) <- E1_bind_filtered$Seq
E1_bind_filtered$Seq <- NULL

normFactors <- matrix(E1_bind_sizeFactor, nrow=ncol(E1_bind_filtered), ncol=nrow(E1_bind_filtered), dimnames=list(1:ncol(E1_bind_filtered),1: nrow(E1_bind_filtered))) %>% t()

normFactors_adj <- matrix(E1_bind_sizeFactor_diploid_adj, nrow=ncol(E1_bind_filtered), ncol=nrow(E1_bind_filtered), dimnames=list(1:ncol(E1_bind_filtered),1: nrow(E1_bind_filtered))) %>% t()

assign(out_cts, E1_bind_filtered, envir = .GlobalEnv)
assign(out_norm_mat, normFactors, envir = .GlobalEnv)
#assign(paste0(out_norm_mat,"_adj"), normFactors_adj, envir = .GlobalEnv)
}

##DESEQ2
run_deseq2 <- function(sg,exon,sample) {
#sample is the column number, so if indlude d21, it is 15, if till d15, it is 1:12
coldata_path <- paste0("./E",exon, "_sg",sg, "_anno.txt")
prefix <- paste0("E",exon,"_sg_",sg)
cts_name<- paste0("Exon",exon,"sg",sg,"_bind_filtered")
norm_name<- paste0("Exon",exon,"sg",sg,"_normFactors")
out_rate_name<- paste0("Exon",exon,"sg",sg,"_rate_raw")
out_rate_name2<- paste0("Exon",exon,"sg",sg,"_rate_selected")
vep1=paste0("Exon",exon,"sg",sg,"_annotation")

coldata <- read.csv(coldata_path ,sep="\t",row.names=1)
cts <- as.matrix(get(cts_name))

cts <- cts[, rownames(coldata)]
normFactors <- get(norm_name)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)
normalizationFactors(dds) <- normFactors[,sample]
dds <- DESeq(dds)
rld <- rlog(dds)
z_score <- assay(rld) %>% as.matrix() %>% t() %>% scale() %>% t() %>% as.data.frame()
colnames(z_score) <- paste0(colnames(z_score), "_z_score")


table_wald <- degComps(dds, combs = "condition", alpha = 0.05, skip = FALSE, type = "apeglm", pairs = FALSE, fdr = "default")

adj <- (get(vep1) %>% select(c("Seq","Consequence")) %>% filter(Consequence == "synonymous_variant"| Consequence=="intron_variant") %>% left_join(as.data.frame(table_wald[[1]])%>%rownames_to_column(var="Seq")%>%select("Seq","log2FoldChange")))$log2FoldChange %>% median()

rate <- as.data.frame(table_wald[[1]])%>%rownames_to_column(var="Seq") %>% mutate(nochange_rate=adj) %>% mutate(adj_rate=log2FoldChange-nochange_rate) %>% mutate(adj_score=adj_rate/lfcSE) %>% dplyr::rename(raw.rate=log2FoldChange, std.err=lfcSE, raw.pvalue=pvalue, raw.padj=padj)

selected_rate <- get(vep1) %>% select(-c("X")) %>% left_join(rate,by="Seq") %>% left_join(as.data.frame(assay(rld)) %>% appendDataFrameColumns(suffix="_rld") %>% rownames_to_column(var="Seq"), by="Seq") %>% left_join(z_score %>% rownames_to_column(var="Seq"), by="Seq")

write.table(rate,sep="\t",file=paste0(prefix,"_rate_raw.txt"), col.names=NA)
write.table(selected_rate,sep="\t",file=paste0(prefix,"_rate_selected.txt"), col.names=NA)

assign(out_rate_name, rate, envir = .GlobalEnv)
assign(out_rate_name2, selected_rate, envir = .GlobalEnv)

coef_wald <- resultsNames(dds)


#Export rld
write.table(as.data.frame(assay(rld)), sep="\t",file=paste0(prefix,"_rld.txt"), col.names=NA)

#Export normalized read count
write.table(counts(dds,normalized=TRUE), sep="\t",file=paste0(prefix,"_norm_count.txt"), col.names=NA)

#Export size factor
write.table(dds@assays@data@listData %>% as.data.frame(),sep="\t",file=paste0(prefix,"_normalization_table.txt"), col.names=NA)

#Export full table
write.table(dds@rowRanges@elementMetadata@listData %>% as.data.frame() ,sep="\t",file=paste0(prefix,"_disper_table.txt"), col.names=NA)

pdf(paste0(prefix,"_dispersion.pdf"))
plotDispEsts(dds, ylim =c(1e-4,2e1))
dev.off()
  
pdf(paste0(prefix,"_scree.pdf"))
plotPCA.hk(rld,intgroup=c("condition", "type","sample"), returnData=FALSE,pc_1=1, pc_2=2, scree=TRUE)
dev.off()
    
#pdf(paste0(prefix,"_PC",1,"and_PC", 2,"_pca.pdf"))
#plotPCA.hk(rld,intgroup=c("condition", "type","sample"), returnData=FALSE,pc_1=1, pc_2=2)
#dev.off()

sampleDistMatrix <- as.matrix( dist( t( assay(rld) ) ) )
rownames(sampleDistMatrix) <- paste( rld$condition,rld$sample,rld$type, sep="-" )
colnames(sampleDistMatrix) <- NULL

pdf(paste0(prefix,"_heatmap.pdf"))
heatmap.2(sampleDistMatrix, trace="none", col=colorRampPalette(rev(brewer.pal(9, "Blues")) )(255), adjRow = c(1,1))
dev.off()
remove(sampleDistMatrix)

###Scatter plot, checking replicate consistency

pdf(paste0(prefix,"_scatter_matrix.pdf"), width=8, height=8)

print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D4")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))

print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D7")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))

print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D11")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count,Clone 5M condition DDX3X Ex14 sg2" ) + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))

print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D15")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))

dev.off()

out_dds = paste0("Exon",exon,"sg",sg,"_dds")
assign(out_dds, dds, envir = .GlobalEnv)

}

#Run the real thing sg1. Choose the Exon that you want to use.
lapply (c(1,2,4,5,6,7,8,10,11,12,13,14,15,16,17), function(x) {
read_annotation(sg=1,e=x)
load_count(sg=1,exon=x)
#sample is a vector that indicate which row to be include in the raw count file. It might be a bit confusing. Require manual checking for the column included which should be the same as the annotation (coldata) file. sample c(13:15) are the Day 15 sample. Therefore, by setting c(1:12), we exclude the Day 15 samples
run_deseq2(sg=1,exon=x,sample=c(1:12))
})

lapply (c(3), function(x) {
read_annotation(sg=1,e=x)
load_count(sg=1,exon=x)
run_deseq2(sg=1,exon=x,sample=c(2:12))
})

lapply (c(9), function(x) {
read_annotation(sg=1,e=x)
load_count(sg=1,exon=x)
run_deseq2(sg=1,exon=x,sample=c(1,3:12))
})

#Combine rate_selected

full_table1 <- mget(ls(pattern = "*sg1_rate_selected")) %>% bind_rows()

#sg2, Choose the Exon that you want to use.
lapply (c(1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17), function(x) {
read_annotation(sg=2,e=x)
load_count(sg=2,exon=x)
run_deseq2(sg=2,exon=x,sample=c(1:12))
})


lapply (c(9), function(x) {
read_annotation(sg=2,e=x)
load_count(sg=2,exon=x)
run_deseq2(sg=2,exon=x,sample=c(1,4:12))
})

#Combine rate_selected

full_table2 <- mget(ls(pattern = "*sg2_rate_selected")) %>% bind_rows()

full_table <- mget(ls(pattern = "full_table*")) %>% bind_rows() %>% mutate(sgRNA=case_when(str_detect(Old_name,"sg1") ~ "sg1", str_detect(Old_name,"sg2") ~ "sg2"))

write.table(full_table ,sep="\t",file=paste0("sg1and2_rate_maintable.txt"), col.names=NA)


sg1_sg2_collapse <- full_table %>% group_by(SGE_exon_group,sg1_sg2_combined_annotation) %>% summarize(Variant_duplication= n(), Variant_Sources = paste(sgRNA, collapse = ','), PAM_status = paste(PAM_codon, collapse = ','), std_err = paste(std.err, collapse = ','), combine_rate = paste(adj_rate, collapse = ','), combine_raw_rate = paste(raw.rate, collapse = ','), combine_baseMean = paste(baseMean, collapse = ',') ) %>% ungroup() %>% separate(combine_rate, sep=",", c("sg1_rate","sg2_rate")) %>% separate(std_err, sep=",", c("sg1_err","sg2_err")) %>% separate(combine_raw_rate, sep=",", c("sg1_raw_rate","sg2_raw_rate")) %>% separate(combine_baseMean, sep=",", c("sg1_baseMean","sg2_baseMean")) %>% mutate(sg2_err=case_when(Variant_duplication=="1" & Variant_Sources == "sg2" ~ sg1_err, TRUE ~ sg2_err)) %>% mutate(sg1_err=case_when(Variant_duplication=="1" & Variant_Sources == "sg2" ~ "NA", TRUE ~ sg1_err))%>% mutate(sg2_rate=case_when(Variant_duplication=="1" & Variant_Sources == "sg2" ~ sg1_rate, TRUE ~ sg2_rate)) %>% mutate(sg1_rate=case_when(Variant_duplication=="1" & Variant_Sources == "sg2" ~ "NA", TRUE ~ sg1_rate))%>% mutate(sg2_raw_rate=case_when(Variant_duplication=="1" & Variant_Sources == "sg2" ~ sg1_raw_rate, TRUE ~ sg2_raw_rate)) %>% mutate(sg1_raw_rate=case_when(Variant_duplication=="1" & Variant_Sources == "sg2" ~ "NA", TRUE ~ sg1_raw_rate)) %>% mutate(sg2_baseMean=case_when(Variant_duplication=="1" & Variant_Sources == "sg2" ~ sg1_baseMean, TRUE ~ sg2_baseMean)) %>% mutate(sg1_baseMean=case_when(Variant_duplication=="1" & Variant_Sources == "sg2" ~ "NA", TRUE ~ sg1_baseMean)) %>% mutate(sg1_rate = as.numeric(sg1_rate)) %>% mutate(sg2_rate = as.numeric(sg2_rate)) %>% mutate(sg1_err = as.numeric(sg1_err)) %>% mutate(sg2_err = as.numeric(sg2_err)) %>% mutate(sg1_raw_rate = as.numeric(sg1_raw_rate)) %>% mutate(sg2_raw_rate = as.numeric(sg2_raw_rate)) %>% mutate(sg1_baseMean = as.numeric(sg1_baseMean)) %>% mutate(sg2_baseMean = as.numeric(sg2_baseMean))

sg1_sg2_combined2 <-sg1_sg2_collapse %>% mutate(weight1=case_when(str_detect(PAM_status,"sg1PAM")~ 0, TRUE~ 1/(sg1_err)^2)) %>% mutate(weight2= case_when(str_detect(PAM_status,"sg2PAM")~ 0, TRUE~1/(sg2_err)^2)) %>% mutate(sum_of_weigth=rowSums(cbind(weight1,weight2), na.rm=TRUE)) %>% mutate(SE_bind = (sum_of_weigth)^(-0.5)) %>% mutate(weighted_sg1rate= weight1*sg1_rate) %>% mutate(weighted_sg2rate= weight2*sg2_rate) %>% mutate(sum_of_weighed_rate=rowSums(cbind(weighted_sg1rate,weighted_sg2rate), na.rm=TRUE)) %>% mutate(combined_rate=sum_of_weighed_rate/sum_of_weigth) %>% filter(combined_rate != "NaN") %>% mutate(combined_Z=combined_rate/SE_bind) %>% mutate(two_tailed_p= pnorm(abs(combined_Z),lower.tail = FALSE) *2) %>% mutate(BH_FDR = p.adjust(two_tailed_p, method = "BH"))



read_all_anno <- function(e) {
  
  file_VEP1=paste0("./combine/Exon",e,"_all_combine.txt")
  vep1=paste0("combine_annotation_",e)
  f_vep1 <- read.csv(file_VEP1 ,sep="\t",header =TRUE)
  assign(vep1, f_vep1, envir = .GlobalEnv)
  
}

lapply(c(1:17), function(x) read_all_anno(x))

full_annotation <- mget(ls(pattern = "combine_annotation*")) %>% bind_rows()

sg1_sg2_combined_annotated2 <- full_annotation %>% left_join(sg1_sg2_combined2, by=c("SGE_exon_group","sg1_sg2_combined_annotation"))

clean_sg1_sg2_combined_annotated2 <- sg1_sg2_combined_annotated2 %>% filter(BH_FDR != "NA")

#Main output table for the subsquent analysis
write.table(clean_sg1_sg2_combined_annotated2 ,sep="\t",file=paste0("sg1and2_combined2.txt"), col.names=NA)
