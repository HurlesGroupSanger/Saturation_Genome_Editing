library(tidyverse)
library(scales)
library(reshape2)
library(ggpubr)
library(gridExtra)

#Load the DDX3X supplement table (supp 5)
total<-read.csv("~/Downloads/Table S6_SGE_variants_annotation.txt", sep='\t', header=T) 

# A clean-up version code

#Plot Confidence
plot_nuc_conf <- function(df, exon) {
  #Levels
  nuc_level = c("T","G","C","A","Codon_del")
  
  #Filter the SNV only. The SNV has one base in both VCF-Ref and VCF-Alt column. Also, create a new column, "plot position".
  totalsnv <- total %>% filter(SGE_exon_group %in% c(exon)) %>% #Filter a specific exon
  filter(Variant_design_type=="snv") %>% #Filter snv only
  filter(nchar(VCF_Ref)==1) %>% filter(nchar(VCF_Alt)==1) %>% #Make sure the snv is snv
  mutate(VCF_Alt=factor(VCF_Alt, levels = nuc_level)) %>% mutate(plot_position=VCF_Hg38_position) %>% #Manipulation
  select("SGE_oligo_name","Variant_design_type","chrom","VCF_Hg38_position","VCF_Ref","VCF_Alt","Consequence","Variant_category","INTRON","EXON","SGE_exon_group","SGE_functional_classification","SGE_NDD_clinical_classification","NDD_clinical_classification_prob_var_functional","NDD_clinical_classification_prob_var_nonfunctional","plot_position") %>%
    mutate(Summary_Consequence=case_when(str_detect(Consequence,"UTR")~ "UTR",
                                         str_detect(Consequence,"splice_donor")~ "Canonical splice",
                                         str_detect(Consequence,"splice_acceptor")~ "Canonical splice",
                                         str_detect(Consequence,"synonymous")~ "Synonymous",
                                         str_detect(Consequence,"stop_gained")~ "Nonsense",
                                         str_detect(Consequence,"stop")~ "Missense",
                                         str_detect(Consequence,"start")~ "Missense",
                                         str_detect(Consequence,"missense")~ "Missense",
                                         str_detect(Consequence,"splice")~ "Splice region",                           
                                         str_detect(EXON,"-")~ "Intronic",
                                         TRUE ~ "Others")) %>% #Define Summary_Consequence
    mutate(Confidence=case_when(NDD_clinical_classification_prob_var_nonfunctional >= 0.9~ "High",
                                NDD_clinical_classification_prob_var_nonfunctional <= 0.5~ "Low",
                                TRUE ~ "Intermediate")) #Define Confidence
  
  #Filter the inframe only. Replace the VCF Alt to codon_del. Also, create a new column, "plot position".
  totalinframedel <- total %>% filter(SGE_exon_group %in% c(exon)) %>% #Filter a specific exon
    filter(str_detect(Variant_design_type,"inf")) %>% #Filter inframe only
    mutate(VCF_Alt="Codon_del",VCF_Alt=factor(VCF_Alt, levels = nuc_level))%>% #Replace VCF_Alt as Codon_del
    mutate(plot_position=VCF_Hg38_position+2) %>% #After testing, it has to be +2 so that the position align well
    select("SGE_oligo_name","Variant_design_type","chrom","VCF_Hg38_position","VCF_Ref","VCF_Alt","Consequence","Variant_category","INTRON","EXON","SGE_exon_group","SGE_functional_classification","SGE_NDD_clinical_classification","NDD_clinical_classification_prob_var_functional","NDD_clinical_classification_prob_var_nonfunctional","plot_position") %>%
    mutate(Summary_Consequence="Codon_del") %>% #Define the consequence as codon_del
    mutate(Confidence=case_when(NDD_clinical_classification_prob_var_nonfunctional >= 0.9~ "High",
                                NDD_clinical_classification_prob_var_nonfunctional <= 0.5~ "Low",
                                TRUE ~ "Intermediate")) #Define Confidence
  
  #Gather the number of x-axis component to plot
  osnv <-seq(min(totalsnv$plot_position), min(totalsnv$plot_position)+221,1) #Hard-coded, 221 is the nucleotide number in the longest exon, Exon13. The idea here is to make sure all plot as the same amount of square/tile. So that the size are the same across all plot.
  maxosnv <- max(totalsnv$plot_position) #Beyond this point, the tile/square were artificially added and can be remove. Appear as white in the plot.
  
  #Create a data frame that make all combination of each plot position and nuc_level
  snv_combi <- expand.grid(VCF_Alt=nuc_level[c(1:4)],plot_position=osnv)
  
  #Combine this all combination table with the data.
  total_snv_filtered <- totalsnv %>% right_join(snv_combi, by=c("VCF_Alt","plot_position")) %>% mutate(width=0.6) %>% #add the column as width. use to define the width of each square. SNV has a smaller width
    mutate(SGE_functional_classification=case_when(plot_position>maxosnv ~ "NULL", TRUE~SGE_functional_classification )) %>%
    mutate(Summary_Consequence=case_when(plot_position>maxosnv ~ "NULL", TRUE~Summary_Consequence )) %>%
    mutate(Confidence=case_when(plot_position>maxosnv ~ "NULL", TRUE~Confidence )) #Replace the columns to NULL whenvever it is more than maxosnv that is the designed maximum coordinate range
  
  #cdel_combi <- expand.grid(VCF_Alt="Codon_del",plot_position=osnv)
  
  #Create a dataframe that has a full matrix of codon position, including those repeated trinucleotide. A possible bug in the future: if the exon ended with a repeated trinucleotide, the last repeated codon will not be added and have to be added manually.
  inframe_full_pos <- data.frame(plot_position=seq(min(totalinframedel$plot_position), max(totalinframedel$plot_position),3))
  
  inframe_full<- full_join(totalinframedel, inframe_full_pos, by="plot_position") %>% arrange(plot_position) %>% fill( colnames(totalinframedel)[-16], .direction = "down")  %>% #Fill the empty row with the value in the previous row whenever the "plot_position" is the same as the previous row.
    mutate(width=0.6*4) #Define the width of the tile. Tested before and 4 is the best.

  # Combine the SNV, inframedel into a table. Factorize the plot_position make all the position in to 1,2,3,4... it might be easier for manipulate some plot features. (but not used directly here. It provide the idea of the range of x coordinate)
  total_snv_filtered_combi<-bind_rows(total_snv_filtered,inframe_full) %>% mutate(plot_position=as.character(plot_position)) %>% mutate(plot_position=factor(plot_position, levels = unique(plot_position))) %>% mutate(factor= as.numeric(plot_position))
  
  #The text in the exon label to show in the plot
  exon_label= paste0("Exon ",exon)
  
  fill <- c("High" = "red", "Intermediate" = "blue", "Low" = "white", "NULL"="white")
  cols <- c("UTR" = "#00BFC4", "Canonical splice" = "#F8766D", "Splice region" = "#00A9FF","Intronic" = "#CD9600", "Synonymous" = "#C77CFF", "Nonsense" = "#FF61CC","Missense" = "#7CAE00", "NULL"="white", "Codon_del"="black","NULL"="white")  #Added in codon_del as black outline. The artificially added tile as white outline and white fill
  
  #ggplot
  p <- ggplot(total_snv_filtered_combi) +
    geom_tile(aes(x=plot_position,y=VCF_Alt,fill = Confidence,color=Summary_Consequence, width=width), size=3, linewidth=0.3, height=0.6,stat="identity") + #plot the data that present in our dataset
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.4), legend.position = "None") + #Remove the legend position setting to show the legend.
    scale_x_discrete(breaks = levels(total_snv_filtered_combi$plot_position)[c(T, rep(F, 9))], position = "top") + #Break the label for showing a label every 10 bases.
    scale_fill_manual(values = fill) +
    scale_color_manual(values = cols) +
    xlab("") + # No x y axis label
    ylab("") +
    geom_tile(data = subset(total_snv_filtered_combi,  is.na(Confidence)), aes(x=plot_position,y=VCF_Alt,fill = NA), fill = "white", size= 3,width=0.6, height=0.6,linewidth=0.3,color="black") + #plot the data as white fill that absent in our dataset
    coord_fixed(ratio = 2)+ labs(title=exon_label) 
  
  gb <- ggplot_build(p)
  
  p + geom_segment(data=gb$data[[2]]%>% filter(y<=4), aes(x=xmin, xend=xmax, y=ymin, yend=ymax), color="black",size=0.3)  #Add diagonal line for the empty SNV data
}


#Plot kinetic
plot_nuc_kinetic <- function(df, exon) {
  #Levels
  nuc_level = c("T","G","C","A","Codon_del")
  
  #Filter the SNV only. The SNV has one base in both VCF-Ref and VCF-Alt column. Also, create a new column, "plot position".
  totalsnv <- total %>% filter(SGE_exon_group %in% c(exon)) %>% #Filter a specific exon
    filter(Variant_design_type=="snv") %>% #Filter snv only
    filter(nchar(VCF_Ref)==1) %>% filter(nchar(VCF_Alt)==1) %>% #Make sure the snv is snv
    mutate(VCF_Alt=factor(VCF_Alt, levels = nuc_level)) %>% mutate(plot_position=VCF_Hg38_position) %>% #Manipulation
    select("SGE_oligo_name","Variant_design_type","chrom","VCF_Hg38_position","VCF_Ref","VCF_Alt","Consequence","Variant_category","INTRON","EXON","SGE_exon_group","SGE_functional_classification","SGE_NDD_clinical_classification","NDD_clinical_classification_prob_var_functional","NDD_clinical_classification_prob_var_nonfunctional","plot_position") %>%
    mutate(Summary_Consequence=case_when(str_detect(Consequence,"UTR")~ "UTR",
                                         str_detect(Consequence,"splice_donor")~ "Canonical splice",
                                         str_detect(Consequence,"splice_acceptor")~ "Canonical splice",
                                         str_detect(Consequence,"synonymous")~ "Synonymous",
                                         str_detect(Consequence,"stop_gained")~ "Nonsense",
                                         str_detect(Consequence,"stop")~ "Missense",
                                         str_detect(Consequence,"start")~ "Missense",
                                         str_detect(Consequence,"missense")~ "Missense",
                                         str_detect(Consequence,"splice")~ "Splice region",                           
                                         str_detect(EXON,"-")~ "Intronic",
                                         TRUE ~ "Others")) %>% #Define Summary_Consequence
    mutate(Confidence=case_when(NDD_clinical_classification_prob_var_nonfunctional >= 0.9~ "High",
                                NDD_clinical_classification_prob_var_nonfunctional <= 0.5~ "Low",
                                TRUE ~ "Intermediate")) #Define Confidence
  
  #Filter the inframe only. Replace the VCF Alt to codon_del. Also, create a new column, "plot position".
  totalinframedel <- total %>% filter(SGE_exon_group %in% c(exon)) %>% #Filter a specific exon
    filter(str_detect(Variant_design_type,"inf")) %>% #Filter inframe only
    mutate(VCF_Alt="Codon_del",VCF_Alt=factor(VCF_Alt, levels = nuc_level))%>% #Replace VCF_Alt as Codon_del
    mutate(plot_position=VCF_Hg38_position+2) %>% #After testing, it has to be +2 so that the position align well
    select("SGE_oligo_name","Variant_design_type","chrom","VCF_Hg38_position","VCF_Ref","VCF_Alt","Consequence","Variant_category","INTRON","EXON","SGE_exon_group","SGE_functional_classification","SGE_NDD_clinical_classification","NDD_clinical_classification_prob_var_functional","NDD_clinical_classification_prob_var_nonfunctional","plot_position") %>%
    mutate(Summary_Consequence="Codon_del") %>% #Define the consequence as codon_del
    mutate(Confidence=case_when(NDD_clinical_classification_prob_var_nonfunctional >= 0.9~ "High",
                                NDD_clinical_classification_prob_var_nonfunctional <= 0.5~ "Low",
                                TRUE ~ "Intermediate")) #Define Confidence
  
  #Gather the number of x-axis component to plot
  osnv <-seq(min(totalsnv$plot_position), min(totalsnv$plot_position)+221,1) #Hard-coded, 221 is the nucleotide number in the longest exon, Exon13. The idea here is to make sure all plot as the same amount of square/tile. So that the size are the same across all plot.
  maxosnv <- max(totalsnv$plot_position) #Beyond this point, the tile/square were artificially added and can be remove. Appear as white in the plot.
  
  #Create a data frame that make all combination of each plot position and nuc_level
  snv_combi <- expand.grid(VCF_Alt=nuc_level[c(1:4)],plot_position=osnv)
  
  #Combine this all combination table with the data.
  total_snv_filtered <- totalsnv %>% right_join(snv_combi, by=c("VCF_Alt","plot_position")) %>% mutate(width=0.6) %>% #add the column as width. use to define the width of each square. SNV has a smaller width
    mutate(SGE_functional_classification=case_when(plot_position>maxosnv ~ "NULL", TRUE~SGE_functional_classification )) %>%
    mutate(Summary_Consequence=case_when(plot_position>maxosnv ~ "NULL", TRUE~Summary_Consequence )) %>%
    mutate(Confidence=case_when(plot_position>maxosnv ~ "NULL", TRUE~Confidence )) #Replace the columns to NULL whenvever it is more than maxosnv that is the designed maximum coordinate range
  
  #cdel_combi <- expand.grid(VCF_Alt="Codon_del",plot_position=osnv)
  
  #Create a dataframe that has a full matrix of codon position, including those repeated trinucleotide. A possible bug in the future: if the exon ended with a repeated trinucleotide, the last repeated codon will not be added and have to be added manually.
  inframe_full_pos <- data.frame(plot_position=seq(min(totalinframedel$plot_position), max(totalinframedel$plot_position),3))
  
  inframe_full<- full_join(totalinframedel, inframe_full_pos, by="plot_position") %>% arrange(plot_position) %>% fill( colnames(totalinframedel)[-16], .direction = "down")  %>% #Fill the empty row with the value in the previous row whenever the "plot_position" is the same as the previous row.
    mutate(width=0.6*4) #Define the width of the tile. Tested before and 4 is the best.
  
  # Combine the SNV, inframedel into a table. Factorize the plot_position make all the position in to 1,2,3,4... it might be easier for manipulate some plot features. (but not used directly here. It provide the idea of the range of x coordinate)
  total_snv_filtered_combi<-bind_rows(total_snv_filtered,inframe_full) %>% mutate(plot_position=as.character(plot_position)) %>% mutate(plot_position=factor(plot_position, levels = unique(plot_position))) %>% mutate(factor= as.numeric(plot_position))
  
  #The text in the exon label to show in the plot
  exon_label= paste0("Exon ",exon)
  
  #fill <- c("fast depleting" = "#d55e00", "slow depleting" = "#0072b2", "enriched" = "#f0e442", "unchanged" = "white","NULL"="white")
  fill <- c("fast depleting" = "#1F008A", "slow depleting" = "#4A89F7", "enriched" = "#FF9400", "unchanged" = "white","NULL"="white")
  cols <- c("UTR" = "#00BFC4", "Canonical splice" = "#F8766D", "Splice region" = "#00A9FF","Intronic" = "#CD9600", "Synonymous" = "#C77CFF", "Nonsense" = "#FF61CC","Missense" = "#7CAE00", "NULL"="white", "Codon_del"="black","NULL"="white")  #Added in codon_del as black outline. The artificially added tile as white outline and white fill
  
  #ggplot
  p <- ggplot(total_snv_filtered_combi) +
    geom_tile(aes(x=plot_position,y=VCF_Alt,fill = SGE_functional_classification,color=Summary_Consequence, width=width), size=3, linewidth=0.3, height=0.6,stat="identity") + #plot the data that present in our dataset
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.4), legend.position = "None") + #Remove the legend position setting to show the legend.
    scale_x_discrete(breaks = levels(total_snv_filtered_combi$plot_position)[c(T, rep(F, 9))], position = "top") + #Break the label for showing a label every 10 bases.
    scale_fill_manual(values = fill) +
    scale_color_manual(values = cols) +
    xlab("") + # No x y axis label
    ylab("") +
    geom_tile(data = subset(total_snv_filtered_combi,  is.na(SGE_functional_classification)), aes(x=plot_position,y=VCF_Alt,fill = NA), fill = "white", size= 3,width=0.6, height=0.6,linewidth=0.3,color="black") + #plot the data as white fill that absent in our dataset
    coord_fixed(ratio = 2)+ labs(title=exon_label) 
  
  gb <- ggplot_build(p)
  
  p + geom_segment(data=gb$data[[2]]%>% filter(y<=4), aes(x=xmin, xend=xmax, y=ymin, yend=ymax), color="black",size=0.3)  #Add diagonal line for the empty SNV data
}

#A better way to put the plot. All exon in one plot
pdf("test.pdf",13,30)
grid.arrange(grobs = lapply(c(1:17), function (x) {plot_nuc_conf(total,x)}), ncol = 1)
dev.off()

pdf("test2.pdf",13,30)
grid.arrange(grobs = lapply(c(1:17), function (x) {plot_nuc_kinetic(total,x)}), ncol = 1)
dev.off()


#Or, plot each exon in a seperate plot
pdf("viz_codon_fixed.pdf",13,6)
lapply(c(1:17), function (x) {plot_nuc_conf(total,x)})
dev.off()
