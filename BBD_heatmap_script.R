#!/usr/bin/env Rscript

library(pheatmap)
library(RColorBrewer)
library(dplyr)

###this is the variable for the width of the pdf, i have set it at 10 but depending on the numer of genes in the pathway you might have to make it bigger
width_pdf=10

###filename to read in
#HM_file<-"NER.output.tsv"
HM_files<-list.files(pattern="\\.tsv$")
list_file<-read.csv("phylo_order.csv", header=T, row.names=1)
###value to cap all data
max_val=1200
for (HM_file in HM_files){
  ###create output file name
  
  output_filename<-paste0("test/",HM_file,".heatmap.pdf")
  
  ###read in table
  HM_data<-read.table(HM_file, header=T, row.names = 1)
  
  ###extract main HM data
  HM_data_HM<-HM_data[c(3:nrow(HM_data)),]
  
  ####tidy_data
  HM_data_HM[HM_data_HM>max_val]<-max_val
  HM_data_HM[HM_data_HM==0]<-NA
  
  ###left fjopin to preserve order
  HM_data_HM_t<-t(HM_data_HM)
  
  # Assuming `list_file` and `HM_data_HM` are your data frames
  # Convert the rownames of `HM_data_HM` to a column
  HM_data_HM_df <- as.data.frame(HM_data_HM_t)  # If `HM_data_HM` is a matrix, convert it to a data frame
  HM_data_HM_df$rownames <- rownames(HM_data_HM_df)
  
  # Similarly, ensure `list_file` has a column to join by (assuming rownames should be used)
  list_file$rownames <- rownames(list_file)
  
  # Perform the left join using the `rownames` column
  HM_data_HM_ord <- left_join(list_file, HM_data_HM_df, by = "rownames")
  
  # Optionally remove the rownames column if no longer needed
  HM_data_HM_ord$rownames <- NULL
  
  ###extract Lifestage column and code
  
  lifestage_data<-HM_data[1,]
  lifestage_data_t<-t(lifestage_data)
  lifestage_data_t_coded<-data.frame(Lifestyle=ifelse(test=lifestage_data_t==1,yes="Free living", no=ifelse(test = lifestage_data_t==2,yes="Extracellular", no=ifelse(test=lifestage_data_t==3, yes = "Cytoplasmic", no="Intranuclear"))))
  
  ### extract mitochondrial, plastid status data
  MPS_data<-data.frame(list_file[,2])
  colnames(MPS_data)[1] <- "Mitocondrial_plastid_status"
 #MPS_data_t<-t(MPS_data)
 # MPS_data_coded<-data.frame(Mitocondrial_plastid_status=ifelse(test=MPS_data==1,yes="Plastid and ATP-producing mitochondria", no="No ATP-producing mitochondria"))
  
  annotation_data<-cbind(lifestage_data_t_coded,MPS_data)
  
  ###set Heatmap colours
  hmcol<-colorRampPalette(brewer.pal(3,"BuPu"))(256)
  
  my_colour = list( "Lifestyle" = c( "Free living" = "red", "Extracellular" = "darkgreen", "Cytoplasmic" = "deepskyblue3", "Intranuclear" = "purple") ,
                  "Mitocondrial_plastid_status" = c("Plastid and ATP-producing mitochondria" = "#FBA475FF", "No ATP-producing mitochondria" = "#4C84A3FF"))
  
  ###produce heatmap
  pdf(output_filename, width=width_pdf)
  pheatmap(t(HM_data_HM), cluster_rows= F, cluster_cols = F, color = hmcol, annotation_row = annotation_data , annotation_colors = my_colour, fontsize =6, cellwidth = 10, na_col = "white")
  dev.off()
}
