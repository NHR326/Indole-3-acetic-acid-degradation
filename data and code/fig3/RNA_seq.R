#ExampleSE9
#######################################################################
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
#######################################################################
RSEM_out = "SE9"
group_info = "SE9_group.csv"
group_category = c("IAA","Glc")
#######################################################################
file_names = list.files(RSEM_out)
for (i in 1:length(file_names)) {
  table = read.table(paste(RSEM_out,file_names[i],sep = "/"),
                   sep = "\t",
                   header = T,
                   quote = "",
                   stringsAsFactors = FALSE) %>%    setNames(c('gene_id','transcript_id','length','effective_length','expected_count','TPM','FPKM'))
if(i<2){RSEM.count=table[,c('transcript_id','expected_count')]} else {
    RSEM.count = dplyr::full_join(RSEM.count, 
                                  table[,c('transcript_id','expected_count')], 
                                  by = c('transcript_id'))
  }
}
colnames(RSEM.count)[2:ncol(RSEM.count)] = gsub(".genes.results","",file_names)
RSEM.count = tibble::column_to_rownames(RSEM.count,"transcript_id")
RSEM_count <- round(RSEM.count)
#######################################################################
Group = read.csv(group_info)
Group$group = factor(Group$group, levels = group_category)
table(Group$group)
RSEM_count = RSEM_count[, Group$sample]
#######################################################################
dds = DESeqDataSetFromMatrix(countData = RSEM_count, colData = Group, design = ~group)
dds = dds[rownames(counts(dds)) > 1, ] 
dds = estimateSizeFactors(dds) 
dds = DESeq(dds)
Contrast = c("group", group_category)
res = results(dds, contrast = Contrast)
DEG = as.data.frame(res)
data1 = cbind(Gene = rownames(DEG),DEG)
rownames(data1) = NULL
data2 = cbind(Gene = rownames(RSEM_count), RSEM_count)
rownames(data2) = NULL
Final_result = merge(data2,data1, by.x = "Gene")
write.csv(Final_result,file = "SE9.csv")
######################################################################
#Print
print_data = read.csv("SE9.csv", header=TRUE, sep=",")
print_data = na.omit(print_data)
FC = 2
PValue = 0.05
#######################################################################
print_data$sig[(-1*log10(print_data$padj) < -1*log10(PValue)|print_data$padj=="NA")|(print_data$log2FoldChange < log2(FC))& print_data$log2FoldChange > -log2(FC)] <- "NotSig"
print_data$sig[-1*log10(print_data$padj) >= -1*log10(PValue) & print_data$log2FoldChange >= log2(FC)] <- "Up"
print_data$sig[-1*log10(print_data$padj) >= -1*log10(PValue) & print_data$log2FoldChange <= -log2(FC)] <- "Down"
#######################################################################
print_data$label=ifelse(print_data$Marker == 1, as.character(print_data$Gene), '')
p1 = ggplot(print_data,aes(print_data$log2FoldChange, -1*log10(print_data$padj))) +    
  geom_point(aes(color = sig)) +                           
  labs(title="SE9",                                
      x="log[2](FC)", 
       y="-log[10](PValue)") + 
  scale_color_manual(values = c("#0570B0", "#B3B3B3", "#EF3B2C")) + 
 geom_hline(yintercept=-log10(PValue),linetype=2)+        
  geom_vline(xintercept=c(-log2(FC),log2(FC)),linetype=2)+ 
  geom_text_repel(aes(x = print_data$log2FoldChange,                   
                      y = -1*log10(print_data$padj),          
                      label=label),                      
                  max.overlaps = 10000,                    
                  size=3,                                
                  box.padding=unit(0.5,'lines'),           
                  point.padding=unit(0.1, 'lines'), 
                  segment.color='black',                   
                  show.legend=FALSE) +
 theme_classic()
ggsave(filename = "SE9.pdf",p1,width=8,heigh=10)
