setwd("D:/A.clarkii/clarkii_skin_featureCounts")

# read and merge featurecounts output
f2 <- read.csv("F2_skin_ctrl.csv",header=T)
f3 <- read.csv("F3_skin_ctrl.csv",header = T)
f4 <- read.csv("F4_skin_ctrl.csv", header = T)
f5 <- read.csv("F5_skin_interaction.csv", header = T)
f6 <- read.csv("F6_skin_interaction.csv", header = T)
f7 <- read.csv("F7_skin_interaction.csv", header = T)
f9 <- read.csv("F9_skin_interaction.csv", header = T)

df_list <- list(f2, f3, f4, f5, f6, f7, f9)
merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
write.csv(merge,file="skin_cts_data.csv",row.names=F)

#****************************************************************************
#Get normalized counts across all samples
library(DESeq2)
setwd("D:/A.clarkii/clarkii_skin_featureCounts")

cts <- read.csv("skin_cts_data.csv",row.names="Geneid")
coldata <- read.csv("skin_metadata.csv",row.names=1)

#convert decimal to int
cts1=round(cts)

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata) %in% colnames(cts1)) || all(colnames(cts1) %in% rownames(coldata))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts1) == rownames(coldata)))

# all variables in the design formula should be converted to factors
coldata$Condition <- factor(coldata$Condition)
coldata$Sex <- factor(coldata$Sex)

dds <- DESeqDataSetFromMatrix(countData = cts1, colData = coldata, design = ~Sex + Condition)

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# get normalised counts
dds <- estimateSizeFactors(dds)
normalized_cts <- counts(dds,normalized =TRUE)
write.csv(normalized_cts,file="skin_normalized_cts.csv",row.names=T)

#****************************************************************************
#PCA from rlog transformed counts
# link to guide: https://tavareshugo.github.io/data-carpentry-rnaseq/03_rnaseq_pca.html
library(DESeq2)
library(tidyverse)

setwd("D:/A.clarkii/clarkii_skin_featureCounts")
cts <- read.csv("skin_cts_data.csv",row.names="Geneid")
coldata <- read.csv("skin_metadata.csv",row.names=1)

#convert decimal to int
cts1=round(cts)

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata) %in% colnames(cts1)) || all(colnames(cts1) %in% rownames(coldata))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts1) == rownames(coldata)))

# all variables in the design formula should be converted to factors
coldata$Condition <- factor(coldata$Condition)
coldata$Sex <- factor(coldata$Sex)

dds <- DESeqDataSetFromMatrix(countData = cts1, colData = coldata, design = ~Sex + Condition)

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

#The rlog function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size. 
#The transformation is useful when checking for outliers. 
#Note that neither rlog transformation nor the VST are used by the differential expression estimation in DESeq, which always occurs on the raw count data, through generalized linear modeling which incorporates knowledge of the variance-mean dependence. The rlog transformation and VST are offered as separate functionality which can be used for visualization, clustering or other machine learning tasks. 

rlogcounts <- assay(rlog(dds, blind=TRUE)) #blind = T will make the log transformation blind to the sample info)

#Use the rlog transformed counts to do a PCA
#To use prcomp, the samples should be the rows and the variables (in this case genes) should be the columns. So we need to transpose the data
transpose <- t(rlogcounts)

sample_pca <- prcomp(transpose,center=T, scale=T)#Often it is a good idea to standardize the variables before doing the PCA. This is often done by centering the data on the mean and then scaling it by dividing by the standard deviation. This ensures that the PCA is not too influenced by genes with higher absolute expression. By default, the prcomp() function does the centering but not the scaling. 

pc_eigenvalues <- sample_pca$sdev^2 # to extract the variance explained by each PC from our sample_pca object.

#pc_eigenvalues is a vector. So to make a plot with ggplot2, we need to transform it to a dataframe/tibble

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

# The pc_eigenvalues table can now be used to produce a Scree Plot, which shows the fraction of total variance explained by each principal component.
# The plot will show both the variance explained by individual PCs as well as the cumulative variance, using a type of visualisation known as a pareto chart
#A Pareto chart is a type of chart that contains both bars and a line graph, where individual values are represented in descending order by bars, and the cumulative total is represented by the line. 

p <- pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")
p
# The plot shows how successive PCs explain less and less of the total variance in the original data. Also note that in this case 7 components are enough to virtually explain all of the variance in the dataset. This makes sense in this case since we only have 7 biological samples.

p <- p + theme_bw() + theme(axis.text=element_text(size=26), axis.title=element_text(size=30,face = "bold"))
p

setwd("./figures")
ggsave("PCA_variance_for_each_PC.png", dpi = 300)

#Visualising samples on PC space
library(ggpubr) 
setwd("D:/A.clarkii/clarkii_skin_featureCounts")
#extract the PC score for each sample from the prcomp object
coldata <- read.csv("skin_metadata.csv") 
pc_scores <- sample_pca$x ## The PC scores are stored in the "x" value of the prcomp object

#The pc_scores object is of class matrix, so we need to first convert it to a data.frame/tibble for ggplot2.
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "Sample")

# print the result
pc_scores

# Finally plot the PCA
p <- pc_scores %>% 
  # join with "sample_info" table 
  full_join(coldata, by = "Sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = Condition)) +
  geom_point(size=5,aes(shape=Sex)) 

p

pc1_variance <- pc_eigenvalues[1,3]
pc2_variance <- pc_eigenvalues[2,3]

p <- p + xlab(paste0("PC1: ",round(pc1_variance,digits = 2),"%")) +
  ylab(paste0("PC2: ",round(pc2_variance,digits = 2),"%"))
  
p

p <- p + theme_bw() + theme(axis.text=element_text(size=26), axis.title=element_text(size=30,face = "bold"))
p

p <- p + theme(legend.title = element_text(size=30), 
                 legend.text = element_text(size=26))

p

setwd("./figures")
ggsave("PCA_skin.png", dpi = 300)

#Exploring correlation between genes and PCs
#This is to see which genes have the most influence on each PC axis. This information is contained in the variable loadings of the PCA, within the rotation value of the prcomp object.

pc_loadings <- sample_pca$rotation

#The pc_loading object is a matrix, which we can convert to a data.frame/tibble for plotting and data manipulation

pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings

#extract the top 10 genes with highest loading on PC1 and PC2
top_genes <- pc_loadings %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes

#Now, we can use this list of gene names to subset the eigenvalues table:
top_loadings <- pc_loadings %>% 
  filter(gene %in% top_genes)

#plot 
loadings_plot <- ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "black",linewidth=1) +
  geom_text_repel(aes(x = PC1, y = PC2, label = gene),
            size = 8,hjust=2, vjust=3) + theme_bw() +
  theme(axis.text=element_text(size=26), 
        axis.title=element_text(size=30,face = "bold")) +
  xlab("PC1") + ylab ("PC2")

loadings_plot

#Interpretation
#Use the loading plot to identify which variables have the largest effect on each component. 
#Loadings can range from -1 to 1. Loadings close to -1 or 1 indicate that the variable strongly influences the component. 
#Loadings close to 0 indicate that the variable has a weak influence on the component.
#Positive loadings indicate a variable and a principal component are positively correlated. Negative loadings indicate a negative correlation. 

setwd("./figures")
ggsave("PCA_loadings_top10genesForeachPC.png", dpi = 300)

# Combining PC scores and PC loadings to see correlation - biplot
#For example, geneA, with a negative score of -0.00605 (loading score), will be associated with whatever samples are shifted toward the negative end of your PC1 on a bi-plot.

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('PCAtools')

library(PCAtools)
setwd("D:/A.clarkii/clarkii_skin_featureCounts")

cts <- read.csv("skin_cts_data.csv",row.names="Geneid")
coldata <- read.csv("skin_metadata.csv",row.names = 1)

#convert decimal to int
cts1=round(cts)

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata) %in% colnames(cts1)) || all(colnames(cts1) %in% rownames(coldata))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts1) == rownames(coldata)))

# all variables in the design formula should be converted to factors
coldata$Condition <- factor(coldata$Condition)
coldata$Sex <- factor(coldata$Sex)

dds <- DESeqDataSetFromMatrix(countData = cts1, colData = coldata, design = ~Sex + Condition)

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

#The rlog function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size. 
#The transformation is useful when checking for outliers. 
#Note that neither rlog transformation nor the VST are used by the differential expression estimation in DESeq, which always occurs on the raw count data, through generalized linear modeling which incorporates knowledge of the variance-mean dependence. The rlog transformation and VST are offered as separate functionality which can be used for visualization, clustering or other machine learning tasks. 

rlogcounts <- assay(rlog(dds, blind=TRUE)) #blind = T will make the log transformation blind to the sample info)


a <- pca(rlogcounts, metadata = coldata1,center=T,scale=T)

biplot(a, showLoadings = TRUE, colby='Condition', shape='Sex',
       lab = NULL, pointSize = 7, sizeLoadingsNames = 7,
       ntopLoadings = 10,max.overlaps = 30, boxedLoadingsNames = FALSE,
       legendPosition="bottom",legendLabSize = 20,axisLabSize=26,
       legendIconSize = 8,colLegendTitle=NULL,shapeLegendTitle=NULL)

setwd("./figures")
ggsave("PCA_biplot.png", dpi = 300)

#****************************************************************************
#LRT to see the effect of sex and condition on gene expression 
library(DESeq2)

setwd("D:/A.clarkii/clarkii_skin_featureCounts")
cts <- read.csv("skin_cts_data.csv",row.names="Geneid")
coldata <- read.csv("skin_metadata.csv",row.names=1)

#convert decimal to int
cts1=round(cts)

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata) %in% colnames(cts1)) || all(colnames(cts1) %in% rownames(coldata))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts1) == rownames(coldata)))

# all variables in the design formula should be converted to factors
coldata$Condition <- factor(coldata$Condition)
coldata$Sex <- factor(coldata$Sex)

dds <- DESeqDataSetFromMatrix(countData = cts1, colData = coldata, 
                              design = ~Sex + Condition)

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# total expressed genes in brain (passed the rowSums > 0 filter)
n_expressed_skin <- nrow(dds)
n_expressed_skin #22684

# testing for sex effect 
dds_LRT_sex <- DESeq(dds, test = "LRT", reduced = ~ Condition)
res_sex <- results(dds_LRT_sex)
res_sex
summary(res_sex)
results <- data.frame(res_sex)
results <- na.omit(results)

res_sig_padj <- results[results$padj < 0.05, ]

# number significant
n_sig_sex <- sum(res_sex$padj < 0.05, na.rm = TRUE)
n_sig_sex #2

# proportion of expressed genes
round(100 * n_sig_sex / n_expressed_skin, 2)  #0.01%

setwd("./LRT_results/")
write.csv(results,file="LRT_results_sex_effect.csv")
write.csv(res_sig_padj,file="LRT_results_sex_effect_sig.csv")

# testing for condition effect 
dds_LRT_condition <- DESeq(dds, test = "LRT", reduced = ~ Sex)
res_condition <- results(dds_LRT_condition)
res_condition
summary(res_condition)
results <- data.frame(res_condition)
results <- na.omit(results)

res_sig_padj <- results[results$padj<0.05,]

setwd("./LRT_results/")
write.csv(results,file="LRT_results_condition_effect.csv")
write.csv(res_sig_padj,file="LRT_results_condition_effect_sig.csv")

#****************************************************************************
#clustering heatmap to check if samples group by sex
library(DESeq2)
library("pheatmap")
library("dplyr")

setwd("D:/A.clarkii/clarkii_skin_featureCounts")
cts <- read.csv("skin_cts_data.csv",row.names="Geneid")
coldata <- read.csv("skin_metadata.csv",row.names=1)

#convert decimal to int
cts1=round(cts)

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata) %in% colnames(cts1)) || all(colnames(cts1) %in% rownames(coldata))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts1) == rownames(coldata)))


# all variables in the design formula should be converted to factors
coldata$Condition <- factor(coldata$Condition)
coldata$Sex <- factor(coldata$Sex)


dds <- DESeqDataSetFromMatrix(countData = cts1, colData = coldata, design = ~Condition + Sex)

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

#The rlog function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size. 
#The transformation is useful when checking for outliers. 
#Note that neither rlog transformation nor the VST are used by the differential expression estimation in DESeq, which always occurs on the raw count data, through generalized linear modeling which incorporates knowledge of the variance-mean dependence. The rlog transformation and VST are offered as separate functionality which can be used for visualization, clustering or other machine learning tasks. 

rlogcounts <- assay(rlog(dds, blind=TRUE)) #blind = T will make the log transformation blind to the sample info)

cor <- cor(rlogcounts) #Calculate the correlation values between samples 

library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)
heat_colors <- brewer.pal(9, "Greys")

annotation_cols <- list(Condition = c(Control = "#808080", Interaction = "#FEBE00"), 
                        Sex = c(Male = "#6B8E23", Female = "#CC8899"))

p <- pheatmap(cor,color=heat_colors,cellheight=35,cellwidth=40, cluster_rows=TRUE, show_rownames=T,show_colnames=T,legend=T,fontsize = 10,
         cluster_cols=T, annotation=select(coldata,c(Condition,Sex)),annotation_colors = annotation_cols)

setwd("./figures/")
png("clustering_heatmap.png",res = 300,units="cm", width = 35, height = 20)
print(p)
dev.off()

#****************************************************************************
# DE analysis - wald test
library(DESeq2)

setwd("D:/A.clarkii/clarkii_skin_featureCounts")
cts <- read.csv("skin_cts_data.csv",row.names="Geneid")
coldata <- read.csv("skin_metadata.csv",row.names=1)

#convert decimal to int
cts1=round(cts)

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata) %in% colnames(cts1)) || all(colnames(cts1) %in% rownames(coldata))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts1) == rownames(coldata)))

# all variables in the design formula should be converted to factors
coldata$Condition <- factor(coldata$Condition)


dds <- DESeqDataSetFromMatrix(countData = cts1, colData = coldata, 
                              design = ~Condition)

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# Explicitly set the factor levels
dds$condition <- relevel(dds$Condition, ref = "Control")
dds <- DESeq(dds)
resultsNames(dds)

# get the model matrix produced by DESeq2
mod_matrix <- model.matrix(design(dds), colData(dds))
mod_matrix

# calculate the vector of coefficient weights in control
ctrl <- colMeans(mod_matrix[dds$condition == "Control", ])
ctrl

# calculate the vector of coefficient weights in interaction
inter <- colMeans(mod_matrix[dds$condition == "Interaction", ])
inter

# The contrast we are interested in is the difference between sun and shade
inter - ctrl #ctrl should be at the end

# get the results for this contrast
res <- results(dds, alpha = 0.05,test="Wald", contrast = inter - ctrl) #So up/downregulation refers to interaction i.e upregulated in teraction or downregulated in interaction
summary(res)
res

results <- data.frame(res)
results <- na.omit(results)
sig_results <- results[results$padj<0.05, ]

setwd("./DE_analysis_results_waldTest/")
write.csv(results,file="DE_results_all_nosex.csv")
write.csv(sig_results,file="DE_results_sig_nosex.csv")

############ test DE by adding sex in the design - this is the final one
setwd("D:/A.clarkii/clarkii_skin_featureCounts")
cts <- read.csv("skin_cts_data.csv",row.names="Geneid")
coldata <- read.csv("skin_metadata.csv",row.names=1)

#convert decimal to int
cts1=round(cts)

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata) %in% colnames(cts1)) || all(colnames(cts1) %in% rownames(coldata))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts1) == rownames(coldata)))

# all variables in the design formula should be converted to factors
coldata$Condition <- factor(coldata$Condition)
coldata$Sex <- factor(coldata$Sex)


dds <- DESeqDataSetFromMatrix(countData = cts1, colData = coldata, 
                              design = ~ Sex + Condition)

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# Explicitly set the factor levels
dds$condition <- relevel(dds$Condition, ref = "Control")
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, alpha = 0.05, contrast = c("Condition","Interaction","Control"), test ="Wald")

summary(res)
res

results <- data.frame(res)
results <- na.omit(results)
sig_results <- results[results$padj<0.05, ]

setwd("./DE_analysis_correct_results_waldtest/")
write.csv(results,file="DE_results_all_includingsex.csv")
write.csv(sig_results,file="DE_results_sig_includingsex.csv")


#****************************************************************************
# edit A.clarkii annotation file to have one value of each braker gene and combine all the NCBI gene IDs corresponding to that Braker gene 
library(dplyr)

setwd("D:/A.clarkii/clarkii_skin_featureCounts/DE_analysis_results_waldTest")

GO <- read.table("A.clarkii_GO_Annotation.txt",sep="\t",header=T)
GO1 <- GO[,c(2,3)]

GO_edited <- GO1 %>%
  group_by_at(vars(-geneID)) %>%
  summarise(geneID = toString(geneID))

omicsbox_table <- read.csv("omicsbox_table.csv",header=T)

merge <- merge(GO_edited,omicsbox_table,by="BrakerGene",all=T)

write.csv(merge1,file="A.clarkii_GO_Annotation_edited.csv",row.names=F)

#****************************************************************************
# merge skin sig DE gene list with GO functional annotation
setwd("D:/A.clarkii/clarkii_skin_featureCounts/DE_analysis_results_waldTest")

genes <- read.csv("DE_results_sig_includingsex.csv")
GO <- read.csv("A.clarkii_GO_Annotation_edited.csv")

merge <- merge(genes,GO,by="BrakerGene",all.x=T)
write.csv(merge,file="DE_results_sig_includingsex_annotation.csv")



#****************************************************************************
#*# Figures

#*# Enrichment figures
library(ggplot2)

setwd("D:/A.clarkii/clarkii_skin_featureCounts/DE_analysis_correct_results_waldtest")

file <- read.csv("enrichment_results_skin.csv")

p <- ggplot(file,aes(y=GO.Name,x=GO.ID)) +
  geom_point(aes(size=Nr.Test,color=FDR)) + 
  scale_size_area(max_size = 5) 
p

p <- p + theme_bw() +
  theme(axis.title = element_blank(),axis.text.x = element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size=10)) + 
  theme(legend.title = element_text(size=12)) +
  theme(legend.box = "horizontal") +
  labs(size='Gene No.')

p  

setwd("./figures/")
ggsave("enriched_dot_plot_skin.png",dpi = 300)


#*# Fig to show logFC, pval and up/down regulation for genes in each functional group in separate panels
setwd("D:/A.clarkii/clarkii_skin_featureCounts/DE_analysis_correct_results_waldtest")

library(ggplot2)
library(ggh4x)
library(ggrepel)

file <- read.csv("DE_results_sig_includingsex.csv")

file1 <- file[c(1:25),]

design <- "AB#
           DEC"
# design function helps control location of panels - the hashtag is for space

ggplot(file1,aes(x=Gene_symbol_NCBI,y=log2FoldChange,color=padj)) + 
  geom_point(size=3) + geom_hline(yintercept=0) + theme_classic() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  facet_manual(vars(Functional.category),design = design) + 
  theme(panel.background = element_rect(fill = 'white', color = 'black')) +
  theme(axis.text=element_text(size=18),axis.title.x = element_text(size=22), axis.title.y=element_text(size=22)) + 
  ylab("log2 Fold Change") + xlab("Gene ID") + 
  theme(legend.text=element_text(size=18)) +
  theme(legend.title=element_text(size=22)) +
  theme(strip.text = element_text(size = 20)) +
  geom_text_repel(aes(label=Gene_symbol_NCBI))

setwd("./figures/")
ggsave("skin_DE_genes_discussion.png", dpi=300)

# ============================================
# Test overlap between sig sex effect genes from LRT and the sig DE genes from Wald test
# ============================================

lrt <- read.csv("LRT_results/LRT_results_sex_effect_sig.csv",row.names="Geneid")
sex_genes <- rownames(lrt)

de_genes <- read.csv("DE_analysis_correct_results_waldtest/DE_results_sig_includingsex.csv",row.names="gene_id")
skin_de <- row.names(de_genes)

overlap <- intersect(sex_genes, skin_de)
