setwd("D:/A.clarkii/clarkii_brain_DE/featurects_output")

#read featurecounts output
F1_BS <- read.table("F1-BSfeaturecounts_output.txt",header=T) # Read the counts file obtained from featurecounts
F1_BS_1 <- F1_BS[,c(1,7)] # Select the columns of interest - for me it was the 1st and 7th column which is gene ID and the column with the raw counts

### Repeat for all samples

F2_BS <- read.table("F2-BSfeaturecounts_output.txt",header=T)
F2_BS_1 <- F2_BS[,c(1,7)]
F3_BS <- read.table("F3-BSfeaturecounts_output.txt",header=T)
F3_BS_1 <- F3_BS[,c(1,7)]
F4_BS <- read.table("F4-BSfeaturecounts_output.txt",header=T)
F4_BS_1 <- F4_BS[,c(1,7)]
F5_BS <- read.table("F5-BSfeaturecounts_output.txt",header=T)
F5_BS_1 <- F5_BS[,c(1,7)]
F6_BS <- read.table("F6-BSfeaturecounts_output.txt",header=T)
F6_BS_1 <- F6_BS[,c(1,7)]
F7_BS <- read.table("F7-BSfeaturecounts_output.txt",header=T)
F7_BS_1 <- F7_BS[,c(1,7)]
F8_BS <- read.table("F8-BSfeaturecounts_output.txt",header=T)
F8_BS_1 <- F8_BS[,c(1,7)]
F9_BS <- read.table("F9-BSfeaturecounts_output.txt",header=T)
F9_BS_1 <- F9_BS[,c(1,7)]

F1_DE <- read.table("F1-DEfeaturecounts_output.txt",header=T)
F1_DE_1 <- F1_DE[,c(1,7)]
F2_DE <- read.table("F2-DEfeaturecounts_output.txt",header=T)
F2_DE_1 <- F2_DE[,c(1,7)]
F3_DE <- read.table("F3-DEfeaturecounts_output.txt",header=T)
F3_DE_1 <- F3_DE[,c(1,7)]
F4_DE <- read.table("F4-DEfeaturecounts_output.txt",header=T)
F4_DE_1 <- F4_DE[,c(1,7)]
F5_DE <- read.table("F5-DEfeaturecounts_output.txt",header=T)
F5_DE_1 <- F5_DE[,c(1,7)]
F6_DE <- read.table("F6-DEfeaturecounts_output.txt",header=T)
F6_DE_1 <- F6_DE[,c(1,7)]
F7_DE <- read.table("F7-DEfeaturecounts_output.txt",header=T)
F7_DE_1 <- F7_DE[,c(1,7)]
F8_DE <- read.table("F8-DEfeaturecounts_output.txt",header=T)
F8_DE_1 <- F8_DE[,c(1,7)]
F9_DE <- read.table("F9-DEfeaturecounts_output.txt",header=T)
F9_DE_1 <- F9_DE[,c(1,7)]

F1_OT <- read.table("F1-OTfeaturecounts_output.txt",header=T)
F1_OT_1 <- F1_OT[,c(1,7)]
F2_OT <- read.table("F2-OTfeaturecounts_output.txt",header=T)
F2_OT_1 <- F2_OT[,c(1,7)]
F3_OT <- read.table("F3-OTfeaturecounts_output.txt",header=T)
F3_OT_1 <- F3_OT[,c(1,7)]
F4_OT <- read.table("F4-OTfeaturecounts_output.txt",header=T)
F4_OT_1 <- F4_OT[,c(1,7)]
F5_OT <- read.table("F5-OTfeaturecounts_output.txt",header=T)
F5_OT_1 <- F5_OT[,c(1,7)]
F6_OT <- read.table("F6-OTfeaturecounts_output.txt",header=T)
F6_OT_1 <- F6_OT[,c(1,7)]
F7_OT <- read.table("F7-OTfeaturecounts_output.txt",header=T)
F7_OT_1 <- F7_OT[,c(1,7)]
F8_OT <- read.table("F8-OTfeaturecounts_output.txt",header=T)
F8_OT_1 <- F8_OT[,c(1,7)]
F9_OT <- read.table("F9-OTfeaturecounts_output.txt",header=T)
F9_OT_1 <- F9_OT[,c(1,7)]

F1_TE <- read.table("F1-TEfeaturecounts_output.txt",header=T)
F1_TE_1 <- F1_TE[,c(1,7)]
F2_TE <- read.table("F2-TEfeaturecounts_output.txt",header=T)
F2_TE_1 <- F2_TE[,c(1,7)]
F3_TE <- read.table("F3-TEfeaturecounts_output.txt",header=T)
F3_TE_1 <- F3_TE[,c(1,7)]
F4_TE <- read.table("F4-TEfeaturecounts_output.txt",header=T)
F4_TE_1 <- F4_TE[,c(1,7)]
F5_TE <- read.table("F5-TEfeaturecounts_output.txt",header=T)
F5_TE_1 <- F5_TE[,c(1,7)]
F6_TE <- read.table("F6-TEfeaturecounts_output.txt",header=T)
F6_TE_1 <- F6_TE[,c(1,7)]
F7_TE <- read.table("F7-TEfeaturecounts_output.txt",header=T)
F7_TE_1 <- F7_TE[,c(1,7)]
F8_TE <- read.table("F8-TEfeaturecounts_output.txt",header=T)
F8_TE_1 <- F8_TE[,c(1,7)]
F9_TE <- read.table("F9-TEfeaturecounts_output.txt",header=T)
F9_TE_1 <- F9_TE[,c(1,7)]

F2_CB <- read.table("F2-CBfeaturecounts_output.txt",header=T)
F2_CB_1 <- F2_CB[,c(1,7)]
F3_CB <- read.table("F3-CBfeaturecounts_output.txt",header=T)
F3_CB_1 <- F3_CB[,c(1,7)]
F4_CB <- read.table("F4-CBfeaturecounts_output.txt",header=T)
F4_CB_1 <- F4_CB[,c(1,7)]
F5_CB <- read.table("F5-CBfeaturecounts_output.txt",header=T)
F5_CB_1 <- F5_CB[,c(1,7)]
F6_CB <- read.table("F6-CBfeaturecounts_output.txt",header=T)
F6_CB_1 <- F6_CB[,c(1,7)]
F7_CB <- read.table("F7-CBfeaturecounts_output.txt",header=T)
F7_CB_1 <- F7_CB[,c(1,7)]
F8_CB <- read.table("F8-CBfeaturecounts_output.txt",header=T)
F8_CB_1 <- F8_CB[,c(1,7)]
F9_CB <- read.table("F9-CBfeaturecounts_output.txt",header=T)
F9_CB_1 <- F9_CB[,c(1,7)]

df_list <- list(F1_BS_1,F1_DE_1,F1_OT_1,F1_TE_1,F2_BS_1,F2_DE_1,F2_OT_1,F2_TE_1,F2_CB_1,
                F3_BS_1,F3_DE_1,F3_OT_1,F3_TE_1,F3_CB_1,F4_BS_1,F4_DE_1,F4_OT_1,F4_TE_1,F4_CB_1,
                F5_BS_1,F5_DE_1,F5_OT_1,F5_TE_1,F5_CB_1,F6_BS_1,F6_DE_1,F6_OT_1,F6_TE_1,F6_CB_1,
                F7_BS_1,F7_DE_1,F7_OT_1,F7_TE_1,F7_CB_1,F8_BS_1,F8_DE_1,F8_OT_1,F8_TE_1,F8_CB_1,
                F9_BS_1,F9_DE_1,F9_OT_1,F9_TE_1,F9_CB_1)

merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
write.csv(merge,file="clarkii_brain_cts.csv",row.names=F)


#****************************************************************************
#PCA from rlog transformed counts
# link to guide: https://tavareshugo.github.io/data-carpentry-rnaseq/03_rnaseq_pca.html
library(DESeq2)
library(tidyverse)

setwd("D:/A.clarkii/clarkii_brain_DE")
cts <- read.csv("clarkii_brain_cts.csv",row.names="Geneid")
coldata <- read.csv("brain_metadata.csv",row.names=1)

#convert decimal to int
cts1=round(cts)

cts2 <- cts1[,-40]
coldata1 <- coldata[-40,]


# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata1) %in% colnames(cts2)) || all(colnames(cts2) %in% rownames(coldata1))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts2) == rownames(coldata1)))

# all variables in the design formula should be converted to factors
coldata1$Condition <- factor(coldata1$Condition)
coldata1$Sex <- factor(coldata1$Sex)
coldata1$Brain_region <- factor(coldata1$Brain_region)


dds <- DESeqDataSetFromMatrix(countData = cts2, colData = coldata1, design = ~Sex + Brain_region + Condition)

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
ggsave("PCA_variance_for_each_PC_remvedF9BS.png", dpi = 300,units="cm",width=60,height=40)

#Visualising samples on PC space
library(ggpubr) 
setwd("D:/A.clarkii/clarkii_brain_DE")
#extract the PC score for each sample from the prcomp object
coldata <- read.csv("brain_metadata.csv")
coldata1 <- coldata[-40,]

pc_scores <- sample_pca$x ## The PC scores are stored in the "x" value of the prcomp object

#The pc_scores object is of class matrix, so we need to first convert it to a data.frame/tibble for ggplot2.
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "Sample")

# print the result
pc_scores

library(ggrepel)
# Finally plot the PCA
p <- pc_scores %>% 
  # join with "sample_info" table 
  full_join(coldata1, by = "Sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = Brain_region,label=Sample,group=Brain_region)) +
  geom_point(size=5,aes(shape=interaction(Condition,Sex))) + scale_shape_manual(values = c(1, 16, 2, 17)) +
  geom_text_repel(size=5,max.overlaps = 15,min.segment.length = 0) +
  stat_conf_ellipse(bary=F, level = 0.95,npoint = 1000, geom = "polygon",alpha=0.3, aes (fill=Brain_region))
   
#library(scales) - to see what default colors are used
#show_col(hue_pal()(5))

p

#no labels
p <- pc_scores %>% 
  # join with "sample_info" table 
  full_join(coldata, by = "Sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = Brain_region,group=Brain_region)) +
  geom_point(size=5,aes(shape=interaction(Condition,Sex))) + scale_shape_manual(values = c(1, 16, 2, 17)) +
  stat_conf_ellipse(bary=F, level = 0.95,npoint = 1000, geom = "polygon",alpha=0.3, aes (fill=Brain_region))


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
ggsave("PCA_brain_all_factors_noLabels_removedF9BS.png", dpi = 300)

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
ggsave("PCA_loadings_top10genesForeachPC_RemovedF9BS.png", dpi = 300)

# Combining PC scores and PC loadings to see correlation - biplot
#For example, geneA, with a negative score of -0.00605 (loading score), will be associated with whatever samples are shifted toward the negative end of your PC1 on a bi-plot.

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('PCAtools')

library(PCAtools)
setwd("D:/A.clarkii/clarkii_brain_DE")

cts <- read.csv("clarkii_brain_cts.csv",row.names="Geneid")
coldata <- read.csv("brain_metadata.csv",row.names = 1)

#convert decimal to int
cts1=round(cts)

cts2 <- cts1[,-40]
coldata1 <- coldata[-40,]


# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata1) %in% colnames(cts2)) || all(colnames(cts2) %in% rownames(coldata1))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts2) == rownames(coldata1)))

# all variables in the design formula should be converted to factors
coldata1$Condition <- factor(coldata1$Condition)
coldata1$Sex <- factor(coldata1$Sex)
coldata1$Brain_region <- factor(coldata1$Brain_region)



dds <- DESeqDataSetFromMatrix(countData = cts2, colData = coldata1, design = ~Sex + Condition + Brain_region)

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

#The rlog function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size. 
#The transformation is useful when checking for outliers. 
#Note that neither rlog transformation nor the VST are used by the differential expression estimation in DESeq, which always occurs on the raw count data, through generalized linear modeling which incorporates knowledge of the variance-mean dependence. The rlog transformation and VST are offered as separate functionality which can be used for visualization, clustering or other machine learning tasks. 

rlogcounts <- assay(rlog(dds, blind=TRUE)) #blind = T will make the log transformation blind to the sample info)


a <- pca(rlogcounts, metadata = coldata1,center=T,scale=T)

biplot(a, showLoadings = TRUE, colby='Brain_region', shape='Sex',
       lab = NULL, pointSize = 5, sizeLoadingsNames = 7,ellipse=F,
       ntopLoadings = 10,max.overlaps = 30, boxedLoadingsNames = FALSE,
       legendPosition="right",legendLabSize = 20,axisLabSize=26,
       legendIconSize = 8,colLegendTitle=NULL,shapeLegendTitle=NULL)

setwd("./figures")
ggsave("PCA_biplot_removedF9BS.png", dpi = 300)


#****************************************************************************
#clustering heatmap to check if samples group by sex
library(DESeq2)
library("pheatmap")
library("dplyr")

setwd("D:/A.clarkii/clarkii_brain_DE")

cts <- read.csv("clarkii_brain_cts.csv",row.names="Geneid")
coldata <- read.csv("brain_metadata.csv",row.names = 1)

#convert decimal to int
cts1=round(cts)

cts2 <- data.frame(cts1[,c(-40)])
coldata1 <- coldata[c(-40),]

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata1) %in% colnames(cts2)) || all(colnames(cts2) %in% rownames(coldata1))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts2) == rownames(coldata1)))


# all variables in the design formula should be converted to factors
coldata1$Condition <- factor(coldata1$Condition)
coldata1$Sex <- factor(coldata1$Sex)
coldata1$Brain_region <- factor(coldata1$Brain_region)


dds <- DESeqDataSetFromMatrix(countData = cts2, colData = coldata1, design = ~Condition + Sex + Brain_region)

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
                        Sex = c(Male = "#6B8E23", Female = "#CC8899"),
                        Brain_region = c("Brain stem"    = "#5BC8E8",
                                         "Cerebellum"    = "#1FC4B0",
                                         "Diencephalon"  = "#8E9CE0",
                                         "Optic tectum"  = "#A0A018",
                                         "Telencephalon" = "#F26FD0"))


p <- pheatmap(cor,color=heat_colors, cellheight=35,cellwidth=40,cluster_rows=TRUE, show_rownames=T,show_colnames=T,legend=T,fontsize = 30,
              cluster_cols=T, annotation=select(coldata1,c(Condition,Sex,Brain_region)),,annotation_colors = annotation_cols)

setwd("./figures/")
png("clustering_heatmap_removedF9BS.png",res = 300,units="cm", width = 82, height = 63)
print(p)
dev.off()


#****************************************************************************
#LRT to see the effect of sex, brain_region and condition on gene expression 
library(DESeq2)

setwd("D:/A.clarkii/clarkii_brain_DE")

cts <- read.csv("clarkii_brain_cts.csv",row.names="Geneid")
coldata <- read.csv("brain_metadata.csv",row.names = 1)

#convert decimal to int
cts1=round(cts)

cts2 <- data.frame(cts1[,c(-40)])
coldata1 <- coldata[c(-40),]

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata1) %in% colnames(cts2)) || all(colnames(cts2) %in% rownames(coldata1))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts2) == rownames(coldata1)))


# all variables in the design formula should be converted to factors
coldata1$Condition <- factor(coldata1$Condition)
coldata1$Sex <- factor(coldata1$Sex)
coldata1$Brain_region <- factor(coldata1$Brain_region)


dds <- DESeqDataSetFromMatrix(countData = cts2, colData = coldata1, 
                              design = ~Sex + Condition + Brain_region)

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# total expressed genes in brain (passed the rowSums > 0 filter)
n_expressed_brain <- nrow(dds)
n_expressed_brain #24907

# testing for sex effect 
dds_LRT_sex <- DESeq(dds, test = "LRT", reduced = ~ Condition + Brain_region)
res_sex <- results(dds_LRT_sex)
res_sex
summary(res_sex)
results <- data.frame(res_sex)
results <- na.omit(results)

res_sig_padj <- results[results$padj < 0.05, ]

# number significant
n_sig_sex <- sum(res_sex$padj < 0.05, na.rm = TRUE)
n_sig_sex #365

# proportion of expressed genes
round(100 * n_sig_sex / n_expressed_brain, 1)  

setwd("./LRT_results/")
write.csv(results,file="LRT_results_sex_effect.csv")
write.csv(res_sig_padj,file="LRT_results_sex_effect_sig.csv")

# testing for condition effect 
dds_LRT_condition <- DESeq(dds, test = "LRT", reduced = ~ Sex + Brain_region)
res_condition <- results(dds_LRT_condition)
res_condition
summary(res_condition)
results <- data.frame(res_condition)
results <- na.omit(results)

res_sig_padj <- results[results$padj<0.05,]

setwd("./LRT_results/")
write.csv(results,file="LRT_results_condition_effect.csv")
write.csv(res_sig_padj,file="LRT_results_condition_effect_sig.csv")

# testing for Brain_region effect 
dds_LRT_BrainRegion <- DESeq(dds, test = "LRT", reduced = ~ Sex + Condition)
res_BrainRegion <- results(dds_LRT_BrainRegion)
res_BrainRegion
summary(res_BrainRegion)
results <- data.frame(res_BrainRegion)
results <- na.omit(results)

res_sig_padj <- results[results$padj<0.05,]

setwd("./LRT_results/")
write.csv(results,file="LRT_results_BrainRegion_effect.csv")
write.csv(res_sig_padj,file="LRT_results_BrainRegion_effect_sig.csv")

#****************************************************************************
# DE analysis - wald test
library(DESeq2)

setwd("D:/A.clarkii/clarkii_brain_DE")

cts <- read.csv("clarkii_brain_cts.csv",row.names="Geneid")
coldata <- read.csv("brain_metadata.csv",row.names = 1)

#convert decimal to int
cts1=round(cts)

cts2 <- data.frame(cts1[,c(-40)])
coldata1 <- coldata[c(-40),]

#### BS
cts_BS <- cts2[,c(1,5,10,15,20,25,30,35)]
coldata_BS <- coldata1[c(1,5,10,15,20,25,30,35),]

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata_BS) %in% colnames(cts_BS)) || all(colnames(cts_BS) %in% rownames(coldata_BS))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts_BS) == rownames(coldata_BS)))

dds_BS <- DESeqDataSetFromMatrix(countData = cts_BS, colData = coldata_BS,design = ~ Sex + Condition) 
levels(dds_BS$Condition)    
levels(dds_BS$Sex)

dds_BS$condition <- relevel(dds_BS$Condition, ref = "Control")


# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds_BS)) > 0
dds_BS <- dds_BS[keep,]

dds_BS <- DESeq(dds_BS)
resultsNames(dds_BS)

resBS <- results(dds_BS, alpha = 0.05, contrast = c("Condition","Interaction","Control"), test ="Wald")#So up/downregulation refers to interaction i.e upregulated in interaction or downregulated in interaction
summary(resBS)
resBS

resultsBS <- data.frame(resBS)
resultsBS <- na.omit(resultsBS)
sig_results_BS <- resultsBS[resultsBS$padj<0.05, ]

setwd("./DE_results_wald_test/")
write.csv(resultsBS,file="DE_results_all_BS.csv")
write.csv(sig_results_BS,file="DE_results_sig_BS.csv")

#### OT
cts_OT <- cts2[,c(3,7,12,17,22,27,32,37,41)]
coldata_OT <- coldata1[c(3,7,12,17,22,27,32,37,41),]

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata_OT) %in% colnames(cts_OT)) || all(colnames(cts_OT) %in% rownames(coldata_OT))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts_OT) == rownames(coldata_OT)))

dds_OT <- DESeqDataSetFromMatrix(countData = cts_OT, colData = coldata_OT,design = ~ Sex + Condition) 
levels(dds_OT$Condition)    
levels(dds_OT$Sex)

dds_OT$condition <- relevel(dds_OT$Condition, ref = "Control")


# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds_OT)) > 0
dds_OT <- dds_OT[keep,]

dds_OT <- DESeq(dds_OT)
resultsNames(dds_OT)

resOT <- results(dds_OT, alpha = 0.05, contrast = c("Condition","Interaction","Control"), test ="Wald")#So up/downregulation refers to interaction i.e upregulated in interaction or downregulated in interaction
summary(resOT)
resOT

resultsOT <- data.frame(resOT)
resultsOT <- na.omit(resultsOT)
sig_results_OT <- resultsOT[resultsOT$padj<0.05, ]

setwd("./DE_results_wald_test/")
write.csv(resultsOT,file="DE_results_all_OT.csv")
write.csv(sig_results_OT,file="DE_results_sig_OT.csv")


#### TE
cts_TE <- cts2[,c(4,8,13,18,23,28,33,38,42)]
coldata_TE <- coldata1[c(4,8,13,18,23,28,33,38,42),]

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata_TE) %in% colnames(cts_TE)) || all(colnames(cts_TE) %in% rownames(coldata_TE))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts_TE) == rownames(coldata_TE)))

dds_TE <- DESeqDataSetFromMatrix(countData = cts_TE, colData = coldata_TE,design = ~ Sex + Condition) 
levels(dds_TE$Condition)    
levels(dds_TE$Sex)

dds_TE$condition <- relevel(dds_TE$Condition, ref = "Control")


# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds_TE)) > 0
dds_TE <- dds_TE[keep,]

dds_TE <- DESeq(dds_TE)
resultsNames(dds_TE)

resTE <- results(dds_TE, alpha = 0.05, contrast = c("Condition","Interaction","Control"), test ="Wald")#So up/downregulation refers to interaction i.e upregulated in interaction or downregulated in interaction
summary(resTE)
resTE

resultsTE <- data.frame(resTE)
resultsTE <- na.omit(resultsTE)
sig_results_TE <- resultsTE[resultsTE$padj<0.05, ]

setwd("./DE_results_wald_test/")
write.csv(resultsTE,file="DE_results_all_TE.csv")
write.csv(sig_results_TE,file="DE_results_sig_TE.csv")


#### CB
cts_CB <- cts2[,c(9,14,19,24,29,34,39,43)]
coldata_CB <- coldata1[c(9,14,19,24,29,34,39,43),]

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata_CB) %in% colnames(cts_CB)) || all(colnames(cts_CB) %in% rownames(coldata_CB))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts_CB) == rownames(coldata_CB)))

dds_CB <- DESeqDataSetFromMatrix(countData = cts_CB, colData = coldata_CB,design = ~ Sex + Condition) 
levels(dds_CB$Condition)    
levels(dds_CB$Sex)

dds_CB$condition <- relevel(dds_CB$Condition, ref = "Control")


# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds_CB)) > 0
dds_CB <- dds_CB[keep,]

dds_CB <- DESeq(dds_CB)
resultsNames(dds_CB)

resCB <- results(dds_CB, alpha = 0.05, contrast = c("Condition","Interaction","Control"), test ="Wald")#So up/downregulation refers to interaction i.e upregulated in interaction or downregulated in interaction
summary(resCB)
resCB

resultsCB <- data.frame(resCB)
resultsCB <- na.omit(resultsCB)
sig_results_CB <- resultsCB[resultsCB$padj<0.05, ]

setwd("./DE_results_wald_test/")
write.csv(resultsCB,file="DE_results_all_CB.csv")
write.csv(sig_results_CB,file="DE_results_sig_CB.csv")


#### DE
cts_DE <- cts2[,c(2,6,11,16,21,26,31,36,40)]
coldata_DE <- coldata1[c(2,6,11,16,21,26,31,36,40),]

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata_DE) %in% colnames(cts_DE)) || all(colnames(cts_DE) %in% rownames(coldata_DE))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts_DE) == rownames(coldata_DE)))

dds_DE <- DESeqDataSetFromMatrix(countData = cts_DE, colData = coldata_DE,design = ~ Sex + Condition) 
levels(dds_DE$Condition)    
levels(dds_DE$Sex)

dds_DE$condition <- relevel(dds_DE$Condition, ref = "Control")


# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds_DE)) > 0
dds_DE <- dds_DE[keep,]

dds_DE <- DESeq(dds_DE)
resultsNames(dds_DE)

resDE <- results(dds_DE, alpha = 0.05, contrast = c("Condition","Interaction","Control"), test ="Wald")#So up/downregulation refers to interaction i.e upregulated in interaction or downregulated in interaction
summary(resDE)
resDE

resultsDE <- data.frame(resDE)
resultsDE <- na.omit(resultsDE)
sig_results_DE <- resultsDE[resultsDE$padj<0.05, ]

setwd("./DE_results_wald_test/")
write.csv(resultsDE,file="DE_results_all_DE.csv")
write.csv(sig_results_DE,file="DE_results_sig_DE.csv")


#****************************************************************************
# merge  DE gene list with GO functional annotation
setwd("D:/A.clarkii/clarkii_brain_DE/DE_results_wald_test/")

GO <- read.csv("A.clarkii_GO_Annotation_edited.csv")

BS <- read.csv("DE_results_sig_final_BS.csv")
merge <- merge(BS,GO,by="Geneid",all.x=T)
write.csv(merge,file="DE_results_sig_final_BS_annotation.csv")

CB <- read.csv("DE_results_sig_final_CB.csv")
merge <- merge(CB,GO,by="Geneid",all.x=T)
write.csv(merge,file="DE_results_sig_final_CB_annotation.csv")

DE <- read.csv("DE_results_sig_final_DE.csv")
merge <- merge(DE,GO,by="Geneid",all.x=T)
write.csv(merge,file="DE_results_sig_final_DE_annotation.csv")

TE <- read.csv("DE_results_sig_final_TE.csv")
merge <- merge(TE,GO,by="Geneid",all.x=T)
write.csv(merge,file="DE_results_sig_final_TE_annotation.csv")

OT <- read.csv("DE_results_sig_final_OT.csv")
merge <- merge(OT,GO,by="Geneid",all.x=T)
write.csv(merge,file="DE_results_sig_final_OT_annotation.csv")

############################################## Venn Diagrams  
# http://www.interactivenn.net/index2.html



##### merge individual brain regions to get a combined mega file

setwd("D:/A.clarkii/clarkii_brain_DE/DE_results_wald_test/final_files")

BS <- read.csv("DE_results_sig_final_BS_annotation.csv")
CB <- read.csv("DE_results_sig_final_CB_annotation.csv")
DE <- read.csv("DE_results_sig_final_DE_annotation.csv")
OT <- read.csv("DE_results_sig_final_OT_annotation.csv")
TE <- read.csv("DE_results_sig_final_TE_annotation.csv")

df_list <- list(BS,CB,DE,OT,TE)

merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
write.csv(merge,file="clarkii_DE_final_Allcombined.csv",row.names=F)

############################################## Heatmap 
library(ggplot2)
library(ggbeeswarm) # for beeswarm plot
library(ggrepel)

setwd("D:/A.clarkii/clarkii_brain_DE/DE_results_wald_test/final_files")

file <- read.csv("file_for_box.csv")
file1 <- na.omit(file)

ggplot(file1,aes(x=Brain_Region,y=log2FoldChange,color=Brain_Region)) + 
  geom_quasirandom(size=4) + theme_bw() + theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=18)) + geom_text_repel(data=file1,aes(label=Gene_symbol_NCBI),position=position_quasirandom(),size=4) +
  theme(axis.text=element_text(size=18),axis.title.x = element_blank(), axis.title.y=element_text(size=22)) +
  scale_fill_manual(values=c("#F8766D","A3A500","00BF7D","00B0F6","E76BF3")) +
  geom_hline(yintercept=0, linetype="dashed",color="black",size=1)

setwd("D:/A.clarkii/clarkii_brain_DE/figures")
ggsave(file="brain_regions_beeswarm.png",dpi=300,width=20, height=20,units="cm")

# ============================================
# Test overlap between sig sex effect genes from LRT and the sig DE genes from Wald test
# ============================================

setwd("D:/A.clarkii/clarkii_brain_DE")

lrt_sex <- read.csv("LRT_results/LRT_results_sex_effect_sig.csv",row.names=1)
sex <- rownames(lrt_sex)

de_genes <- read.csv("DE_results_wald_test/final_files/clarkii_DE_final_Allcombined.csv",row.names=1)
brain_de <- rownames(de_genes)

overlap <- intersect(sex, brain_de)
