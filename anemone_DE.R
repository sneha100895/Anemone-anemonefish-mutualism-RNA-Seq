setwd("D:/A.clarkii/anemone_DE/RSEM")

#Read gene.results from RSEM for each sample and merge to get a dataframe of combined cts
NFA1 <- read.table("NFA1RNA.genes.results",header=T)
NFA1_1 <- NFA1[,c(1,5)]

NFA2 <- read.table("NFA2RNA.genes.results",header=T)
NFA2_1 <- NFA2[,c(1,5)]

NFA3 <- read.table("NFA3RNA.genes.results",header=T)
NFA3_1 <- NFA3[,c(1,5)]

NFA4 <- read.table("NFA4RNA.genes.results",header=T)
NFA4_1 <- NFA4[,c(1,5)]

NFA5 <- read.table("NFA5RNA.genes.results",header=T)
NFA5_1 <- NFA5[,c(1,5)]

F5T <- read.table("F5TRNA.genes.results",header=T)
F5T_1 <- F5T[,c(1,5)]

F6T <- read.table("F6TRNA.genes.results",header=T)
F6T_1 <- F6T[,c(1,5)]

F7T <- read.table("F7TRNA.genes.results",header=T)
F7T_1 <- F7T[,c(1,5)]

F8T <- read.table("F8TRNA.genes.results",header=T)
F8T_1 <- F8T[,c(1,5)]

F9T <- read.table("F9TRNA.genes.results",header=T)
F9T_1 <- F9T[,c(1,5)]

df_list <- list(NFA1_1, NFA2_1, NFA3_1, NFA4_1, NFA5_1, F5T_1, F6T_1, F7T_1, 
                F8T_1, F9T_1)
merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
write.csv(merge,file="anemone_cts.csv",row.names=F)


#****************************************************************************
#PCA from rlog transformed counts
# link to guide: https://tavareshugo.github.io/data-carpentry-rnaseq/03_rnaseq_pca.html
library(DESeq2)
library(tidyverse)

setwd("D:/A.clarkii/anemone_DE")
cts <- read.csv("anemone_cts.csv",row.names="gene_id")
coldata <- read.csv("anemone_metadata.csv",row.names=1)

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata) %in% colnames(cts)) || all(colnames(cts) %in% rownames(coldata))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts) == rownames(coldata)))

# all variables in the design formula should be converted to factors
coldata$Condition <- factor(coldata$Condition)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

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
setwd("D:/A.clarkii/anemone_DE")
#extract the PC score for each sample from the prcomp object
coldata <- read.csv("anemone_metadata.csv") 
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
  ggplot(aes(x = PC1, y = PC2, colour = Condition,label=Sample)) +
  geom_point(size=5) + geom_text(hjust=0.9, vjust=1.6,size=4)

p

#np label
p <- pc_scores %>% 
  # join with "sample_info" table 
  full_join(coldata, by = "Sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = Condition)) +
  geom_point(size=5) 
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

p <- p + stat_conf_ellipse(bary=F, level = 0.95,npoint = 1000, geom = "polygon",alpha=0.1, aes (fill=Condition))
p

setwd("./figures")
ggsave("PCA_no_label_anemone.png", dpi = 300)

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

library(ggrepel)
#plot 
loadings_plot <- ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "black",linewidth=1) +
  geom_text_repel(aes(x = PC1, y = PC2, label = gene),
                  size = 8,hjust=0.9, vjust=1,max.overlaps = 15) + theme_bw() +
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
ggsave("PCA_loadings_top10genesForeachPC.png", dpi = 300,height=40, width=40,units = "cm")

# Combining PC scores and PC loadings to see correlation - biplot
#For example, geneA, with a negative score of -0.00605 (loading score), will be associated with whatever samples are shifted toward the negative end of your PC1 on a bi-plot.

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('PCAtools')

library(PCAtools)
setwd("D:/A.clarkii/anemone_DE")

cts <- read.csv("anemone_cts.csv",row.names="gene_id")
coldata <- read.csv("anemone_metadata.csv",row.names = 1)

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata) %in% colnames(cts)) || all(colnames(cts) %in% rownames(coldata))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts) == rownames(coldata)))

# all variables in the design formula should be converted to factors
coldata$Condition <- factor(coldata$Condition)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

#The rlog function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size. 
#The transformation is useful when checking for outliers. 
#Note that neither rlog transformation nor the VST are used by the differential expression estimation in DESeq, which always occurs on the raw count data, through generalized linear modeling which incorporates knowledge of the variance-mean dependence. The rlog transformation and VST are offered as separate functionality which can be used for visualization, clustering or other machine learning tasks. 

rlogcounts <- assay(rlog(dds, blind=TRUE)) #blind = T will make the log transformation blind to the sample info)


a <- pca(rlogcounts, metadata = coldata,center=T,scale=T)

biplot(a, showLoadings = T, colby='Condition', ellipse=T,
       lab = NULL, pointSize = 7, sizeLoadingsNames = 7,
       ntopLoadings = 10,max.overlaps = 20, boxedLoadingsNames = FALSE,
       legendPosition="bottom",legendLabSize = 20,axisLabSize=26,
       legendIconSize = 8,colLegendTitle=NULL,shapeLegendTitle=NULL,
       xlim=c(-300,300),ylim=c(-400,400))

setwd("./figures")
ggsave("PCA_biplot.png", dpi = 300,units="cm",height=80,width=70)


#****************************************************************************
#clustering heatmap 
library(DESeq2)
library("pheatmap")
library("dplyr")

setwd("D:/A.clarkii/anemone_DE")
cts <- read.csv("anemone_cts.csv",row.names="gene_id")
coldata <- read.csv("anemone_metadata.csv",row.names=1)

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata) %in% colnames(cts)) || all(colnames(cts) %in% rownames(coldata))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts) == rownames(coldata)))


# all variables in the design formula should be converted to factors
coldata$Condition <- factor(coldata$Condition)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~Condition)

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

annotation_cols <- list(Condition = c(Control = "#808080", Interaction = "#FEBE00")) 
                        

p <- pheatmap(cor,color=heat_colors,cellheight=35,cellwidth=40, cluster_rows=TRUE, show_rownames=T,show_colnames=T,legend=T,fontsize = 20,
              cluster_cols=T, annotation=select(coldata,c(Condition)),annotation_colors = annotation_cols)

setwd("./figures/")
png("clustering_heatmap.png",res = 300,units="cm", width = 30, height = 20)
print(p)
dev.off()

#****************************************************************************
# DE analysis - wald test
library(DESeq2)

setwd("D:/A.clarkii/anemone_DE")
cts <- read.csv("anemone_cts.csv",row.names="gene_id")
coldata <- read.csv("anemone_metadata.csv",row.names=1)

# check if all colnames of count matrix are in the rownames of sample info; if false exit the script
(all(rownames(coldata) %in% colnames(cts)) || all(colnames(cts) %in% rownames(coldata))) 

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts) == rownames(coldata)))

# all variables in the design formula should be converted to factors
coldata$Condition <- factor(coldata$Condition)


dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, 
                              design = ~Condition)

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# Explicitly set the factor levels
dds$condition <- relevel(dds$Condition, ref = "Control")
dds <- DESeq(dds)
counts <- counts(dds, normalized=TRUE)
write.csv(counts, file="normalised_cts.csv", row.names=T)

resultsNames(dds)

# get the results for interaction vs ctrl
res <- results(dds, alpha = 0.05,test="Wald", contrast = c("Condition","Interaction","Control")) #So up/downregulation refers to interaction i.e upregulated in teraction or downregulated in interaction
summary(res)
res

results <- data.frame(res)
results <- na.omit(results)
sig_results <- results[results$padj<0.05, ]

setwd("./DE_analysis_results_waldTest/")
write.csv(results,file="DE_results_anemone.csv")
write.csv(sig_results,file="DE_results_anemone_sig.csv")

#****************************************************************************
# Add annotation info to the sig DE genes 

setwd("D:/A.clarkii/anemone_DE/DE_analysis_results_waldTest/")

genes <- read.csv("DE_results_anemone_sig_final.csv")

annot <- read.csv("Annotation_GO.csv")

merge <- merge(genes,annot,by="SeqName",all.x=T)

write.csv(merge,file="DE_results_anemone_sig_final_annot.csv",
          row.names=F)

#****************************************************************************
#figure for GO functions

library(ggplot2)
library(ggh4x)

setwd("D:/A.clarkii/anemone_DE/enrichment")

file <- read.csv("summary_functional_classification_for_fig.csv")

de <- read.csv("D:/A.clarkii/anemone_DE/DE_analysis_results_waldTest/DE_results_anemone_sig_final_annot.csv")

file1 <- merge(file,de,by="SeqName",all.x=T)

file1 <- file1[order(file1$Functional.Classification), ]

write.csv(file1,file="summary_functional_classification_for_fig_withDEseqRes.csv",row.names=F)

design <- "ADEF
           GHIJ
           KBC#"
           
# design function helps control location of panels - the hashtag is for space

library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)
heat_colors <- brewer.pal(9, "Greys")

p <- ggplot(file1,aes(x=SeqName,y=log2FoldChange,color=padj)) + 
  geom_point(size=3) + theme_classic() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  facet_manual(vars(Functional.Classification), design = design) + 
  theme(panel.background = element_rect(fill = 'white', color = 'black')) +
  theme(axis.text=element_text(size=18),axis.title.x = element_text(size=22), axis.title.y=element_text(size=22)) + 
  ylab("log2 Fold Change") + xlab("Gene ID") + 
  theme(legend.text=element_text(size=18)) +
  theme(legend.title=element_text(size=22)) +
  theme(strip.text = element_text(size = 20))

p
setwd("D:/A.clarkii/anemone_DE/figures")
ggsave("GO_functions_edited.png", dpi=300,width=50,height=25, units="cm")



