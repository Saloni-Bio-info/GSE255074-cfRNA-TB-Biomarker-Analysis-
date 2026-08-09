#GSE255074

#-------loading the packages--------
library(limma)
library(edgeR)
library(GEOquery)
library(biomaRt)
library(ggplot2)

#for microarray data- we extract expression data, 
#for RNAseq data - we download and load it manually

#-------counts file-------
#read.csv for csv file
#read.delim for tsv or txt file or gz (compressed) file
counts<-read.delim(file.choose(),
                   row.names = 1,
                   check.names = FALSE)

#--------metadata file--------
#obtain the sample metadata
gse<-getGEO("GSE255074", GSEMatrix = TRUE)
gset<-gse[[1]]

#extract the phenotype data
metadata<-pData(gset)

head(metadata)
dim(metadata)
colnames(metadata)

#-------count matrix-------
#genes should be rows and samples should be columns in counts file
dim(counts)
head(counts)
#if not, then transpose 
#counts<-t(counts)

#check for missing values
sum(is.na(counts))

class(counts)
str(counts)

#removing ensembl version numbers (genes name in counts)
rownames(counts)<-gsub("\\..*", "", rownames(counts))
View(counts)

#remove low counts genes - by calculating the total count for each gene across all samples.
counts<- counts[rowSums(counts)>10, ]

#Check metadata
View(metadata)

#table of respective sample ID (title) column in metadata - inspect the values and frequency in title column 
table(metadata$title)

#ensure samples order matches in counts and metadata
all(colnames(counts)==metadata$title)

#if FALSE=order is different
#then, create a column named sample in metadata containing clean, matching title to counts file
metadata$Sample <- sub(".*patient: ", "", metadata$title)

#ensure samples order matches in counts and metadata, should return TRUE
all(colnames(counts)==metadata$title)

#Create a group variable (group is the column of metadata which contains control vs disease information)
group <- factor(metadata$characteristics_ch1.33)
table(group)

#value 2 confirms there are two conditions to compare. 
length(unique(group))

#make an object dge as DGEList will containing everything needed for Differential Gene Expression analysis. 
dge <- DGEList(counts = counts, group = group)

#filter low expression genes - as RNA-seq datasets contain many genes with extremely low counts.
keep <- filterByExpr(dge, group = group)
sum(keep)
#keep only the better expressed genes as dge object for Differential Gene Expression analysis.
dge <- dge[keep, , keep.lib.sizes = FALSE]
dim(dge)

# ==========================
# Boxplot BEFORE normalization
# ==========================

boxplot(
  cpm(dge, log = TRUE),
  las = 2,
  main = "Before TMM Normalization"
)

#calculate TMM normalization factors
dge <- calcNormFactors(dge)
dge$samples

# ==========================
# Boxplot AFTER normalization
# ==========================

boxplot(
  cpm(dge, log = TRUE, normalized.lib.sizes = TRUE),
  las = 2,
  main = "After TMM Normalization"
)

#craete a design matrix as design matrix is the instruction sheet that tells the statistical model -
#that These samples are in Group A, and these samples are in Group B (TB positive, negative as 0 and 1).
design <- model.matrix(~0 + group)
design

#voom transformation - with quality weights - it converts counts to logCPM, assign precision weight to address heteroscedasticity.
#lmFit() assumes consistent variability
vm <- voomWithQualityWeights(dge, design, plot = TRUE)

# ------- Mean-variance trend-------Which genes are noisy?
# Low-expression genes are noisier.
# Red curve estimates the variance trend

# -------Sample-specific weights-------Which samples should I trust more?
# Higher bar = more reliable sample
# Lower bar = less reliable sample


#ERROR- syntactically valid names in R so we do this :
colnames(design) <- c("TB_negative", "TB_positive")

#lmFit() fits a linear regression model to every gene
fit <- lmFit(vm, design)

#tells limma what are we comparing 
contrast.matrix <- makeContrasts(
  TBvsControl = TB_positive - TB_negative,
  levels = design
)

contrast.matrix
fit2 <- contrasts.fit(fit, contrast.matrix)

#statistics- Compute moderated t-statistics and p-values
#fit2 is basically fit1 + the comparision that we want
fit2 <- eBayes(fit2)

#lmFit() builds the model
#contrasts.fit() chooses the comparison
#eBayes() tests whether that comparison is statistically significant

#extract all the results in a table
results <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)

dim(results)
head(results)

#filtering significant genes
deg <- results[
  abs(results$logFC) > 1 &
    results$adj.P.Val < 0.05,
]

dim(deg)
nrow(deg)

up <- deg[deg$logFC > 1, ]
down <- deg[deg$logFC < -1, ]

nrow(up)
nrow(down)

#saving the files
write.csv(results, "All_DEGs.csv", row.names = TRUE)

write.csv(deg, "Significant_DEGs.csv", row.names = TRUE)

write.csv(up, "Upregulated_Genes.csv", row.names = TRUE)

write.csv(down, "Downregulated_Genes.csv", row.names = TRUE)

#check
list.files()

#visualisation
# -------volcano plot-------
# Create a copy of results for plotting
volcano <- results

# Classify genes
volcano$Significance <- "Not Significant"

volcano$Significance[
  volcano$adj.P.Val < 0.05 &
    volcano$logFC > 1
] <- "Upregulated"

volcano$Significance[
  volcano$adj.P.Val < 0.05 &
    volcano$logFC < -1
] <- "Downregulated"

# Check number of genes in each category
table(volcano$Significance)

# Draw volcano plot
ggplot(
  volcano,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Significance
  )
) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )
#save the image 
ggsave("Volcano_Plot.png", width = 8, height = 6, dpi = 300)


# Gene Set Enrichment Analysis (GSEA)
#we used limma, so we'll use the logFC column from All_DEGs.csv to rank genes for GSEA
#loading the libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

#checking the names of gene
colnames(results)
head(rownames(results))

#GSEA using clusterProfiler cannot directly use Ensembl IDs for GO analysis
#We first need to convert them to Entrez IDs using the org.Hs.eg.db annotation package

# Convert Ensembl IDs to Entrez IDs
gene_conversion <- bitr(
  rownames(results),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
#inspecting the conversion table
head(gene_conversion)

# Add Ensembl IDs as a column
results$ENSEMBL <- rownames(results)

# Merge DEG results with converted gene IDs
gsea_data <- merge(gene_conversion, results, by = "ENSEMBL")

#checking the result
head(gsea_data)

#creating ranked gene list

gene_list <- gsea_data$logFC
names(gene_list) <- gsea_data$ENTREZID

gene_list <- sort(gene_list, decreasing = TRUE)

# Remove duplicate Entrez IDs
gene_list <- gene_list[!duplicated(names(gene_list))]

#GSEA using GO
gsea_result <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

#making the dot plot
dotplot(gsea_result, showCategory = 20)

#saving the results
gsea_table <- as.data.frame(gsea_result)

write.csv(
  gsea_table,
  "GSEA_GO_BP_Results.csv",
  row.names = FALSE
)

#save the plot
png("GSEA_DotPlot.png", width = 1800, height = 1500, res = 300)

dotplot(gsea_result, showCategory = 15)

dev.off()
