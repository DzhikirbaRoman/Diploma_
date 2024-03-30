library(Rtsne)
#BiocManager::install("affy")
library(affy)
#BiocManager::install("ArrayExpress")
library(ArrayExpress)
# library(arrayQualityMetrics)
library(sva)
library(clusterProfiler)
#BiocManager::install("limma")
library(limma)
library(stringr)
library(RColorBrewer)
library(factoextra)
library(stats)
library(usedist)
library(easypackages)
library(data.table)
library(openxlsx)
library(pheatmap)
library(DOSE)
library(enrichplot)
library(ggupset)
#install.packages("xlsx")
library("xlsx")
#install.packages("readxl")
library("readxl")


getwd()
setwd('/Users/RTIntelektFBT/Desktop/Roman_Project/try_from_almostZERO')
getwd()

# DEGS analysis
sub_exprs = read.table("sub_exprs.tsv", sep="\t")
#View(sub_exprs)
sub_pdata = read.table("sub_sub_pdata.tsv", sep="\t", header = TRUE)
#sub_pdata = read.csv("sub_pdata.csv", sep=",")
View(sub_pdata)
#sub_pdata[,"etimated_sex_unimodality_p"]
#dim(sub_exprs)
#dim(sub_pdata)
table(sub_pdata[,c( "trim_term")])
as.factor(sub_pdata$trim_term)
sub_pdata$trim_term


fit_mod = model.matrix(~ 0 + as.factor(sub_pdata$trim_term), data=sub_pdata)
#colnames(fit_mod)
colnames(fit_mod) = c('i','ii')
fit <- lmFit(sub_exprs, fit_mod)  # fit each probeset to model
contrast.matrix <- makeContrasts(ii-i, levels=fit_mod)
fit2 <- contrasts.fit(fit, contrast.matrix)
efit <- eBayes(fit2)        # empirical Bayes adjustment

t = topTable(efit, number = nrow(sub_exprs))
nrow(t)
#d = t
#View(d)
#d = d[d$adj.P.Val<0.05,]
#d = d[which(abs(d$logFC)>1),]
#nrow(d)
diff <- t
#nrow(diff)
#rownames(diff)
#View(diff)

library(org.Hs.eg.db)
#keys(org.Hs.eg.db)
#View(sub_exprs)

sel = AnnotationDbi::select(org.Hs.eg.db, rownames(diff), columns = c("ENTREZID","GENENAME"), keytype="SYMBOL")
sel = sel[!is.na(sel$ENTREZID),]
#View(sel)

sel2 = AnnotationDbi::select(org.Hs.eg.db, rownames(sub_exprs), columns =  c("ENTREZID","GENENAME"), keytype="SYMBOL")
sel2 = sel2[!is.na(sel2$ENTREZID),]

length(unique(rownames(sub_exprs)))==length(rownames(sub_exprs))
length(sel2$ENTREZID) == length(unique(sel2$ENTREZID))

length(sel2$ENTREZID)
length(unique(sel2$ENTREZID))


NA %in% sel2$ENTREZID
NA %in% rownames(sub_exprs)
sel2$ENTREZID

sub_exprs <- sub_exprs[sel2$SYMBOL,]
#View(sub_exprs)

nrow(sub_exprs)== nrow(sel2)

!(FALSE %in% (rownames(sub_exprs) == sel2$SYMBOL))
nrow(diff)
#sel2$SYMBOL

# it is true
#rownames(sub_exprs) = sel2$SYMBOL 

#View

#write.csv(sub_exprs, file.path(datapath, "sub_exprs_tmp.csv"))

#View(sel)
#View(diff)
nrow(sel)
nrow(diff)

merged_sel = merge(sel,diff, by.x = "SYMBOL", by.y = "row.names")


nrow(merged_sel)
nrow(sub_exprs)

#View(merged_sel)

length(levels(as.factor(sub_pdata$trim_term)))



difexp_exprs11 = sub_exprs[which(rownames(sub_exprs) %in% merged_sel$SYMBOL),]
nrow(difexp_exprs11)

rownames(merged_sel)
length(rownames(merged_sel))
length(merged_sel$SYMBOL)

rownames(sub_exprs)
length(rownames(sub_exprs))

#View(merged_sel)

missing_rows <- setdiff(rownames(sub_exprs),merged_sel$SYMBOL)
missing_rows
missing_rows <- setdiff(merged_sel$SYMBOL,rownames(sub_exprs))
missing_rows


merged_sel[merged_sel$SYMBOL %in% c("MMD2","HBD","TEC"),]

merged_sel[merged_sel$ENTREZID == 100187828 ,]$SYMBOL <- "HBD.1"
merged_sel[merged_sel$ENTREZID == 100505381 ,]$SYMBOL <- "MMD2.1"
merged_sel[merged_sel$ENTREZID == 100124696 ,]$SYMBOL <- "TEC.1"

sub_exprs[rownames(sub_exprs) %in% c("MMD2","HBD","TEC","MMD2.1","HBD.1","TEC.1"),]



missing_rows <- setdiff(rownames(sub_exprs),merged_sel$SYMBOL)
missing_rows
missing_rows <- setdiff(merged_sel$SYMBOL,rownames(sub_exprs))
missing_rows


#View(sub_exprs)
trim = levels(as.factor(sub_pdata$trim_term))[1]
trim
difexp_exprs = sub_exprs[which(rownames(sub_exprs) %in% merged_sel$SYMBOL),]
#View(difexp_exprs)
#View(merged_sel)
#merged_sel1 <- merged_sel[merged_sel$SYMBOL %in% rownames(difexp_exprs),]
merged_sel <- merged_sel[match(rownames(difexp_exprs),merged_sel$SYMBOL),]

#!(FALSE %in% (merged_sel$SYMBOL == rownames(difexp_exprs)))
trim_sample_names = sub_pdata[sub_pdata$trim_term==trim,]$arraydatafile_exprscolumnnames
#trim_sample_names

difexp_exprs_trim = difexp_exprs[,which(colnames(difexp_exprs) %in% trim_sample_names)]

difexp_exprs_trim
rowMeans(difexp_exprs_trim)
merged_sel[,paste(trim,"Average",sep = " ")] = rowMeans(difexp_exprs_trim)
#View(merged_sel)


trim = levels(as.factor(sub_pdata$trim_term))[2]
trim
#  trim
#  sub_exprs[which(rownames(sub_exprs) %in% merged_sel$SYMBOL),]
difexp_exprs = sub_exprs[which(rownames(sub_exprs) %in% merged_sel$SYMBOL),]
#View(difexp_exprs)
#View(merged_sel)

#  setdiff(rownames(difexp_exprs), merged_sel$SYMBOL)
#  merged_sel[is.na(merged_sel$SYMBOL),]
#  nrow(difexp_exprs)
#  nrow(merged_sel)

!(FALSE %in% (merged_sel$SYMBOL == rownames(difexp_exprs)))
trim_sample_names = sub_pdata[sub_pdata$trim_term==trim,]$arraydatafile_exprscolumnnames

difexp_exprs_trim = difexp_exprs[,which(colnames(difexp_exprs) %in% trim_sample_names)]
merged_sel[,paste(trim,"Average",sep = " ")] = rowMeans(difexp_exprs_trim)
View(merged_sel)

View(sub_exprs)
### DONEEE !!!!!
write.xlsx(merged_sel,"difexp_1-2.xlsx")
write.xlsx(sub_exprs, "sub_exprs_mapped_1-2.xlsx")
#View(difexp_1_2)
# Set the first column as row names
#rownames(difexp_1_2) <- difexp_1_2[[1]]




#
#### HEATMAP Creation

library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)

# sub_exprs_mapped_1-2.xlsx is sub_expr from file, after I processed it 
expression <- read_excel("sub_exprs_mapped_1-2.xlsx", sheet = "Sheet1")
expression
# 1_2_metadata.tsv is sub_sub_P-data from file
metadata <- read.table("1_2_metadata.tsv", header = TRUE, sep = "\t")
# 1_2_DEGS is difexp_1-2.xlsx from file and with filtering non signficant
#DEGs <- read.table("1_2_DEGs.tsv", header = TRUE, sep = "\t")
#DEGs <- DEGs[order(DEGs$logFC, decreasing = TRUE),]

DEGs <- read_excel("difexp_1-2.xlsx", sheet = "Sheet1")
#df <- read.csv(paste0(in_path, 'expression2.csv'), row.names = 1)
# Remove the first column from the data frame
DEGs <- DEGs[-1]
# filter DEGS
DEGs <- DEGs[DEGs$logFC > 0 & DEGs$adj.P.Val < 0.05 | DEGs$logFC < 0 & DEGs$adj.P.Val < 0.05, ]

DEGs <- DEGs[order(DEGs$logFC, decreasing = TRUE),]
#View(DEGs)

# Assuming 'DEGs' is your original dataframe
# Get the first 50 rows
first_50 <- DEGs[1:50, ]

# Get the last 50 rows
last_50 <- DEGs[(nrow(DEGs) - 49):nrow(DEGs), ]

# Combine the first and last 50 rows into a new dataframe
combined_df <- rbind(first_50, last_50)
combined_df
# Now 'combined_df' contains the first 50 and last 50 rows of 'DEGs'


#expression
#str(expression)
expression <- as.data.frame(expression)
rownames(expression) <- expression[, 1]
expression <- expression[-1]
#View(expression)

## PCA
pca <- prcomp(
  t(expression), # transpose our data frame to obtain PC scores for samples, not genes
  scale = TRUE # we want the data scaled to have unit variance for each gene
)
#head(pca$x[, 1:5])
pca_summary <- summary(pca)
#pca_summary
# Now access the importance information for the first 5 PCs
#pca_summary$importance[, 1:5]

# Make the first two PCs into a data frame for plotting with `ggplot2`
pca_df <- data.frame(pca$x[, 1:2]) %>%
  # Turn samples IDs stored as row names into a column
  tibble::rownames_to_column("arraydatafile_exprscolumnnames") %>%
  # Bring only the variables that we want from the metadata into this data frame
  # here we are going to join by `refinebio_accession_code` values
  dplyr::inner_join(
    dplyr::select(metadata, arraydatafile_exprscolumnnames,Gestational.Age.Appr ,Combined.Fetus.Sex ,trim_term),
    by = "arraydatafile_exprscolumnnames"
  )
#rm(pca_df)
#str(pca_df)
#pca_df1
metadata$Gestational.Age.Appr
#pca_summary$importance[2, 1:2]

PC1_2_pVar <- unname(c(pca_summary$importance[2, 1:2]))
#PC1_2_pVar
# Make a scatterplot using `ggplot2` functionality
pca_plot <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = trim_term, # label points with different colors for each `subgroup`
    shape = Combined.Fetus.Sex
  )
) +
  geom_point(size = 4) + # Plot individual points to make a scatterplot
  geom_text(aes(label = Gestational.Age.Appr),vjust = -0.6 ) +
  theme_classic() +# Format as a classic-looking plot with no gridlines
  labs(x = paste("PC1 (", PC1_2_pVar[1]*100, "%)", sep = ""),
       y = paste("PC2 (", PC1_2_pVar[2]*100, "%)", sep = ""))
# Print out the plot here
pca_plot
###



scaled_counts<- t(apply(expression, 1, scale)) #center and scale each column (Z-score) then transpose
scaled_counts[combined_df$SYMBOL,]
#View(scaled_counts)

scaled_counts
colnames(scaled_counts) <- colnames(expression)

  
heatmap <- HeatmapAnnotation(
  Condition = metadata$trim_term,  # Use the 'condition' column as annotation data
  col = list(Trimester = c("First Trimester" = "green", "Second Trimester" = "blue"))  # Define colors for annotation levels
  ) %>%
    Heatmap(
      scaled_counts[combined_df$SYMBOL,],
      cluster_rows = TRUE,
      column_labels = colnames(scaled_counts),
      name = "Z-score",
      cluster_columns = TRUE,
      top_annotation = .,
      show_column_names = FALSE,
      column_title = "First and second trimestr comparison"
    )
  
heatmap



heatmap <- HeatmapAnnotation(
  Trimestr = metadata$trim_term,  # Use the 'condition' column as annotation data
  Sex = metadata$Combined.Fetus.Sex,
  col = list(Trimester = c("First Trimester" = "green", "Second Trimester" = "blue"))  # Define colors for annotation levels
) %>%
  Heatmap(
    scaled_counts[combined_df$SYMBOL,],
    cluster_rows = TRUE,
    column_labels = colnames(scaled_counts),
    name = "Z-score",
    cluster_columns = TRUE,
    top_annotation = .,
    show_column_names = FALSE,
    column_title = "First and second trimestr comparison"
  )
heatmap




# Gestational.Age
heatmap <- HeatmapAnnotation(
  Condition = metadata$Gestational.Age.Appr,  # Use the 'condition' column as annotation data
  col = list(Trimester = c("First Trimester" = "green", "Second Trimester" = "blue"))  # Define colors for annotation levels
) %>%
  Heatmap(
    scaled_counts[combined_df$SYMBOL,],
    cluster_rows = TRUE,
    column_labels = colnames(scaled_counts),
    name = "Z-score",
    cluster_columns = TRUE,
    top_annotation = .,
    show_column_names = FALSE
  )
heatmap


# boys heatmap
metadata_boys <- metadata[metadata$Combined.Fetus.Sex == "Male",]
metadata_boys_samples <- c(metadata_boys$arraydatafile_exprscolumnnames)
metadata_boys_samples
expression
colnames(metadata_boys_samples)
expression_boys <- expression[colnames(expression) %in% metadata_boys_samples]
#expression_boys

#expression
#View(expression)
scaled_counts<- t(apply(expression_boys, 1, scale)) #center and scale each column (Z-score) then transpose
scaled_counts[combined_df$SYMBOL,]
#View(scaled_counts)

scaled_counts
colnames(scaled_counts) <- colnames(expression_boys)


heatmap <- HeatmapAnnotation(
  Condition = metadata_boys$trim_term,  # Use the 'condition' column as annotation data
  col = list(Trimester = c("First Trimester" = "green", "Second Trimester" = "blue"))  # Define colors for annotation levels
) %>%
  Heatmap(
    scaled_counts[combined_df$SYMBOL,],
    cluster_rows = TRUE,
    column_labels = colnames(scaled_counts),
    name = "Z-score",
    cluster_columns = TRUE,
    top_annotation = .,
    show_column_names = FALSE
  )

heatmap


metadata_girls <- metadata[metadata$Combined.Fetus.Sex == "Female",]
metadata_girls_samples <- c(metadata_girls$arraydatafile_exprscolumnnames)
metadata_girls_samples
expression
colnames(metadata_girls_samples)
expression_girls <- expression[colnames(expression) %in% metadata_girls_samples]
#expression_boys

#expression
#View(expression)
scaled_counts<- t(apply(expression_girls, 1, scale)) #center and scale each column (Z-score) then transpose
scaled_counts[combined_df$SYMBOL,]
#View(scaled_counts)

scaled_counts
colnames(scaled_counts) <- colnames(expression_girls)


heatmap <- HeatmapAnnotation(
  Condition = metadata_girls$trim_term,  # Use the 'condition' column as annotation data
  col = list(Trimester = c("First Trimester" = "green", "Second Trimester" = "blue"))  # Define colors for annotation levels
) %>%
  Heatmap(
    scaled_counts[combined_df$SYMBOL,],
    cluster_rows = TRUE,
    column_labels = colnames(scaled_counts),
    name = "Z-score",
    cluster_columns = TRUE,
    top_annotation = .,
    show_column_names = FALSE
  )

heatmap




###
#### Overrepresentation analysis ORA 
#
# Read in the data
df <- read_excel("difexp_1-2.xlsx", sheet = "Sheet1")
#df <- read.csv(paste0(in_path, 'expression2.csv'), row.names = 1)

# Remove the first column from the data frame
df <- df[-1]
#View(df)


original_gene_list <- df$logFC
names(original_gene_list) <- df$SYMBOL
original_gene_list
# omit any NA values 
gene_list <- na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

#gene_list
#df_2_4_7_8 <- df[,c(1,4,7,8)]
#View(df_2_4_7_8)

# Annotate according to differential expression
# I can use different paramaters, lile logFC> 0 or logFC >1
df <- df %>% mutate(diffexpressed = case_when(
  logFC > 1 & adj.P.Val < 0.05 ~ 'All',
  logFC < -1 & adj.P.Val < 0.05 ~ 'All',
  adj.P.Val > 0.05 ~ 'NO'
))
df


# Get the genes that are present in your dataframe
genes_in_data <- df$SYMBOL
genes_in_data
# Read in the .gmt file
file <- "BP_subset_GO_Hs.symbols.gmt"
pwl2 <- read.gmt(file) 
View(pwl2)
# Subset to the genes that are present in our dataset
pwl2 <- pwl2[pwl2$gene %in% genes_in_data,] 
View(pwl2)
# Save the filtered background gene set
filename <- 'GO_BP.RDS'
saveRDS(pwl2, filename)

##### SQUIDTIP! If you want to parse several .gmt files at once, you can use a loop:
#
#gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
#for (file in gmt_files){
#  file <- gmt_files[1]
#  pwl2 <- read.gmt(file) 
#  pwl2 <- pwl2[pwl2$gene %in% genes_in_data,]
#  filename <- paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RDS', sep = '')
#  saveRDS(pwl2, filename)
#}
#####
# Remove non-significant genes
#View(df_DifExrp)
df <- df[df$diffexpressed != 'NO', ]
# Substitute names so they are annotated nicely in the heatmap later
#df_DifExrp$diffexpressed <- gsub('DOWN', 'Healthy', gsub('UP', 'Severe', df_DifExrp$diffexpressed))
View(df)
unique(df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
deg_results_list <- split(df, df$diffexpressed)
#View(deg_results_list)

## Run ClusterProfiler -----------------------------------------------
# Settings
out_path <- "/Users/RTIntelektFBT/Desktop/Roman_Project/try_from_almostZERO"
name_of_comparison <- '2-1_trimestrs' # for our filename
background_genes <- 'Biological_processes_GO' # for our filename
bg_genes <- readRDS('GO_BP.RDS') # read in the background genes
#View(bg_genes)
padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways
filename <- paste0(out_path, 'clusterProfiler/', name_of_comparison, '_', background_genes) # filename of our PEA results
##### SQUIDTIP! An option to read in your background genes by only defining your 'background_genes' variable
#if(background_genes == 'KEGG'){
#  bg_genes <- readRDS(paste0(bg_path, 'kegg.RDS'))
#} else if(background_genes == 'reactome'){
#  bg_genes <- readRDS(paste0(bg_path, 'reactome.RDS'))
#} else if(background_genes == 'go.bp'){
#  bg_genes <- readRDS(paste0(bg_path, 'go.bp.RDS'))
#} else {
#  stop('Invalid background genes. Select one of the following: KEGG, Reactome, GO, or add new pwl to function')
#}
####
#View(deg_results_list)
#View(bg_genes)
unique(bg_genes$gene)
deg_results_list[1]



# Run clusterProfiler on each sub-dataframe
View(deg_results_list)
res <- lapply(names(deg_results_list),
              function(x) enricher(gene = deg_results_list[[x]]$SYMBOL,
                                   TERM2GENE = bg_genes))
View(res)
names(res) <- names(deg_results_list)
#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)

res_df <- do.call(rbind, res_df)
#head(res_df)
View(res_df)


res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', rownames(res_df)))

View(res_df)

# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
target_pws
length(res_df$ID)
length(target_pws)
res_df <- res_df[res_df$ID %in% target_pws, ]
View(res_df)

# Save clusterprofiler results
write.csv(res_df, "2-1_trim_resclusterp.csv", row.names = FALSE)

#df_res_df_sign <- res_df_sign[res_df_sign$Description == "GOBP_MITOCHONDRIAL_TRANSLATION",]
#View(df_res_df_sign)


### For visualisation

#install.packages('pheatmap')
#install.packages("DOSE")
#install.packages("enrichplot")
#install.packages("ggupset")
library(pheatmap)
library(DOSE)
library(enrichplot)
library(ggupset)

# Read in the data
res_df <- read.csv('2-1_trim_resclusterp.csv')
View(res_df)
#res_df_sign_i <- res_df_sign[,1:11]
#View(res_df_sign_i)

bg_genes <- readRDS('GO_BP.RDS')
View(bg_genes)
#bg_genes$term <- as.character(bg_genes$term)

# Convert clusterProfiler object to a new "enrichResult" object
# Select only upregulated genes in Severe
#View(res_df)
# I wanna try all genes, not only UP regulated
res_df <- res_df %>% filter(diffexpressed == 'All') %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 

#res_df <- res_df %>% dplyr::select(!c('minuslog10padj', 'diffexpressed')) 

# the next is already done
#rownames(res_df) <- res_df$ID
# For visualisation purposes, let's shorten the pathway names
res_df$Description <- gsub('(H|h)iv', 'HIV', 
                           gsub('pd 1', 'PD-1',
                                gsub('ecm', 'ECM', 
                                     gsub('(I|i)nterleukin', 'IL', 
                                          gsub('(R|r)na', 'RNA', 
                                               gsub('(D|d)na', 'DNA',
                                                    gsub(' i ', ' I ', 
                                                         gsub('(A|a)tp ', 'ATP ', 
                                                              gsub('(N|n)adh ', 'NADH ', 
                                                                   gsub('(N|n)ad ', 'NAD ',
                                                                        gsub('t cell', 'T cell',
                                                                             gsub('b cell', 'B cell',
                                                                                  gsub('built from .*', ' (...)',
                                                                                       gsub('mhc', 'MHC',
                                                                                            gsub('mhc class i', 'MHC I', 
                                                                                                 gsub('mhc class ii', 'MHC II', 
                                                                                                      stringr::str_to_sentence(
                                                                                                        gsub('_', ' ',  
                                                                                                             gsub('GOBP_|KEGG_|REACTOME_', '', res_df$Description)))))))))))))))))))



#View(genes_in_data)
#is.character(genes_in_data)
rownames(res_df) <- res_df$ID 
View(res_df)

enrichres <- new("enrichResult",
                 readable = FALSE,
                 result = res_df,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 organism = "human",
                 ontology = "BP",
                 gene = df$SYMBOL,
                 keytype = "SYMBOL",
                 universe = unique(bg_genes$gene),
                 gene2Symbol = character(0),
                 geneSets = bg_genes)

class(enrichres)

# Barplot
barplot_20 <- barplot(enrichres, showCategory = 20) 
barplot_20
mutate(enrichres, qscore = -log(p.adjust, base = 10)) %>% 
  barplot(x = "qscore")
barplot_path <- "/enrichment/barplot.png"
ggsave(barplot_path,barplot_20)

# Dotplot
dotplot(enrichres, showCategory = 15) + ggtitle("Up-Down_regulated")

# Cnetplot
cnetplot(enrichres)
cnetplot(enrichres, categorySize="pvalue",foldChange = gene_list, showCategory=8)

goplot(enrichres)

View(enrichres)

# Heatplot
heatplot(enrichres, showCategory = 5)

# Treeplot
enrichres2 <- pairwise_termsim(enrichres) # calculate pairwise similarities of the enriched terms using Jaccard’s similarity index
#treeplot(enrichres2)

# Enrichment map 
emapplot(enrichres2)

# Upsetplot
upsetplot(enrichres)






#
##
### GSEA gen set enrichment analysis
##
#BiocManager::install("fgsea")
library("fgsea")
# ----------------------
# GSEA tutorial
# ----------------------
# Setting up environment ===================================================
# Clean environment
#rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
#gc() # free up memory and report the memory usage
#options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
# Set seed
#set.seed(123456)
# Set project library
.libPaths('C:/Users/RTIntelektFBT/Desktop/Roman_Project/try_from_almostZERO')
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(fgsea)
# Set relevant paths
list.files()
in_path <- "Datasets/"
out_path <- "PEA/Results/"

getwd()

setwd('/Users/RTIntelektFBT/Desktop/Roman_Project/try_from_almostZERO')
getwd()

## Set relevant paths
#list.files()
#in_path <- "Datasets/"
#out_path <- "PEA/Results/"
#bg_path <- "PEA/Background_genes/"

# Functions ===================================================
## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: prepare_gmt --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}

# Analysis ====================================================

## 1. Read in data -----------------------------------------------------------
#df <- read.table("1_2_DEGs.tsv", header = TRUE, sep = "\t")

DEGs <- read_excel("difexp_1-2.xlsx", sheet = "Sheet1")
#df <- read.csv(paste0(in_path, 'expression2.csv'), row.names = 1)
# Remove the first column from the data frame
DEGs <- DEGs[-1]
#View(DEGs)

head(DEGs)
# Get all the genes in your dataset and assign them to my_genes 
# THE SAME AS IN THE TUTORIAL SQUID
df_1_4_7_8 <- DEGs[,c(1,4,7,8)]
#my_genes <- df_2_4_7_8$SYMBOL

# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# Copy the .gmt file to your folder, in my case, its 'PEA/Background_genes/'
# Then read in the .gmt file
gmt_file <- "BP_subset_GO_Hs.symbols.gmt"

## 2. Prepare background genes -----------------------------------------------

# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df_1_4_7_8$SYMBOL
bg_path <- "C:/Users/RTIntelektFBT/Desktop/Roman_Project/try_from_almostZERO"
gmt_files <- list.files(path = bg_path,pattern = '.gmt', full.names = TRUE)
gmt_files
bg_genes <- prepare_gmt(gmt_files, my_genes, savefile = FALSE)
#View(bg_genes)
head(df_1_4_7_8)
# Prepare your ranked list of genes
rankings <- sign(df_1_4_7_8$logFC)*(-log10(df_1_4_7_8$P.Value)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- df_1_4_7_8$SYMBOL # genes as names
#View(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
head(rankings)
plot(rankings)

###
#max(rankings)
#min(rankings)
# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
#max_ranking <- max(rankings[is.finite(rankings)])
#min_ranking <- min(rankings[is.finite(rankings)])
#rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
#rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
#rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
###

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

head(GSEAres)
View(GSEAres)
GSEAres_ordered <- GSEAres[order(padj), ]
View(GSEAres_ordered)

## 6. Check results ------------------------------------------------------
# Top 6 enriched pathways (ordered by p-val)
head(GSEAres[order(pval), ])

sum(GSEAres[, padj < 0.01])
sum(GSEAres[, pval < 0.01])


number_of_top_pathways_up = 10
number_of_top_pathways_down = 10
topPathwaysUp <- GSEAres[ES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(padj), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

#pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()

# Select only independent pathways, removing redundancies/similar pathways
collapsedPathways <- collapsePathways(GSEAres[order(padj)][padj < 0.05], bg_genes, rankings)
mainPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
#pdf(file = paste0('GSEA/Selected_pathways/', paste0(filename, background_genes, '_gsea_mainpathways.pdf')), width = 20, height = 15)
plotGseaTable(bg_genes[mainPathways], rankings, GSEAres, gseaParam = 0.5)
#dev.off()

#If you’d like to export the tables, just uncomment the 2 lines above. You can also export to .png, or other formats:
#png(file = paste0(filename, ‘_gsea_mainpathways.png’), width = 1500, height = 800)
#plotGseaTable(bg_genes[mainPathways], rankings, GSEAres, gseaParam = 0.5)





# plot the most significantly enriched pathway
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)


bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]]

#View(GSEAres_ordered)
#head(GSEAres_ordered,1)
#head(GSEAres[order(padj), ], 2)$pathway
#GSEAres_ordered[2,]

bg_genes[[GSEAres_ordered[2,]$pathway]]
View(bg_genes)
# plot the 2 most significantly enriched pathway
plotEnrichment(bg_genes[[GSEAres_ordered[2,]$pathway]],
               rankings) + 
  labs(title = GSEAres_ordered[2,]$pathway)

# plot the 2 most significantly enriched pathway

plotEnrichment(bg_genes[[GSEAres_ordered[2,]$pathway]],
               rankings) + 
  labs(title = GSEAres_ordered[2,]$pathway)

plotEnrichment(bg_genes[[GSEAres_ordered[142,]$pathway]],
               rankings) + 
  labs(title = GSEAres_ordered[142,]$pathway)


bg_genes[[GSEAres_ordered[1,]$pathway]]
GSEAres_ordered[1,]$pathway
bg_genes[1]

class(bg_genes)
bg_genes[[1]]

GSEAres_ordered$pathway[1:5]

GSEAres_ordered[1:5]$pathway

plotGseaTable(
  pathways = bg_genes[GSEAres_ordered$pathway[1:5]],
  stats = rankings,
  fgseaRes = GSEAres_ordered,
  gseaParam = 1,
  colwidths = c(5, 3, 0.8, 1.2, 1.2),
  pathwayLabelStyle = NULL,
  headerLabelStyle = NULL,
  valueStyle = NULL,
  axisLabelStyle = NULL,
  render = NULL
)


#install.packages("Seurat")
library(Seurat)
#ps <- plotCoregulationProfileReduction(bg_genes[topPathways], expression, reduction="tsne")
#ps


is.data.frame(GSEAres_ordered)
class(GSEAres_ordered)

GSEAres_ordered_df <- as.data.frame(GSEAres_ordered)
dotplot(GSEAres_ordered, showCategory=10) 


class(GSEAres_ordered)

emapplot(GSEAres_ordered_df, showCategory = 10)

ridgeplot(GSEAres_ordered_df) + labs(x = "enrichment distribution")


expression

gesecaRes <- geseca(pathways = bg_genes, E = expression, minSize = 15, maxSize = 500)
head(gesecaRes, 10)
View(gesecaRes)
plotGesecaTable(gesecaRes |> head(10), bg_genes, expression[,sub_pdata[sub_pdata$trim_term == "First Trimester",]$arraydatafile_exprscolumnnames])
plotGesecaTable(gesecaRes |> head(10), bg_genes, expression[,sub_pdata[sub_pdata$trim_term == "Second Trimester",]$arraydatafile_exprscolumnnames])
plotGesecaTable(gesecaRes |> head(10), bg_genes, expression)


topPathways[1]
bg_genes[topPathways[1]]
expression
colnames(expression)
class(bg_genes[topPathways[1]])

ordered_sup_pdata <- sub_pdata[,c("Gestational.Age.Appr","arraydatafile_exprscolumnnames")] 
ordered_sup_pdata <- ordered_sup_pdata[order(sub_pdata$Gestational.Age.Appr),]
ordered_sup_pdata[order(sub_pdata$Gestational.Age.Appr),]$arraydatafile_exprscolumnnames
ordered_sup_pdata
expression

expression_ordered <- expression[,ordered_sup_pdata[order(sub_pdata$Gestational.Age.Appr),]$arraydatafile_exprscolumnnames]
expression_ordered

plotGesecaTable(gesecaRes |> head(10), bg_genes, expression_ordered)


plotCoregulationProfile(pathway = bg_genes[[topPathways[1]]], 
                        E = expression_ordered,titles = colnames(expression_ordered),
                        conditions = ordered_sup_pdata$Gestational.Age.Appr)






plotCoregulationProfile(pathway = bg_genes[[topPathways[1]]], 
                        E = expression,titles = colnames(expression),
                        conditions = sub_pdata$Gestational.Age.Appr)



plotCoregulationProfile(pathway = bg_genes[[topPathways[1]]], 
                        E = expression[,sub_pdata[sub_pdata$trim_term == "First Trimester",]$arraydatafile_exprscolumnnames],
                        )

plotCoregulationProfile(pathway = bg_genes[[topPathways[1]]], 
                        E =  expression[,sub_pdata[sub_pdata$trim_term == "Second Trimester",]$arraydatafile_exprscolumnnames],
                        )

plotCoregulationProfile(pathway = bg_genes[["GOBP_MYELOID_LEUKOCYTE_ACTIVATION"]], 
                        E = expression[,sub_pdata[sub_pdata$trim_term == "First Trimester",]$arraydatafile_exprscolumnnames],
)

plotCoregulationProfile(pathway = bg_genes[["GOBP_MYELOID_LEUKOCYTE_ACTIVATION"]], 
                        E =  expression[,sub_pdata[sub_pdata$trim_term == "Second Trimester",]$arraydatafile_exprscolumnnames],
)


bg_genes[["GOBP_MYELOID_LEUKOCYTE_ACTIVATION"]]


bg_genes[["GOBP_CHROMOSOME_SEGREGATION"]]

plotCoregulationProfile(pathway = bg_genes[["GOBP_CHROMOSOME_SEGREGATION"]], 
                        E = expression[,sub_pdata[sub_pdata$trim_term == "First Trimester",]$arraydatafile_exprscolumnnames],
)

plotCoregulationProfile(pathway = bg_genes[["GOBP_CHROMOSOME_SEGREGATION"]], 
                        E =  expression[,sub_pdata[sub_pdata$trim_term == "Second Trimester",]$arraydatafile_exprscolumnnames],
)

GOBP_CHROMOSOME_SEGREGATION

#sub_pdata[sub_pdata$trim_term == "First Trimester",]$secondaryaccession
#sub_pdata$trim_term
#sub_pdata

colnames(expression)

is.factor(sub_pdata$trim_term)
sub_pdata$Gestational.Age.Appr









#de_genes <- read_excel("difexp_1_2.xlsx")
library(ggrepel)

de_genes <- df
View(de_genes)
de_genes$diffexpressed <- "Not Significant"
de_genes$diffexpressed[de_genes$logFC > 1 & de_genes$adj.P.Val < 0.05] <- "Up-regulated"
de_genes$diffexpressed[de_genes$logFC < -1 & de_genes$adj.P.Val < 0.05] <- "Down-regulated"
my_colors <- c("black","blue","red","green","purple")
names(my_colors) <- c("Not Significant","Down-regulated","Up-regulated","Down-reg_top20","Up-reg_top20")
de_genes$delabel<- NA

de_genes$pi_value <- abs(-log10(de_genes$adj.P.Val)*de_genes$logFC)
#View(de_genes)

filtered_genes <- de_genes %>%
  filter(diffexpressed %in% c("Up-regulated","Down-regulated")) %>%
  arrange(desc(pi_value))
View(de_genes)
threshold <- filtered_genes %>% slice(20) %>% pull(pi_value)
View(filtered_genes)
threshold
de_genes$diffexpressed[de_genes$pi_value >= threshold & de_genes$diffexpressed %in% c("Up-regulated")] <- "Up-reg_top20"
de_genes$diffexpressed[de_genes$pi_value >= threshold & de_genes$diffexpressed %in% c("Down-regulated")] <- "Down-reg_top20"
de_genes <- de_genes[order(de_genes$pi_value, decreasing = TRUE),]
de_genes

de_genes$delabel[de_genes$diffexpressed %in% c("Up-reg_top20","Down-reg_top20")] <- de_genes$SYMBOL

View(de_genes)

plot1 <- ggplot(data = de_genes, aes(x = logFC, y = -log10(adj.P.Val), col = diffexpressed, label = delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = my_colors) +
  geom_vline(xintercept = c(-0.072,0.072,1,-1),col="red",linetype = 2) +
  geom_hline(yintercept = -log10(0.05246),col="red",linetype = 2) +
  theme(text = element_text(size = 20)) 
plot1
ggsave("Nice_Voclano_Plot.png",plot=plot1,width=10, height=8,path="C:\Users\RTIntelektFBT\Desktop\Volcano_Plot")














# Assuming df is your data frame
df <- data.frame(
  Name = c("Alice", "Bob", "Charlie"),
  Age = c(25, 30, 22),
  Score = c(90, 85, 92)
)

# Extracting the "Name" column as a vector
name_vector <- df[["Name"]]
name_vector
name_vector <- df[,"Name"]
# Extracting the "Name" column using the $ operator
name_vector <- df$Name



###
plotEnrichment(bg_genes[['REACTOME_FORMATION_OF_FIBRIN_CLOT_CLOTTING_CASCADE']],
               rankings) + 
  labs(title = 'Reactome pathway: Formation of fibrin clot / clotting cascade') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  scale_y_continuous('Enrichment score (ES)') +
  geom_line(col = 'purple', linewidth = 2) 
###







