##############################################################################################
##############################################################################################
############################## Bioinformatics in Biomedicine #################################
################################## Pratical Lecture 9 #######################################
############ Identification of transcriptome alterations in disease ##########################
##############################################################################################
##############################################################################################


##############################################################################################
# Load Transcriptome data (RNA-seq) inside the RData object
##############################################################################################

table <- read.delim("C:/Users/beatr/Desktop/3º ano - 2022,2023/Projeto/Dados/pituitary.readcounts.tab", row.names = 1)
table <- table[ , colnames(table) != "Transcript.ID"]
table <- table[ , colnames(table) != "Length"]
dados <- as.matrix(table)
class(dados) <- "numeric"


##############################################################################################
# Create R object for limma function 
##############################################################################################

library(limma); library(edgeR)
data <- DGEList(dados)

##############################################################################################
# Identify the differentially expressed genes between high and low samples
##############################################################################################

#Create design matrix with coefficients to be used in linear models
foxo4_levels <- read.table("C:/Users/beatr/Desktop/3º ano - 2022,2023/Projeto/Dados/foxo4agepituitary.tab", sep = "\t", row.names = 1)
#foxo4_levels <- foxo1_levels["ENSG00000184481.17",-c(1,2)]
sampleType <- as.vector(ifelse(as.numeric(foxo4_levels >= 60), "Old", "Young"))
design <- model.matrix(~ 0 + sampleType); 
colnames(design) <- gsub("sampleType", "", colnames(design))
rownames(design) <- rownames(foxo4_levels)
data <- data[, rownames(foxo4_levels)]

# Filter low expressed genes
keep <- filterByExpr(data, design)
data <- data[keep,]


#Normalize data
data_norm <- voom(data, design)


## Fit a linear model per gene
data_model <- lmFit(data_norm, design)

# Defines the comparisons of interest
comparisons <- makeContrasts(Young-Old, levels=design)
data_comparisons <- contrasts.fit(data_model, contrasts = comparisons)


# Calculate levels of significance. Uses an empirical Bayes algorithm
# to shrink the gene-specific variances towards the average variance
# across all genes (Moderated T-test)

data_final <- eBayes(data_comparisons)

#Get differentially expressed genes
allgenes <- topTable(data_final, n="inf", adjust.method="fdr")
DEG <- topTable(data_final, n="inf", adjust.method="fdr",p.value = 0.05)


#Get differentially expressed genes (apply logFC cut-off)
#DEG <- topTable(data_final, n="inf", adjust.method="fdr",p.value = 0.05, lfc = 2) 
table(sign(DEG[,1]))


#Write DEG in a table (to open after in excel)
write.table(DEG, "DEG_foxo1_testis.csv", sep=";", quote=F, col.names=NA)
write.table(allgenes, "allgenes_foxo4_pituitary.csv", sep=";", quote=F, col.names=NA)

