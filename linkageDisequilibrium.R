# Packages
# Install packages that are necessary:
install.packages("MatrixEQTL")
install.packages("plyr")
install.packages("rlist") 
install.packages("genetics")

# Load them to your Rsession
library(MatrixEQTL)
library(plyr)
library(rlist)
library(genetics)

# MatrixEQTL with ANOVA model causes vector memory problem in Rstudio
# Follows these steps:
# 1) Open terminal
# 2) cd ~
# 3) touch .Renviron
# 4) open .Renviron
# 5) Write R_MAX_VSIZE=700Gb in the first line of .txt and save
# 6) Restart R and start eQTL analysis

# Load gene_expression, genotype, probe_loc and snp_loc data 
# gene_names: (shape: R x 1)
gene_names = read.csv("AddNeuroMed or HBTRC/eqtl_genes.csv", 
                      header = TRUE, 
                      sep = ",", 
                      row.names = 1)
gene_names <- as.matrix(gene_names)

# subject_names: (shape: S x 1)
subject_names = read.csv("AddNeuroMed or HBTRC/eqtl_subjects.csv", 
                         header = TRUE, 
                         sep = ",", 
                         row.names = 1)
subject_names <- as.matrix(subject_names)

# snp_names: (shape: N x 1)
snp_names = read.csv("AddNeuroMed or HBTRC/eqtl_snps.csv", 
                     header = TRUE, 
                     sep = ",", 
                     row.names = 1)
snp_names <- as.matrix(snp_names)

# genes_value: Matrix with gene expression levels (shape: R x S)
genes_value = read.csv("AddNeuroMed or HBTRC/eqtl_expression.csv", 
                       header = TRUE, 
                       sep = ",", 
                       quote = "")

# snps_value: Matrix with genotypes (shape: N x S)
snps_value = read.csv("AddNeuroMed or HBTRC/eqtl_genotype.csv", 
                      header = TRUE, 
                      sep = ",", 
                      quote = "")

# probePos: Table with gene_id, chromosome number, start location, end location (shape: R x 4)
probePos = read.csv("AddNeuroMed or HBTRC/eqtl_probe_loc.csv", 
                    header = TRUE, 
                    sep = ",", 
                    quote = "")

# snpPos: Table with snp_names, chromosome number, and position on that chromosome (shape: N x 3)
snpPos = read.csv("AddNeuroMed or HBTRC/eqtl_snp_loc.csv", 
                  header = TRUE, 
                  sep = ",", 
                  quote = "")

# snpToGene: Table with snps and corresponding genes - 1st col: snp_names, 2nd col: gene_ids (shape: N x 2)
snpToGene = read.csv("AddNeuroMed or HBTRC/eqtl_snp_gene.csv", 
                     header = TRUE, 
                     sep = ",", 
                     quote = "")

# Prepare data; change row_names, and col_names, delete column
rownames(genes_value) <- genes_value[['id']]
genes_value <- subset(genes_value, select = -c(1))
colnames(genes_value) <- subject_names

rownames(snps_value) <- snp_names
snps_value <- subset(snps_value, select = -c(1))
colnames(snps_value) <- subject_names

rownames(snpPos) <- snp_names
snpToGene <- subset(snpToGene, select = -c(1))

# Create Sliced matrix for eQTL analysis
genes <- SlicedData$new()
genes$CreateFromMatrix(as.matrix(genes_value))

snps <- SlicedData$new()
snps$CreateFromMatrix(as.matrix(snps_value))

all(colnames(snps) == colnames(genes)) # check

# eQTL analysis
eQTL <- Matrix_eQTL_main(snps, genes,
                         useModel = modelANOVA,
                         output_file_name=NULL,
                         output_file_name.cis=NULL,
                         pvOutputThreshold.cis=1e-3, 
                         snpspos=as.data.frame(snpPos[,1:3]),
                         genepos=as.data.frame(probePos[1:4]))

# Minor allele frequency (MAF)
mafs = apply(as.matrix(snps_value),1,mean)/2
mafs = data.frame(snps=names(mafs), maf = mafs)

# Process cis-eQTL analysis
cis_eqtl_res = eQTL$cis$eqtls
cis_eqtls = merge(cis_eqtl_res, mafs, by="snps") # Add mafs information to the cis-analysis

cis_snps = cis_eqtl_res[['snps']] # Retrive resulting snps from cis-analysis
cis_snps <- as.matrix(cis_snps)

## Retrieve corresponding genes
cis_genes <- list()
for (snp in cis_snps){
  index = which(snp == snpToGene[[1]])
  cis_genes <- c(cis_genes, snpToGene[[2]][index])
}
cis_genes <- matrix(unlist(cis_genes))
cis_snps_gene <- cbind(cis_snps, cis_genes)
colnames(cis_snps_gene) <- c('snps', 'genes')

## Merge with corresponding genes 
cis_eqtls_all <- cbind(cis_eqtls, cis_snps_gene[, "genes"])
names(cis_eqtls_all)[names(cis_eqtls_all) == "cis_snps_gene[, \"genes\"]"] <- "genes_of_snps"
cis_eqtls_all <- cis_eqtls_all[c(1, 7, 2, 3, 4, 5, 6)]


# LINKAGE DISEQUILIBRUM
# Load genotype, snps name, alleles of snps
# snp_names: (shape: N x 1)
snp_names = read.csv("/Users/gulkamisli/Google Drive /Oxford/MSc Dissertation 2019/Tensorflow Model/eQTL Analysis/eqtl_snps.csv", 
                     header = TRUE, 
                     sep = ",", 
                     row.names = 1)
snp_names <- as.matrix(snp_names)

# alleles: (shape: N x 2)
alleles = read.csv("/Users/gulkamisli/Google Drive /Oxford/MSc Dissertation 2019/Tensorflow Model/eQTL Analysis/ANM/eqtl_snp_alleles.csv", 
                   header = TRUE, 
                   sep = ",")
alleles <- as.matrix(alleles)
alleles <- subset(alleles, select = -c(1))
rownames(alleles) <- snp_names

# snps_value: Matrix with genotypes (shape: N x S)
genotypes = read.csv("/Users/gulkamisli/Google Drive /Oxford/MSc Dissertation 2019/Tensorflow Model/eQTL Analysis/ANM/eqtl_ANM_genotype.csv", 
                     header = TRUE, 
                     sep = ",", 
                     quote = "")

rownames(genotypes) <- snp_names
genotypes <- subset(genotypes, select = -c(1))
colnames(genotypes) <- subject_names

# Store genotypes as list with each row of genotypes matrix as numeric vector
genotypes_list <- lapply(seq_len(nrow(genotypes)), function(i) as.numeric(genotypes[i,]) )

# Create another list to store genotypes as alleles (i.e. "A/T", "A/A", )
link_genotypes <- vector("list", length(snp_names))

# Function to convert genotypes (0, 1, 2) to allele pairs ("A/T", "A/A", "T/T")
link_genotypes <- lapply(seq_len(nrow(genotypes)), function(i) 
  as.genotype.allele.count(genotypes_list[[i]], alleles = c(alleles[i,"a0"], alleles[i,"a1"])))

# Change the names to snp_names to keep track them LD
names(link_genotypes) <- snp_names

asGenotypes <- sapply(seq_len(nrow(genotypes)), function(i) link_genotypes[[i]])
colnames(asGenotypes) <- snp_names

# Find eQTL snps before generating genotypes and having LD analysis
cisSnpIndex <- sapply(seq_len(length(cis_snps)), function(i) match(cis_snps[i], snp_names))
selectedAsGenotypes <- asGenotypes[,cisSnpIndex]

# Make genotypes from cis_snps
selectedCisMakeGenotypes <- makeGenotypes(selectedAsGenotypes)

# Linkage disequilibrium analysis for selected variants
LD_cis_ANM <- LD(selectedCisMakeGenotypes)

# Create a variable to keep track of LD snps
ld_snps <- colnames(LD_cis_ANM$`R^2`)
ld_snps <- as.matrix(ld_snps)

# Loading any of the files is as follows:
load(file = "AddNeuroMed/eQTL analysis/linkageDisequilibrum_ANM_2.Rdata") # for AddNeuroMed
load(file = "HBTRC/eQTL analysis/linkageDisequilibrum_ANM_2.Rdata") # for HBTRC
