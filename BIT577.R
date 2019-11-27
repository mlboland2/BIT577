rm(list = ls())

#load required packages/install if necessary
library(phyloseq)
library(readxl)
library(ggplot2) 
library(DESeq2)

#import data; can do using import dataset tool in "Environment"; description of each table in project_workflow txt file
#change file location to match yours if not importing from "Environment"
otu_mat <- read_excel("G:/My Drive/BIT577_data/OTU_table.xlsx")
tax_mat <- read_excel("G:/My Drive/BIT577_data/taxonomy.xlsx")
samples_df <- read_excel("G:/My Drive/BIT577_data/metadata.xlsx")

#set row names
row.names(otu_mat) <- otu_mat$otu
#remove first column since phyloseq matrix has to be numeric
otu_mat <- otu_mat[,-1]

row.names(tax_mat) <- tax_mat$OTU
tax_mat <- tax_mat[,-1]

row.names(samples_df) <- samples_df$Sample

#convert to matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#set variables to factors
samples_df$TreatmentGroup <- as.factor(samples_df$TreatmentGroup)
samples_df$Depth <- as.factor(samples_df$Depth)
samples_df$Mortality <- as.factor(samples_df$Mortality)
samples_df$Ageweeks <- as.factor(samples_df$Ageweeks)
samples_df$Temperature <- as.factor(samples_df$Temperature)

#convert to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

#top10OTUS
TopNOTUs <- names(sort(taxa_sums(OTU), TRUE)[1:10])
top10 <- prune_taxa(TopNOTUs, OTU)

#transform to relative abundance
relative <- transform_sample_counts(OTU, function(x) x/sum(x))
#filter for variance greater than 0.000001
relative_filtered <- filter_taxa(relative, function(x) var(x) > 1e-05, TRUE)

#combine into phyloseq object
physeq <- phyloseq(relative_filtered, TAX, samples)

#abundnace plot by depth
abundance_by_depth_plot <- plot_bar(subset_taxa(physeq, class != ""), x="Depth", fill = "class")

#view plot
abundance_by_depth_plot

#abundance plot by sample
abundance_by_sample_plot <- plot_bar(subset_taxa(physeq, class !=""), x="Sample", fill = "class")
abundance_by_sample_plot

#Jaccard plot
dist <- phyloseq::distance(physeq, "jaccard")
imds <- ordinate(physeq, "MDS", distance=dist)
jaccard_plot <- plot_ordination(physeq, imds, shape="TreatmentGroup", color="Depth")
jaccard_plot <- jaccard_plot + theme_bw() + geom_point(size = 5)
jaccard_plot

#Bray
dist2 <- phyloseq::distance(physeq, "bray")
imds2 <- ordinate(physeq, "MDS", distance=dist2)
bray_plot <- plot_ordination(physeq, imds2, shape="TreatmentGroup", color="Depth")
bray_plot <- bray_plot + theme_bw() + geom_point(size = 5)
bray_plot

#revert back to absolute abundance (differential abundance will not work with relative abundances)
physeq2 <- phyloseq(OTU, TAX, samples)

##differential abundance
diffa <- phyloseq_to_deseq2(physeq2, ~ Depth)
diffa <- DESeq(diffa, test="Wald", fitType = "parametric")
res = results(diffa, cooksCutoff = FALSE)
alpha= 0.01
sigtab = res[which(res$padj < alpha),]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq2)[rownames(sigtab),], "matrix"))
head(sigtab)

#by compost
diffa <- phyloseq_to_deseq2(physeq2, ~ TreatmentGroup)
diffa <- DESeq(diffa, test="Wald", fitType = "parametric")
res = results(diffa, cooksCutoff = FALSE)
alpha= 0.01
sigtab = res[which(res$padj < alpha),]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq2)[rownames(sigtab),], "matrix"))
head(sigtab)



