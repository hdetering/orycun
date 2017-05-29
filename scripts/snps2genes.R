# connect to Ensembl's rabbit genes database
library("biomaRt")
library(LDsnpR)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
ensembl <- useDataset("ocuniculus_gene_ensembl", mart=ensembl)

# get all protein-coding genes
genes <- data.frame(
  getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","strand","external_gene_name","description"), 
        filters=c("biotype"), 
        values=list("protein_coding"), 
        mart=ensembl)
)

# get chromosome <-> ID mapping
# (order according to ENSEMBL's OryCun2.0 reference genome)
chroms <- read.delim("meta/oryCun2.chr.csv", header=F, col.names=c("chromosome_id","chromosome_name","len"), strip.white=T)
genes <- merge(genes, chroms[,c("chromosome_id","chromosome_name")])

# only genes on autosomes are relevant
genes <- subset(genes, !(chromosome_name %in% c("MT", "X", "Y")))
genes <- genes[order(genes$chromosome_name, genes$ensembl_gene_id), ]
write.table(genes[,c(8,1,2,3,4,5,6,7)], "orycun_genes_ensembl.csv", row.names=F, sep='\t', quote=F)
write.table(genes[,c(8,1,2,3,4,5)], "genes.csv", row.names=F, sep='\t', quote=F)

# export mapping for use with LDsnpR
write.table(genes[,c("ensembl_gene_id","chromosome_id","start_position","end_position")], 
            "genes.ldsnpr.csv", row.names=F, sep='\t', quote=F)

# filter and format SNP data
# issues:
#  - a few coordinates reported by STACKS have errors (bp positions of >4e+9) -> need to be filtered out
#  - SNP ids need to be numeric -> reformat: "<tag_id>_<offset>" to "<tag_id><offset>"
snps <- read.delim("meta/snps.ldsnpr.csv", stringsAsFactors=F)
snps.filt <- subset(snps, pos<4e+9)
snps.filt$snp <- sapply(snps.filt$snp, 
                              function(x) {
                                a <- as.numeric(unlist(strsplit(x,'_')))
                                100*a[1]+a[2]
                              })
write.table(snps.filt, "snps.ldsnpr.f.csv", sep='\t', row.names=F, quote=F)

###############################################################
# Perform the actual assignment of SNPs to genes using LDnspR #
###############################################################

# flank.genes.* used to be 40000, but we changed it to zero so only linked SNPs are considered
# !IMPORTANT! Adapt call to load.ld in LDsnpR:::compute.weights:22
# !IMPORTANT! --> add parameter: ', rsprefix = ""'

snp.ld.analysis(snpdata.url = "snps.ldsnpr.f.csv", 
                genome.url = "genes.ldsnpr.csv",
                include.ld.data = F,
                scoring.function = "min",
                flank.genes.left = 40000,
                flank.genes.right = 40000,
                generate.plink.set = F,
                outfile.path = "LDsnpRout.40k.min.csv")

snp.ld.analysis(snpdata.url = "snps.ldsnpr.f.csv", 
                genome.url = "genes.ldsnpr.csv",
                include.ld.data = F,
                scoring.function = "get.snps",
                flank.genes.left = 40000,
                flank.genes.right = 40000,
                generate.plink.set = F,
                outfile.path = "LDsnpRout.40k.snps.csv")


################################
#    output postprocessing     #
################################

# combine results in terminal:
#join -1 3 -2 4 <(tail -n+2 genes.csv | sort -k3) <(tail -n+6 LDsnpRout.40k.min.csv | cut -f3-7 | sort -k4) | tr ' ' '\t' > genes.scores.40k.csv

g.scores <- read.delim("genes.scores.40k.csv", header=F, stringsAsFactors=F,
                     col.names=c("ensembl_gene_id","chr_id","chr_name","gene_start","gene_end","strand",
                                 "window_start","window_end","window_size","score"))
# apply Benjamini-Hochberg correction for multiple testing to p values
g.scores$p.fdr <- p.adjust(g.scores$score, method='fdr')

# does p-value magnitude correlate with gene window size?
cor(g.scores$window_size, -log10(g.scores$score), use="complete")
cor(g.scores$window_size, -log10(g.scores$p.fdr), use="complete")

g.snps <- read.delim("LDsnpRout.40k.snps.csv",
                     col.names=c("row","space","start","end","width","names","score.data"), 
                     comment.char="#", stringsAsFactors=F)

# fix LDsnpR output (snp with id "400000" is written as "4e+05"...)
system("sed -i 's/4e+05/400000/g' LDsnpRout.40k.snps.csv")
# determine SNP providing minimum p-value for each gene
system(paste("./scripts/get_min_snps.py results/snp.scores.geno.csv", 
             "LDsnpRout.40k.min.csv", 
             "LDsnpRout.40k.snps.csv",
             "> LDsnpRout.40k.gs.csv"))
# load for each gene the SNP having the minimum p-value (output of "get_min_snps.py")
g.min_snps <- read.delim("LDsnpRout.40k.gs.csv", header=F, col.names=c("ensembl_gene_id", "min_snp"))
g.scores.snps <- merge(g.scores, g.min_snps)

# add SNP attributes
snps <- read.delim("results/snp.scores.geno.csv", header=F, stringsAsFactors=F,
                   col.names=c("snp_id","chr","pos","strand1","A","B","test","aff","unaff","p","empty"))
g.scores.snps <- merge(g.scores.snps, snps[,c("snp_id","pos","A","B","test","aff","unaff")], 
                       by.x="min_snp", by.y="snp_id", all.x=T, sort=F)
# restore original ordering of rows and columns
g.scores.snps <- g.scores.snps[order(g.scores.snps$ensembl_gene_id),]
g.scores.snps <- g.scores.snps[,c(2:12,1,13:18)]


# export results
write.table(g.scores.snps, "genes_scores_snps.40k.csv", sep='\t', row.names=F)
