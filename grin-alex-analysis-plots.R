
##############################################
# Some initial set-up

rm(list=ls());              # start with clean slate in R
options(stringsAsFactors=F) # turn off the most annonying default in R

library(openxlsx)
#install.packages("circlize")
library(circlize)
library(ComplexHeatmap)
####################################
# Load clinical, lesion and expression data

ddir<-"C:/Users/aelsayed/Box/BookChapterApril2022/ForGit/"

datfl<-paste0(ddir, "example-data-rounded.xlsx")
clin<-read.xlsx(datfl,  sheet = "clinical")
lsn<-read.xlsx(datfl,  sheet = "lesions")
RNA <- read.xlsx(datfl,  sheet = "expression")
geneann <- read.xlsx(datfl,  sheet = "annotation")

################################################
# Load GRIN-ALEX library:
source(paste0(ddir, "GRIN-ALEX-library.R"))

#########################################################
# Run GRIN using hg19 gene annotation data
##################################################

# To retrieve chromosome size data for GRCh37 (hg19) genome build from chr.info txt file available on UCSC genome browser
chr.size=get.chrom.length("Human_GRCh37")

GRIN.results=grin.stats(lsn, geneann, chr.size)

write.grin.xlsx(GRIN.results, "TALL2017-GRIN-result.xlsx")

######################################
# GRIN Plots:
####################################

# 1) To generate a whole-genome lesion plot
#######################################################
genomewide.plot=genomewide.lsn.plot(GRIN.results) # this function use the list of GRIN.results

#################################################
# 2) To generate an oncoprint of top genes:
##########################################

grin.results=as.data.frame(GRIN.results$gene.hits)
gene.lsn.data=GRIN.results$gene.lsn.data
ensembl.annotation=cbind.data.frame(geneann$gene, geneann$gene.name)
colnames(ensembl.annotation)=c("gene", "gene.name")

# prepare oncoprint for genes with q2 n.subjects<0.01 (Second Order Constellation GRIN Test)
oncoprint.genes = grin.results[grin.results$q2.nsubj < 0.01, ] # specify cutoff for selected genes 
oncoprint.genes=oncoprint.genes[2] # extract gene name of selected genes
oncoprint.genes=unlist(oncoprint.genes)
oncoprint.genes=as.vector(oncoprint.genes)

# Prepare matrix for selected genes with each row as a gene and each column is a patient
oncoprint.mtx=grin.oncoprint.mtx(oncoprint.genes, gene.lsn.data, ensembl.annotation)

##############################
# Pass the matrix of selected genes to oncoprint function:
#######################################################

# assign colour for each mutation type
col = c("gain" = "red", "loss" = "blue", "mutation" = "olivedrab", "fusion" = "gold")

# specify for each mutation type if the bar should take the whole or part of the cell, this is important if the same patient has gain and mutation for example  
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.00005, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big red
  gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["gain"], col = NA))
  },
  # big blue
  loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["loss"], col = NA))
  },
  # small green
  mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "pt"), h*0.4, 
              gp = gpar(fill = col["mutation"], col = NA))
  },
  # small gold
  fusion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "pt"), h*0.4, 
              gp = gpar(fill = col["fusion"], col = NA))
  }
)


column_title = "GRIN OncoPrint q2.nsubj < 0.01" # optional

heatmap_legend_param = list(title = "Lesion Category", at = c("gain", "loss", "mutation", "fusion"), 
                            labels = c("Gain", "Loss", "Mutation", "Fusion"))

# use oncoprint function from complexheatmap library to plot the oncoprint
oncoPrint(oncoprint.mtx,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)

#################################################
# N.B) instead of using genes with q2<0.01, we can pass a list of known driver genes and generate an oncoprint using same steps listed above
################################

oncoprint.genes=as.vector(c("ENSG00000148400", "ENSG00000171862",
                            "ENSG00000156531", "ENSG00000162367", "ENSG00000139687",
                            "ENSG00000138795", "ENSG00000184937", "ENSG00000127152", "ENSG00000135363",
                            "ENSG00000118513", "ENSG00000102974", "ENSG00000139083", "ENSG00000097007",
                            "ENSG00000164438", "ENSG00000159216", "ENSG00000133703", "ENSG00000006468",
                            "ENSG00000171843", "ENSG00000157554", "ENSG00000123473", "ENSG00000236311",
                            "ENSG00000096968", "ENSG00000105639"))

oncoprint.genes=unlist(oncoprint.genes)
oncoprint.genes=as.vector(oncoprint.genes)

oncoprint.mtx=grin.oncoprint.mtx(oncoprint.genes, gene.lsn.data, ensembl.annotation)

##############################
# 3) Generate regional plots for top genes (one gene per page)-plot all lesion types based on gene coordinates
##########################################

# this function will plot top affected genes by each lesion type specifying q.max=0.10
pdf("GRIN-Plots-top50genes.pdf",width = 8,height = 6, onefile = TRUE) # open a PDF to add lesion data plots, one gene per page
top.grin.gene.plots(GRIN.results)
dev.off()


##############################
# 4) Generate plots for one gene
##########################################
one.gene.plot=grin.gene.plot(GRIN.results, "LYL1")


###################################################################
####################################################################
# To run The associate lesions with expression (ALEX) analysis
##############################################################################
##############################################################################

gene.lsn=prep.gene.lsn.data(lsn, geneann)    # Prepare gene and lesion data for later computations
gene.lsn.overlap= find.gene.lsn.overlaps(gene.lsn)  # Use the results of prep.gene.lsn.data to find lesion-gene overlaps

lsn.type.matrix=prep.lsn.type.matrix(gene.lsn.overlap) # prepare lesion type matrix (JUST ONE ROW FOR EACH GENE)

lsn.grp.mtx=as.data.frame(lsn.type.matrix)
id.hits = colnames(lsn.grp.mtx)

# prepare expression data for ALEX analysis
expr.mtx=RNA
rownames(expr.mtx)=expr.mtx[,1] ## Abdel; take the first column with gene names outside the expression matrix
expr.mtx=expr.mtx[,-c(1,2)]

id.expr=colnames(expr.mtx)
lsn.grp.mtx=lsn.grp.mtx[ ,(names(lsn.grp.mtx) %in% id.expr)]

id.lsn.final=colnames(lsn.grp.mtx)
expr.mtx=expr.mtx[ ,(names(expr.mtx) %in% id.lsn.final)]

# To keep only genes with expression value > 0 in at least 5 patients
zero.expr=rowSums(expr.mtx==0)
expr.grp.keep=(ncol(expr.mtx)-zero.expr)>=5 
expr.mtx=expr.mtx[expr.grp.keep,]

# keep only gene that has sum of expression >1
expr.sum=rowSums(expr.mtx)
head(expr.sum)
expr.keep=(expr.sum)>1 
expr.mtx=expr.mtx[expr.keep,]  

############################
# To exclude any gene that has no lesions in at least 5 patients

n.none=rowSums(lsn.grp.mtx== "none")
lsn.grp.keep=(ncol(lsn.grp.mtx)-n.none)>=5
lsn.grp.mtx=lsn.grp.mtx[lsn.grp.keep,] 

# To keep only shared set of genes between expression and lesion matrices
lsn.genes=rownames(lsn.grp.mtx)
expr.mtx=expr.mtx[rownames(expr.mtx) %in% lsn.genes,] 

expr.genes=rownames(expr.mtx)
lsn.grp.mtx=lsn.grp.mtx[rownames(lsn.grp.mtx) %in% expr.genes,]

## order lsn matrix by gene ID first then colnames for patient IDs alphabitically
lsn.grp.mtx=lsn.grp.mtx[order(rownames(lsn.grp.mtx)),]
lsn.grp.mtx=lsn.grp.mtx[,order(colnames(lsn.grp.mtx))]

## order expr matrix by gene ID first then colnames for patient IDs alphabitically
expr.mtx=expr.mtx[order(rownames(expr.mtx)),]
expr.mtx=expr.mtx[,order(colnames(expr.mtx))]

# To make sure that both lesion and expression matrices have same patients and genes order
all(rownames(lsn.grp.mtx)==rownames(expr.mtx))
all(colnames(lsn.grp.mtx)==colnames(expr.mtx))


# To generate row.mtch file
lsn.ensembl.id=rownames(lsn.grp.mtx)
expr.ensembl.id=rownames(expr.mtx)

row.mtch=cbind.data.frame(expr.row=expr.ensembl.id,
                          hit.row=lsn.ensembl.id)


###################################
## To run Kruskal-wallis test to evaluate gene expression by lesion groups:

expr.matrix= as.matrix(expr.mtx)
lsn.matrix=as.matrix(lsn.grp.mtx)

KW.res=KW.hit.express(expr.matrix,
                      lsn.matrix,
                      row.mtch)

##################
# To prepare and print KW test results

colnames(KW.res)=c("expr.row", "gene", "p.KW")
kw.res.annotated=merge(geneann,KW.res,by="gene", all.y=TRUE)

# Compute FDR with the Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat)
kw.res.annotated$FDR.q=pc06.fdr(kw.res.annotated$p.KW)


################################
# A) To count number of patients with and without lesions

N.pts.without.lsn=rowSums(lsn.matrix == "none")
N.pts.gain=rowSums(lsn.matrix == "gain")  
N.pts.loss=rowSums(lsn.matrix == "loss")   
N.pts.mut=rowSums(lsn.matrix == "mutation")   
N.pts.fus=rowSums(lsn.matrix == "fusion")   
N.pts.multiple.lsns=rowSums(lsn.matrix == "multiple") 

# B) to get mean expression level for each gene in patients with and without lesions

mean.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, mean)

mean.pts.fus=mean.by.lsn[,1]
mean.pts.gain=mean.by.lsn[,2]
mean.pts.loss=mean.by.lsn[,3]
mean.pts.multiple=mean.by.lsn[,4]
mean.pts.mut=mean.by.lsn[,5]
mean.pts.none=mean.by.lsn[,6]

# C) to get median expression level for each gene in patients with and without lesions

median.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, median)

median.pts.fus=median.by.lsn[,1]
median.pts.gain=median.by.lsn[,2]
median.pts.loss=median.by.lsn[,3]
median.pts.multiple=median.by.lsn[,4]
median.pts.mut=median.by.lsn[,5]
median.pts.none=median.by.lsn[,6]

# D) to get stdev of expression level for each gene in patients with and without lesions 

sd.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, sd)

sd.pts.fus=sd.by.lsn[,1]
sd.pts.gain=sd.by.lsn[,2]
sd.pts.loss=sd.by.lsn[,3]
sd.pts.multiple=sd.by.lsn[,4]
sd.pts.mut=sd.by.lsn[,5]
sd.pts.none=sd.by.lsn[,6]

# E) to write all the results in one excel sheet

kw.res.annotated.final = cbind(kw.res.annotated, N.pts.gain, N.pts.without.lsn,
                               N.pts.loss, N.pts.fus, N.pts.mut, N.pts.multiple.lsns,
                               mean.pts.gain, mean.pts.none, mean.pts.loss, mean.pts.fus,
                               mean.pts.mut, mean.pts.multiple, median.pts.gain, median.pts.none,
                               median.pts.loss, median.pts.fus,
                               median.pts.mut, median.pts.multiple, sd.pts.gain, sd.pts.none,
                               sd.pts.loss, sd.pts.fus, sd.pts.mut, sd.pts.multiple)

write.csv(kw.res.annotated.final, "KW-results-TALL2017-lsn-Expr-association-median.csv")

##############################################################
# To extract genes with q<0.05 for boxplots:

selected.genes=kw.res.annotated[kw.res.annotated$FDR.q<0.05,]
selected.IDS=selected.genes$gene
selected.lsns=as.data.frame(lsn.grp.mtx[rownames(lsn.grp.mtx) %in% selected.IDS,])
selected.lsns=t(selected.lsns)
selected.expr=as.data.frame(expr.mtx[rownames(expr.mtx) %in% selected.IDS,])
selected.expr=t(selected.expr)

all(rownames(selected.lsns)==rownames(selected.expr))
all(colnames(selected.lsns)==colnames(selected.expr))

gene.annotation <- gene.data %>% mutate_all(na_if,"")
gene.annotation$gene.name <- ifelse(is.na(gene.annotation$gene.name),gene.annotation$gene, gene.annotation$gene.name) # add ensembl.ID to gene.name column if empty

pdf("TALL2017-KW-q0.05-boxplots-Final.pdf",width = 8,height = 5, onefile = TRUE) # open a PDF to add boxplots, one gene per page

n=ncol(selected.lsns)

library(ggplot2)
library(forcats)
library(EnvStats)

for (i in 1:n)
{
  ensembl.id=colnames(selected.lsns)[i]
  ensembl.df=as.data.frame(ensembl.id)
  colnames(ensembl.df) <- ("gene")
  genes.df.merged=merge(gene.annotation,ensembl.df,by="gene", all.y=TRUE)
  gene.names=as.character(genes.df.merged$gene.name)
  lsn.box=unlist(selected.lsns[,i])
  expr.box=unlist(selected.expr[,i])
  df=as.data.frame(cbind(lsn.box, expr.box))
  df$expr.box=as.numeric(df$expr.box)
  {
    p=ggplot(df, aes(x = fct_reorder(lsn.box, expr.box, .desc =TRUE), y = expr.box)) +
      geom_boxplot(aes(fill = fct_reorder(lsn.box, expr.box,  .desc =TRUE))) + 
      geom_jitter(position=position_jitter(0.02)) + # Increasing the number from 0.02 increase the space between the dots but add some extra dots if number of patients is small
      xlab(paste0(gene.names, '_lsn'))+
      ylab(paste0(gene.names, '_Expression'))+ 
      theme_bw(base_size = 12) +
      scale_fill_discrete(guide = guide_legend(title = "Lesion")) + 
      theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"))+
      theme(legend.text=element_text(size=8))+
      #theme(legend.position="none")+
      stat_n_text()
    print(p)
  }
}

dev.off()

