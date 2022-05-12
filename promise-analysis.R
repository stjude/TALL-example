#####################################
# Obtain packages

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("PROMISE")

#install.packages("openxlsx")

library(PROMISE)
library(openxlsx)

#################################
# Specify local directory

ddir<-"C:/Users/xcao12/Box Sync/BookChapterApril2022/ForGit/"
source(paste0(ddir, "pc06-fdr.R"))

datfl<-paste0(ddir, "chr9-example-data.xlsx")
clin<-read.xlsx(datfl,  sheet = "clinical")
RNA <- read.xlsx(datfl,  sheet = "expression")
geneann <- read.xlsx(datfl,  sheet = "annotation")

############################
# Function for data merging

mtx.phe.mtch<-function(mtx,      # gene matrix 
                       phe,      # phenotype data frame
                       nam=NULL) # add feature names as first row of gene matrix
{
  mtx.mtch<-is.element(dimnames(mtx)[[2]], dimnames(phe)[[1]])
  phe.mtch<-is.element(dimnames(phe)[[1]], dimnames(mtx)[[2]])
  
  mtx<-mtx[, mtx.mtch]
  phe<-phe[phe.mtch,]
  
  mtx.ord<-order(dimnames(mtx)[[2]])
  phe.ord<-order(dimnames(phe)[[1]])
  
  mtx<-mtx[, mtx.ord]
  phe<-phe[phe.ord,]
  if (!is.null(nam))
  { mtx<-cbind.data.frame(probeid=dimnames(mtx)[[1]], mtx)}
  
  res<-list(mtx=mtx, phe=phe)
  return(res)
}

###############################
# Data preparation

RNA<-RNA[!is.na(RNA$ensembl.ID), ]
rownames(RNA)<-RNA$ensembl.ID
rownames(clin)<-clin$ID

clin$EFScensor[clin$First_Event%in%c("None", "Censored")]<-0
clin$EFScensor[clin$First_Event%in%c("Death","Progression", "Relapse" , "Second Malignant Neoplasm")]<-1

clin$OScensor[clin$Vital_Status%in%"Alive"]<-0
clin$OScensor[clin$Vital_Status%in%"Dead"]<-1

RNAdat<-mtx.phe.mtch(RNA, clin)
nprb<-nrow(RNA)

RNA_mtx<-as.matrix(RNAdat$mtx)
RNA_clin<-RNAdat$phe

#Generate the exprSet and Geneset
metadata<-data.frame(labelDescription=c(dimnames(RNA_clin)[[2]]),
                     row.names=c(dimnames(RNA_clin)[[2]]))
phenoData<-new("AnnotatedDataFrame", data=RNA_clin, varMetadat=metadata)                     
exprset<-new("ExpressionSet", exprs=RNA_mtx, phenoData = phenoData)

prdat<-cbind.data.frame(stat.coef=c(1,1,1), 
                        stat.func=c("spearman.rstat", "jung.rstat", "jung.rstat"), 
                        endpt.vars=c("MRD29", "Event_Days,EFScensor", "OS_Days,OScensor"))
rownames(prdat)<-c("MRD29", "EFS", "OS")

################################
# Performs the analysis

test1 <- PROMISE(exprSet=exprset,
                 geneSet=NULL, 
                 promise.pattern=prdat, 
                 strat.var=NULL,
                 proj0=FALSE,
                 nbperm=TRUE, 
                 max.ntail=10,
                 seed=13, 
                 nperms=1000)

generes <- test1$generes
perm.p.index <- grep("perm.p", names(generes))

qq <- apply(generes[, perm.p.index], 2, function(prb){return(pc06.fdr(prb))})
colnames(qq) <- gsub("perm.p", "perm.q", colnames(qq))
generes2<-cbind.data.frame(generes, qq)
#generes2<-generesf

generesf<-merge(geneann, generes2, by.x=1, by.y=1)
rownames(generesf)<-NULL
generesf<-generesf[order(generesf$"PROMISE.perm.q"), ]

write.xlsx(generesf, paste0(ddir, "PROMISE-results.xlsx"), rowNames = FALSE)
