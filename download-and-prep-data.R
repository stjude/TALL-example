##############################################
# Some initial set-up

rm(list=ls());              # start with clean slate in R
options(stringsAsFactors=F) # turn off the most annonying default in R


##############################
# Specify Local Directory for Preparing Data

local.dir="D:/Genomic_Projects2_Abdelrahman/Book-chapter/"

########################################
# Obtain needed packages

#install.packages("readxl")
#install.packages("writexl")

library(readxl)
library(writexl)

########################################################
# Download supplementary data for Liu et al (2017)

ng.supp.link="https://static-content.springer.com/esm/art%3A10.1038%2Fng.3909/MediaObjects/41588_2017_BFng3909_MOESM2_ESM.xlsx"
ng.supp.file=paste0(local.dir,basename(ng.supp.link))

download.file(ng.supp.link,
              ng.supp.file,
              mode="wb")

########################################
# Download the detailed clinical data from TARGET

target.clin.link="https://target-data.nci.nih.gov/Public/ALL/clinical/Phase2/harmonized/TARGET_ALL_ClinicalData_Phase_II_Validation_20211118.xlsx"
target.clin.file=paste0(local.dir,basename(target.clin.link))

download.file(target.clin.link,
              target.clin.file,
              mode="wb")

############################################
# Prepare the clinical data

# Read the TARGET clinical data
target.clin.data=read_xlsx(target.clin.file,sheet=1) 
target.clin.data=as.data.frame(target.clin.data)     

# Read the clinical data in the supplementary materials
supp.clin.data=read_xlsx(ng.supp.file,sheet=1)       
supp.clin.data=as.data.frame(supp.clin.data)         

# Clean up the USI patient ID in the TARGET data
target.clin.data$USI=gsub("TARGET-10-","",
                          target.clin.data$`TARGET USI`,
                          fixed=T)

# Merge the two clinical data sets
comb.clin.data=merge(target.clin.data,
                     supp.clin.data,
                     by="USI")

# clean up the column names of the combined data set
colnames(comb.clin.data)=gsub(".x","_TARGET",
                              colnames(comb.clin.data),fixed=T)

colnames(comb.clin.data)=gsub(".y","_supplement",
                              colnames(comb.clin.data),fixed=T)

colnames(comb.clin.data)=gsub("\r","",
                              colnames(comb.clin.data),fixed=T)

colnames(comb.clin.data)=gsub("\n","",
                              colnames(comb.clin.data),fixed=T)

colnames(comb.clin.data)=gsub(" ","_",
                              colnames(comb.clin.data),fixed=T)

###########################################
# Double check the data consistency across the two sources

table(comb.clin.data$Gender_TARGET,
      comb.clin.data$Gender_supplement)

table(comb.clin.data$Race_TARGET,
      comb.clin.data$Race_supplement)



##################################################
# Get names of spreadsheets in supplementary data file

supp.sheets=excel_sheets(ng.supp.file)
supp.sheets


######################################################
# Read the RNAseq expression data from the supplementary data file

RNAseq.data=read_xlsx(ng.supp.file,
                      sheet="Table S5 RNAseq FPKM")
RNAseq.data=as.data.frame(RNAseq.data)
rownames(RNAseq.data)=RNAseq.data[,1]
RNAseq.data=RNAseq.data[,-1]
RNAseq.data=as.matrix(RNAseq.data)
RNAseq.data=log2(RNAseq.data+1)

##################################################
# Read the sequence mutation data

SeqMut.data=read_xlsx(ng.supp.file,
                      sheet="Table S8 Sequence mutations")
SeqMut.data=as.data.frame(SeqMut.data)

####################################################
# Read the fusion data

fusion.data=read_xlsx(ng.supp.file,
                      sheet="Table S11 fusions")
fusion.data=as.data.frame(fusion.data)

###################################################
# Read the Copy Number Abnormality data from the supplementary data file

CNA.data=read_xlsx(ng.supp.file,
                   sheet="Table S13 CNA")
CNA.data=as.data.frame(CNA.data)

#################################################
# Prepare combined lesion data for GRIN analysis

SeqMut.lsns=cbind.data.frame(ID=SeqMut.data$sample,
                             chrom=gsub("chr","",SeqMut.data$chromosome), 
                             loc.start=SeqMut.data$start,
                             loc.end=SeqMut.data$start,
                             lsn.type="mutation")

fusion.lsnsA=cbind.data.frame(ID=fusion.data$sample,
                              chrom=gsub("chr","",fusion.data$chr_a),
                              loc.start=fusion.data$position_a,
                              loc.end=fusion.data$position_a,
                              lsn.type="fusion")

fusion.lsnsB=cbind.data.frame(ID=fusion.data$sample,
                              chrom=gsub("chr","",fusion.data$chr_b),
                              loc.start=fusion.data$position_b,
                              loc.end=fusion.data$position_b,
                              lsn.type="fusion")

fusion.lsns=rbind.data.frame(fusion.lsnsA,
                             fusion.lsnsB)

CNA.lsns=cbind.data.frame(ID=CNA.data$Case,
                          chrom=CNA.data$Chromosome,
                          loc.start=CNA.data$Start,
                          loc.end=CNA.data$End,
                          lsn.type="copy.number")

CNA.lsns$lsn.type[CNA.data$log2_Ratio<(-0.2)]="loss"
CNA.lsns$lsn.type[CNA.data$log2_Ratio>(+0.2)]="gain"

CNA.lsns$chrom=as.character(CNA.lsns$chrom)
CNA.lsns$chrom=gsub("23","X",CNA.lsns$chrom)
CNA.lsns$chrom=gsub("24","Y",CNA.lsns$chrom)

lsn.data=rbind.data.frame(SeqMut.lsns,
                          fusion.lsns,
                          CNA.lsns)


##############################################
# Find and correct a few typos in SJTALL IDs

RNAseq.clms=colnames(RNAseq.data)
clin.RNAseq.IDs=comb.clin.data$RNAseq_id_D

not.in.clin=setdiff(RNAseq.clms,clin.RNAseq.IDs)
not.in.RNA=setdiff(clin.RNAseq.IDs,RNAseq.clms)

not.in.clin
not.in.RNA

comb.clin.data$RNAseq_id_D=gsub("SJTALL022433_D2",
                                "SJTALL022433_D1",
                                comb.clin.data$RNAseq_id_D)

comb.clin.data$RNAseq_id_D=gsub("SJTALL171_E",
                                "SJTALL171_D",
                                comb.clin.data$RNAseq_id_D)

RNAseq.clms=colnames(RNAseq.data)
clin.RNAseq.IDs=comb.clin.data$RNAseq_id_D

not.in.clin=setdiff(RNAseq.clms,clin.RNAseq.IDs)
not.in.RNA=setdiff(clin.RNAseq.IDs,RNAseq.clms)

comb.clin.data$RNAseq_id_D=gsub(not.in.RNA,"",
                                comb.clin.data$RNAseq_id_D)


###############################################
# Replace SJTALL IDs with USI IDs in RNAseq data

clin.IDs=comb.clin.data[,c("RNAseq_id_D","USI")]
RNA.IDs=cbind.data.frame(RNAseq_id_D=colnames(RNAseq.data),
                         clm.index=1:ncol(RNAseq.data))

mtch.IDs=merge(clin.IDs,RNA.IDs,by=1)

colnames(RNAseq.data)[mtch.IDs$clm.index]=mtch.IDs$USI


##########################################
# Simplify clinical data for example analysis

ex.clin.data=comb.clin.data[,c("USI","Gender_TARGET",
                               "Race_TARGET","Ethnicity",
                               "Age_at_Diagnosis_in_Days",
                               "Year_of_Diagnosis",
                               "WBC_at_Diagnosis",
                               "MRD_Day_29",
                               "Event_Free_Survival_Time_in_Days",
                               "First_Event",
                               "Overall_Survival_Time_in_Days",
                               "Vital_Status")]

colnames(ex.clin.data)=c("ID","Sex","Race","Ethnicity",
                         "Age_Days","Year_Dx",
                         "WBC","MRD29",
                         "Event_Days","First_Event",
                         "OS_Days","Vital_Status")

##############################################
# Write prepared data in an Excel file

RNAseq.data=cbind.data.frame(gene.id=rownames(RNAseq.data),
                             RNAseq.data)

enesembl.annotation.file=("D:/Genomic_Projects2_Abdelrahman/Book-chapter/Dataset-Final/ensembl.annotation.csv")
enesembl.annotation=read.csv(enesembl.annotation.file)

RNAseq.data.final=merge(enesembl.annotation,RNAseq.data,by="gene.id", all.y=TRUE)
RNAseq.data.final=RNAseq.data.final[,-1]

nice.data=list(clinical=ex.clin.data,
               lesions=lsn.data,
               expression=RNAseq.data.final)

write_xlsx(nice.data,
           paste0(local.dir,"TALL-dataset.xlsx"))


clin=ex.clin.data
lsns=lsn.data
RNA=RNAseq.data.final

save(clin,lsns,RNA,
     file=paste0(local.dir,"TALL-dataset.Rdata"))











