
#################################################
# GRIN2.0-ALEX library
#######################################################

# get.chrom.length: This function will retrieve chromosome size data from chr.info
# txt file available on the UCSC genome browser based on the specified genome assembly.

get.chrom.length=function(genome.assembly)   # one of four genome assemblies that include "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39" and "Mouse_HGCm38" can be specified
{
  # retrieve chromosome size data for GRCh38 (hg38) genome build
  if (genome.assembly=="Human_GRCh38")
  {
    chr.size.hg38= circlize::read.chromInfo(species = "hg38")
    chr.size.hg38=as.data.frame(chr.size.hg38)
    chr.size=cbind.data.frame(chrom=chr.size.hg38$chromosome,
                              size=chr.size.hg38$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))

    return(chr.size)
  }

  # retrieve chromosome size data for GRCh37 (hg19) genome build
  if (genome.assembly=="Human_GRCh37")
  {
    chr.size.hg19= circlize::read.chromInfo(species = "hg19")
    chr.size.hg19=as.data.frame(chr.size.hg19)
    chr.size=cbind.data.frame(chrom=chr.size.hg19$chromosome,
                              size=chr.size.hg19$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))

    return(chr.size)
  }

  # retrieve chromosome size data for Mouse_HGCm39 (mm39) genome build
  if (genome.assembly=="Mouse_HGCm39")
  {
    chr.size.mm39= circlize::read.chromInfo(species = "mm39")
    chr.size.mm39=as.data.frame(chr.size.mm39)
    chr.size=cbind.data.frame(chrom=chr.size.mm39$chromosome,
                              size=chr.size.mm39$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))

    return(chr.size)
  }

  # retrieve chromosome size data for Mouse_HGCm38 (mm10) genome build
  if (genome.assembly=="Mouse_HGCm38")
  {
    chr.size.mm38= circlize::read.chromInfo(species = "mm10")
    chr.size.mm38=as.data.frame(chr.size.mm38)
    chr.size=cbind.data.frame(chrom=chr.size.mm38$chromosome,
                              size=chr.size.mm38$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))

    return(chr.size)
  }
}

################################################################################
# get.ensembl.annotation: This function can be used to directly retrieve gene and
# regulatory features annotation data from ensembl biomaRt database based on the
# specified genome assembly

get.ensembl.annotation=function(genome.assembly)  # one of four genome assemblies that include "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39" and "Mouse_HGCm38" can be specified

{
  # retrieve gene annotation data for human_GRCh38 (hg38) genome assembly from ensembl biomaRt version 104
  if (genome.assembly=="Human_GRCh38")
  {
    ensembl_GRCh38 = useEnsembl(biomart="genes",
                                dataset="hsapiens_gene_ensembl",  # specify dataset for homosapiens
                                version = "104")                  # specifying version is critical to get stable search results. If we did not specify version, query will extract data from the most updated version
    chromosomes = c(1:22, "X", "Y")                               # specify chromosomes of interst
    hg38_gene_annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                               'start_position','end_position',
                                               'description','external_gene_name',
                                               'gene_biotype', 'strand', 'band'),
                                  filters = 'chromosome_name',  values = chromosomes,
                                  mart = ensembl_GRCh38)         # attributes specify data to be retreived from biomaRt database

    gene=hg38_gene_annotation[,1]
    chrom=hg38_gene_annotation[,2]
    loc.start=hg38_gene_annotation[,3]
    loc.end=hg38_gene_annotation[,4]
    description=hg38_gene_annotation[,5]
    gene.name=hg38_gene_annotation[,6]
    biotype=hg38_gene_annotation[,7]
    chrom.strand=hg38_gene_annotation[,8]
    chrom.band=hg38_gene_annotation[,9]

    gene.data=cbind.data.frame(gene=gene,  # we should change "gene" in the package functions to "Ensembl_ID"
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)

    # To retrieve data for regulatory features mapped to GRCh38
    regulatory.hg38 = useEnsembl(biomart="regulation",
                                 dataset="hsapiens_regulatory_feature", # regulatory features includes (promoters, promoter flanking regions, enhancers, CTCF binding sites, TF binding sites and open chromatin regions)
                                 version = '104')                 # specify version for stable results
    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site",
                          "TF binding site", "Open chromatin")
    hg38.regulatory= getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                        "chromosome_start","chromosome_end",
                                        "feature_type_description"),
                           filters =c("regulatory_feature_type_name", "chromosome_name")
                           , values =list(regulatory.features,chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                           mart =regulatory.hg38)

    gene=hg38.regulatory[,1]
    chrom=hg38.regulatory[,2]
    loc.start=hg38.regulatory[,3]
    loc.end=hg38.regulatory[,4]
    description=hg38.regulatory[,5]

    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)


    res=list(gene.annotation=gene.data,        # gene.data represent annotated genes
             reg.annotation=regulation.data)   # regulation.data represent predicted regulatory features
    return(res)
  }

  # retrieve gene annotation data for human_GRCh37 (hg19) genome assembly from ensembl biomaRt version 75
  if (genome.assembly=="Human_GRCh37")
  {
    ensemblGRCh37 <- useEnsembl(biomart = 'ensembl',
                                dataset = 'hsapiens_gene_ensembl',
                                version = '75')
    chromosomes = c(1:22, "X", "Y")
    hg19_gene_annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                               'start_position','end_position',
                                               'description','external_gene_id',
                                               'gene_biotype', 'strand', 'band'),
                                  filters = 'chromosome_name',  values = chromosomes,
                                  mart = ensemblGRCh37)

    gene=hg19_gene_annotation[,1]
    chrom=hg19_gene_annotation[,2]
    loc.start=hg19_gene_annotation[,3]
    loc.end=hg19_gene_annotation[,4]
    description=hg19_gene_annotation[,5]
    gene.name=hg19_gene_annotation[,6]
    biotype=hg19_gene_annotation[,7]
    chrom.strand=hg19_gene_annotation[,8]
    chrom.band=hg19_gene_annotation[,9]

    gene.data=cbind.data.frame(gene=gene,
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)

    # To retrieve data for regulatory features mapped to GRCh37
    regulatory.hg19 = useEnsembl(biomart="regulation",
                                 dataset="hsapiens_regulatory_feature",
                                 version = 'GRCh37')
    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site",
                          "TF binding site", "Open chromatin")
    hg19.regulatory= getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                        "chromosome_start","chromosome_end",
                                        "feature_type_description"),
                           filters =c("regulatory_feature_type_name", "chromosome_name")
                           , values =list(regulatory.features,chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                           mart =regulatory.hg19)

    gene=hg19.regulatory[,1]
    chrom=hg19.regulatory[,2]
    loc.start=hg19.regulatory[,3]
    loc.end=hg19.regulatory[,4]
    description=hg19.regulatory[,5]

    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)

    res=list(gene.annotation=gene.data,
             reg.annotation=regulation.data)

    return(res)
  }

  # retrieve gene annotation data for Mouse_HGCm39 genome assembly from ensembl biomaRt version 104
  if (genome.assembly=="Mouse_HGCm39")
  {
    ensembl.HGCm39 = useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl",
                                version = "104")
    chromosomes = c(1:19, "X", "Y")
    HGCm39.gene_annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                                 'start_position','end_position',
                                                 'description','external_gene_name',
                                                 'gene_biotype', 'strand', 'band'),
                                    filters = 'chromosome_name',  values = chromosomes,
                                    mart = ensembl.HGCm39)

    gene=HGCm39.gene_annotation[,1]
    chrom=HGCm39.gene_annotation[,2]
    loc.start=HGCm39.gene_annotation[,3]
    loc.end=HGCm39.gene_annotation[,4]
    description=HGCm39.gene_annotation[,5]
    gene.name=HGCm39.gene_annotation[,6]
    biotype=HGCm39.gene_annotation[,7]
    chrom.strand=HGCm39.gene_annotation[,8]
    chrom.band=HGCm39.gene_annotation[,9]

    gene.data=cbind.data.frame(gene=gene,
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)

    # To retrieve data for mouse_regulatory features mapped to mm39
    regulatory.mm39 = useEnsembl(biomart="regulation",
                                 dataset="mmusculus_regulatory_feature",
                                 version = '104')
    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site",
                          "TF binding site", "Open chromatin")
    mm39.regulatory= getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                        "chromosome_start","chromosome_end",
                                        "feature_type_description"),
                           filters =c("regulatory_feature_type_name", "chromosome_name")
                           , values =list(regulatory.features,chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                           mart =regulatory.mm39)

    gene=mm39.regulatory[,1]
    chrom=mm39.regulatory[,2]
    loc.start=mm39.regulatory[,3]
    loc.end=mm39.regulatory[,4]
    description=mm39.regulatory[,5]

    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)

    res=list(gene.annotation=gene.data,
             reg.annotation=regulation.data)
    return(res)
  }

  # retrieve gene annotation data for Mouse_HGCm38 from ensembl biomaRt version 102
  if (genome.assembly=="Mouse_HGCm38")
  {

    ensemblGRCm38 <- useEnsembl(biomart = 'genes',
                                dataset = 'mmusculus_gene_ensembl',
                                version = '102')
    chromosomes = c(1:19, "X", "Y")
    Mouse_mm10.gene_annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                                     'start_position','end_position',
                                                     'description','external_gene_name',
                                                     'gene_biotype', 'strand', 'band'),
                                        filters = 'chromosome_name',  values = chromosomes,
                                        mart = ensemblGRCm38)

    gene=Mouse_mm10.gene_annotation[,1]
    chrom=Mouse_mm10.gene_annotation[,2]
    loc.start=Mouse_mm10.gene_annotation[,3]
    loc.end=Mouse_mm10.gene_annotation[,4]
    description=Mouse_mm10.gene_annotation[,5]
    gene.name=Mouse_mm10.gene_annotation[,6]
    biotype=Mouse_mm10.gene_annotation[,7]
    chrom.strand=Mouse_mm10.gene_annotation[,8]
    chrom.band=Mouse_mm10.gene_annotation[,9]

    gene.data=cbind.data.frame(gene=gene,
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)

    # To retrieve data for mouse_regulatory features mapped to mm10
    regulatory.mm10 = useEnsembl(biomart="regulation",
                                 dataset="mmusculus_regulatory_feature",
                                 version = '102')
    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site",
                          "TF binding site", "Open chromatin")
    mm10.regulatory= getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                        "chromosome_start","chromosome_end",
                                        "feature_type_description"),
                           filters =c("regulatory_feature_type_name", "chromosome_name")
                           , values =list(regulatory.features,chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                           mart =regulatory.mm10)

    gene=mm10.regulatory[,1]
    chrom=mm10.regulatory[,2]
    loc.start=mm10.regulatory[,3]
    loc.end=mm10.regulatory[,4]
    description=mm10.regulatory[,5]

    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)

    res=list(gene.annotation=gene.data,
             reg.annotation=regulation.data)

    return(res)
  }
}

#################################################################################
# order.index.gene.data: This function order and index gene data by chromosome,
# start location, and end location.

order.index.gene.data=function(gene.data)
{
  g=nrow(gene.data)
  gene.ord=order(gene.data[,"chrom"],
                 gene.data[,"loc.start"],
                 gene.data[,"loc.end"])

  gene.data=gene.data[gene.ord,]
  new.chrom=which(gene.data[-1,"chrom"]!=gene.data[-g,"chrom"])
  chr.start=c(1,new.chrom+1)
  chr.end=c(new.chrom,g)
  gene.index=cbind.data.frame(chrom=gene.data[chr.start,"chrom"],
                              row.start=chr.start,
                              row.end=chr.end)
  gene.data$gene.row=1:g

  res=list(gene.data=gene.data,
           gene.index=gene.index)
  return(res)
}

#################################################################################
# order.index.lsn.data: This function order and index lesion data by type,
# chromosome, and subject.

order.index.lsn.data=function(lsn.data)

{
  l=nrow(lsn.data)
  lsn.ord=order(lsn.data[,"lsn.type"],
                lsn.data[,"chrom"],
                lsn.data[,"ID"])
  lsn.data=lsn.data[lsn.ord,]

  lsn.chng=which((lsn.data[-1,"lsn.type"]!=lsn.data[-l,"lsn.type"])|
                   (lsn.data[-1,"chrom"]!=lsn.data[-l,"chrom"])|
                   (lsn.data[-1,"ID"]!=lsn.data[-l,"ID"]))
  lsn.start=c(1,lsn.chng+1)
  lsn.end=c(lsn.chng,l)
  lsn.index=cbind.data.frame(lsn.type=lsn.data[lsn.start,"lsn.type"],
                             chrom=lsn.data[lsn.start,"chrom"],
                             ID=lsn.data[lsn.start,"ID"],
                             row.start=lsn.start,
                             row.end=lsn.end)
  lsn.data$lsn.row=1:l

  res=list(lsn.data=lsn.data,
           lsn.index=lsn.index)
  return(res)
}

##################################################################################
# prep.gene.lsn.data:This function prepare gene and lesion data for later computations.

prep.gene.lsn.data=function(lsn.data,       # lesion data: ID, chrom, loc.start, loc.end, lsn.type
                            gene.data,      # gene locus data: gene, chrom, loc.start, loc.end
                            mess.freq=10)   # message frequency: display message every mess.freq^{th} lesion block
{
  saf=options()$stringsAsFactors
  options(stringsAsFactors=F)

  # order lesion data by type, chromosome, and subject
  lsn.dset=order.index.lsn.data(lsn.data)
  lsn.data=lsn.dset$lsn.data
  lsn.index=lsn.dset$lsn.index

  # order and index gene locus data by chromosome and position
  gene.dset=order.index.gene.data(gene.data)
  gene.data=gene.dset$gene.data
  gene.index=gene.dset$gene.index

  # Extract some basic information
  g=nrow(gene.data) # number of genes
  l=nrow(lsn.data)  # number of lesions

  # Create gene position data
  message(paste0("Formatting gene position data for counting: ",date()))
  gene.pos.data=rbind.data.frame(cbind.data.frame(ID="",  # gene start data
                                                  lsn.type="",
                                                  lsn.row=NA,
                                                  gene=gene.data[,"gene"],
                                                  gene.row=gene.data[,"gene.row"],
                                                  chrom=gene.data[,"chrom"],
                                                  pos=gene.data[,"loc.start"],
                                                  cty=1),
                                 cbind.data.frame(ID="", # gene end data
                                                  lsn.type="",
                                                  lsn.row=NA,
                                                  gene=gene.data[,"gene"],
                                                  gene.row=gene.data[,"gene.row"],
                                                  chrom=gene.data[,"chrom"],
                                                  pos=gene.data[,"loc.end"],
                                                  cty=4)
  )
  # order gene position data
  ord=order(gene.pos.data[,"chrom"],
            gene.pos.data[,"pos"],
            gene.pos.data[,"cty"])
  gene.pos.data=gene.pos.data[ord,]

  # Create lesion position data with one row for each edge of each lesion
  message(paste0("Formatting lesion position data for counting:  ",date()))
  lsn.pos.data=rbind.data.frame(cbind.data.frame(ID=lsn.data[,"ID"],
                                                 lsn.type=lsn.data[,"lsn.type"],
                                                 lsn.row=lsn.data[,"lsn.row"],
                                                 gene="",
                                                 gene.row=NA,
                                                 chrom=lsn.data[,"chrom"],
                                                 pos=lsn.data[,"loc.start"],
                                                 cty=2),
                                cbind.data.frame(ID=lsn.data[,"ID"],
                                                 lsn.type=lsn.data[,"lsn.type"],
                                                 lsn.row=lsn.data[,"lsn.row"],
                                                 gene="",
                                                 gene.row=NA,
                                                 chrom=lsn.data[,"chrom"],
                                                 pos=lsn.data[,"loc.end"],
                                                 cty=3))
  # order lesion position data
  ord=order(lsn.pos.data[,"chrom"],
            lsn.pos.data[,"pos"],
            lsn.pos.data[,"cty"])
  lsn.pos.data=lsn.pos.data[ord,]

  # Combine gene & lesion data
  message(paste0("Combining formatted gene and lesion postion data: ",date()))
  gene.lsn.data=rbind.data.frame(gene.pos.data,
                                 lsn.pos.data)

  # Order and index gene & lesion data
  ord=order(gene.lsn.data[,"chrom"],
            gene.lsn.data[,"pos"],
            gene.lsn.data[,"cty"])
  gene.lsn.data=gene.lsn.data[ord,]
  m=nrow(gene.lsn.data)
  gene.lsn.data[,"glp.row"]=1:m

  # compute vector to order gene.lsn.data by lsn.row and gene.row
  ord=order(gene.lsn.data[,"lsn.row"],
            gene.lsn.data[,"gene.row"],
            gene.lsn.data[,"cty"])

  # use that vector to add gene.lsn.data row.start and row.end indices to lsn.data
  lsn.pos=gene.lsn.data[ord[1:(2*l)],]
  lsn.data[,"glp.row.start"]=lsn.pos[2*(1:l)-1,"glp.row"]
  lsn.data[,"glp.row.end"]=lsn.pos[2*(1:l),"glp.row"]


  # use that vector to add gene.lsn.data row.start and row.end indices to gene.data
  gene.pos=gene.lsn.data[ord[-(1:(2*l))],]
  gene.data[,"glp.row.start"]=gene.pos[2*(1:g)-1,"glp.row"]
  gene.data[,"glp.row.end"]=gene.pos[2*(1:g),"glp.row"]

  # Double-check table pointers from lsn.data and gene.data to gene.lsn.data

  message(paste0("Verifying structure of combined gene and lesion data: ",date()))
  glp.gene.start=gene.lsn.data[gene.data$glp.row.start,c("gene","chrom","pos")]
  colnames(glp.gene.start)=c("gene","chrom","loc.start")
  ok.glp.gene.start=all(glp.gene.start==gene.data[,c("gene","chrom","loc.start")])

  glp.gene.end=gene.lsn.data[gene.data$glp.row.end,c("gene","chrom","pos")]
  colnames(glp.gene.end)=c("gene","chrom","loc.end")
  ok.glp.gene.end=all(glp.gene.end==gene.data[,c("gene","chrom","loc.end")])

  glp.lsn.start=gene.lsn.data[lsn.data$glp.row.start,c("ID","chrom","pos","lsn.type")]
  colnames(glp.lsn.start)=c("ID","chrom","loc.start","lsn.type")
  ok.glp.lsn.start=all(glp.lsn.start==lsn.data[,c("ID","chrom","loc.start","lsn.type")])

  glp.lsn.end=gene.lsn.data[lsn.data$glp.row.end,c("ID","chrom","pos","lsn.type")]
  colnames(glp.lsn.end)=c("ID","chrom","loc.end","lsn.type")
  ok.glp.lsn.end=all(glp.lsn.end==lsn.data[,c("ID","chrom","loc.end","lsn.type")])

  # Double-check table pointers from gene.lsn.data to gene.data and lsn.data
  glp.gene.start=gene.lsn.data[gene.lsn.data$cty==1,c("gene.row","gene","chrom","pos")]
  ok.gene.start=all(glp.gene.start[,c("gene","chrom","pos")]==gene.data[glp.gene.start$gene.row,c("gene","chrom","loc.start")])

  glp.gene.end=gene.lsn.data[gene.lsn.data$cty==4,c("gene.row","gene","chrom","pos")]
  ok.gene.end=all(glp.gene.end[,c("gene","chrom","pos")]==gene.data[glp.gene.end$gene.row,c("gene","chrom","loc.end")])

  glp.lsn.start=gene.lsn.data[gene.lsn.data$cty==2,c("lsn.row","ID","chrom","pos","lsn.type")]
  ok.lsn.start=all(glp.lsn.start[,c("ID","chrom","pos","lsn.type")]==lsn.data[glp.lsn.start$lsn.row,c("ID","chrom","loc.start","lsn.type")])

  glp.lsn.end=gene.lsn.data[gene.lsn.data$cty==3,c("lsn.row","ID","chrom","pos","lsn.type")]
  ok.lsn.end=all(glp.lsn.end[,c("ID","chrom","pos","lsn.type")]==lsn.data[glp.lsn.end$lsn.row,c("ID","chrom","loc.end","lsn.type")])

  all.ok=all(c(ok.glp.gene.start,ok.glp.gene.end,
               ok.glp.lsn.start,ok.glp.lsn.end,
               ok.gene.start,ok.gene.end,
               ok.lsn.start,ok.lsn.end))

  if (!all.ok)
    stop("Error in constructing and indexing combined lesion and gene data.")

  message(paste0("Verified correct construction and indexing of combined lesion and gene data: ",date()))

  return(list(lsn.data=lsn.data, # Input lesion data
              gene.data=gene.data, # Input gene annotation data
              gene.lsn.data=gene.lsn.data, # Data.frame ordered by gene and lesions start position. Gene start position is coded as 1 in the cty column and gene end position is coded as 4. Lesion start position is coded as 2 in the cty column and lesion end position is coded as 3.
              gene.index=gene.index, # Data.frame that shows ordered row start and row end for each chromosome in the gene.lsn.data table
              lsn.index=lsn.index)) # Data.frame that shows row start and row end for each lesion in the gene.lsn.data table

}

#################################################################################
# find.gene.lsn.overlaps: Use the results of prep.gene.lsn.data function to find
# lesion-gene overlaps.

find.gene.lsn.overlaps=function(gl.data) # A list of five data.frames that represent the output results of the prep.gene.lsn.data function

{
  gene.data=gl.data$gene.data
  lsn.data=gl.data$lsn.data
  lsn.index=gl.data$lsn.index
  gene.index=gl.data$gene.index
  gene.lsn.data=gl.data$gene.lsn.data
  m=nrow(gene.lsn.data)

  message(paste0("Scanning through combined lesion and gene data to find gene-lesion overlaps: ",date()))


  gene.row.mtch=NULL  # initialize vector for rows of gene data matched to rows of lesion data
  lsn.row.mtch=NULL   # initialize vector for rows of lesion data matched to rows of gene data
  current.genes=NULL  # initialize vector of genes overlapping this point of the scan
  current.lsns=NULL   # initialize vector of lesions overlapping this point of the scan
  for (i in 1:m)      # loop over rows of gene.lsn.data
  {
    # enter a gene
    if (gene.lsn.data$cty[i]==1)
    {
      # add this gene to the set of current genes
      current.genes=c(current.genes,
                      gene.lsn.data$gene.row[i])

      # match this gene to set of current lesions
      lsn.row.mtch=c(lsn.row.mtch,current.lsns)          # add current lesions to lsn.row.mtch
      gene.row.mtch=c(gene.row.mtch,                     # add this gene for each current lesion
                      rep(gene.lsn.data$gene.row[i],
                          length(current.lsns)))

    }

    # exit a gene
    if (gene.lsn.data$cty[i]==4)
    {
      # drop this gene from the set of current genes
      current.genes=setdiff(current.genes,
                            gene.lsn.data$gene.row[i])
    }


    if (gene.lsn.data$cty[i]==2)       # enter a lesion
    {
      lsn.row.mtch=c(lsn.row.mtch,
                     rep(gene.lsn.data$lsn.row[i],length(current.genes)))
      gene.row.mtch=c(gene.row.mtch,current.genes)
      current.lsns=c(current.lsns,
                     gene.lsn.data$lsn.row[i])
    }

    if (gene.lsn.data$cty[i]==3)
    {
      current.lsns=setdiff(current.lsns,
                           gene.lsn.data$lsn.row[i])
    }
  }

  message(paste0("Completed scan of combined gene-lesion data: ",date()))

  # Generate the gene-lesion hit data
  gene.lsn.hits=cbind.data.frame(gene.data[gene.row.mtch,
                                           c("gene.row","gene","chrom","loc.start","loc.end")],
                                 lsn.data[lsn.row.mtch,
                                          c("lsn.row","ID","chrom","loc.start","loc.end","lsn.type")])

  colnames(gene.lsn.hits)=c("gene.row","gene","gene.chrom","gene.loc.start","gene.loc.end",
                            "lsn.row","ID","lsn.chrom","lsn.loc.start","lsn.loc.end","lsn.type")

  res=list(lsn.data=lsn.data, # Input lesion data
           gene.data=gene.data, # Input gene annotation data
           gene.lsn.data=gene.lsn.data, # Data.frame ordered by gene and lesions start position. Gene start position is coded as 1 in the cty column and gene end position is coded as 4. Lesion start position is coded as 2 in the cty column and lesion end position is coded as 3
           gene.lsn.hits=gene.lsn.hits, # Each row represent a gene overlapped by a certain lesion. Gene column shows the overlapped gene and ID column has the patient ID
           gene.index=gene.index, # Data.frame that shows row start and row end for each chromosome in the gene.lsn.data table
           lsn.index=lsn.index) # Data.frame that shows row start and row end for each lesion in the gene.lsn.data table

  return(res)

}

##################################################################################
# count.hits: The function computes the number of hits affecting each gene by lesion
# type. It also compute the number of subjects with a hit in each annotated gene.

count.hits=function(ov.data) # result of find.gene.lsn.overlaps function
{
  lsn.data=ov.data$lsn.data
  lsn.index=ov.data$lsn.index
  gene.lsn.hits=ov.data$gene.lsn.hits
  gene.lsn.data=ov.data$gene.lsn.data
  gene.data=ov.data$gene.data
  gene.index=ov.data$gene.index

  g=nrow(gene.data)

  # Compute the number of hits matrix
  lsn.types=sort(unique(lsn.index[,"lsn.type"]))
  k=length(lsn.types)

  nhit.mtx=matrix(0,g,k)
  colnames(nhit.mtx)=lsn.types

  nhit.tbl=table(gene.lsn.hits$gene.row,
                 gene.lsn.hits$lsn.type)
  nhit.rows=as.numeric(rownames(nhit.tbl))

  for (i in 1:ncol(nhit.tbl))
    nhit.mtx[nhit.rows,colnames(nhit.tbl)[i]]=nhit.tbl[,i]

  # Compute the matrix of the number of subjects with a hit
  gene.subj.type=paste0(gene.lsn.hits$gene.row,"_",
                        gene.lsn.hits$ID,"_",
                        gene.lsn.hits$lsn.type)
  dup.gene.subj.type=duplicated(gene.subj.type)

  subj.gene.hits=gene.lsn.hits[!dup.gene.subj.type,]

  nsubj.mtx=matrix(0,g,k)
  colnames(nsubj.mtx)=lsn.types

  nsubj.tbl=table(subj.gene.hits$gene.row,
                  subj.gene.hits$lsn.type)

  nsubj.rows=as.numeric(rownames(nsubj.tbl))

  for (i in 1:ncol(nsubj.tbl))
    nsubj.mtx[nsubj.rows,colnames(nsubj.tbl)[i]]=nsubj.tbl[,i]

  res=list(lsn.data=lsn.data,  # Input lesion data
           lsn.index=lsn.index, # Data.frame that shows row start and row end for each lesion in the gene.lsn.data table
           gene.data=gene.data, # Input gene annotation data
           gene.index=gene.index, # Data.frame that shows ordered row start and row end for each chromosome in the gene.lsn.data table
           nhit.mtx=nhit.mtx, # A data matrix with number of hits in each gene by lesion type
           nsubj.mtx=nsubj.mtx, # A data matrix with number of affected subjects by lesion type
           gene.lsn.data=gene.lsn.hits, # Each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene and ID column has the patient ID
           glp.data=gene.lsn.data) # Data.frame ordered by gene and lesions start position. Gene start position is coded as 1 in the cty column and gene end position is coded as 4. Lesion start position is coded as 2 in the cty column and lesion end position is coded as 3

  return(res)

}

#################################################################################
# Following four functions row.prob.subj.hit, row.bern.conv, p.order and pc06.fdr
# will be used by prob.hits function to compute p-values for the probablity of each gene
# to be affected by one or a constellation of multiple types of lesions

# row.prob.subj.hit: Compute the probability that a subject has a hit for each gene.

row.prob.subj.hit=function(P,   # matrix of lesion hit probabilities, rows for genes, columns for lesion types
                           IDs) # vector of subject IDs for each lesion, length must equal ncol(P)
{
  if (length(IDs)!=ncol(P))
    stop("length(IDs) must equal ncol(P).")

  g=nrow(P)
  l=ncol(P)

  ord=order(IDs)
  IDs=IDs[ord]
  P=matrix(P[,ord],g,l)
  l=length(IDs)

  new.ID=which(IDs[-1]!=IDs[-l])
  ID.start=c(1,new.ID+1)
  ID.end=c(new.ID,l)

  n=length(ID.start)
  m=nrow(P)
  Pr=matrix(NA,m,n)
  for (i in 1:n)
  {
    pr.mtx=matrix(P[,ID.start[i]:ID.end[i]],m,ID.end[i]-ID.start[i]+1)
    Pr[,i]=1-exp(rowSums(log(1-pr.mtx)))
  }

  return(Pr)
}

################################################################################
# row.bern.conv: This function Compute a convolution of Bernoullis for each row of
# a Bernoulli success probability matrix.

row.bern.conv=function(P, #matrix of lesion hit probabilities, rows for genes, columns for lesion types
                       max.x=NULL) # Maximum number of subjects or maximum number of hits

{
  m=nrow(P)
  n=ncol(P)
  if (is.null(max.x))
    max.x=(ncol(P))
  Pr=matrix(0,m,max.x+1)
  Pr[,1]=1

  for (i in 1:n)
  {
    P1=Pr*P[,i]
    P0=Pr*(1-P[,i])
    Pr0=P0
    Pr1=cbind(0,P1)
    Pr1[,max.x+1]=Pr1[,max.x+1]+Pr1[,max.x+2]
    Pr=Pr0+Pr1[,-(max.x+2)]
  }
  rs.Pr=rowSums(Pr)
  Pr=Pr/rs.Pr
  return(Pr)
}

#################################################################################
# p.order: This function compute ordered p-values and prepare the data for the
# constellation test which evaluates if the gene is affected by multiple types of lesions.

p.order=function(P) # matrix of lesion hit probabilities, rows for genes, columns for lesion types
{
  k=ncol(P)
  p.mtx=apply(P,1,sort,na.last=T)
  p.mtx=t(p.mtx)
  #n.pvals=rowSums(!is.na(p.mtx))

  res=p.mtx
  for (i in 1:k)
  {
    res[,i]=pbeta(p.mtx[,i],i,k-i+1)
  }
  return(res)
}

###############################################################################
# pc06.fdr: Compute FDR with the Pounds & Cheng (2006) estimator of
# the proportion of tests with a true null (pi.hat)
# Reference: https://pubmed.ncbi.nlm.nih.gov/16777905/

pc06.fdr=function(p)

{
  pi.hat=min(1,2*mean(p,na.rm=T))
  q=pi.hat*p.adjust(p,method="fdr")
  return(q)
}


#################################################################################
# prob.hits: The function evaluates the probablity of a gene to be affected by one
# or a constellation of multiple types of lesions.

prob.hits=function(hit.cnt, # Output results of the count.hits function
                   chr.size=NULL) # A data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs. It should has two columns named "chrom" with the chromosome number and "size" for the size of the chromosome in base pairs.
{
  if (is.null(chr.size))
    chr.size=impute.chrom.size(hit.cnt$lsn.data,
                               hit.cnt$gene.data)

  ###################
  # order and index gene.lsn.data
  ord=order(hit.cnt$gene.lsn.data$lsn.type,
            hit.cnt$gene.lsn.data$lsn.chrom,
            hit.cnt$gene.lsn.data$gene.row,
            hit.cnt$gene.lsn.data$ID)
  hit.cnt$gene.lsn.data=hit.cnt$gene.lsn.data[ord,]
  m=nrow(hit.cnt$gene.lsn.data)

  new.sect=which((hit.cnt$gene.lsn.data$gene.chrom[-1]!=hit.cnt$gene.lsn.data$gene.chrom[-m])|
                   (hit.cnt$gene.lsn.data$lsn.type[-1]!=hit.cnt$gene.lsn.data$lsn.type[-m])|
                   (hit.cnt$gene.lsn.data$gene.row[-1]!=hit.cnt$gene.lsn.data$gene.row[-m]))
  sect.start=c(1,new.sect+1)
  sect.end=c(new.sect,m)
  gene.lsn.index=cbind.data.frame(lsn.type=hit.cnt$gene.lsn.data$lsn.type[sect.start],
                                  chrom=hit.cnt$gene.lsn.data$gene.chrom[sect.start],
                                  gene.row=hit.cnt$gene.lsn.data$gene.row[sect.start],
                                  row.start=sect.start,
                                  row.end=sect.end,
                                  n.lsns=sect.end-sect.start+1)
  k=nrow(gene.lsn.index)

  new.chr=which(gene.lsn.index$chrom[-1]!=gene.lsn.index$chrom[-k])
  chr.start=c(1,new.chr+1)
  chr.end=c(new.chr,k)
  gene.lsn.chr.index=cbind.data.frame(lsn.type=gene.lsn.index$lsn.type[chr.start],
                                      chrom=gene.lsn.index$chrom[chr.start],
                                      row.start=chr.start,
                                      row.end=chr.end,
                                      n.rows=chr.end-chr.start+1)

  nr.li=nrow(hit.cnt$lsn.index)
  new.chr=which((hit.cnt$lsn.index$lsn.type[-1]!=hit.cnt$lsn.index$lsn.type[-nr.li])|
                  (hit.cnt$lsn.index$chrom[-1]!=hit.cnt$lsn.index$chrom[-nr.li]))
  chr.start=c(1,new.chr+1)
  chr.end=c(new.chr,nr.li)
  lsn.chr.index=cbind.data.frame(lsn.type=hit.cnt$lsn.index$lsn.type[chr.start],
                                 chrom=hit.cnt$lsn.index$chrom[chr.start],
                                 row.start=chr.start,
                                 row.end=chr.end)

  b=nrow(gene.lsn.chr.index)
  g=nrow(hit.cnt$nhit.mtx)
  nlt=ncol(hit.cnt$nhit.mtx)

  p.nsubj=p.nhit=matrix(1,g,nlt)
  colnames(p.nsubj)=colnames(p.nhit)=colnames(hit.cnt$nhit.mtx)

  for (i in 1:b)
  {
    # find rows for affected genes
    gli.start.row=gene.lsn.chr.index$row.start[i]
    gli.end.row=gene.lsn.chr.index$row.end[i]
    gld.start.row=gene.lsn.index$row.start[gli.start.row]
    gld.end.row=gene.lsn.index$row.end[gli.end.row]
    gld.rows=gld.start.row:gld.end.row
    gene.rows=unique(hit.cnt$gene.lsn.data$gene.row[gld.rows])
    n.genes=length(gene.rows)

    # find rows for lesions of this type on this chromosomes
    lsn.chr.mtch=which((lsn.chr.index$lsn.type==gene.lsn.chr.index$lsn.type[i])&
                         (lsn.chr.index$chrom==gene.lsn.chr.index$chrom[i]))
    lsn.index.start.row=lsn.chr.index$row.start[lsn.chr.mtch]
    lsn.index.end.row=lsn.chr.index$row.end[lsn.chr.mtch]
    lsn.start.row=hit.cnt$lsn.index$row.start[lsn.index.start.row]
    lsn.end.row=hit.cnt$lsn.index$row.end[lsn.index.end.row]
    lsn.rows=lsn.start.row:lsn.end.row
    n.lsns=length(lsn.rows)
    lsn.type=hit.cnt$lsn.data$lsn.type[lsn.start.row]

    message(paste0("Computing p-values for ",
                   n.genes," gene(s) on chromosome ",
                   gene.lsn.chr.index$chrom[i],
                   " affected by ",
                   n.lsns," ",lsn.type,
                   " (data block ",i," of ",b,"): ",date()))

    # find chromosome size
    chr.mtch=which(gene.lsn.chr.index$chrom[i]==chr.size$chrom)
    chrom.size=chr.size$size[chr.mtch]

    # obtain gene sizes, lesion sizes, and gene hit probabilities
    lsn.size=hit.cnt$lsn.data$loc.end[lsn.rows]-hit.cnt$lsn.data$loc.start[lsn.rows]+1
    gene.size=hit.cnt$gene.data$loc.end[gene.rows]-hit.cnt$gene.data$loc.start[gene.rows]+1

    log.pr=log(rep(lsn.size,each=n.genes)+rep(gene.size,times=n.lsns))-log(chrom.size)
    pr.gene.hit=matrix(exp(log.pr),n.genes,n.lsns)
    pr.gene.hit[pr.gene.hit>1]=1

    lsn.subj.IDs=hit.cnt$lsn.data$ID[lsn.rows]
    pr.subj=row.prob.subj.hit(pr.gene.hit,lsn.subj.IDs)

    max.nsubj=max(hit.cnt$nsubj.mtx[gene.rows,lsn.type])
    max.nhit=max(hit.cnt$nhit.mtx[gene.rows,lsn.type])

    pr.nhit=row.bern.conv(pr.gene.hit,max.nhit)
    pr.nsubj=row.bern.conv(pr.subj,max.nsubj)

    for (j in 1:n.genes)
    {
      nsubj=hit.cnt$nsubj.mtx[gene.rows[j],lsn.type]
      nhit=hit.cnt$nhit.mtx[gene.rows[j],lsn.type]
      p.nsubj[gene.rows[j],lsn.type]=sum(pr.nsubj[j,(nsubj+1):(max.nsubj+1)])
      p.nhit[gene.rows[j],lsn.type]=sum(pr.nhit[j,(nhit+1):(max.nhit+1)])
    }
  }

  rownames(p.nhit)=rownames(hit.cnt$nhit.mtx)
  rownames(p.nsubj)=rownames(hit.cnt$nsubj.mtx)
  colnames(hit.cnt$nhit.mtx)=paste0("nhit.",colnames(hit.cnt$nhit.mtx))
  colnames(hit.cnt$nsubj.mtx)=paste0("nsubj.",colnames(hit.cnt$nsubj.mtx))
  colnames(p.nhit)=paste0("p.",colnames(hit.cnt$nhit.mtx))
  colnames(p.nsubj)=paste0("p.",colnames(hit.cnt$nsubj.mtx))

  # Compute q-values
  message(paste0("Computing q-values: ",date()))
  q.nhit=p.nhit
  q.nsubj=p.nsubj
  for (i in 1:ncol(q.nhit))
  {
    pi.hat=min(1,2*mean(p.nhit[,i],na.rm=T))
    q.nhit[,i]=pi.hat*p.adjust(p.nhit[,i],method="fdr")
    pi.hat=min(1,2*mean(p.nsubj[,i],na.rm=T))
    q.nsubj[,i]=pi.hat*p.adjust(p.nsubj[,i],method="fdr")
  }

  colnames(q.nhit)=paste0("q.",colnames(hit.cnt$nhit.mtx))
  colnames(q.nsubj)=paste0("q.",colnames(hit.cnt$nsubj.mtx))

  # Now get ordered p-values
  message(paste0("Computing p-values for number of lesion types affecting genes: ",date()))
  p.ord.nhit=p.order(p.nhit)
  colnames(p.ord.nhit)=paste0("p",1:ncol(p.nhit),".nhit")

  p.ord.nsubj=p.order(p.nsubj)
  colnames(p.ord.nsubj)=paste0("p",1:ncol(p.nsubj),".nsubj")

  # q-values of ordered p-values
  q.ord.nhit=p.ord.nhit
  q.ord.nsubj=p.ord.nsubj
  message(paste0("Computing q-values for number of lesion types affecting genes: ",date()))
  for (i in 1:ncol(p.ord.nhit))
  {
    pi.hat=min(1,2*mean(p.ord.nhit[,i],na.rm=T))
    q.ord.nhit[,i]=pi.hat*p.adjust(p.ord.nhit[,i],method="fdr")
    pi.hat=min(1,2*mean(p.ord.nsubj[,i],na.rm=T))
    q.ord.nsubj[,i]=pi.hat*p.adjust(p.ord.nsubj[,i],method="fdr")
  }
  colnames(q.ord.nsubj)=paste0("q",1:ncol(p.nsubj),".nsubj")
  colnames(q.ord.nhit)=paste0("q",1:ncol(p.nhit),".nhit")

  gd.clms=setdiff(colnames(hit.cnt$gene.data),c("glp.row.start","glp.row.end"))
  lsn.clms=setdiff(colnames(hit.cnt$lsn.data),c("glp.row.start","glp.row.end"))

  gd.clms=c("gene.row",setdiff(gd.clms,"gene.row"))
  lsn.clms=c("lsn.row",setdiff(lsn.clms,"lsn.row"))

  gene.res=cbind.data.frame(hit.cnt$gene.data[,gd.clms],
                            hit.cnt$nsubj.mtx,
                            p.nsubj,
                            q.nsubj,
                            p.ord.nsubj,
                            q.ord.nsubj,
                            hit.cnt$nhit.mtx,
                            p.nhit,
                            q.nhit,
                            p.ord.nhit,
                            q.ord.nhit)

  res=list(gene.hits=gene.res, # A data table of GRIN results that includes gene annotation, number of subjects and number of hits affecting each locus. in addition, p and FDR adjusted q-values showing the probability of each locus being affected by one or a constellation of multiple types of lesions are also included in the GRIN results.
           lsn.data=hit.cnt$lsn.data[,lsn.clms], # Input lesion data
           gene.data=hit.cnt$gene.data[,gd.clms], # Input gene data
           gene.lsn.data=hit.cnt$gene.lsn.data, # Each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene and ID column has the patient ID
           chr.size=chr.size, # A data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs
           gene.index=hit.cnt$gene.index, # Data.frame that shows row start and row end for each chromosome in the gene.lsn.data table
           lsn.index=hit.cnt$lsn.index) # Data.frame that shows row start and row end for each lesion in the gene.lsn.data table

  return(res)
}

##############################################################################
# grin.stats: The function run the Genomic Random Interval (GRIN) analysis to
# determine whether a certain locus has an abundance of lesions that is
# statistically significant.

grin.stats=function(lsn.data,            # data.frame with columns ID (subject identifier), chrom, loc.start, loc.end, lsn.type
                    gene.data=NULL,      # data.frame with columns gene, chrom, loc.start, loc.end
                    chr.size=NULL,       # data.frame with columns chrom and size
                    genome.version=NULL) # character string with genome version. Currently, four genome assemblies are accepted including "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39", and "Mouse_HGCm38"

{
  if (is.null(genome.version)&&(is.null(gene.data)||is.null(chr.size)))
  {
    genome.version=select.list(c("Human_GRCh38",
                                 "Human_GRCh37",
                                 "Mouse_HGCm39",
                                 "Mouse_HGCm38"))
  }

  if (is.character(genome.version))
  {
    ensembl.data=get.ensembl.annotation(genome.version)
    if (is.null(gene.data))
    {
      gene.data=ensembl.data$gene.annotation
    }
    chrom.size=get.chrom.length(genome.version)
    if (is.null(chr.size))
    {
      chr.size=chrom.size
    }
  }

  prep.data=prep.gene.lsn.data(lsn.data,
                               gene.data)
  find.overlap=find.gene.lsn.overlaps(prep.data)
  hit.cnt=count.hits(find.overlap)
  hit.pvals=prob.hits(hit.cnt,
                      chr.size)
  return(hit.pvals)
}

##################################################################################
# write.grin.xlsx: This function Write GRIN results to an excel file with multiple
# tabs that include GRIN results, lesion data, gene annotation data, chromosome size,
# gene-lesion overlap and methods paragraph.

write.grin.xlsx=function(grin.result, # output of the grin.stats function
                         output.file) # output file name ".xlsx"
{
  saf=options()$stringsAsFactors
  options(stringsAsFactors = F)

  rpt.res=grin.result[c("gene.hits",
                        "gene.lsn.data",
                        "lsn.data",
                        "gene.data",
                        "chr.size")]

  ################################
  # Contents of the data sheets
  sheet.int=cbind.data.frame(sheet.name=c("gene.hits",
                                          "gene.lsn.data",
                                          "lsn.data",
                                          "gene.data",
                                          "chr.size"),
                             col.name="entire sheet",
                             meaning=c("GRIN statistical results",
                                       "gene-lesion overlaps",
                                       "input lesion data",
                                       "input gene location data",
                                       "input chromosome size data"))

  ######################################
  # Interpretation of gene hit stats
  gh.cols=colnames(rpt.res[["gene.hits"]])
  genehit.int=cbind.data.frame(sheet.name="gene.hits",
                               col.name=gh.cols,
                               meaning=gh.cols)
  rownames(genehit.int)=gh.cols
  rpt.clms=c("gene.row","gene","loc.start","loc.end")
  rpt.defs=c("gene data row index",
             "gene name",
             "locus of left edge of gene",
             "locus of right edge of gene")

  rpt.indx=rep(NA,length(rpt.clms))
  for (i in 1:length(rpt.clms))
    rpt.indx[i]=which(genehit.int$col.name==rpt.clms[i])

  genehit.int$meaning[rpt.indx]=rpt.defs

  nsubj.clms=which(substring(gh.cols,1,5)=="nsubj")
  genehit.int$meaning[nsubj.clms]=paste0("Number of subjects with a ",
                                         substring(gh.cols[nsubj.clms],7)," lesion ",
                                         genehit.int$meaning[nsubj.clms])

  p.nsubj.clms=which(substring(gh.cols,1,7)=="p.nsubj")
  genehit.int$meaning[p.nsubj.clms]=paste0("p-value for the number of subjects with a ",
                                           substring(gh.cols[p.nsubj.clms],9)," lesion ",
                                           genehit.int$meaning[p.nsubj.clms])

  q.nsubj.clms=which(substring(gh.cols,1,7)=="q.nsubj")
  genehit.int$meaning[q.nsubj.clms]=paste0("FDR estimate for the number of subjects with a ",
                                           substring(gh.cols[q.nsubj.clms],9)," lesion ",
                                           genehit.int$meaning[q.nsubj.clms])

  nhit.clms=which(substring(gh.cols,1,4)=="nhit")
  genehit.int$meaning[nhit.clms]=paste0("Number of ",
                                        substring(gh.cols[nhit.clms],6)," lesions ",
                                        genehit.int$meaning[nhit.clms])

  p.nhit.clms=which(substring(gh.cols,1,6)=="p.nhit")
  genehit.int$meaning[p.nhit.clms]=paste0("p-value for the number of ",
                                          substring(gh.cols[p.nhit.clms],8)," lesions ",
                                          genehit.int$meaning[p.nhit.clms])

  q.nhit.clms=which(substring(gh.cols,1,6)=="q.nhit")
  genehit.int$meaning[q.nhit.clms]=paste0("FDR estimate for the number of ",
                                          substring(gh.cols[q.nhit.clms],8)," lesions ",
                                          genehit.int$meaning[q.nhit.clms])


  p1.nsubj.clms=paste0("p",1:length(nsubj.clms),".nsubj")
  genehit.int[p1.nsubj.clms,"meaning"]=paste0("p-value for the number of subjects ",
                                              "with any ",1:length(nsubj.clms)," type(s) of lesion ",
                                              "overlapping the gene locus")

  q1.nsubj.clms=paste0("q",1:length(nsubj.clms),".nsubj")
  genehit.int[q1.nsubj.clms,"meaning"]=paste0("FDR estimate for the number of subjects ",
                                              "with any ",1:length(nsubj.clms)," type(s) of lesion ",
                                              "overlapping the gene locus")

  p1.nhit.clms=paste0("p",1:length(nhit.clms),".nhit")
  genehit.int[p1.nhit.clms,"meaning"]=paste0("p-value for the number of ",
                                             "any ",1:length(nhit.clms)," type(s) of lesion ",
                                             "overlapping the gene locus")

  q1.nhit.clms=paste0("q",1:length(nhit.clms),".nhit")
  genehit.int[q1.nhit.clms,"meaning"]=paste0("FDR estimate for the number of",
                                             " any ",1:length(nhit.clms)," type(s) of lesion ",
                                             "overlapping the gene locus")

  #######################################
  # Lesion data interpretation
  lsn.clms=colnames(rpt.res[["lsn.data"]])
  lsn.int=cbind.data.frame(sheet.name="lsn.data",
                           col.name=lsn.clms,
                           meaning=lsn.clms)

  rpt.clms=c("ID","chrom","loc.start","loc.end","lsn.type")
  int.clms=c(which(lsn.int[,"col.name"]=="ID"),
             which(lsn.int[,"col.name"]=="chrom"),
             which(lsn.int[,"col.name"]=="loc.start"),
             which(lsn.int[,"col.name"]=="loc.end"),
             which(lsn.int[,"col.name"]=="lsn.type"))
  lsn.int[int.clms,"meaning"]=c("Input Subject Identifier",
                                "Input Chromosome",
                                "Input Gene Locus Left Edge",
                                "Input Gene Locus Right Edge",
                                "Input Lesion Type")

  ############################################
  # Gene data interpretation
  gene.clms=colnames(rpt.res[["gene.data"]])
  gene.int=cbind.data.frame(sheet.name="gene.data",
                            col.name=gene.clms,
                            meaning="Echoed from Input")
  rpt.clms=c("gene","chrom","loc.start","loc.end",
             "glp.row.start","glp.row.end")
  int.clms=c(which(gene.int[,"col.name"]=="gene"),
             which(gene.int[,"col.name"]=="chrom"),
             which(gene.int[,"col.name"]=="loc.start"),
             which(gene.int[,"col.name"]=="loc.end"))
  gene.int[int.clms,"meaning"]=c("Input Gene Locus Name",
                                 "Input Gene Locus Chromosome",
                                 "Input Gene Locus Left Edge",
                                 "Input Gene Locus Right Edge")

  ###############################################
  # Chromosome size data interpretation
  chr.clms=colnames(rpt.res[["chr.size"]])
  chr.int=cbind.data.frame(sheet.name="chr.size",
                           col.name=chr.clms,
                           meaning="Echoed from Input")
  int.clms=c(which(chr.int[,"col.name"]=="chrom"),
             which(chr.int[,"col.name"]=="size"))
  chr.int[int.clms,"meaning"]=c("Input Chromosome",
                                "Input Chromosome Size")


  rpt.res$interpretation=rbind.data.frame(sheet.int,
                                          genehit.int,
                                          lsn.int,
                                          gene.int,
                                          chr.int)


  #########################################
  # Write the methods paragraph
  lsn.types=sort(unique(rpt.res[["lsn.data"]][,"lsn.type"]))
  mp=c("The genomic random interval model [ref 1] was used to evaluate ",
       "the statistical significance of the number of subjects with ",
       paste0("each type of lesion (",
              paste0(lsn.types,collapse=","),")"),
       "in each gene.  For each type of lesion, robust false discovery ",
       "estimates were computed from p-values using Storey's q-value [ref 2] ",
       "with the Pounds-Cheng estimator of the proportion of hypothesis ",
       "tests with a true null [ref 3].  Additionally, p-values for the ",
       paste0("number of subjects with any 1 to ",length(lsn.types)," types of lesions "),
       "were computed using the beta distribution derived for order statistics ",
       "of the uniform distribution [ref 4].",
       "",
       "REFERENCES",
       "[ref 1] Pounds S, et al (2013)  A genomic random interval model for statistical analysis of genomic lesion data.  Bioinformatics, 2013 Sep 1;29(17):2088-95 (PMID: 23842812).",
       "[ref 2] Storey J (2002).  A direct approach to false discovery rates.  Journal of the Royal Statistical Society Series B.  64(3): 479-498.  (doi.org/10.1111/1467-9868.00346).",
       "[ref 3] Pounds S and Cheng C (2005)  Robust estimation of the false discovery rate.  Bioinformatics 22(16): 1979-87.  (PMID: 16777905).  ",
       "[ref 4] Casella G and Berger RL (1990)  Statistical Inference.  Wadsworth & Brooks/Cole: Pacific Grove, California.  Example 5.5.1.")

  rpt.res$methods.paragraph=cbind.data.frame(methods.paragraph=mp)

  m=length(rpt.res)
  for (i in 1:m)
  {
    rpt.res[[i]]=as.data.frame(rpt.res[[i]])
  }

  write_xlsx(rpt.res,output.file)

  options(stringsAsFactors=saf)

  return(invisible())

}

############################################################################
# prep.binary.lsn.mtx: The function can be used to prepare a lesion matrix with each gene
# affected by certain lesion type as a row and each patient as a column.

prep.binary.lsn.mtx=function(ov.data,   # result of find.gene.lsn.overlaps function
                             min.ngrp=0)     # omit rows with fewer than min.ngrp hit or fewer than min.ngrp not hit

{
  gene.lsn=ov.data$gene.lsn.hits

  # Order and index data by gene and lesion type
  gene.lsn.type=paste0(gene.lsn$gene,"_",gene.lsn$lsn.type)
  ord=order(gene.lsn$gene,gene.lsn$lsn.type)
  gene.lsn=gene.lsn[ord,]
  m=nrow(gene.lsn)
  new.gene.lsn=which((gene.lsn$gene[-1]!=gene.lsn$gene[-m])|
                       (gene.lsn$lsn.type[-1]!=gene.lsn$lsn.type[-m]))

  row.start=c(1,new.gene.lsn+1)
  row.end=c(new.gene.lsn,m)
  gene.lsn.indx=gene.lsn$gene.lsn[row.start]
  gene.lsn.index=cbind.data.frame(gene=gene.lsn$gene[row.start],
                                  lsn.type=gene.lsn$lsn.type[row.start],
                                  row.start=row.start,
                                  row.end=row.end)
  k=nrow(gene.lsn.index)
  uniq.ID=unique(gene.lsn$ID)
  n=length(uniq.ID)
  gene.lsn.mtx=matrix(0,k,n)
  colnames(gene.lsn.mtx)=as.character(uniq.ID)
  rownames(gene.lsn.mtx)=paste0(gene.lsn.indx$gene,"_",
                                gene.lsn.index$lsn.type)

  for (i in 1:k)
  {
    rows=(gene.lsn.index$row.start[i]:gene.lsn.index$row.end[i])
    ids=as.character(gene.lsn$ID[rows])
    gene.lsn.mtx[i,ids]=1
  }

  n.hit=rowSums(gene.lsn.mtx)
  min.n=pmin(n.hit,n-n.hit)

  keep.row=which(min.n>=min.ngrp)
  gene.lsn.mtx=matrix(gene.lsn.mtx[keep.row,],
                      length(keep.row),n)
  colnames(gene.lsn.mtx)=as.character(uniq.ID)
  rownames(gene.lsn.mtx)=paste0(gene.lsn.index$gene,"_",
                                gene.lsn.index$lsn.type)[keep.row]


  return(gene.lsn.mtx)
}


#############################################################################
# prep.lsn.type.matrix: The function can be used to prepare a lesion matrix with
# each gene as a row and each patient as a column.

prep.lsn.type.matrix=function(ov.data,
                              min.ngrp=0)
{
  gene.lsn=ov.data$gene.lsn.hits

  # order and index gene-lesion overlaps by subject ID and gene
  ord=order(gene.lsn$ID,
            gene.lsn$gene)
  gene.lsn=gene.lsn[ord,]
  m=nrow(gene.lsn)
  new.block=which((gene.lsn$ID[-1]!=gene.lsn$ID[-m])|
                    (gene.lsn$gene[-1]!=gene.lsn$gene[-m]))
  row.start=c(1,new.block+1)
  row.end=c(new.block,m)
  ID.gene.index=cbind.data.frame(ID=gene.lsn$ID[row.start],
                                 gene=gene.lsn$gene[row.start],
                                 row.start=row.start,
                                 row.end=row.end)

  # initialize the result matrix
  lsn.types=sort(unique(gene.lsn$lsn.type))
  uniq.genes=unique(ID.gene.index$gene)
  uniq.IDs=unique(ID.gene.index$ID)

  grp.mtx=matrix("none",
                 length(uniq.genes),
                 length(uniq.IDs))
  rownames(grp.mtx)=uniq.genes
  colnames(grp.mtx)=uniq.IDs

  # fill in the matrix
  n.index=nrow(ID.gene.index)
  for (i in 1:n.index)
  {
    gene.lsn.rows=(ID.gene.index$row.start[i]:ID.gene.index$row.end[i])
    block.ID=ID.gene.index$ID[i]
    block.gene=ID.gene.index$gene[i]
    block.lsns=gene.lsn$lsn.type[gene.lsn.rows]
    block.lsns=unique(block.lsns)
    if (length(block.lsns)==1) grp.mtx[block.gene,block.ID]=block.lsns
    else grp.mtx[block.gene,block.ID]="multiple"
  }

  n.none=rowSums(grp.mtx== "none")
  lsn.grp.keep=(ncol(grp.mtx)-n.none)>=min.ngrp
  grp.mtx=grp.mtx[lsn.grp.keep,]

  return(grp.mtx)
}

#############################################################################
# default.grin.colors: This function specifies default colors for each lesion type
# in GRIN plots.

default.grin.colors=function(lsn.types) # Unique lesion types as specified in the lesion data file
{
  message("Computationally assigning lesion type colors for GRIN plots.")
  uniq.types=sort(unique(lsn.types))
  n.types=length(uniq.types)
  default.colors=c("gold","red","blue",
                   "olivedrab","purple",
                   "cyan", "brown",
                   "orange","steelblue")
  if (length(n.types)>length(default.colors))
    stop(paste0("Too many lesion types for default grin colors; please assign colors manually."))

  res=default.colors[1:n.types]
  names(res)=uniq.types
  return(res)

}

##################################################################################
# compute.gw.coordinates: This function Compute plotting coordinates necessary for
# the genome-wide lesion plot.

compute.gw.coordinates=function(grin.res,  # GRIN results (output of the grin.stats function)
                                scl=1000000) # length of chromosome units in base pairs

{
  # Compute new coordinates for chromosomes
  cum.size=cumsum(grin.res$chr.size$size/scl)
  n.chr=nrow(grin.res$chr.size)
  grin.res$chr.size$x.start=c(0,cum.size[-n.chr])
  grin.res$chr.size$x.end=cum.size

  grin.res$gene.data$x.start=NA
  grin.res$gene.data$x.end=NA
  grin.res$gene.hits$x.start=NA
  grin.res$gene.hits$x.end=NA
  grin.res$lsn.data$x.start=NA
  grin.res$lsn.data$x.end=NA

  gene.data.ord=order(grin.res$gene.data$gene.row)
  grin.res$gene.data=grin.res$gene.data[gene.data.ord,]

  gene.hits.ord=order(grin.res$gene.hits$gene.row)
  grin.res$gene.hits=grin.res$gene.hits[gene.hits.ord,]

  lsn.ord=order(grin.res$lsn.data$lsn.row)
  grin.res$lsn.data=grin.res$lsn.data[lsn.ord,]


  ngi=nrow(grin.res$gene.index)
  for (i in 1:ngi)
  {
    chr.mtch=which(grin.res$gene.index$chrom[i]==grin.res$chr.size$chrom)
    chr.start=grin.res$chr.size$x.start[chr.mtch]

    gene.rows=(grin.res$gene.index$row.start[i]:grin.res$gene.index$row.end[i])
    grin.res$gene.data$x.start[gene.rows]=chr.start+grin.res$gene.data$loc.start[gene.rows]/scl
    grin.res$gene.data$x.end[gene.rows]=chr.start+grin.res$gene.data$loc.end[gene.rows]/scl

    grin.res$gene.hits$x.start[gene.rows]=chr.start+grin.res$gene.hits$loc.start[gene.rows]/scl
    grin.res$gene.hits$x.end[gene.rows]=chr.start+grin.res$gene.hits$loc.end[gene.rows]/scl
  }

  nli=nrow(grin.res$lsn.index)
  for (i in 1:nli)
  {
    chr.mtch=which(grin.res$lsn.index$chrom[i]==grin.res$chr.size$chrom)
    chr.start=grin.res$chr.size$x.start[chr.mtch]

    lsn.rows=(grin.res$lsn.index$row.start[i]:grin.res$lsn.index$row.end[i])
    grin.res$lsn.data$x.start[lsn.rows]=chr.start+grin.res$lsn.data$loc.start[lsn.rows]/scl
    grin.res$lsn.data$x.end[lsn.rows]=chr.start+grin.res$lsn.data$loc.end[lsn.rows]/scl
  }

  return(grin.res) # modified GRIN results to allow adding genome-wide coordinates
}

##################################################################################
# genomewide.lsn.plot:This function create genomewide lesion plot (all chromosomes).

genomewide.lsn.plot=function(grin.res, # GRIN results (output of the grin.stats function)
                             lsn.colors=NULL, # Lesion colors
                             max.log10q=50) # Maximum log10 q value to be added to the plot

{
  if (!is.element("x.start",colnames(grin.res$lsn.data)))
    grin.res=compute.gw.coordinates(grin.res)

  grin.res$lsn.data$x.ID=as.numeric(as.factor(grin.res$lsn.data$ID))
  n=max(grin.res$lsn.data$x.ID)
  n.chr=nrow(grin.res$chr.size)

  # set up plotting region
  plot(c(-0.2,1.2)*n,c(0,-1.1*grin.res$chr.size$x.end[n.chr]),
       type="n",axes=F,xlab="",ylab="")

  # background colors for chromosomes
  rect(0,-grin.res$chr.size$x.start,
       n,-grin.res$chr.size$x.end,
       col=c("lightgray","gray"),
       border=NA)

  rect(-0.075*n,-grin.res$chr.size$x.start,
       -0.2*n,-grin.res$chr.size$x.end,
       col=c("lightgray","gray"),
       border=NA)

  rect(1.075*n,-grin.res$chr.size$x.start,
       1.2*n,-grin.res$chr.size$x.end,
       col=c("lightgray","gray"),
       border=NA)

  lsn.types=sort(unique(grin.res$lsn.index$lsn.type))

  if (is.null(lsn.colors))
  {
    lsn.colors=default.grin.colors(lsn.types)
  }
  grin.res$lsn.data$lsn.colors=lsn.colors[grin.res$lsn.data$lsn.type]

  grin.res$lsn.data$lsn.size=grin.res$lsn.data$x.end-grin.res$lsn.data$x.start
  ord=order(grin.res$lsn.data$lsn.size,decreasing=T)

  rect(grin.res$lsn.data$x.ID-1,
       -grin.res$lsn.data$x.start,
       grin.res$lsn.data$x.ID,
       -grin.res$lsn.data$x.end,
       col=grin.res$lsn.data$lsn.colors,
       border=grin.res$lsn.data$lsn.colors)

  text(c(n,0)[(1:n.chr)%%2+1],
       pos=c(4,2)[(1:n.chr)%%2+1],
       -(grin.res$chr.size$x.start+grin.res$chr.size$x.end)/2,
       grin.res$chr.size$chrom,
       cex=0.5)

  legend(n/2,-1.025*grin.res$chr.size$x.end[n.chr],
         fill=lsn.colors,cex=0.75,
         legend=names(lsn.colors),
         xjust=0.5,
         ncol=length(lsn.colors),
         border=NA,bty="n")

  nsubj.mtx=unlist(grin.res$gene.hits[,paste0("nsubj.",lsn.types)])
  qval.mtx=unlist(grin.res$gene.hits[,paste0("q.nsubj.",lsn.types)])
  nsubj.data=cbind.data.frame(gene=grin.res$gene.hits$gene,
                              x.start=grin.res$gene.hits$x.start,
                              x.end=grin.res$gene.hits$x.end,
                              nsubj=nsubj.mtx,
                              log10q=-log10(qval.mtx),
                              lsn.type=rep(lsn.types,each=nrow(grin.res$gene.hits)))
  nsubj.data=nsubj.data[nsubj.data$nsubj>0,]
  nsubj.data$lsn.colors=lsn.colors[nsubj.data$lsn.type]

  ord=order(nsubj.data$nsubj,decreasing=T)
  nsubj.data=nsubj.data[ord,]

  segments(1.075*n,
           -(nsubj.data$x.start+nsubj.data$x.end)/2,
           1.075*n+0.125*nsubj.data$nsubj/max(nsubj.data$nsubj)*n,
           col=nsubj.data$lsn.colors)

  nsubj.data$log10q[nsubj.data$log10q>max.log10q]=max.log10q

  ord=order(nsubj.data$log10q,decreasing=T)
  nsubj.data=nsubj.data[ord,]

  segments(-0.075*n,
           -(nsubj.data$x.start+nsubj.data$x.end)/2,
           -0.075*n-0.125*nsubj.data$log10q/max(nsubj.data$log10q)*n,
           col=nsubj.data$lsn.colors)

  text(-(0.075+0.20)*n/2,0,
       "-log10(q)",cex=0.75,
       pos=3)

  text((1.075+1.2)*n/2,0,
       "Subjects",cex=0.75,
       pos=3)

  text(c(-0.075,1.075)*n,
       -grin.res$chr.size$x.end[n.chr],
       0,cex=0.75,pos=1)

  text(c(-0.2,1.2)*n,
       -grin.res$chr.size$x.end[n.chr],
       c(max(nsubj.data$log10q),
         max(nsubj.data$nsubj)),
       cex=0.75,pos=1)
}

###################################################################################
# grin.gene.plot: This function plot lesion data and add GRIN results
# (number of subjects with each lesion type, -log10p and -log10q values) for one gene
# (regional gene plot).

grin.gene.plot=function(grin.res, # GRIN results (output of the grin.stats function)
                        gene, # gene name
                        lsn.clrs=NULL, # Specified colors per lesion type
                        expand=0.3) # Controls ratio of the gene locus (start and end position) to the whole plot

{
  # Find the requested gene
  if (length(gene)!=1)
    stop("Exactly one gene must be specified!")

  gene.data=grin.res[["gene.data"]]
  gene.mtch=which(gene.data[,"gene.name"]==gene)
  if (length(gene.mtch)==0)
    stop(paste0(gene," not found in gene.data."))

  if (length(gene.mtch)>1)
    stop(paste0("Multiple matches of ",gene," found in gene data."))

  gene.chr=gene.data[gene.mtch,"chrom"]
  gene.start=gene.data[gene.mtch,"loc.start"]
  gene.end=gene.data[gene.mtch,"loc.end"]
  gene.size=(gene.end-gene.start)+1

  # Find lesions on the same chromosome as the gene
  lsn.dset=order.index.lsn.data(grin.res[["lsn.data"]])
  lsn.data=lsn.dset$lsn.data

  lsn.types=unique(lsn.data$lsn.type)
  if (is.null(lsn.clrs))
    lsn.clrs=default.grin.colors(lsn.types)

  lsn.index=lsn.dset$lsn.index
  lsn.ind.mtch=which(lsn.index$chrom==gene.chr)
  if (length(lsn.ind.mtch)==0)
    stop(paste0("No lesions overlap ",gene,"."))

  lsn.chr.rows=NULL
  for (i in lsn.ind.mtch)
  {
    blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
    lsn.chr.rows=c(lsn.chr.rows,blk.rows)
  }
  lsn.chr.rows=unlist(lsn.chr.rows)

  lsn.chr.data=lsn.data[lsn.chr.rows,]

  if (any(lsn.chr.data$chrom!=gene.chr))
    stop(paste0("Error in finding lesions on same chromosome as gene ",gene,"."))

  # Find lesions at overlap the gene
  ov.rows=which((lsn.chr.data$loc.start<=gene.end)&(lsn.chr.data$loc.end>=gene.start))
  if (length(ov.rows)==0)
    stop(paste0("No lesions overlap gene ",gene,"."))

  lsn.gene=lsn.chr.data[ov.rows,]

  # define plotting data
  x.start=gene.start-expand*gene.size
  x.end=gene.end+expand*gene.size

  lsn.gene$subj.num=as.numeric(as.factor(lsn.gene$ID))
  lsn.gene$lsn.clr=lsn.clrs[lsn.gene$lsn.type]
  lsn.gene$type.num=as.numeric(as.factor(lsn.gene$lsn.type))

  n.type=max(lsn.gene$type.num)
  n.subj=max(lsn.gene$subj.num)

  lsn.gene$y0=-lsn.gene$subj.num+(lsn.gene$type.num-1)/n.type
  lsn.gene$y1=-lsn.gene$subj.num+lsn.gene$type.num/n.type

  lsn.gene$x0=pmax(lsn.gene$loc.start,x.start)
  lsn.gene$x1=pmin(lsn.gene$loc.end,x.end)


  plot(c(x.start-0.20*(x.end-x.start),x.end),
       c(+0.1,-1.3)*n.subj,type="n",
       main="",
       xlab="",
       ylab="",axes=F)

  rect(x.start,
       -(1:n.subj),
       x.end,
       -(1:n.subj)+1,
       col=c("snow","gainsboro")[1+(1:n.subj)%%2],
       border=NA)

  segments(c(gene.start,gene.end),
           rep(-n.subj,2),
           c(gene.start,gene.end),
           rep(0,2),
           col="darkgray")

  rect(lsn.gene$x0,
       lsn.gene$y0,
       lsn.gene$x1,
       lsn.gene$y1,
       col=lsn.gene$lsn.clr,
       border=lsn.gene$lsn.clr)

  segments(c(gene.start,gene.end),
           rep(-n.subj,2),
           c(gene.start,gene.end),
           rep(0,2),
           col="darkgray",lty=2)

  text(x.start,-lsn.gene$subj.num+0.5,
       lsn.gene$ID,pos=2,cex=0.5)
  text(c(gene.start,gene.end),
       -n.subj,
       c(gene.start,gene.end),
       pos=1)
  text((gene.start+gene.end)/2,
       0,gene,pos=3,cex=1.5)

  lgd=legend((x.start+x.end)/2,-1.10*n.subj,
             fill=lsn.clrs,
             legend=names(lsn.clrs),
             ncol=length(lsn.clrs),
             xjust=0.5,border=NA,
             cex=0.75,bty="n")

  text(lgd$text$x[1]-0.05*diff(range(lgd$text$x)),
       -c(1.20,1.25,1.30)*n.subj,
       c("n","-log10p","-log10q"),pos=2)

  gene.stats=grin.res[["gene.hits"]]
  stat.mtch=which(gene.stats$gene.name==gene)
  gene.stats=gene.stats[stat.mtch,]

  text(lgd$text$x,-1.20*n.subj,
       gene.stats[,paste0("nsubj.",names(lsn.clrs))],
       cex=0.75)
  text(lgd$text$x,-1.25*n.subj,
       round(-log10(gene.stats[,paste0("p.nsubj.",names(lsn.clrs))]),2),
       cex=0.75)
  text(lgd$text$x,-1.30*n.subj,
       round(-log10(gene.stats[,paste0("q.nsubj.",names(lsn.clrs))]),2),
       cex=0.75)
}

############################################################################
# top.grin.gene.plots: This function plot lesion data and add GRIN results
# (number of subjects, -log10p and -log10q values) for top genes affected by each
# lesion type in the GRIN analysis (regional gene plots).

top.grin.gene.plots=function(grin.res,                      # GRIN results (output of the grin.stats function)
                             lsn.clrs=NULL,                 # vector of color names
                             q.max=0.10,                    # q-value threshold for plots
                             max.plots.per.type=25,         # maximum number of plots per lesion type
                             max.plots.overall=250)         # overall maximum number of genes with plots
{
  lsn.types=unique(grin.res$lsn.data$lsn.type)
  if (is.null(lsn.clrs))
  {
    lsn.clrs=default.grin.colors(lsn.types)
  }

  gene.stats=grin.res[["gene.hits"]]
  gene.stats=gene.stats[gene.stats$biotype=="protein_coding",]
  gene.stats=gene.stats[!duplicated(gene.stats[ , "gene.name"]),]
  max.plots.per.type=min(max.plots.per.type,
                         nrow(gene.stats))

  top.genes=NULL
  p.clms=c(paste0("p.nsubj.",lsn.types),
           paste0("p",1:length(lsn.types),".nsubj"))
  q.clms=c(paste0("q.nsubj.",lsn.types),
           paste0("q",1:length(lsn.types),".nsubj"))
  for (i in 1:length(p.clms))
  {
    ord=order(gene.stats[,p.clms[i]])
    gene.stats=gene.stats[ord,]
    keep=which(gene.stats[,q.clms[i]]<q.max)
    if (length(keep)>0)
    {
      keep=keep[keep<max.plots.per.type]
      top.genes=c(top.genes,gene.stats[keep,"gene"])
    }
  }
  top.genes=unique(top.genes)
  if (length(top.genes)==0)
    stop("No genes meet selection criteria for plotting.")

  top.gene.rows=which(is.element(gene.stats$gene,top.genes))
  top.gene.stats=gene.stats[top.gene.rows,]
  ord.crit=rowSums(log10(top.gene.stats[,paste0("q.nsubj.",lsn.types)]))
  ord=order(ord.crit)
  top.gene.stats=top.gene.stats[ord,]
  ntop=nrow(top.gene.stats)
  top.gene.stats=top.gene.stats[1:min(ntop,max.plots.overall),]
  top.genes=top.gene.stats$gene.name


  for (i in 1:length(top.genes))
  {
    grin.gene.plot(grin.res,
                   top.genes[i],
                   lsn.clrs=lsn.clrs)
  }
}

#########################################################################
# grin.oncoprint.mtx: This function can be used to prepare the lesion matrix
# that the user can pass to the oncoprint function from complexheat map library
# to geneare an OncoPrint for a selcted list of genes.

grin.oncoprint.mtx=function(grin.res, # list of grin results
                            oncoprint.genes) # vector of ensembl IDs for the selected list of genes

{
  selected=unlist(oncoprint.genes)
  selected=as.vector(selected)
  selected.genes= grin.res$gene.lsn.data[grin.res$gene.lsn.data$gene %in% selected,]
  selected.genes=selected.genes[,c(2,7,11)]  # extract patient IDs and lsn type for each gene in the selected genes list
  row.data=paste(selected.genes[,1],
                 selected.genes[,2],
                 selected.genes[,3],
                 sep="_")
  dup.data=duplicated(row.data)
  select.genes=selected.genes[!dup.data,]

  ord=order(select.genes$gene,
            select.genes$ID,
            select.genes$lsn.type)
  select.genes=select.genes[ord,]

  uniq.genes=unique(select.genes$gene)
  uniq.subj=unique(select.genes$ID)
  n.genes=length(uniq.genes)
  n.subj=length(uniq.subj)
  mtx=matrix("",n.genes,n.subj)   # create a matrix with each gene as a row
  colnames(mtx)=uniq.subj
  rownames(mtx)=uniq.genes

  k=nrow(select.genes)
  for (i in 1:k)
  {
    subj.id=select.genes[i,"ID"]
    gene.id=select.genes[i,"gene"]
    mtx[gene.id,subj.id]=paste0(mtx[gene.id,subj.id],
                                select.genes[i,"lsn.type"],";")
  }
  mtx=as.data.frame(mtx)
  mtx<-tibble::rownames_to_column(mtx, "ensembl.ID")

  gene.annotation= grin.res$gene.data
  ensembl.annotation=cbind(gene.annotation$gene, gene.annotation$gene.name)
  colnames(ensembl.annotation)=c("ensembl.ID", "gene.name")

  mtx.final=merge(ensembl.annotation,mtx,by="ensembl.ID", all.y=TRUE)  # add gene name
  mtx.final=mtx.final[,-1]
  rownames(mtx.final)=mtx.final[,1]
  mtx.final=mtx.final[,-1]

  return(mtx.final)

}

#########################################################################################
# Associate Lesions with Expression (ALEX) library
########################################################################################

####################################################################
# alex.prep.lsn.expr: This function prepares lesion and expression data matrices
# for the KW.hit.express function that runs kruskal-Wallis test (it only keeps genes
# with both lesion and expression data with rows ordered by gene name and columns
# ordered by patient's ID)

alex.prep.lsn.expr=function(expr.mtx, # normalized log2 transformed expression data with genes in rows and subjects in columns (first column is ensembl IDs)
                            lsn.data, # lesion data in GRIN compatible format
                            gene.annotation, # gene annotation data
                            min.pts.expr=NULL, # minimum number of patients with expression level of a gene > 0
                            min.pts.lsn=NULL)  # minimum number of patients with any type of lesions in a certain gene

{
  gene.lsn=prep.gene.lsn.data(lsn.data, gene.annotation)    # Prepare gene and lesion data for later computations
  gene.lsn.overlap= find.gene.lsn.overlaps(gene.lsn)  # Use the results of prep.gene.lsn.data to find lesion-gene overlaps

  message(paste0("Preparing gene-lesion matrix: ",date()))
  lsn.type.matrix=prep.lsn.type.matrix(gene.lsn.overlap) # prepare lesion type matrix (JUST ONE ROW FOR EACH GENE)
  lsn.grp.mtx=as.data.frame(lsn.type.matrix)
  id.hits = colnames(lsn.grp.mtx)

  # prepare expression data for ALEX analysis
  expr.mtx=expr.mtx
  rownames(expr.mtx)=expr.mtx[,1]
  expr.mtx=expr.mtx[,-1]

  id.expr=colnames(expr.mtx)
  lsn.grp.mtx=lsn.grp.mtx[ ,(names(lsn.grp.mtx) %in% id.expr)]

  id.lsn.final=colnames(lsn.grp.mtx)
  expr.mtx=expr.mtx[ ,(names(expr.mtx) %in% id.lsn.final)]

  if (is.numeric(min.pts.expr))
  {
    # To keep only genes with expression value > 0 in at least the number of patients specified in min.pts.expr
    zero.expr=rowSums(expr.mtx==0)
    expr.grp.keep=(ncol(expr.mtx)-zero.expr)>=min.pts.expr
    expr.mtx=expr.mtx[expr.grp.keep,]
  }
  ############################
  # To exclude any gene that has no lesions in at least the number of patients specified in min.pts.lsn
  if (is.numeric(min.pts.lsn))
  {
    n.none=rowSums(lsn.grp.mtx== "none")
    lsn.grp.keep=(ncol(lsn.grp.mtx)-n.none)>=min.pts.lsn
    lsn.grp.mtx=lsn.grp.mtx[lsn.grp.keep,]
  }
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
  check.rownames=all(rownames(lsn.grp.mtx)==rownames(expr.mtx))
  if (check.rownames==FALSE)
    stop("Gene-lesion matrix rownames must match expression matrix rownames.")

  check.colnames=all(colnames(lsn.grp.mtx)==colnames(expr.mtx))
  if (check.colnames==FALSE)
    stop("Gene-lesion matrix patient IDs must match expression matrix patient IDs.")

  # To generate row.mtch file
  lsn.ensembl.id=rownames(lsn.grp.mtx)
  expr.ensembl.id=rownames(expr.mtx)

  row.mtch=cbind.data.frame(expr.row=expr.ensembl.id,
                            hit.row=lsn.ensembl.id)


  res=list(alex.expr=expr.mtx,      # expression data (data table with gene ensembl IDs as row names and each column is a patient)
           alex.lsn=lsn.grp.mtx,    # overlapped gene lesion data table with gene ensembl IDs as row names and each column is a patient
           alex.row.mtch=row.mtch)  # ensembl ID for one row of the alex.lsn table to be paired with one row of the expression data for the Kruskal-Wallis test

  return(res)

}

##############################################################################################
# KW.hit.express: This function uses Kruskal-Wallis test to evaluate association
# between lesion groups and expression level of the same corresponding gene.

KW.hit.express=function(alex.data,          # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data ready for KW test, "alex.lsn" with overlapped gene lesion data and row.mtch)
                        gene.annotation,    # gene annotation data
                        min.grp.size=NULL)  # minimum group size to perform the test (there should be at least two groups with number of patients > min.grp.size)

{
  expr.matrix=as.matrix(alex.data$alex.expr)
  lsn.matrix=as.matrix(alex.data$alex.lsn)
  row.mtch=alex.data$alex.row.mtch

  message(paste0("Computing KW P-value: ",date()))
  p=apply(row.mtch,1,one.KW.pvalue,
          expr.mtx=expr.matrix,
          hit.grps=lsn.matrix,
          min.grp.size=min.grp.size)

  message(paste0("Preparing results table: ",date()))
  kw.res=cbind(row.mtch,p.KW=p)
  # To prepare and print KW test results

  colnames(kw.res)=c("expr.row", "gene", "p.KW")
  kw.res.annotated=merge(gene.annotation,kw.res,by="gene", all.y=TRUE)

  # Compute FDR with the Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat)
  #kw.res.annotated$q.KW=pc06.fdr(kw.res.annotated$p.KW)
  kw.res.annotated$q.KW=p.adjust(kw.res.annotated$p.KW,method="fdr")

  # To compute number of subjects, mean, median and stdev by lesion type
  unique.grps=sort(unique(lsn.matrix))
  unique.grps=sort(unique(unique.grps))
  count.by.lesion <- sapply(unique.grps,function(x)rowSums(lsn.matrix==x)) #count occurrences of x in each row
  colnames(count.by.lesion) = paste(colnames(count.by.lesion),"n.subjects",sep="_")
  count.by.lesion=as.matrix(count.by.lesion)

  mean.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, mean)
  colnames(mean.by.lsn) = paste(colnames(mean.by.lsn),"mean",sep="_")
  mean.by.lsn=as.matrix(mean.by.lsn)

  median.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, median)
  colnames(median.by.lsn) = paste(colnames(median.by.lsn),"median",sep="_")
  median.by.lsn=as.matrix(median.by.lsn)

  sd.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, sd)
  colnames(sd.by.lsn) = paste(colnames(sd.by.lsn),"sd",sep="_")
  sd.by.lsn=as.matrix(sd.by.lsn)

  kw.res.final=cbind(kw.res.annotated, count.by.lesion, mean.by.lsn, median.by.lsn, sd.by.lsn)
  kw.res.final=kw.res.final[!duplicated(kw.res.final[ , "gene.name"]),]
  res=kw.res.final

  return(res)

}

########################################
# one.KW.pvalue: Function Computes the KW p-value for one row of the lesion data matrix
# paired with one row of the expression data matrix

one.KW.pvalue=function(one.row.mtch,    # one row of the row.mtch matrix
                       expr.mtx,        # expression data matrix (alex.expr, output of the alex.prep.lsn.expr function)
                       hit.grps,        # lesion matrix (alex.lsn, output of the alex.prep.lsn.expr function)
                       min.grp.size=NULL) # minimum group size to perform the test (there should be at least two groups with number of patients > min.grp.size)

{
  # extract data for the test
  expr.row=one.row.mtch["expr.row"]
  hit.row=one.row.mtch["hit.row"]
  y=expr.mtx[expr.row,]
  grps=hit.grps[hit.row,]

  # identify data to include in the test
  grps=define.grps(grps,min.grp.size)
  exc=(grps=="EXCLUDE")
  grps=grps[!exc]
  y=y[!exc]

  # return NA if there is no test to perform
  uniq.grps=unique(grps)
  if (length(uniq.grps)<2) return(NA)

  # perform the KW test
  kw.res=kruskal.test(y~grps)

  # return the KW p-value
  return(kw.res$p.value)
}


#######################################################
# define.grps: Function defines the lesion groups to be included in the Kruskal-Wallis test.

define.grps=function(grps,
                     min.grp.size=NULL)
{
  grp.tbl=table(grps)
  small.grp=(grp.tbl<min.grp.size)
  exclude.grp=grps%in%(names(grp.tbl)[small.grp])
  grps[exclude.grp]="EXCLUDE"
  return(grps)
}

###############################################
# stat.by.group: Function to compute stats for row summary (mean, median and sd)

stat.by.group=function(x,              # vector of data
                       g,              # vector of group labels corresponding to x
                       stat,           # stats to be computed
                       all.grps=NULL,  # vector of all the possible group labels
                       ...)            # additional arguments to stat

{
  if (is.null(all.grps))               # if all.grps isn't specified, define it from g
    all.grps=sort(unique(g))

  all.grps=sort(unique(all.grps))      # ensure grps don't appear twice and are in the same order for each application

  res=rep(NA,length(all.grps))         # initialize the result vector
  names(res)=all.grps                  # name the elements of the result vector

  for (i in all.grps)                  # loop over the groups
  {
    in.grp=(g%in%i)
    if (any(in.grp))
      res[i]=stat(x[in.grp],...)
  }

  return(res)

}

###################################################
row.stats.by.group=function(X,
                            G,
                            stat,
                            ...)

{
  if (any(dim(X)!=dim(G)))
    stop("X and G are of incompatible dimensions.  X and G must have the same dimensions.")

  if (any(colnames(X)!=colnames(G)))
    stop("X and G must have column names matched.")

  all.grps=sort(unique(G))
  all.grps=sort(unique(all.grps))

  k=length(all.grps)

  m=nrow(X)

  res=matrix(NA,m,k)
  colnames(res)=all.grps
  rownames(res)=rownames(X)

  for (i in 1:m)
  {
    res[i,]=stat.by.group(X[i,],G[i,],stat,all.grps,...)
  }

  return(res)
}

##################################################################################################
# alex.boxplots: Function return box plots for expression data by lesion groups for selected
# number of genes based on a specified q-value of the kruskal-wallis test

alex.boxplots=function(alex.data,          # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data ready for KW test, "alex.lsn" with overlapped gene lesion data and row.mtch)
                       alex.kw.results,    # ALEX Kruskal-Wallis test results (output of the KW.hit.express function)
                       q,                  # minimum q value for a gene to be included in output PDF file of box plots
                       gene.annotation)    # gene annotation data

{
  # To extract genes below specified q.value for boxplots:
  KW.results=alex.kw.results[!is.na(alex.kw.results$q.KW),]
  selected.genes=KW.results[KW.results$q.KW<q,]

  selected.IDS=selected.genes$gene
  selected.lsns=alex.data$alex.lsn[rownames(alex.data$alex.lsn) %in% selected.IDS,]
  selected.lsns=t(selected.lsns)
  selected.expr=alex.data$alex.expr[rownames(alex.data$alex.expr) %in% selected.IDS,]
  selected.expr=t(selected.expr)

  all(rownames(selected.lsns)==rownames(selected.expr))
  all(colnames(selected.lsns)==colnames(selected.expr))

  gene.annotation=gene.annotation %>% mutate_all(na_if,"")
  gene.annotation$gene.name <- ifelse(is.na(gene.annotation$gene.name),gene.annotation$gene, gene.annotation$gene.name) # add ensembl.ID to gene.name column if empty

  n=ncol(selected.lsns)

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
        theme(legend.text=element_text(size=8))
      print(p)
    }
  }
}


###########################################
# alex.waterfall.prep: Function prepares lesion and expression data of a selected gene for the
# alex.waterfall.plot function

alex.waterfall.prep=function(alex.data,             # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data ready for KW test, "alex.lsn" with overlapped gene lesion data and row.mtch)
                             alex.kw.results,       # ALEX Kruskal-Wallis test results (output of the KW.hit.express function)
                             gene,                  # gene name or ensembl.ID
                             lsn.data)              # lesion data in a GRIN compatible format


{
  # Find the row of the integrated results with results for the gene

  int.row=which((alex.kw.results$gene==gene)|
                  (alex.kw.results$gene.name==gene))

  if (length(int.row)!=1) stop("gene must match exactly ONE row of alex.kw.results")

  ens.ID=alex.kw.results$gene[int.row]
  gene.ID=alex.kw.results$gene.name[int.row]

  gene.lsn.exp=as.data.frame(colnames(alex.data$alex.lsn))
  colnames(gene.lsn.exp)="ID"
  # Add the lesion data for this gene to the gene.lsn.expical data
  rownames(gene.lsn.exp)=gene.lsn.exp$ID
  gene.lsn=paste0(gene.ID,".lsn")
  lsn.row=which(rownames(alex.data$alex.lsn)==ens.ID)
  gene.lsn.exp[,gene.lsn]=""
  lsn.ids=intersect(gene.lsn.exp$ID,colnames(alex.data$alex.lsn))
  gene.lsn.exp[lsn.ids,gene.lsn]=unlist(alex.data$alex.lsn[lsn.row,lsn.ids])


  # Add the expression data for this gene to the gene.lsn.exp data
  gene.rna=paste0(gene.ID,".RNA")
  gene.lsn.exp[,gene.rna]=NA
  rna.ids=intersect(gene.lsn.exp$ID,colnames(alex.data$alex.expr))
  rna.row=which(rownames(alex.data$alex.expr)==ens.ID)
  gene.lsn.exp[rna.ids,gene.rna]=unlist(alex.data$alex.expr[rna.row,rna.ids])

  # Find lesions that overlap the gene
  lsn.mtch=(lsn.data$chrom==alex.kw.results$chrom[int.row])&
    (lsn.data$loc.end>=alex.kw.results$loc.start[int.row])&
    (lsn.data$loc.start<=alex.kw.results$loc.end[int.row])

  lsns=lsn.data[lsn.mtch,]

  res=list(gene.lsn.exp=gene.lsn.exp,
           lsns=lsns,
           stats=alex.kw.results[int.row,],
           gene.ID=gene.ID)
}

###################################################################
# alex.waterfall.plot: Function return a waterfall plot for expression data by lesion groups of
# a selected gene.

alex.waterfall.plot=function(waterfall.prep,
                             lsn.data,         # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data ready for KW test, "alex.lsn" with overlapped gene lesion data and row.mtch)
                             lsn.clrs=NULL,
                             ex.clrs=rgb(c(0,0,1,1),
                                     rep(0,4),
                                     c(1,1,0,0),
                                     alpha=c(1,0.5,0.5,1)),
                                     delta=0.5)

{
  unique.grps=sort(unique(lsn.data$lsn.type))
  if (is.null(lsn.clrs))
  {
    lsn.grps.clr=default.grin.colors(unique.grps)
    common.grps.clr=c(none="gray", multiple="black")
    lsn.clrs=c(lsn.grps.clr, common.grps.clr)
  }

  gene.ID=waterfall.prep$gene.ID
  gene.lsn.exp=waterfall.prep$gene.lsn.exp

  lsn.clm=paste0(gene.ID,".lsn")
  rna.clm=paste0(gene.ID,".RNA")

  gene.lsn.exp.ord=order(gene.lsn.exp[,lsn.clm],
                         gene.lsn.exp[,rna.clm])

  gene.lsn.exp=gene.lsn.exp[gene.lsn.exp.ord,]

  pt.num=1:nrow(gene.lsn.exp)
  names(pt.num)=gene.lsn.exp$ID


  ####################################
  # Set up plotting region
  plot(c(-1.1,+1.5),
       c(0.1,-1.1)*nrow(gene.lsn.exp),
       type="n",axes=F,
       xlab="",ylab="")

  ###################################
  # DNA lesion plot has x in -1.1 to -0.1

  # gene locus
  loc.rng=c(waterfall.prep$stats$loc.start,
            waterfall.prep$stats$loc.end)
  loc.lng=diff(loc.rng)

  pos.rng=loc.rng+c(-1,1)*delta*loc.lng
  waterfall.prep$lsns$x.start=(waterfall.prep$lsns$loc.start-pos.rng[1])/(diff(pos.rng))-1.1
  waterfall.prep$lsns$x.end=(waterfall.prep$lsns$loc.end-pos.rng[1])/(diff(pos.rng))-1.1

  x.locus=(loc.rng-pos.rng[1])/diff(pos.rng)-1.1

  # background for lesion plot
  rect(-1.1,-nrow(gene.lsn.exp),
       -0.1,0,
       col=lsn.clrs["none"],
       border=lsn.clrs["none"])

  waterfall.prep$lsns$size=waterfall.prep$lsns$loc.end-
    waterfall.prep$lsns$loc.start+1
  ord=rev(order(waterfall.prep$lsns$size))
  waterfall.prep$lsns=waterfall.prep$lsns[ord,]
  rect(pmax(waterfall.prep$lsns$x.start,-1.1),
       -pt.num[waterfall.prep$lsns$ID],
       pmin(waterfall.prep$lsns$x.end,-0.1),
       -pt.num[waterfall.prep$lsns$ID]+1,
       col=lsn.clrs[waterfall.prep$lsns$lsn.type],
       border=lsn.clrs[waterfall.prep$lsns$lsn.type])

  segments(x.locus,-nrow(gene.lsn.exp),
           x.locus,0,col="white",lty=3)

  text(x.locus,-1.05*nrow(gene.lsn.exp),
       loc.rng,cex=0.75)
  text(-0.55,+0.05*nrow(gene.lsn.exp),
       paste0(gene.ID," DNA Lesions"))

  ##############################################
  # RNA expression plot has x in +0.1 to +1.1


  rna.lbls=pretty(gene.lsn.exp[,rna.clm])
  ok.lbls=(rna.lbls>min(gene.lsn.exp[,rna.clm],na.rm=T))&
    (rna.lbls<max(gene.lsn.exp[,rna.clm],na.rm=T))
  rna.lbls=rna.lbls[ok.lbls]

  gene.lsn.exp$x.rna=(gene.lsn.exp[,rna.clm]-min(gene.lsn.exp[,rna.clm],na.rm=T))/diff(range(gene.lsn.exp[,rna.clm],na.rm=T))
  gene.lsn.exp$x.rna=gene.lsn.exp$x.rna+0.1

  x.rna.lbls=(rna.lbls-min(gene.lsn.exp[,rna.clm],na.rm=T))/diff(range(gene.lsn.exp[,rna.clm],na.rm=T))+0.1

  n.mdn=function(x)
  {
    mdn=median(x,na.rm=T)
    n=sum(!is.na(x))
    res=unlist(c(n=n,mdn=mdn))
    return(res)
  }

  lsn.mdn=aggregate(x=gene.lsn.exp$x.rna,
                    by=list(lsn=gene.lsn.exp[,lsn.clm]),
                    FUN=n.mdn)
  rownames(lsn.mdn$x)=lsn.mdn$lsn

  x.mdn=lsn.mdn$x[gene.lsn.exp[,lsn.clm],"mdn"]

  text(x.rna.lbls,-1.05*nrow(gene.lsn.exp),
       rna.lbls,cex=0.75)
  segments(x.rna.lbls,0,
           x.rna.lbls,-nrow(gene.lsn.exp),
           lty=3)

  rect(pmin(gene.lsn.exp$x.rna,x.mdn),-(1:nrow(gene.lsn.exp)),
       pmax(gene.lsn.exp$x.rna,x.mdn),-(1:nrow(gene.lsn.exp))+1,
       col=lsn.clrs[gene.lsn.exp[,lsn.clm]],
       border=NA)

  segments(x.mdn,-(1:nrow(gene.lsn.exp)),
           x.mdn,-(1:nrow(gene.lsn.exp))+1,
           col=lsn.clrs[gene.lsn.exp[,lsn.clm]])

  text(0.55,0.05*nrow(gene.lsn.exp),
       paste0(gene.ID," RNA Expression"))

  #############################
  # Add legend

  lsn.inc=names(lsn.clrs)%in%gene.lsn.exp[,lsn.clm]
  legend(1.15,-0.25*nrow(gene.lsn.exp),
         fill=unlist(lsn.clrs[lsn.inc]),
         legend=names(lsn.clrs[lsn.inc]),
         cex=0.75,border=NA,bty="n")

}

#############################################################
# top.alex.waterfall.plots: Function can be used to generate waterfall plots
# for top significant genes in the KW results table:

top.alex.waterfall.plots=function(out.dir,            # path to the folder to add waterfall plots
                                  alex.data,          # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data ready for KW test, "alex.lsn" with overlapped gene lesion data and row.mtch)
                                  alex.kw.results,    # ALEX Kruskal-Wallis test results (output of the KW.hit.express function)
                                  q,                  # minimum q value for a gene to be included in output PDF f box plots
                                  lsn.data)           # gene annotation data


{
  # To extract genes below specified q.value for waterfall plots:
  KW.results=alex.kw.results[!is.na(alex.kw.results$q.KW),]
  selected.genes=KW.results[KW.results$q.KW<q,]
  selected.genes=selected.genes[!is.na(selected.genes$gene.name),]
  selected.genes=selected.genes[!duplicated(selected.genes[ , "gene.name"]),]
  top.genes=selected.genes$gene.name

  for (i in 1:length(top.genes))
  {
    pdf(paste0(out.dir,top.genes[i],".pdf"),
        height=8,width=10)
    temp=alex.waterfall.prep(alex.data, alex.kw.results,top.genes[i],lsn.data)
    plot.gene=alex.waterfall.plot(temp, lsn.data)
    dev.off()
  }
}

##################################################################################
#

# compute distance for lesions
dist.lsn=function(lsn.mtx) # subjects in columns, genes in rows
{
  n=ncol(lsn.mtx)
  res=matrix(NA,n,n)
  colnames(res)=rownames(res)=colnames(lsn.mtx)
  diag(res)=0
  for (i in 1:(n-1))
    for (j in (i+1):n)
    {
      res[i,j]=res[j,i]=sum(lsn.mtx[,i]!=lsn.mtx[,j])
    }
  return(res)
}



##############################################################
# alex.pathway: Function separately compute the distance between subjects in th dataset based on lesion and
# expression data of genes specified by the user to be part of the pathway of interst and return
# two panels of clusterd subjects based on the computed distances.
# Function also return a matrix of lesion and expression data for each subject ordered based
# on the hiearchial clustering analysis (same order of lesion and expression panels of the figure)

alex.pathway=function(alex.data,          # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data ready for KW test, "alex.lsn" with overlapped gene lesion data and row.mtch)
                      lsn.data,           # lesion data in a GRIN compatible format
                      pathways,           # data.table with three columns "gene.name" with gene symbols, "ensembl.id" with gene ensembl ID and "pathway" that has the pathway name
                      selected.pathway)   # pathway of interest

{
  expr=alex.data$alex.expr
  lsn=alex.data$alex.lsn
  selected.genes=pathways[pathways$pathway==selected.pathway,]
  pathway.ensembl=selected.genes$ensembl.id

  path.expr=expr[(rownames(expr)%in%pathway.ensembl),]
  path.lsn=lsn[(rownames(lsn)%in%pathway.ensembl),]

  path.expr=as.matrix(path.expr)
  path.lsn=as.matrix(path.lsn)

  unique.grps=sort(unique(lsn.data$lsn.type))
  lsn.grps.clr=default.grin.colors(unique.grps)
  common.grps.clr=c(none="gray", multiple="black")
  lsn.clrs=c(lsn.grps.clr, common.grps.clr)

  lsn.grps=names(lsn.clrs)
  clrs=as.character(lsn.clrs)

  path.lsn.clr=path.lsn
  # to replace lesion groups with colors in the lsn matrix
  path.lsn.clr[path.lsn.clr %in% lsn.grps] <- clrs[match(path.lsn.clr, lsn.grps, nomatch = 0)]

  path.genes=selected.genes$ensembl.id
  path.gene.names=selected.genes$gene.name
  names(path.gene.names)=path.genes

  # to change ensembl IDs with gene name
  rownames(path.expr)=path.gene.names[rownames(path.expr)]
  rownames(path.lsn)=path.gene.names[rownames(path.lsn)]

  dist.lsn.genes=dist.lsn(t(path.lsn))
  hcl.expr.genes=hclust(dist(path.expr),"complete")
  hcl.lsn.genes=hclust(as.dist(dist.lsn.genes),"complete")

  # Explore clustering
  path.ex.dist=dist(t(path.expr),"euclidean")
  path.lsn.dist=dist.lsn(path.lsn)
  path.ex.dist=as.matrix(path.ex.dist)
  path.lsn.dist=as.matrix(path.lsn.dist)

  path.ex.hcl=hclust(as.dist(path.ex.dist),method="ward.D2")
  path.ls.hcl=hclust(as.dist(path.lsn.dist),method="ward.D2")
  sum.none=rowSums(path.lsn.clr=="gray")
  subj.ord=path.ls.hcl$order
  subj.labels=path.ls.hcl$labels

  #####################################
  # side-by-side heatmap
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(c(0,+1.1*ncol(path.lsn.clr)),
       c(0,-(nrow(path.lsn.clr)+3+nrow(path.expr))),
       type="n",xlab="",ylab="",
       axes=F)

  for (i in 1:nrow(path.lsn.clr))
    rect(0:(ncol(path.lsn.clr)-1),-(i-1),
         1:ncol(path.lsn.clr),-i,
         col=path.lsn.clr[hcl.lsn.genes$order[i],subj.ord],
         border=NA)
  text(ncol(path.lsn),-(1:nrow(path.lsn))+0.5,
       rownames(path.expr)[hcl.lsn.genes$order],
       pos=4,cex=0.75)

  # Add legend
  lsn.inc=names(lsn.clrs)%in%path.lsn
  legend("topright", inset=c(-0.25,0.01),
         fill=unlist(lsn.clrs[lsn.inc]),
         legend=names(lsn.clrs[lsn.inc]),
         cex=0.72,border=NA,bty="n")

  for (i in 1:nrow(path.expr))
  {
    y=path.expr[hcl.lsn.genes$order[i],]
    z=(y-mean(y))/sd(y)
    clr=rgb((z>0),0,(z<0),alpha=sqrt(1-exp(-abs(z))))
    rect(0:(ncol(path.lsn.clr)-1),-nrow(path.lsn.clr)-3-(i-1),
         1:ncol(path.lsn.clr),-nrow(path.lsn.clr)-3-i,
         col=clr[subj.ord],
         border=NA)
  }

  legend("bottomright", inset=c(-0.25,0.2),
         legend=c("Z.expr<0",0,"Z.expr>0"),
         fill = c("blue", "white", "red"),
         cex=0.72,border=NA,bty="n")

  text(ncol(path.expr),-nrow(path.lsn.clr)-3-1:nrow(path.expr)+0.5,
       rownames(path.expr)[hcl.lsn.genes$order],cex=0.75,pos=4)

  pts.labels=as.data.frame(subj.labels)
  pts.labels<-tibble::rownames_to_column(pts.labels, "subj.ord")
  pts.order=as.data.frame(subj.ord)
  pts.order$index=1:nrow(pts.order)
  ordered.subj.final=merge(pts.labels,pts.order,by="subj.ord", all.y=TRUE)
  ordered.subj.final=ordered.subj.final[order(ordered.subj.final$index),]
  path.lsn.df=as.data.frame(path.lsn)
  ordered.lsn.data=path.lsn.df[ordered.subj.final$subj.labels]
  rownames(ordered.lsn.data) = paste(rownames(ordered.lsn.data),"_lsn")
  path.expr.df=as.data.frame(path.expr)
  ordered.expr.data=path.expr.df[ordered.subj.final$subj.labels]
  rownames(ordered.expr.data) = paste(rownames(ordered.expr.data),"_expr")
  ordered.path.data=rbind(ordered.lsn.data, ordered.expr.data)

  return(ordered.path.data)

}

###########################################################################
# Distance matrix covariance

dist.cov=function(dist.mtx1,
                  dist.mtx2,
                  pre.cent1=F,
                  pre.cent2=F)

{
  n=unique(c(dim(dist.mtx1),dim(dist.mtx2)))
  if (length(n)>1)
    stop("dist.mtx must be a square matrix.")

  if (!pre.cent1)  dist.mtx1=U.center(dist.mtx1)
  if (!pre.cent2)  dist.mtx2=U.center(dist.mtx2)

  prd=dist.mtx1*dist.mtx2
  res=(sum(prd)-sum(diag(prd)))/(n*(n-3))
  return(res)
}


###########################################################
# alex.gsda: Function run gene-set distance analysis (GSDA) for association between lesion and expression data in
# different gene-sets obtained from Molecular Signature Database (MSigDB).

alex.gsda=function(alex.data,          # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data ready for KW test, "alex.lsn" with overlapped gene lesion data and row.mtch)
                   gene.data,          # gene annotation data
                   big.gset=250,         # largest number of genes per gene-set
                   small.gset=2,       # smallest number of genes per gene-set
                   mess.freq=100)      # print a progress message every mess.freq gene-sets

{
  expr=as.matrix(alex.data$alex.expr)
  lsn=as.matrix(alex.data$alex.lsn)
  ensembl.annotation=gene.data

  # obtain gene-set data
  message(paste0("Obtaining gene-set data from Molecular Signature Database (MSigDB): ",date()))
  msigdb.data=msigdbr(species = "Homo sapiens")
  msigdb.data=as.data.frame(msigdb.data)

  #######################################
  # Prepare GSDA vset data

  vset.data=cbind.data.frame(vID=msigdb.data$ensembl_gene,
                             Vgene.name=msigdb.data$human_gene_symbol,
                             vset=msigdb.data$gs_name)


  # Keep only gene-sets with reasonable number of genes that range between numbers specified by the user
  keep.vset=(vset.data$vID%in%rownames(expr))&(vset.data$vID%in%rownames(lsn))
  vset.data=vset.data[keep.vset,]
  vset.tbl=table(vset.data$vset)
  big.vset=which(vset.tbl>big.gset)
  small.vset=which(vset.tbl<small.gset)
  drop.vset=(vset.data$vset%in%(names(vset.tbl)[c(big.vset,small.vset)]))
  vset.data=vset.data[!drop.vset,]
  ord=order(vset.data$vset)
  vset.data=vset.data[ord,]
  gset.data=vset.data

  # order and index gset.data
  message(paste0("Ordering and indexing gene-set data: ",date()))
  ord.gset=order(gset.data$vset,
                 gset.data$vID)
  gset.data=gset.data[ord.gset,]
  m=nrow(gset.data)
  new.gset=which(gset.data$vset[-1]!=gset.data$vset[-m])
  row.start=c(1,new.gset+1)
  row.end=c(new.gset,m)
  gset.index=cbind.data.frame(gset=gset.data$vset[row.start],
                              row.start=row.start,
                              row.end=row.end)
  n.gset=nrow(gset.index)

  # Prepare lesion and expression data matrices
  lsn=t(lsn) # transpose lesion data
  expr=t(expr)   # transpose expression data

  # Compute v for degrees of freedom
  n=nrow(lsn)
  v=n*(n-3)/2-1
  # Compute distance correlation p-values for each gene set (lsn and expression matrices for each gene-set)
  p.gset=rep(NA,n.gset)
  gset.genes=rep("",n.gset)
  gset.gene.names=rep("",n.gset)
  gset.count=rep("",n.gset)
  r.gset=rep(NA,n.gset)
  t.stat.gset=rep(NA,n.gset)

  for (i in 1:n.gset)
  {
    if (i%%mess.freq==1)
      message(paste("Computing distance correlation for gene-set ",
                    i," of ",n.gset,": ",date()))
    gset.rows=gset.index$row.start[i]:gset.index$row.end[i]
    gene.clms=unique(gset.data$vID[gset.rows])
    genes.df=as.data.frame(gene.clms)
    colnames(genes.df) <- ("gene")
    genes.df.merged=merge(ensembl.annotation,genes.df,by="gene", all.y=TRUE)
    gene.names=as.character(genes.df.merged$gene.name)
    gset.genes[i]=paste(gene.clms,collapse=" | ")
    gset.gene.names[i]=paste(gene.names,collapse=" | ")
    gset.count[i]=length(gene.clms)

    lsn.geneset=(lsn[,gene.clms])
    lsn.geneset=t(lsn.geneset) # transpose lesion matrix so that dist.lsn function can compute distance between subjects not genes
    lsn.dist=dist.lsn(lsn.geneset)
    lsn.dist=as.matrix(lsn.dist)
    lsn.dist=U.center(lsn.dist)
    lsn.dcov=dist.cov(lsn.dist,lsn.dist,T,T)

    Y=expr[,gene.clms]
    Y.dist=dist(Y,"euclidean") # Compute euclidean distance for expression data for all genes in a specific gene-set
    Y.dist=as.matrix(Y.dist)
    Y.dist=U.center(Y.dist)
    expr.dcov=dist.cov(Y.dist,Y.dist,T,T)

    dcov=dist.cov(Y.dist,lsn.dist,T,T)
    dcor=dcov/sqrt(lsn.dcov*expr.dcov)
    r.gset[i]=dcor
    dt=sqrt(v)*(dcor/sqrt(1-dcor*dcor))
    t.stat.gset[i]=dt
    p.gset[i]=2*pt(-abs(dt),v)
  }

  q=p.adjust(p.gset,"fdr")
  #q=q*min(1,2*mean(p.gset), na.rm = TRUE)

  res=cbind.data.frame(gset=gset.index$gset,
                       r=r.gset,
                       t.stat=t.stat.gset,
                       p=p.gset,
                       q=q,
                       n.genes=gset.count,
                       genes=gset.genes,
                       gene.name=gset.gene.names)

  ord=order(res$p)

  expr.data=t(expr)
  expr.data=as.data.frame(expr.data)
  expr.data=cbind(id=rownames(expr.data),
                  expr.data)
  gset.data=as.data.frame(gset.data)
  lsn.mtx=t(lsn)
  lsn.mtx=as.data.frame(lsn.mtx)
  lsn.mtx=cbind.data.frame(gene.lsn.type=rownames(lsn.mtx),
                           lsn.mtx)

  opts=cbind.data.frame(
                        min.gset.size=small.gset,
                        max.gset.size=big.gset,
                        message.frequency=mess.freq)

  methods.paragraph=c("Lesion data were compared to gene location data to generate",
                      "a binary (0 or 1) matrix indicating whether each subject had a lesion of",
                      "each type in each gene.  This matrix had one row for each type of lesion ",
                      "and each row of gene location data.  Analysis was then limited to ",
                      paste0("for each of the ",nrow(gset.index)," gene-sets with between ",paste0(big.gset),
                      paste0("and ",small.gset)," genes appearing in the expression matrix, "),
                      "a distance correlation test (Cao and Pounds 2021) was used to evaluate the association of ",
                      "the Euclidean distance between subjects' expression profiles with ",
                      "the Manhattan distance between subjects' lesion matrix profiles.  To address multiple-testing, we computed ",
                      "q-values (Storey 2002) using a moment-based estimate of the proportion of tests with a true null (Pounds and Cheng 2006).",
                      "",
                      "REFERENCES","",
                      "Cao X and Pounds S (2021) Gene-set distance analysis (GSDA): a powerful tool for gene-set association analysis.",
                      "                          BMC Bioinformatics. 2021 Apr 21;22(1):207. doi: 10.1186/s12859-021-04110-x. PMID: 33882829 ",
                      "",
                      "Pounds S and Cheng C (2006)  Robust estimation of the false discovery rate.",
                      "                             Bioinformatics. 2006 Aug 15;22(16):1979-87.  doi: 10.1093/bioinformatics/btl328.  PMID: 16777905",
                      "",
                      "Storey J (2002)  A direct approach to false discovery rates.",
                      "                 Journal of the Royal Statistical Society Series B.  64(3): 479-498.  https://doi.org/10.1111/1467-9868.00346")


  methods.paragraph=as.data.frame(methods.paragraph)

  full.res=list(lsn.gsda=res,
                expr.data=expr.data,
                gene.data=ensembl.annotation,
                lsn.mtx=lsn.mtx,
                opts=opts,
                methods.paragraph=methods.paragraph)

  return(full.res)

}


##################################################################################
# grin.assoc.lsn.outome: function first version
# run association analysis between the binary lesion matrix (output of prep.binary.lsn.mtx function)
# and clinical outcomes including MRD, EFS and OS

grin.assoc.lsn.outcome=function(lsn.mtx,            # output of prep.binary.lsn.mtx function
                               clin.data,          # clinical data table. First column denoted as "ID" and each row is a patient
                               annotation.data,    # gene annotation data file
                               mrd=NULL,           # MRD data denoted as "MRD.binary"
                               efs=NULL,           # user should provide two columns denoted as "efs.time" and "efs.censor"
                               efs.censor=NULL,
                               os=NULL,            # user should provide two columns denoted as "os.time" and "os.censor"
                               os.censor=NULL,
                               covariate=NULL)    # covariates that the model will adjust for

{
  lsn.mtx=t(lsn.mtx)
  lsn.df=as.data.frame(lsn.mtx)

  clin.data=clin.data
  clin.ids=clin.data$ID

  lsn.df=lsn.df[(rownames(lsn.df) %in% clin.ids),]
  lsn.ids=rownames(lsn.df)
  clin.data=clin.data[clin.data$ID %in% lsn.ids,]
  lsn.df=lsn.df[order(rownames(lsn.df)),]
  clin.data=clin.data[order(clin.data$ID),]
  # To make sure that both lesion and expression matrices have same patients and genes order
  check.rownames=all(rownames(lsn.df)==clin.data$ID)
  if (check.rownames==FALSE)
    stop("Gene-lesion matrix rownames must match patient IDs in the clinical data.")

  merged.data=cbind(lsn.df, clin.data)

  # Run logistic regression models for association with selected variables
  if (is.character(mrd))
    mrd.data=merged.data[!is.na(merged.data$MRD.binary),]

  message(paste0("Running logistic regression models for specified variable: ",date()))

  mrd.lsn.clms=mrd.data[ ,grepl("ENSG", names(mrd.data))]
  mrd=mrd.data$MRD.binary
  mrd.data.final=cbind(mrd,mrd.lsn.clms)
  mrd.positive= mrd.data.final[mrd.data.final$mrd==1,]
  mrd.negative= mrd.data.final[mrd.data.final$mrd==0,]

  mrd.positive=mrd.positive[,-1]
  mrd.negative=mrd.negative[,-1]

  mrd.positive.with.lsn=as.data.frame(colSums(mrd.positive== 1))
  mrd.positive.without.lsn=as.data.frame(colSums(mrd.positive== 0))
  mrd.negative.with.lsn=as.data.frame(colSums(mrd.negative== 1))
  mrd.negative.without.lsn=as.data.frame(colSums(mrd.negative== 0))

  mrd.count=cbind(mrd.positive.with.lsn, mrd.positive.without.lsn,
                  mrd.negative.with.lsn, mrd.negative.without.lsn)

  colnames(mrd.count)=c("mrd.positive.with.lsn", "mrd.positive.without.lsn",
                        "mrd.negative.with.lsn", "mrd.negative.without.lsn")


  if (is.null(covariate)) {
    log.fit<-lapply(mrd.lsn.clms,function(x) glm(MRD.binary  ~ x ,data = mrd.data, family = "binomial"))
    ## To extract coefficients
    coeff = lapply(log.fit,function(f) summary(f)$coefficients[,1])

    ## Number of columns should be 4 if we are adjusting for risk groups but for univariate models it should be 2 (column for intercept and column for model coefficient):
    coeff.mat <- do.call(rbind,lapply(coeff,matrix,ncol=2,byrow=TRUE))
    mod_coeff = lapply(coeff, head, 2)
    #write.csv(mod_coeff1,file ="NOGO_MRD1_COEFF.csv") ## Write name of gene or SNP in addition to intercept and variable coefficient to a CSV file

    ## To extract model error:
    error = lapply(log.fit,function(f) summary(f)$coefficients[,2])
    error.mat <- do.call(rbind,lapply(error,matrix,ncol=2,byrow=TRUE))

    ## To calculate odds ratio and 95% CI:
    ## N.B) Here, we are using the standard errors not the robust SE
    odds.ratio= exp(coeff.mat[,2])
    odds.mat=as.matrix(odds.ratio)
    lower.confint = exp(coeff.mat[,2] - 1.96*error.mat[,2])
    lower.mat=as.matrix(lower.confint)
    higher.confint = exp(coeff.mat[,2] + 1.96*error.mat[,2])
    higher.mat=as.matrix(higher.confint)

    ## To extract P-value:
    Pvalue<-lapply(log.fit,function(f) summary(f)$coefficients[,4])
    pvalue.mat=do.call(rbind,lapply(Pvalue,matrix,ncol=2,byrow=TRUE))
    pvalue.mrd= pvalue.mat[,2]
    pi.hat=min(1,2*mean(pvalue.mrd,na.rm=T))
    q.mrd=pi.hat*p.adjust(pvalue.mrd,method="fdr")
    results.logistic=cbind(odds.mat, lower.mat, higher.mat, pvalue.mrd, q.mrd)
    logistic.lsn=as.data.frame(colnames(lsn.mtx))
    colnames(logistic.lsn)="Gene_lsn"
    results.logistic.final=cbind(logistic.lsn, results.logistic)
    colnames(results.logistic.final)=c("Gene_lsn", "MRD.odds.ratio", "MRD.lower95", "MRD.upper95", "MRD.p-value", "MRD.q-value")
    MRD.results.count.added=cbind(results.logistic.final, mrd.count)

  } else {

    # MRD with covariates
    covariates=covariate
    covariates=mrd.data[,covariates]
    log.fit.adj<-lapply(mrd.lsn.clms,function(x) glm(data = mrd.data, MRD.binary  ~ x + covariates , family = "binomial"))
    ## To extract coefficients
    coeff.adj = lapply(log.fit.adj,function(f) summary(f)$coefficients[,1])

    ## Number of columns should be 4 if we are adjusting for risk groups but for univariate models it should be 2 (column for intercept and column for model coefficient):
    coeff.mat.adj <- do.call(rbind,lapply(coeff.adj,matrix,byrow=TRUE))
    num.grps=length(coeff.mat.adj)/ncol(lsn.mtx)
    coeff.mat.adj <- do.call(rbind,lapply(coeff.adj,matrix,ncol=num.grps,byrow=TRUE))
    #mod_coeff = lapply(coeff, head, 2)
    #write.csv(mod_coeff1,file ="NOGO_MRD1_COEFF.csv") ## Write name of gene or SNP in addition to intercept and variable coefficient to a CSV file

    ## To extract model error:
    error.adj = lapply(log.fit.adj,function(f) summary(f)$coefficients[,2])
    error.mat.adj <- do.call(rbind,lapply(error.adj,matrix,ncol=num.grps,byrow=TRUE))

    ## To calculate odds ratio and 95% CI:
    ## N.B) Here, we are using the standard errors not the robust SE
    odds.ratio.adj= exp(coeff.mat.adj[,2])
    odds.mat.adj=as.matrix(odds.ratio.adj)
    lower.confint.adj = exp(coeff.mat.adj[,2] - 1.96*error.mat.adj[,2])
    lower.mat.adj=as.matrix(lower.confint.adj)
    higher.confint.adj = exp(coeff.mat.adj[,2] + 1.96*error.mat.adj[,2])
    higher.mat.adj=as.matrix(higher.confint.adj)

    ## To extract P-value:
    Pvalue.adj<-lapply(log.fit.adj,function(f) summary(f)$coefficients[,4])
    pvalue.mat.adj=do.call(rbind,lapply(Pvalue.adj,matrix,ncol=num.grps,byrow=TRUE))
    pvalue.mrd.adj= pvalue.mat.adj[,2]
    pi.hat=min(1,2*mean(pvalue.mrd.adj,na.rm=T))
    q.mrd.adj=pi.hat*p.adjust(pvalue.mrd.adj,method="fdr")

    results.logistic.adj=cbind(odds.mat.adj, lower.mat.adj, higher.mat.adj, pvalue.mrd.adj, q.mrd.adj)

    logistic.lsn=as.data.frame(colnames(lsn.mtx))
    colnames(logistic.lsn)="Gene_lsn"
    results.logistic.final.adj=cbind(logistic.lsn, results.logistic.adj)
    colnames(results.logistic.final.adj)=c("Gene_lsn", "MRD.odds.ratio.adj", "MRD.lower95.adj", "MRD.upper95.adj", "MRD.p-value.adj", "MRD.q-value.adj")
    MRD.results.count.added.adj=cbind(results.logistic.final.adj, mrd.count)

  }

  # Run cox proportional hazard models for association with EFS
  if (is.character(efs))
    efs.data=merged.data[!is.na(merged.data$efs.time),]

  message(paste0("Running COXPH models for association with EFS: ",date()))

  efs.lsn.clms=efs.data[ ,grepl("ENSG", names(efs.data))]
  efs.time=efs
  efs.censor=efs.censor

  censor.efs=efs.data$efs.censor
  efs.data.final=cbind(censor.efs,efs.lsn.clms)
  efs.event= efs.data.final[efs.data.final$censor.efs==1,]
  efs.without.event= efs.data.final[efs.data.final$censor.efs==0,]

  efs.event=efs.event[,-1]
  efs.without.event=efs.without.event[,-1]

  efs.event.with.lsn=as.data.frame(colSums(efs.event== 1))
  efs.event.without.lsn=as.data.frame(colSums(efs.event== 0))
  efs.no.event.with.lsn=as.data.frame(colSums(efs.without.event== 1))
  efs.no.event.without.lsn=as.data.frame(colSums(efs.without.event== 0))

  efs.count=cbind(efs.event.with.lsn, efs.event.without.lsn,
                  efs.no.event.with.lsn, efs.no.event.without.lsn)

  colnames(efs.count)=c("efs.event.with.lsn", "efs.event.without.lsn",
                        "efs.no.event.with.lsn", "efs.no.event.without.lsn")

  if (is.null(covariate)) {
    EFSfit<-lapply(efs.lsn.clms,function(x) coxph(Surv(efs.time, efs.censor) ~ x , data=efs.data))

    ## To extract coefficients:
    EFS.Coeff<-lapply(EFSfit,function(f) summary(f)$coefficients[,1])
    ## Number of columns depend on how many variables are in the model, for example SNP+adjusting for risk group (should have 3 columns SNP+STANDARD+HIGH.RISK)
    EFS.coeff.mat <- do.call(rbind,lapply(EFS.Coeff,matrix,ncol=1,byrow=TRUE))

    ## To extract hazard.ratio:
    EFS.hazard<-lapply(EFSfit,function(f) summary(f)$coefficients[,2])
    EFS.hazard.mat = as.matrix(EFS.hazard)

    ## To extract model error:
    EFS.error = lapply(EFSfit,function(f) summary(f)$coefficients[,3])
    EFS.error.mat <- do.call(rbind,lapply(EFS.error,matrix,ncol=1,byrow=TRUE))

    ## To calculate lower95% CI:
    EFS.lower95 = exp(EFS.coeff.mat[,1] - 1.96*EFS.error.mat[,1])
    Efs.lower.mat = as.matrix(EFS.lower95)

    ## To calculate upper95% CI:
    EFS.upper95 = exp(EFS.coeff.mat[,1] + 1.96*EFS.error.mat[,1])
    EFS.upper.mat = as.matrix(EFS.upper95)

    ## To extract p-value:
    EFS.Pvalue <-lapply(EFSfit,function(f) summary(f)$coefficients[,5])
    EFS.Pvalue.mat = as.matrix(EFS.Pvalue)
    EFS.pvalue.mat=do.call(rbind,lapply(EFS.Pvalue,matrix,ncol=1,byrow=TRUE))
    pvalue.efs= EFS.pvalue.mat[,1]
    pi.hat=min(1,2*mean(pvalue.efs,na.rm=T))
    q.efs=pi.hat*p.adjust(pvalue.efs,method="fdr")
    ## To prepare final EFS results:
    EFS.results=cbind(EFS.hazard.mat, Efs.lower.mat, EFS.upper.mat, pvalue.efs, q.efs)
    efs.lsn=as.data.frame(colnames(lsn.mtx))
    colnames(efs.lsn)="Gene_lsn"
    results.efs.final=cbind(efs.lsn, EFS.results)
    colnames(results.efs.final)=c("Gene_lsn", "EFS.HR", "EFS.lower95", "EFS.upper95", "EFS.p-value", "EFS.q-value.adj")
    EFS.results.count.added=cbind(results.efs.final, efs.count)

  } else {

    # EFS with covariates
    covariates=covariate
    covariates=efs.data[,covariates]

    EFSfit.adj<-lapply(efs.lsn.clms,function(x) coxph(Surv(efs.time, efs.censor) ~ x+covariates , data=efs.data))

    ## To extract coefficients:
    EFS.Coeff.adj<-lapply(EFSfit.adj,function(f) summary(f)$coefficients[,1])
    ## Number of columns depend on how many variables are in the model, for example SNP+adjusting for risk group (should have 3 columns SNP+STANDARD+HIGH.RISK)
    EFS.coeff.mat.adj <- do.call(rbind,lapply(EFS.Coeff.adj,matrix,byrow=TRUE))
    num.grps=length(EFS.coeff.mat.adj)/ncol(lsn.mtx)
    EFS.coeff.mat.adj <- do.call(rbind,lapply(EFS.Coeff.adj,matrix,ncol=num.grps,byrow=TRUE))

    ## To extract hazard.ratio:
    EFS.hazard.adj<-lapply(EFSfit.adj,function(f) summary(f)$coefficients[,2])
    EFS.hazard.mat.adj <- do.call(rbind,lapply(EFS.hazard.adj,matrix,ncol=num.grps,byrow=TRUE))


    ## To extract model error:
    EFS.error.adj = lapply(EFSfit.adj,function(f) summary(f)$coefficients[,3])
    EFS.error.mat.adj <- do.call(rbind,lapply(EFS.error.adj,matrix,ncol=num.grps,byrow=TRUE))

    ## To calculate lower95% CI:
    EFS.lower95.adj = exp(EFS.coeff.mat.adj[,1] - 1.96*EFS.error.mat.adj[,1])
    Efs.lower.mat.adj = as.matrix(EFS.lower95.adj)

    ## To calculate upper95% CI:
    EFS.upper95.adj = exp(EFS.coeff.mat.adj[,1] + 1.96*EFS.error.mat.adj[,1])
    EFS.upper.mat.adj = as.matrix(EFS.upper95.adj)

    ## To extract p-value:
    EFS.Pvalue.adj <-lapply(EFSfit.adj,function(f) summary(f)$coefficients[,5])
    EFS.Pvalue.mat.adj = as.matrix(EFS.Pvalue.adj)
    EFS.pvalue.mat.adj=do.call(rbind,lapply(EFS.Pvalue.adj,matrix,ncol=num.grps,byrow=TRUE))
    pvalue.efs.adj=EFS.pvalue.mat.adj[,1]
    pi.hat=min(1,2*mean(pvalue.efs.adj,na.rm=T))
    q.efs.adj=pi.hat*p.adjust(pvalue.efs.adj,method="fdr")

    ## To prepare final EFS results:
    EFS.results.adj=cbind(EFS.hazard.mat.adj[,1], Efs.lower.mat.adj, EFS.upper.mat.adj, pvalue.efs.adj, q.efs.adj)

    efs.lsn=as.data.frame(colnames(lsn.mtx))
    colnames(efs.lsn)="Gene_lsn"
    results.efs.final.adj=cbind(efs.lsn, EFS.results.adj)
    colnames(results.efs.final.adj)=c("Gene_lsn", "EFS.HR.adj", "EFS.lower95.adj", "EFS.upper95.adj", "EFS.p-value.adj", "EFS.q-value.adj")
    EFS.results.count.added.adj=cbind(results.efs.final.adj, efs.count)

  }

  # Run cox proportional hazard models for association with OS
  if (is.character(os))
    os.data=merged.data[!is.na(merged.data$os.time),]

  message(paste0("Running COXPH models for association with OS: ",date()))

  os.lsn.clms=os.data[ ,grepl("ENSG", names(os.data))]
  os.time=os
  os.censor=os.censor
  censor.os=os.data$os.censor
  os.data.final=cbind(censor.os,os.lsn.clms)
  os.event= os.data.final[os.data.final$censor.os==1,]
  os.without.event= os.data.final[os.data.final$censor.os==0,]

  os.event=os.event[,-1]
  os.without.event=os.without.event[,-1]

  os.event.with.lsn=as.data.frame(colSums(os.event== 1))
  os.event.without.lsn=as.data.frame(colSums(os.event== 0))
  os.no.event.with.lsn=as.data.frame(colSums(os.without.event== 1))
  os.no.event.without.lsn=as.data.frame(colSums(os.without.event== 0))

  os.count=cbind(os.event.with.lsn, os.event.without.lsn,
                 os.no.event.with.lsn, os.no.event.without.lsn)

  colnames(os.count)=c("os.event.with.lsn", "os.event.without.lsn",
                       "os.no.event.with.lsn", "os.no.event.without.lsn")

  if (is.null(covariate)) {
    OSfit<-lapply(os.lsn.clms,function(x) coxph(Surv(os.time, os.censor) ~ x , data=os.data))

    ## To extract coefficients:
    OS.Coeff<-lapply(OSfit,function(f) summary(f)$coefficients[,1])
    ## Number of columns depend on how many variables are in the model, for example SNP+adjusting for risk group (should have 3 columns SNP+STANDARD+HIGH.RISK)
    OS.coeff.mat <- do.call(rbind,lapply(OS.Coeff,matrix,ncol=1,byrow=TRUE))

    ## To extract hazard.ratio:
    OS.hazard<-lapply(OSfit,function(f) summary(f)$coefficients[,2])
    OS.hazard.mat = as.matrix(OS.hazard)

    ## To extract model error:
    OS.error = lapply(OSfit,function(f) summary(f)$coefficients[,3])
    OS.error.mat <- do.call(rbind,lapply(OS.error,matrix,ncol=1,byrow=TRUE))

    ## To calculate lower95% CI:
    OS.lower95 = exp(OS.coeff.mat[,1] - 1.96*OS.error.mat[,1])
    Os.lower.mat = as.matrix(OS.lower95)

    ## To calculate upper95% CI:
    OS.upper95 = exp(OS.coeff.mat[,1] + 1.96*OS.error.mat[,1])
    OS.upper.mat = as.matrix(OS.upper95)

    ## To extract p-value:
    OS.Pvalue <-lapply(OSfit,function(f) summary(f)$coefficients[,5])
    OS.Pvalue.mat = as.matrix(OS.Pvalue)
    OS.pvalue.mat=do.call(rbind,lapply(OS.Pvalue,matrix,ncol=1,byrow=TRUE))
    pvalue.os=OS.pvalue.mat[,1]
    pi.hat=min(1,2*mean(pvalue.os,na.rm=T))
    q.os=pi.hat*p.adjust(pvalue.os,method="fdr")
    ## To prepare final EFS results:
    OS.results=cbind(OS.hazard.mat, Os.lower.mat, OS.upper.mat, pvalue.os, q.os)
    os.lsn=as.data.frame(colnames(lsn.mtx))
    colnames(os.lsn)="Gene_lsn"
    results.os.final=cbind(os.lsn, OS.results)
    colnames(results.os.final)=c("Gene_lsn",  "OS.HR", "OS.lower95", "OS.upper95", "OS.p-value", "OS.q-value")
    OS.results.count.added=cbind(results.os.final, os.count)

  } else {

    # OS with covariates
    covariates=covariate
    covariates=os.data[,covariates]

    OSfit.adj<-lapply(os.lsn.clms,function(x) coxph(Surv(os.time, os.censor) ~ x+covariates , data=os.data))


    ## To extract coefficients:
    OS.Coeff.adj<-lapply(OSfit.adj,function(f) summary(f)$coefficients[,1])
    ## Number of columns depend on how many variables are in the model, for example SNP+adjusting for risk group (should have 3 columns SNP+STANDARD+HIGH.RISK)
    OS.coeff.mat.adj <- do.call(rbind,lapply(OS.Coeff.adj,matrix,byrow=TRUE))
    num.grps=length(OS.coeff.mat.adj)/ncol(lsn.mtx)
    OS.coeff.mat.adj <- do.call(rbind,lapply(OS.Coeff.adj,matrix,ncol=num.grps,byrow=TRUE))

    ## To extract hazard.ratio:
    OS.hazard.adj<-lapply(OSfit.adj,function(f) summary(f)$coefficients[,2])
    OS.hazard.mat.adj <- do.call(rbind,lapply(OS.hazard.adj,matrix,ncol=num.grps,byrow=TRUE))


    ## To extract model error:
    OS.error.adj = lapply(OSfit.adj,function(f) summary(f)$coefficients[,3])
    OS.error.mat.adj <- do.call(rbind,lapply(OS.error.adj,matrix,ncol=num.grps,byrow=TRUE))

    ## To calculate lower95% CI:
    OS.lower95.adj = exp(OS.coeff.mat.adj[,1] - 1.96*OS.error.mat.adj[,1])
    Os.lower.mat.adj = as.matrix(OS.lower95.adj)

    ## To calculate upper95% CI:
    OS.upper95.adj = exp(OS.coeff.mat.adj[,1] + 1.96*OS.error.mat.adj[,1])
    OS.upper.mat.adj = as.matrix(OS.upper95.adj)

    ## To extract p-value:
    OS.Pvalue.adj <-lapply(OSfit.adj,function(f) summary(f)$coefficients[,5])
    OS.Pvalue.mat.adj = as.matrix(OS.Pvalue.adj)
    OS.pvalue.mat.adj=do.call(rbind,lapply(OS.Pvalue.adj,matrix,ncol=num.grps,byrow=TRUE))
    pvalue.os.adj=OS.pvalue.mat.adj[,1]
    pi.hat=min(1,2*mean(pvalue.os.adj,na.rm=T))
    q.os.adj=pi.hat*p.adjust(pvalue.os.adj,method="fdr")

    ## To prepare final EFS results:
    OS.results.adj=cbind(OS.hazard.mat.adj[,1], Os.lower.mat.adj, OS.upper.mat.adj, pvalue.os.adj, q.os.adj)

    os.lsn=as.data.frame(colnames(lsn.mtx))
    colnames(os.lsn)="Gene_lsn"
    results.os.final.adj=cbind(os.lsn, OS.results.adj)
    colnames(results.os.final.adj)=c("Gene_lsn", "OS.HR.adj", "OS.lower95.adj", "OS.upper95.adj", "OS.p-value.adj", "OS.q-value.adj")
    OS.results.count.added.adj=cbind(results.os.final.adj, os.count)

  }

  if (is.null(covariate)) {
    combined.unadjusted.results=cbind(MRD.results.count.added,
                                      EFS.results.count.added,
                                      OS.results.count.added)
    gene.list=str_split_fixed(combined.unadjusted.results$Gene_lsn, "_", 2)
    colnames(gene.list)=c("gene","lsn")
    annotation.merged=merge(annotation.data,gene.list,by="gene", all.y=TRUE)
    combined.unadjusted.results.final=cbind(annotation.merged, combined.unadjusted.results)
    return(combined.unadjusted.results.final)
  } else {
    combined.adjusted.results=cbind(MRD.results.count.added.adj,
                                    EFS.results.count.added.adj,
                                    OS.results.count.added.adj)
    gene.list=str_split_fixed(combined.adjusted.results$Gene_lsn, "_", 2)
    colnames(gene.list)=c("gene","lsn")
    annotation.merged=merge(annotation.data,gene.list,by="gene", all.y=TRUE)
    combined.adjusted.results.final=cbind(annotation.merged, combined.adjusted.results)
    return(combined.adjusted.results.final)
  }
}


########################################################
# Lesion boundaries evaluation
###################################################

grin.lsn.boundaries=function(lsn.data,
                             chrom.size)
{
  lsn.data=lsn.data
  chr.size=chrom.size
  chroms=chr.size$chrom
  lsn.bound = data.frame()
  {
    for (i in chroms) {
      chr=i
      chr.lsns=lsn.data[lsn.data$chrom==chr,]
      ord=order(chr.lsns$loc.start,
                chr.lsns$loc.end)
      uniq.start=unique(chr.lsns$loc.start)
      uniq.end=unique(chr.lsns$loc.end)
      chrom.size=chr.size$size[chr.size$chrom==chr]
      all.start=unique(c(1,uniq.start,uniq.end+1))
      all.end=unique(c(uniq.start-1,uniq.end,chrom.size))
      all.start=sort(all.start)
      all.end=sort(all.end)
      chr.lsn.loci=cbind.data.frame(gene=paste0("chr",chr,"_",all.start,"_",all.end),
                                    chrom=chr,
                                    loc.start=all.start,
                                    loc.end=all.end)
      chr.lsn.bound = chr.lsn.loci
      lsn.bound=rbind(lsn.bound, chr.lsn.loci)

    }
    lsn.bound$diff=(lsn.bound$loc.end - lsn.bound$loc.start)
    lsn.bound=lsn.bound[!lsn.bound$diff<=0,]
    
    return(lsn.bound)

  }
}


############################################################################################

grin.lsn.bound.plot=function(grin.res, # GRIN results (output of the grin.stats function)
                             lsn.colors=NULL, # Lesion colors
                             max.log10q=50) # Maximum log10 q value to be added to the plot

{
  if (!is.element("x.start",colnames(grin.res$lsn.data)))
    grin.res=compute.gw.coordinates(grin.res)

  grin.res$lsn.data$x.ID=as.numeric(as.factor(grin.res$lsn.data$ID))
  n=max(grin.res$lsn.data$x.ID)
  n.chr=nrow(grin.res$chr.size)

  # set up plotting region
  plot(c(-0.2,1.2)*n,c(0,-1.1*grin.res$chr.size$x.end[n.chr]),
       type="n",axes=F,xlab="",ylab="")

  # background colors for chromosomes
  rect(0*n,-grin.res$chr.size$x.start,
       1*n,-grin.res$chr.size$x.end,
       col=c("lightgray","gray"),
       border=NA)

  lsn.types=sort(unique(grin.res$lsn.index$lsn.type))

  if (is.null(lsn.colors))
  {
    lsn.colors=default.grin.colors(lsn.types)
  }
  grin.res$lsn.data$lsn.colors=lsn.colors[grin.res$lsn.data$lsn.type]

  grin.res$lsn.data$lsn.size=grin.res$lsn.data$x.end-grin.res$lsn.data$x.start
  ord=order(grin.res$lsn.data$lsn.size,decreasing=T)

  nsubj.mtx=unlist(grin.res$gene.hits[,paste0("nsubj.",lsn.types)])
  qval.mtx=unlist(grin.res$gene.hits[,paste0("q.nsubj.",lsn.types)])
  nsubj.data=cbind.data.frame(gene=grin.res$gene.hits$gene,
                              x.start=grin.res$gene.hits$x.start,
                              x.end=grin.res$gene.hits$x.end,
                              nsubj=nsubj.mtx,
                              log10q=-log10(qval.mtx),
                              lsn.type=rep(lsn.types,each=nrow(grin.res$gene.hits)))
  nsubj.data=nsubj.data[nsubj.data$nsubj>0,]
  nsubj.data$lsn.colors=lsn.colors[nsubj.data$lsn.type]

  ord=order(nsubj.data$nsubj,decreasing=T)
  nsubj.data=nsubj.data[ord,]

  nsubj.data$log10q[nsubj.data$log10q>max.log10q]=max.log10q

  ord=order(nsubj.data$log10q,decreasing=T)
  nsubj.data=nsubj.data[ord,]

  segments(0*n,
           -(nsubj.data$x.start+nsubj.data$x.end)/2,
           0*n+1*nsubj.data$log10q/max(nsubj.data$log10q)*n,
           col=nsubj.data$lsn.colors)

  text(c(n,0)[(1:n.chr)%%2+1],
       pos=c(4,2)[(1:n.chr)%%2+1],
       -(grin.res$chr.size$x.start+grin.res$chr.size$x.end)/2,
       grin.res$chr.size$chrom,
       cex=0.75)

  legend(n/2,-1.025*grin.res$chr.size$x.end[n.chr],
         fill=lsn.colors,cex=1,
         legend=names(lsn.colors),
         xjust=0.5,
         ncol=length(lsn.colors),
         border=NA,bty="n")

   text(1*n/2,0,
       "-log10(q)",cex=1,
       pos=3)


  text(c(0,1)*n,
       -grin.res$chr.size$x.end[n.chr],
       c(0,max.log10q),
       cex=0.75,pos=1)
}

