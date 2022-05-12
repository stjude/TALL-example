############################
# Supporting packages
library(writexl)
library(circlize)
library(biomaRt)

####################################
## Notes (03/03/2022) by Abdel
# I modified grin.stats function by adding prep.data and find.overlap before running count.hits function in line 64 

################################################
## Notes (03-07-2022)
# I replaced prep.lsn.type.matrix with an older version of the function (return required results without any errors)
# I added three functions (genomewide.lsn.plot, compute.gw.coordinates and default.grin.colors) to prepare genomewide lesion plot (manhattan plot)

###################################################
# Notes (03-08-2022)
# I modified top.grin.gene.plots function to plot only protein coding genes by gene.name not by ensembl.ID. Some biotypes might have so many genes by the same gene name such as (Y_RNA and U6) which stop the code
# I also modified grin.gene.plot function to plot genes by gene name not by ensembl.ID

#####################################
# Notes (8/24/21)
# moved order.index.lsn.data and order.index.gene.data up before
# creation of glp.data to try and remedy mismatch of pointers from
# glp.data to gene.data


#######################################
# interactive GRIN analysis
# rough idea

# grin.interactive=function()
# {
#   lsn.file=file.choose() # user specifies input file
#   lsn.data=read.csv(lsn.file)
#   genome.assembly=select.list("hg19","hg18")
#   gene.data=get.gene.data(genome.assembly)
#   chr.data=get.chrom.size.data(genome.assembly)
#   res=grin.stats(lsn.data,gene.data,chr.size)
#   output.dir=paste0(dirname(lsn.file),"/GRIN/",)
#   write.grin.xlsx(res)
# }

######################################
# Complete GRIN analysis

grin.stats=function(lsn.data,            # data.frame with columns ID (subject identifier), chrom, loc.start, loc.end, lsn.type
                    gene.data=NULL,      # data.frame with columns gene, chrom, loc.start, loc.end
                    chr.size=NULL,       # data.frame with columns chrom and size
                    genome.version=NULL) # character string with genome version
  
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
    ensembl.data=get.ensembl.grin.data(genome.version)
    if (is.null(gene.data))
    {
      gene.data=ensembl.data$gene.data
      gene.data$gene=gene.data$Ensembl_ID
    }
    if (is.null(chr.size))
      chr.size=ensembl.data$chr.size
  }

  prep.data=prep.gene.lsn.data(lsn.data,
                               gene.data)
  find.overlap=find.gene.lsn.overlaps(prep.data)
  hit.cnt=count.hits(find.overlap)
  hit.pvals=prob.hits(hit.cnt,
                      chr.size)
  return(hit.pvals)
}

########################################
# Generate plots of GRIN analysis

top.grin.gene.plots=function(grin.res,                      # result of grin.stats
                             lsn.clrs=NULL,                 # vector of color names with
                             q.max=0.10,                    # q-value threshold for plots
                             max.plots.per.type=25,         # maximum number of plots per lesion type
                             max.plots.overall=250)         # overall maximum number 
{
  lsn.types=unique(grin.res$lsn.data$lsn.type)
  if (is.null(lsn.clrs))
  {
    lsn.clrs=default.grin.colors(lsn.types)
  }
  
  gene.stats=grin.res[["gene.hits"]]
  gene.stats=gene.stats[gene.stats$biotype=="protein_coding",]
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
  #top.genes.sel=top.gene.stats[top.gene.stats$biotype=="protein_coding",]
  top.genes=top.gene.stats$gene.name
  
  
  for (i in 1:length(top.genes))
  {
    grin.gene.plot(grin.res,
                   top.genes[i],
                   lsn.clrs=lsn.clrs)
  }
}


####################################################
# Plot lesion data and GRIN results for one gene

grin.gene.plot=function(grin.res,
                        gene,
                        lsn.clrs=NULL,
                        expand=0.3)
  
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



###############################################################
prob.hits=function(hit.cnt,chr.size=NULL)
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
  
  res=list(gene.hits=gene.res,
           lsn.data=hit.cnt$lsn.data[,lsn.clms],
           gene.data=hit.cnt$gene.data[,gd.clms],
           gene.lsn.data=hit.cnt$gene.lsn.data,
           chr.size=chr.size,
           gene.index=hit.cnt$gene.index,
           lsn.index=hit.cnt$lsn.index)
  
  return(res)
}



###########################################
# Prepare gene and lesion data for later computations

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
 
  return(list(lsn.data=lsn.data,
              gene.data=gene.data,
              gene.lsn.data=gene.lsn.data,
              gene.index=gene.index,
              lsn.index=lsn.index))
   
}
  

#############################################
# Use the results of prep.gene.lsn.data to find lesion-gene overlaps

find.gene.lsn.overlaps=function(gl.data)

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
 
  res=list(lsn.data=lsn.data,
           gene.data=gene.data,
           gene.lsn.data=gene.lsn.data,
           gene.lsn.hits=gene.lsn.hits,
           gene.index=gene.index,
           lsn.index=lsn.index)
  
  return(res)
   
}

#####################################################
# Count hits

count.hits=function(ov.data)
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
  
  res=list(lsn.data=lsn.data,
           lsn.index=lsn.index,
           gene.data=gene.data,
           gene.index=gene.index,
           nhit.mtx=nhit.mtx,
           nsubj.mtx=nsubj.mtx,
           gene.lsn.data=gene.lsn.hits,
           glp.data=gene.lsn.data)
  
  return(res)
  
}

######################################
# Create a binary lesion hit matrix 

prep.binary.lsn.mtx=function(overlap.data,   # result of find.overlaps
                             min.ngrp=0)     # omit rows with fewer than min.ngrp hit or fewer than min.ngrp not hit
  
{
  gene.lsn=overlap.data$gene.lsn.hits
  
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
  
  keep.row=which(min.n>min.ngrp)
  gene.lsn.mtx=matrix(gene.lsn.mtx[keep.row,],
                      length(keep.row),n)
  colnames(gene.lsn.mtx)=as.character(uniq.ID)
  rownames(gene.lsn.mtx)=paste0(gene.lsn.index$gene,"_",
                                gene.lsn.index$lsn.type)[keep.row]
  
  
  return(gene.lsn.mtx)
}


######################################################
# prepare lesion type matrix

prep.lsn.type.matrix=function(overlap.data)
{
  gene.lsn=overlap.data$gene.lsn.hits
  
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
  
  return(grp.mtx)
}


#######################################
# Write GRIN results to an xlsx file

write.grin.xlsx=function(grin.result,
                         output.file) # output of grin.stats
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



###################################################################
# order and index gene data by chromosome, start location, and end location

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



################################################
# order and index lesion data by type, chromosome, and subject

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


#########################################
# compute ordered p-values

p.order=function(P)
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


##############################
# Compute a convolution of Bernoullis for each row of a Bernoulli success probability matrix

row.bern.conv=function(P,
                       max.x=NULL)
  
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

#######################################
# Compute the probability that a subject has a hit for each gene

row.prob.subj.hit=function(P,   # matrix of lesion hit probabilities, rows for genes, columns for lesions
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


##########################################
# Compute convolution for a series of independent Bernoulli random variables
# 

p.bern.conv=function(p,                   # vector of success probabilities
                     x=NULL,              # number of successes
                     tail.side="right")   # compute Pr(X<=x) for "left" or Pr(X>=x) for "right"
{
  # Limit to non-zero probabilities
  p=p[p>0]
  n=length(p)
  if (n==0) return(1)
  
  # Initialize result vector
  res=rep(0,n+1)
  res[1]=1
  
  for (i in 1:n)
  {
    pr1=exp(log(res[1:i])+log(p[i]))
    pr0=exp(log(res[1:i])+log(1-p[i]))
    res[1:(i+1)]=c(pr0,0)+c(0,pr1)
  }

  if (!is.null(x))
  {
    if (tail.side=="right") res=sum(res[(x+1):(n+1)])
    if (tail.side=="left")  res=sum(res[(x+1):(n+1)])
  }
  return(res)
}


#######################################################################

genomewide.lsn.plot=function(grin.res,
                             lsn.colors=NULL,
                             max.log10q=50)
  
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


#######################################
# Compute genome wide plotting coordinates

compute.gw.coordinates=function(grin.res,scl=1000000)
  
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
  
  return(grin.res)
}


####################################
# default colors for GRIN plots

default.grin.colors=function(lsn.types)
{
  message("Computationally assigning lesion type colors for GRIN plots.")
  uniq.types=sort(unique(lsn.types))
  n.types=length(uniq.types)
  default.colors=c("gold","red","blue",
                   "olivedrab","purple",
                   "brown","cyan",
                   "orange","steelblue",
                   "black")
  if (length(n.types)>length(default.colors))
    stop(paste0("Too many lesion types for default grin colors; please assign colors manually."))
  
  res=default.colors[1:n.types]
  names(res)=uniq.types
  return(res)
  
}

#################################################
# To prepare the matrix that we will pass to the oncoprint function from complexheat map library

grin.oncoprint.mtx=function(oncoprint.genes, # vector of oncoprint genes ensembl IDs
                            gene.lsn.data,   # data frame of mapped lesions for each patients (one component of GRIN.results list)
                            ensembl.annotation) # data frame with two columns "gene" for ensembl.ID and "gene.name" with gene symbols
  
{
  selected.genes= gene.lsn.data[gene.lsn.data$gene %in% oncoprint.genes,]
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
  mtx<-tibble::rownames_to_column(mtx, "gene")
  
  mtx.final=merge(ensembl.annotation,mtx,by="gene", all.y=TRUE)  # add gene name 
  mtx.final=mtx.final[,-1]
  rownames(mtx.final)=mtx.final[,1]
  mtx.final=mtx.final[,-1]
  
  return(mtx.final)
  
}


#####################
# get.ensembl.annotation function can be used to retrieve gene and regulatory features annotation data from ensembl database based on the specified genome assembly
###############################################

get.ensembl.annotation=function(genome.assembly)
  
{
  if (genome.assembly=="Human_GRCh38")
  {
    ensembl_GRCh38 = useEnsembl(biomart="genes",                  # retrieve data for human_GRCh38 from ensembl biomaRt version 104
                                dataset="hsapiens_gene_ensembl",  # specifying mirror="useast" is faster than not specifying any mirrors but we can not use mirror argument if we specified version (connection will be made to the main ensembl site)
                                version = "104")                  # specifying version is very important to extract data for GRCH38. If we did not specify version, query will extract data from the most updated version
    chromosomes = c(1:22, "X", "Y")
    hg38_gene_annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                               'start_position','end_position', 
                                               'description','external_gene_name', 
                                               'gene_biotype', 'strand', 'band'), 
                                  filters = 'chromosome_name',  values = chromosomes, 
                                  mart = ensembl_GRCh38)
    
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
    
    regulatory.hg38 = useEnsembl(biomart="regulation", 
                                 dataset="hsapiens_regulatory_feature",     # To retrieve data for human_regulatory features mapped to GRCh38 
                                 version = '104')
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
    
    
    res=list(gene.annotation=gene.data,
             reg.annotation=regulation.data)
    return(res)
  }
  
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
    
    gene.data=cbind.data.frame(gene=gene,  # we should change "gene" in the package functions to "Ensembl_ID"
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)
    
    regulatory.hg19 = useEnsembl(biomart="regulation", 
                                 dataset="hsapiens_regulatory_feature",     # To retrieve data for human_regulatory features mapped to GRCh37 
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
    
    gene.data=cbind.data.frame(gene=gene,  # we should change "gene" in the package functions to "Ensembl_ID"
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)
    
    regulatory.mm39 = useEnsembl(biomart="regulation", 
                                 dataset="mmusculus_regulatory_feature",     # To retrieve data for mouse_regulatory features mapped to mm39
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
    
    gene.data=cbind.data.frame(gene=gene,  # we should change "gene" in the package functions to "Ensembl_ID"
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)
    
    regulatory.mm10 = useEnsembl(biomart="regulation", 
                                 dataset="mmusculus_regulatory_feature",     # To retrieve data for mouse_regulatory features mapped to mm39
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


#####################
# get.chrom.length function will retrieve chromosome size data from chr.info txt file available on UCSC genome browser based on the specified genome assembly
###############################################

get.chrom.length=function(genome.assembly)
{
  if (genome.assembly=="Human_GRCh38")
  {
    chr.size.hg38= circlize::read.chromInfo(species = "hg38") #retrieve chromosome size data for GRCh38 (hg38) genome build from chr.info txt file available on UCSC genome browser
    chr.size.hg38=as.data.frame(chr.size.hg38)
    chr.size=cbind.data.frame(chrom=chr.size.hg38$chromosome,
                              size=chr.size.hg38$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))
    
    return(chr.size)
  }
  
  if (genome.assembly=="Human_GRCh37")
  {
    chr.size.hg19= circlize::read.chromInfo(species = "hg19")   #retrieve chromosome size data for GRCh37 (hg19) genome build from chr.info txt file available on UCSC genome browser
    chr.size.hg19=as.data.frame(chr.size.hg19)
    chr.size=cbind.data.frame(chrom=chr.size.hg19$chromosome,
                              size=chr.size.hg19$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))
    
    return(chr.size)
  }
  
  if (genome.assembly=="Mouse_HGCm39")
  {
    chr.size.mm39= circlize::read.chromInfo(species = "mm39")  #retrieve chromosome size data for Mouse_HGCm39 (mm39) genome build from chr.info txt file available on UCSC genome browser
    chr.size.mm39=as.data.frame(chr.size.mm39)
    chr.size=cbind.data.frame(chrom=chr.size.mm39$chromosome,
                              size=chr.size.mm39$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))
    
    return(chr.size)
  }
  
  if (genome.assembly=="Mouse_HGCm38")
  {
    chr.size.mm38= circlize::read.chromInfo(species = "mm10") #retrieve chromosome size data for Mouse_HGCm38 (mm10) genome build from chr.info txt file available on UCSC genome browser
    chr.size.mm38=as.data.frame(chr.size.mm38)
    chr.size=cbind.data.frame(chrom=chr.size.mm38$chromosome,
                              size=chr.size.mm38$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))
    
    return(chr.size)
  }
}



####################################################################
# The associate lesions with expression (ALEX) library of functions
##############################################################################

KW.hit.express=function(expr.mtx, # normalized expression matrix with genes in rows and subjects in columns
                        hit.grps, # each row assigns each subject to a hit group (no hits, one lesion or multiple lesions; each row is a subject)
                        row.mtch, # matrix or data.frame with columns "expr.row" and "hit.row" (Abdel: hit.row is a gene in GRIN output)
                        min.grp.size=4) # minimum number of subjects to include in one group
  
{
  p=apply(row.mtch,1,one.KW.pvalue,
          hit.grps=hit.grps,
          expr.mtx=expr.mtx,
          min.grp.size=min.grp.size)
  
  res=cbind(row.mtch,p.KW=p)
  
  return(res)
  
}

########################################
# Compute the KW p-value for one row of the hit group matrix paired
# with one row of the expression data matrix
# (to be used in the apply function)

one.KW.pvalue=function(one.row.mtch,    # one row of the row.mtch matrix
                       expr.mtx,        # expression data matrix
                       hit.grps,        # hit group matrix
                       min.grp.size=4)  # minimum group size to perform the test
  
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

define.grps=function(grps,min.grp.size=4)
{
  grp.tbl=table(grps)
  small.grp=(grp.tbl<=min.grp.size)
  exclude.grp=grps%in%(names(grp.tbl)[small.grp])
  grps[exclude.grp]="EXCLUDE"
  return(grps)
}



################################
# compute Kruskal-Wallis p-value for each row

row.KW.pvalue=function(R,               # data matrix with z-ranked rows
                       grps,            # vector of group labels
                       zrank.done=F)    # indicates whether rows of R have already been z-ranked
  
{
  if (!zrank.done)                    # if not already done, rank each row
    R=row.zrank(R)                  
  
  grps=as.character(grps)             # represent groups as character data
  L=model.matrix(~0+grps)             # get the group label matrix
  n=colSums(L)                        # group-specific sample size (Abdel: number of patients in each group)
  k=ncol(L)                           # number of groups
  for (i in 1:k) L[,i]=L[,i]/n[i]     # create matrix to compute mean ranks
  mean.ranks=R%*%L                    # compute mean ranks for each gene and group
  cs=(mean.ranks^2)%*%n               # compute KW chi-square statistic
  p=pchisq(cs,k-1,lower.tail=F)       # compute KW p-value  # Abdel; KW-p value comparing all groups 
  return(p)                           # return KW p-value
}


########################
# rank and z-transform each row of a data matrix

row.zrank=function(X)
  
{
  n=ncol(X)
  r=nrow(X)
  R=apply(X,1,rank,na.last="keep") # rank each row
  R=matrix(R,n,r)
  mn=colMeans(R,na.rm=T)           # mean rank for each row ## Abdel, this might be rowmeans not colnmeans????????
  std=apply(R,2,sd,na.rm=T)        # stdev of ranks for each row
  R=t(R)                           # transpose the rows
  R=(R-mn)/std                     # z-transform the ranks
  return(R)                        # return the result
}

# Function to compute stats for row summary (mean, median and sd)

stat.by.group=function(x,              # vector of data
                       g,              # vector of group labels corresponding to x
                       stat,           # name of a function to compute a stat
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

###########################################
# Compute FDR with the Pounds & Cheng (2006) estimator of 
# the proportion of tests with a true null (pi.hat)
# Reference: https://pubmed.ncbi.nlm.nih.gov/16777905/

pc06.fdr=function(p)
  
{
  pi.hat=min(1,2*mean(p,na.rm=T))
  q=pi.hat*p.adjust(p,method="fdr")
  return(q)
}
