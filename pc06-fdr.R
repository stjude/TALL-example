


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