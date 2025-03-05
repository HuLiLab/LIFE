ts=proc.time()[3]
options(stringsAsFactors=F)
iffun1=function(g1,g2) {return(g1/g2+g2/g1)}
funstr='iffun1'



args=commandArgs(trailingOnly=T)
if(length(args)==0){
cat("usage: Rscript this_script.R <expression matrix, genes by samples, log2, filtered.rds>\n")
stop()

}
fnrd=args[1]

fntrunc=tools::file_path_sans_ext(fnrd)
cat('loading file ', fnrd, '\n')
mt=readRDS(fnrd)
print(dim(mt))
cat('time used: ', proc.time()[3]-ts, '\n')

nr=nrow(mt)
nc=ncol(mt)
n=0

totalpair=(nr-1+1)*(nr-1)/2

ticksize=ceiling(totalpair/100)


cat('initializing result matrix for total pairs ',totalpair,'...\n')

retmat=matrix(NA, ncol=4, nrow=totalpair)
colnames(retmat)=c('i','j','sd','mean')


cat('time used: ', proc.time()[3]-ts, '\n')
for(i in 1:(nr-1))
  {
  #cat('i=',i,'\n')
  for(j in (i+1):nr)
    {
    n=n+1
    lret2=iffun1(mt[i,],mt[j,])
    sd_t=sd(lret2)
    mean_t=mean(lret2)
    retmat[n,]=c(i,j,sd_t, mean_t)
    if(n%%ticksize==0)
      {
      cat('n=',n,', ', floor(n/totalpair*100),'%  ')
      cat('time used: ', proc.time()[3]-ts, '\n')
      }
    }
  }
fnout=paste0(fntrunc, '.',funstr,'.rds')
#cat('time used: ', proc.time()[3]-ts, '\n')
cat('writing to file ', fnout, '\n')
saveRDS(retmat, file=fnout)
cat('time used: ', proc.time()[3]-ts, '\n')


