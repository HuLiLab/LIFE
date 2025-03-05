library(stringr)
ts=proc.time()[3]
options(stringsAsFactors=F)
iffun1=function(g1,g2) {return(g1/g2+g2/g1)}
iffun2=function(g1,g2) {return(g1*g2)}

args=commandArgs(trailingOnly=T)
fnret=args[1]

fnge=paste0(tools::file_path_sans_ext(tools::file_path_sans_ext(fnret)), '.rdata')

if(grepl('iffun2',fnret)){iffun=iffun2}
fntrunc=tools::file_path_sans_ext(fnret)

cat('loading ',fnge,'\n')
mt=readRDS(fnge)
genes=rownames(mt)
cat('loading ',fnret,'\n')
retmat=readRDS(fnret)
cat('retmat dim:');print(dim(retmat))

cat('iffun set to:');print(iffun)

cat('time used: ', proc.time()[3]-ts, '\n')
cat('range of sd:');print(range(retmat[,'sd']))
cat('range of mean:');print(range(retmat[,'mean']))
cat('time used: ', proc.time()[3]-ts, '\n')
nr=nrow(retmat)
ntopbot=1000

cat("number of rows: ", nr, ", top bot to select:", ntopbot, "\n")
if(nr<ntopbot*2){stop("not enough rows to select, have to be at least 2*ntopbot")}

cat('ordering by sd ...\n')
o=order(retmat[,'sd'], decreasing=F) # takes about 20 min for full result set
cat('time used: ', proc.time()[3]-ts, '\n')

cat('sorting ...\n')
top1000=o[1:ntopbot]
bot1000=o[(nrow(retmat)-ntopbot):nrow(retmat)]

topsel=retmat[top1000,]
botsel=retmat[bot1000,]

topsel=topsel[order(topsel[,'sd']),]
botsel=botsel[order(botsel[,'sd'],decreasing=T),]
cat('range of topsel-sd:');print(range(topsel[,'sd']))
cat('range of botsel-sd:');print(range(botsel[,'sd']))


fnout=paste0(fntrunc, '.top',ntopbot,'.csv')
cat('writing to file ', fnout, '\n')
mout=cbind(topsel, genes[topsel[,'i']], genes[topsel[,'j']])
colnames(mout)=c(colnames(topsel),'gene_i_name','gene_j_name')
write.table(mout, file=fnout, row.names=F, col.names=T, sep=',')

fnout=paste0(fntrunc, '.bot',ntopbot,'.csv')
cat('writing to file ', fnout, '\n')
mout=cbind(botsel, genes[botsel[,'i']], genes[botsel[,'j']])
colnames(mout)=c(colnames(botsel),'gene_i_name','gene_j_name') 
write.table(mout, file=fnout, row.names=F, col.names=T, sep=',')

cat('time used: ', proc.time()[3]-ts, '\n')
#stop()
nplot=50
fnpdf=paste0(fntrunc, '.top',nplot,'.pdf')
cat('writing to file', fnpdf,'\n')
pdf(fnpdf)
for(pos in 1:nplot){
#cat('pos=',pos,'\n')
i=topsel[pos,'i']
j=topsel[pos,'j']
lret=iffun(mt[i,],mt[j,])
#print(sd(lret))
#print(mean(lret))
xl=paste0('gene ',i,' - ',genes[i])
yl=paste0('gene ',j,' - ',genes[j])
if(!is.null(e2s))
  {
  xl=paste0('gene ',i,' - ', e2s[genes[i]],' - ',genes[i])
  yl=paste0('gene ',j,' - ', e2s[genes[j]],' - ',genes[j])
  }

plot(x=mt[i,],y=mt[j,], 
     xlab=xl,
     ylab=yl);

}
dev.off()
fnpdf=paste0(fntrunc, '.bot',nplot,'.pdf')
cat('writing to file', fnpdf,'\n')
pdf(fnpdf)
for(pos in 1:nplot){
#cat('pos=',pos,'\n')
i=botsel[pos,'i']
j=botsel[pos,'j']
lret=iffun(mt[i,],mt[j,])
#print(sd(lret))
#print(mean(lret))
xl=paste0('gene ',i,' - ',genes[i])
yl=paste0('gene ',j,' - ',genes[j])
if(!is.null(e2s))
  {
  xl=paste0('gene ',i,' - ', e2s[genes[i]],' - ',genes[i])
  yl=paste0('gene ',j,' - ', e2s[genes[j]],' - ',genes[j])

  }

plot(x=mt[i,],y=mt[j,], 
     xlab=xl,
     ylab=yl);

}
dev.off()

cat('time used: ', proc.time()[3]-ts, '\n')

