options(stringsAsFactors=F)

ts=proc.time()[3]


args=commandArgs(trailingOnly=T)
if(length(args)==0){
cat("usage: Rscript this_script.R <expression matrix, genes by samples .rds>\n")
stop()

}
fnin=args[1]
print(fnin)

fnot=tools::file_path_sans_ext(fnin)
em=readRDS(fnin)
rsd=genefilter::rowSds(em)

sdmt=rsd
cat('taking log2 before normal test\n')
mt=log2(em+2)
genes=rownames(mt)
  
#shapiro.test for normality for each gene: a small p-value rejects that a gene conforms to normal distribution
cat('shapiro test for normality...\n')
lpv_spr=rep(NA, nrow(mt))
for(i in 1:nrow(mt))
  {
  sdi=sdmt[i]
  if(sdi==0){lpv_spr[i]=0
    }else{
    ret=shapiro.test(mt[i,])
    lpv_spr[i]=ret$p.value
    }
  }
lpv_spr_fdr=p.adjust(lpv_spr, method='fdr')
names(lpv_spr_fdr)=rownames(mt)
cat('time used: ', proc.time()[3]-ts, '\n')
bsel=lpv_spr_fdr>0.05
#bsel=lpv_spr>0.1 # even fewer genes


ems=mt[bsel,]

#hist(ems[1,], main=paste0(tt2, 'n=', ncol(mt), ' ', rownames(ems)[1], ' shapiro fdr=', lpv_spr_fdr[rownames(ems)[1]]))
#hist(em[!bsel,][1,], main=paste0(tt2, 'n=', ncol(mt), ' ', rownames(em[!bsel,])[1], ' shapiro fdr=', lpv_spr_fdr[rownames(em[!bsel,])[1]]))
emsl2p2=ems#log2(ems+2)
nc=ncol(emsl2p2)
#[1] 663
nr=nrow(emsl2p2)
#[1] 54038

fno=paste0(fnot, '.noMT_log2p2_normfdr005.',nc,'.',nr,'.rds')
saveRDS(emsl2p2, file=fno)

cat('end of script, time used: ',proc.time()[3]-ts, '\n')

