#!/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import time
import sys
import os
import re
import glob
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import rpy2.robjects as robjects
r = robjects.r

#np.random.seed(1337)  # for reproducibility
def append_to_mainfn(fnin, str_app):
  (strpath, strfn) = os.path.split(os.path.abspath(fnin))
  (strmainfn, strext) = os.path.splitext(strfn)

  (strmainfn2, strext2) = os.path.splitext(strmainfn)#secondary ext
  if strmainfn2=='' or strext2=='':#secondary split not needed
    pass
    #strmainfn=strmainfn; strext=strext#no operation
  else:#secondary splittinig needed
    strmainfn=strmainfn2
    strext=strext2+strext
  str_r = os.path.join(strpath, strmainfn+ str_app +strext)
  return str_r
#end def append_to_mainfn(fnin, str_app):


def float_to_str(v):
  #make a float value suitable for file name by removing its decimal dot.
  #
  s='%g'%v
  if '.' in s: s=s.replace('.','')
  return s

########################################################################################################

respath='.'
ntop=1000
 
l_fntop=glob.glob(respath+'/*iffun*.normfdr005.top1000.csv')
allrestrunc=os.path.basename(os.path.abspath(respath))+'.normfdr005_top1000'; ntop=1000



print(l_fntop)
all_alldeg=list()
def r_plot_deg(str_title, list_val, fnpdf):
  
  r_val = robjects.IntVector(list_val)
  r.assign("r_val", r_val)
  r.assign("fnpdf", fnpdf)
  r.assign("plttitle", str_title)

  robjects.r('''
        tbl=table(r_val)
        fncsv=paste0(tools::file_path_sans_ext(fnpdf), '.csv')
        mo=cbind(names(tbl), tbl)
        colnames(mo)=c('deg','count')
        cat('writing to file ', fncsv,'\n')
        write.csv(mo, file=fncsv, row.names=F, col.names=T)
        xrng=c(1:max(tbl))
        for(xval in as.character(xrng)){
          if(length(tbl[xval])==0){tbl[xval]=0}
          }
        x=as.numeric(names(tbl))
        tbl=tbl[order(x, decreasing=F)]
        #y=as.numeric(tbl)
        #range(x)
        #range(y)
        cat("plotting to file ", fnpdf, '\n')
        pdf(fnpdf)
        bpret=barplot(height=tbl, ylim=c(0, max(tbl)+10),
              main=plttitle,
              cex.names=0.6)#x tick label size
        text(x = bpret, y = as.numeric(tbl), label = as.numeric(tbl), pos = 3, cex = 0.4, col = "blue")
        dev.off()
        ''')

#end def
def r_plot_hubcutoff_vs_nodescovered(str_title, l_hc, l_tc, fnpdf):
  r_hc=robjects.IntVector(l_hc)
  r_tc=robjects.IntVector(l_tc)

  r.assign("r_hc", r_hc)
  r.assign("r_tc", r_tc)

  r.assign("fnpdf", fnpdf)
  r.assign("plttitle", str_title)

  robjects.r('''
        cat("plotting to file ", fnpdf, '\n')
        pdf(fnpdf)
        #bpret=barplot(height=tbl, ylim=c(0, max(tbl)+10),
        #      main=plttitle,
        #      cex.names=0.6)#x tick label size
        #text(x = bpret, y = as.numeric(tbl), label = as.numeric(tbl), pos = 3, cex = 0.4, col = "blue")
        plot(x=r_hc, y=r_tc, main=plttitle)
        dev.off()
        ''')
#end def r_plot_hubcutoff_vs_nodescovered(str_title, l_hc, l_tc, fnpdf):

hubcutoff=3
minw=1
for fn in l_fntop:
  df_t_a=pd.read_csv(fn, sep=',',na_filter=False)

  g=nx.Graph()#DiGraph()

  fntrunc=os.path.splitext(os.path.basename(fn))[0]


  #df_t=pd.read_csv(fn, sep=',')
  df_t=df_t_a
  df_t=df_t.sort_values(by='sd', ascending=True)                                                                             

  if df_t.shape[0]==ntop: 
    pass
  else:
    fntrunc=fntrunc+'.top%d'%ntop
  #end if
  df_t=df_t.ix[0:ntop-1,]

  nr, nc=df_t.shape
  enf=set()
  for ir in range(nr):
    si=df_t.ix[ir,'gene_i_name'] 
    sj=df_t.ix[ir,'gene_j_name'] 

    sd=df_t.ix[ir,'sd']
    w=1/sd#-np.log10(sd)
    if w<0: w=minw
    g.add_edge(si, sj, val=sd, weight=w)
  #end for ir
  for n in g.nodes:
    nbrs=g.neighbors(n)
    #print(n);print(nbrs);input()
    ws=0
    for nb in nbrs:
#      print(g[n][nb]);input()
      ws+=g[n][nb]['weight']
    #
    g.nodes[n]['weight']=ws
  fngml=fntrunc+'.gml'
  print('writing to file %s'%fngml)
  nx.write_gml(G=g, path=fngml)


  dgv=g.degree()#new version of networkx returns this as DegreeView
  dg=dict()
  for k,v in dgv: dg[k]=v
  alldeg=dg.values()
#  continue

  hubs=[k for k in dg.keys() if dg[k]>=hubcutoff]
  s_nbr=set()
  for hub in hubs:
    nbrs=g.neighbors(hub)
    for v in nbrs: s_nbr.add(v)
  #
  str_hubstat='hubcutoff=%d, %d hubs covered %d genes out of a total of %d genes and %d edges '%(hubcutoff, len(hubs), len(s_nbr), len(g.nodes()), len(g.edges()))
  fnout_t=fntrunc+'.hubstats'
  print(str_hubstat)
  open(fnout_t,'w').write(str_hubstat+'\n')

  #sys.exit()
  #continue
  all_alldeg+=alldeg
  fnpdf=os.path.join(respath, fntrunc+'.deghist.pdf')
  r_plot_deg("total %d unique nodes, %d edges"%(len(g.nodes()), len(g.edges())), alldeg, fnpdf)


  #now iterate through different hubcutoffs and plot total nodes covered by hubs vs. hubcutoff
  uniq_deg=sorted(list(set(alldeg)))
  l_hc=list()#hub cutoffs
  l_tc=list()#total nodes covered
  for hubcutoff_i in [0]+uniq_deg:
    hubs=[k for k in dg.keys() if dg[k]>=hubcutoff_i]
    s_nbr=set()
    for hub in hubs:
      nbrs=[v for v in g.neighbors(hub)]
      for v in nbrs+[hub]: s_nbr.add(v)# plus hubs themselves
    #
    tot_cov=len(s_nbr)
    l_hc.append(hubcutoff_i)
    l_tc.append(tot_cov)
  #end hubcutoff_i
  fnpdf=os.path.join(respath, fntrunc+'.hubcutoff_covered.pdf')
  r_plot_hubcutoff_vs_nodescovered('hubcutoff vs. nodes covered, out of %d nodes'%len(g.nodes()),l_hc,l_tc, fnpdf)
  #break
#end for fn 

