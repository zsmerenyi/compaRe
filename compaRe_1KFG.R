####### Load packages #######
  rm(list=ls())
  suppressMessages(library(readr))
  suppressMessages(library(tidyr))
  suppressMessages(library(stringr))
  suppressMessages(library(ape))
  suppressMessages(library(tibble))
  suppressMessages(library(phangorn))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(parallel))

####### load functitons #######
  source("CFfunctions9G_ENG.R") 
  source("compaRe_functions_210324ENG.R")
  
####### set col scheme #######
  SD<-format(Sys.Date(), format="%y%m%d")
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  set.seed(2);col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  col_vector <- unique(col_vector[c(unique(c(44,45,47,46,26,20,24,54,58,22,19,18,8,56,7,25,4,38,35,34,31,30,29,27,6,2,1,42)),c(1:60)[-which(c(1:60) %in% unique(c(47,46,26,20,24,54,58,22,19,18,8,56,7,25,4,38,35,34,31,30,29,27,6,2,1,42)))])] )
  sz<-function(vmi){length(unique(vmi))}

#######   SETUP   ###########
  fold0<-"genetrees"                                         # folder name where the reconciled gene trees are located
  nproc<-20                                                   # number of CPUs
  nametypex<-"JGI"                                            # prefix / sufix / JGI 
  spectree<-"species.tree.nwk"                                # species tree   
  
#######   INTRO   ###########
# from the fold0 reconciled gene trees are replaced to the folder clusters and renamed
  if(dir.exists("clusters")){}else{dir.create("clusters")} 
  rawreconcil<-fold0                                      # folder of reconciled gene trees
  lr<-list.files(path = rawreconcil, full.names = T)
  file.copy(lr,paste0("./clusters"))                      # copy of genetrees
  oldNames<-sapply(str_split(lr,"/"), "[", 2)             # 2nd part of the full location
  newNames0<-str_replace(oldNames,"mcl_","Cluster")       # replace mcl_ to Cluster
  newNames0<-sapply(str_split(newNames0,"\\."), "[", 1)   # exclude extension
  newNames<-str_replace(newNames0,"MTSH_finSH",".tree")   # replace the second part of the name
  file.rename(from = file.path("clusters", oldNames), to = file.path("clusters", newNames))
  dictnames<-tibble(old=oldNames,new=newNames)
  write.table(dictnames, file=paste0('dictnames_','.csv'), quote=FALSE, sep=',', col.names = T, row.names = F)

####### set folders ###########
  fold1<-getwd()
  fold2<-paste0("RESpl",format(Sys.time(), format="%y%m%d")) #%H%M
  setwd(fold1);if(dir.exists(fold2)){}else{dir.create(fold2)} 
  clusters<-"clusters" 

####### Species Tree ####### 
  spt<-LadderizeTree(read.nexus(spectree))
  sptp<-spt
  sptp$edge.length<-NULL
  File2<-paste0("SpeciesTreeNodes",".pdf")
  pdf(file = File2, width=16,height=40, pointsize = 18, family = "sans", bg = "white")
  plot.phylo(sptp,cex=1.0,label.offset = 15.5)
  tiplabels(cex=0.79,col = "red",frame="none",adj = c(-0.3,0.5))
  nodelabels(cex=0.69,col = "red",frame="none")
  dev.off()

  nodes<-(length(spt$tip.label)+1):(length(spt$tip.label)+spt$Nnode)
  b=124
  dsz<-sapply(nodes,function(b) length(Descendants(spt,b)[[1]]))
  dst<-tibble(node=nodes,desc=dsz)
  dst<-dplyr::arrange(dst,desc(desc))
  gerinc<-dst[which(dst$desc>dst[1,"desc"][[1]]/3),]
  sptp<-spt
  sptp$edge.length<-NULL
  File2<-paste0("simpSpeciesTreeNodes",".pdf")
  pdf(file = File2, width=8,height=14, pointsize = 18, family = "sans", bg = "white")
  plot.phylo(sptp,cex=0.4,label.offset = 5.5,edge.color="lightgrey")
  tiplabels(cex=0.5,col = "red",frame="none",adj = c(-0.3,0.5))
  nodelabels(node=gerinc$node, cex=0.6,col = "red",frame="none")
  dev.off()

####### Open gene trees ####### 
setwd(fold1)
lf<-list.files(pattern = "Cluster", path = paste0("./",clusters))  
lf2<-lf

####### RUNNING #######  
st<-Sys.time() 
allDLT<-mclapply(lf2,FUN=compaRe2, mc.cores=nproc,spt=spt,col_vector=col_vector,
                 clusters=clusters,fold2=fold2,fold1=fold1,nametype = nametypex,abrazol=F) 
Sys.time()-st
save(allDLT,file=paste0("allDLT",SD,".saved"))

####### VISUALISATION of one cluster #######  
mapping<-function(x,csize=1.2){ 
  load(x)
    EK<-as_tibble(et);EK$cols<-"#000000"
    for(co in unique(OGC$cols)){
      EK[which.edge(bt,OGC[which(OGC$cols==co),"value"][[1]]),"cols"]<-co
    }
    bt$edge.length<-NULL
    fasize<-30*log10(meret)
    csokk<-1
    File2<-paste0("./",fold2,"/Fig_ell_",aneve,".pdf")
    pdf(file = File2, width=16,height=fasize, pointsize = 18, family = "sans", bg = "white")
    plot.phylo(bt,cex=0.8,edge.width = 6, label.offset = 10, edge.color = EK$cols)
    tiplabels(OGC[match(bt$tip.label, OGC$value),"OG"][[1]],adj = c(-1,0), bg=OGC[match(bt$tip.label, OGC$value),"cols"][[1]],cex=1.3*csokk)
    if(length(snfak4)>0){nodelabels(pch=16,node=as.numeric(names(snfak4)), frame="none", col="lightblue",cex=4*csokk)}
    if(length(vagfak4)>0){nodelabels(pch=16,node=as.numeric(names(vagfak4)), frame="none", col="yellow",cex=4*csokk)}
    if(length(dupnode)>0){nodelabels(dupnode,node=dupnode,bg=NULL, frame="none", col="red",cex=1.3*csokk)}
    if(length(specnode)>0){nodelabels(specnode,node=specnode, frame="none", col="green3",cex=1.3*csokk)}
    dev.off()
  
  tree<-spt;tree$edge.length <- NULL;size<-1
  sizecorr<-log10(length(tree$tip.label))/csize
  for(h in 4:7){
    mit<-DLT[,h][[1]]
    mi<-colnames(DLT)[h]
    
          while(F){
          File2<-paste0("./",fold2,"/",aneve,"_",mi,".png")
          png (file=File2,width=12,height=20,units="in",res=300)
          plot.phylo(tree, main=paste0(aneve," ",mi), label.offset = 0.4, node.depth = 2, show.tip.label = T, cex = sizecorr, edge.width=2)
          nodelabels(pch = 16, col =c("blue","red")[(mit<0*1)+1], node = DLT$treenode, cex = abs(mit) * size)
          nodelabels(text = mit, node = DLT$treenode, adj = c(0.5,0.5), bg = NULL, frame = "none", cex = sizecorr, col="white")
          dev.off()
          }
    File2<-paste0("./",fold2,"/",aneve,"_",mi,".pdf")
    pdf(file = File2, width=16,height=20, pointsize = 18, family = "sans", bg = "white")
    plot.phylo(tree, main=paste0(aneve," ",mi), label.offset = 0.4, node.depth = 2, show.tip.label = T, cex = sizecorr, edge.width=2)
    nodelabels(pch = 16, col =c("blue","red")[(mit<0*1)+1], node = DLT$treenode, cex = abs(mit) * size)
    nodelabels(text = mit, node = DLT$treenode, adj = c(0.5,0.5), bg = NULL, frame = "none", cex = sizecorr, col="white")
    dev.off()
  }
} 

#EXAMPLE:
mapping(paste0("./",fold2,"/Cluster2004.saved"),csize = 1.5)


####### mergeing of OGC files #######  
lfuzo<-function(x,ncpu){
  NL<-x #<-resout
  ll<-length(x)
  tag<-round(ll/ncpu,0)
  a=1
  LL<-list()
  for(i in 1:ncpu){
    b=a+(tag-1)
    LL[[i]]<-c(i,a,b)
    a=b+1
  }
  x<-LL[[1]]
  rag<-function(x){
    aktl<-NL[x[2]:x[3]]
    rag<-plyr::ldply(aktl,rbind)
    return(rag)
  }
  er<-mclapply(LL,FUN=rag,mc.cores = ncpu)
  errag<-plyr::ldply(er,rbind)
  errag<-tibble::as_tibble(errag)
  return(errag) 
}        

OFL<-list.files(path=fold2,pattern = "OGC",full.names = T)
LOF<-lapply(OFL,function(v) read_csv(v))
OGT<-lfuzo(LOF,10)
write.table(OGT, file=paste0('OGCT_ALL_',SD,'.csv'), quote=FALSE, sep=',', col.names = T, row.names = F)


# cls ellenorzes-------------------------------------------------------------------------------------------------
library(dplyr);library(tidyr)
cls<-read_csv("cls_RAG4_1FKG_1215.csv")
OGT<-read_csv("OGCT_ALL_220125.csv")
OGT$m_cl<-gsub("Cluster","mcl_",OGT$cluster)
COGT<-bind_cols(cls[,c(1,3)],OGT[match(cls$protID,OGT$value),])
COGT<-drop_na(COGT,cluster)
#COGT$m_cl<-gsub("Cluster","mcl_",COGT$cluster)
which(COGT$mcl!=COGT$m_cl) # nincs ahol elterne!!
cls[which(duplicated(cls$protID)),] #duplikalt feherje a clsben nincs!
setdiff(OGT$value,cls$protID) # nincs olyan feherje az OGCT-ben ami clsben ne lenne
dupok<-OGT[which(duplicated(OGT$value)),]
ell<-OGT[which(OGT$value %in% dupok$value),]
ell<-bind_cols(ell,cls[match(ell$value,cls$protID),])
#write.table(ell, file=paste0('ell','.csv'), quote=FALSE, sep=',', col.names = T, row.names = F)

sok<-setdiff(OGT$m_cl,cls$mcl) # ezek a clusterek bentmaradtak cls-hez kepest
ezmind<-OGT[which(OGT$m_cl %in% sok),]
ellc<-cls[which(cls$protID %in% ezmind$value),]
intersect(ellc$mcl,ezmind$m_cl)
setdiff(ezmind$value,ellc$protID)
torolni<-unique(ezmind$cluster)
write.table(torolni, file=paste0('torolni','.csv'), quote=FALSE, sep=',', col.names = T, row.names = F)
#----------------------------------------------------------------------------------------------------------------

####### mergeing of DLT files #######  
dltk<-list.files(path = fold2, pattern = "DLT", full.names = T)
st<-Sys.time() 
dltl<-lapply(dltk,function(x) read_csv(x,col_names = T))
Sys.time()-st
save(dltl,file=paste0("compare_1KFG",SD,".saved"))
DLT<-do.call(rbind,dltl)
nev<-"F1K"
DLTV<-DLT 
write.table(DLTV, file=paste0('DLTV_',nev,'.csv'), quote=FALSE, sep=',', col.names = T, row.names = F)

ir <- DLTV %>% group_by(Cluster)
cldb<-group_split(ir)

####### adding birth and extinction values #######  
clDeathBirth<-function(i){
  DLTakt<-i
  DLTakt$clmrca<-0
  DLTakt$fulloss<-0
  if(max(DLTakt$Orthogroup_gains)==0){
    clmrcanode<-min(DLTakt[which(DLTakt$Terminal_gains>0),"treenode"])
  }else{
    clmrcanode<-min(DLTakt[which(DLTakt$Orthogroup_gains>0),"treenode"]) # Gene family birth
  }
  DLTakt[DLTakt$treenode==clmrcanode,"clmrca"]<-1
  DLTakt[which(DLTakt$Net_gains<0 & DLTakt$Copynum==0),"fulloss"]<-1     # Gene family extinction
  return(DLTakt)
}
cDB<-mclapply(cldb, FUN=clDeathBirth, mc.cores = 30)
DBT<-lfuzo(cDB,20)
write.table(DBT, file=paste0('DLTgy_',nev,'.csv'), quote=FALSE, sep=',', col.names = T, row.names = F)


####### TOTAL MAPPING #######  
DLTgy<-read_csv("DLTgy_F1K.csv")
DLTV<-DLTgy 
sz(DLTV$Cluster)

DLTV2<-DLTV
setdiff(F1$cluster,DLTV2$Cluster)
sz(F1$cluster)
sz(DLTV2$Cluster)
DLTSUM<-ddply(DLTV2,c("treenode"),summarise,sumLosses=sum(Losses),sumGains=sum(Gains),sumNet_gains=sum(Net_gains),sumCopynum=sum(Copynum),sumclmrca=sum(clmrca),sumfulloss=sum(fulloss))#,sumclmrca=sum(clmrca),sumfulloss=sum(fulloss)
DLTSUM$Cluster<-nev<-"All_20772cl_"
DLTSUM$sumfulloss<-DLTSUM$sumfulloss*(-1)
write.table(DLTSUM, file=paste0('DLTSUM_',nev,'.csv'), quote=FALSE, sep=',', col.names = T, row.names = F)

tblue <- t_col("blue", perc = 80, name = "tblue")
tred <- t_col("red", perc = 80, name = "tred")

csize=2
tree<-spt;tree$edge.length <- NULL;size0<-0.002 #0.0008 #0.003
sizecorr<-log10(length(tree$tip.label))/csize
aneve<-nev

data.frame(colnames(DLTSUM))
h=5
for(h in 2:7){
  akt<-DLTSUM[(length(spt$tip.label)+1):nrow(DLTSUM),]
  mit<-akt[,h]
  if(h==5){size<-size0/5}else{size<-size0}
  if(h==6 | h==7){size<-size0*5}else{size<-size}
  mi<-colnames(DLTSUM)[h]
  addn<-size*10000
  File2<-paste0(aneve,addn,"_",mi,".pdf")
  pdf(file = File2, width=12,height=18, pointsize = 18, family = "sans", bg = "white")
  plot.phylo(tree, edge.color="lightgrey",main=paste0(aneve," ",mi), label.offset = 0.8, node.depth = 2, show.tip.label = T, cex = 0.6*sizecorr, edge.width=2)
  nodelabels(pch = 16, col =c(tblue,tred)[(mit<0*1)+1], node = akt$treenode, cex = abs(mit) * size)
  nodelabels(text = mit, node = akt$treenode, adj = c(0.5,0.5), bg = NULL, frame = "none", cex = 0.9*sizecorr, col="black")
  nodelabels(bg="white",frame="none",col="green3",adj=c(3,0.5), cex=0.1)
  dev.off()
  
  File2<-paste0(aneve,addn,"_",mi,".png")
  png (file=File2,width=8,height=12,units="in",res=600)
  plot.phylo(tree, edge.color="lightgrey",main=paste0(aneve," ",mi), label.offset = 0.8, node.depth = 2, show.tip.label = T, cex = 0.55*sizecorr, edge.width=2)
  nodelabels(pch = 16, col =c(tblue,tred)[(mit<0*1)+1], node = akt$treenode, cex = abs(mit) * size)
  nodelabels(text = mit, node = akt$treenode, adj = c(0.5,0.5), bg = NULL, frame = "none", cex = 0.9*sizecorr, col="black")
  nodelabels(bg="white",frame="none",col="green3",adj=c(3,0.5), cex=0.1)
  dev.off()
}
save(list=ls(),file=paste0("compare_",nev,SD,".saved"))












