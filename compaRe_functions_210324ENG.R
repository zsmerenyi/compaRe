####### FUNCTIONS ###########
sz<-function(vmi){length(unique(vmi))}

mapping<-function(x,csize=1.2,genfa=T){ 
  load(x)
  if(genfa){ # Gene Tree visualisation with coloured Orthogroups
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
  } 
  
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

# counting duplications
Asz<-function(lx,colsz){sum(cfpl3[[lx]][which(cfpl3[[lx]]$tree_nodes %in% Ascphy),colsz])}
Bsz<-function(lx,colsz){sum(cfpl3[[lx]][which(cfpl3[[lx]]$tree_nodes %in% Basphy),colsz])}
Allsz<-function(lx,colsz){sum(cfpl3[[lx]][,colsz])}
Aszrat<-function(lx,colsz){mean(abs(cfpl3[[lx]][which(cfpl3[[lx]]$tree_nodes %in% Ascphy),colsz])/cfpl3[[lx]][which(cfpl3[[lx]]$tree_nodes %in% Ascphy),"edge.length"])}
Bszrat<-function(lx,colsz){mean(abs(cfpl3[[lx]][which(cfpl3[[lx]]$tree_nodes %in% Basphy),colsz])/cfpl3[[lx]][which(cfpl3[[lx]]$tree_nodes %in% Basphy),"edge.length"])}
# colnames(cltn[[1]])  # colnames that we wanted to use for duplication rate calculation
# colsz<-"Losses"      # name of column
# tree<-tree_lad       # ladderized tree
# whichnode=F          # F for all nodes or node names if we want to specify certain nodes
# ratename<-"lossrate" # output column name for rate
# showprevcop=F        # copynum of parental node
Pszrat<-function(akt,colsz,tree,whichnode,showprevcop=F,ratename="rate"){
  akt<-dplyr::arrange(akt,treenode) 
  if(class(whichnode)!="numeric"){
    whichnode<-treenodes<-akt$treenode
    prevnode<-sapply(whichnode, function(k) tree$edge[tree$edge[,2]==k,1]) # parental node
    treeprevnodes <- sapply(prevnode, function(x) which(akt$treenode==x))
    kivesz<-which(is.na(treeprevnodes>0)) # bzc first node has no length
    treeprevnodes[kivesz]<-1
    treeprevnodes<-unlist(treeprevnodes)
  }else{
    treenodes<-which(akt$treenode %in% whichnode) # filter for nodes
    prevnode<-sapply(whichnode, function(k) tree$edge[tree$edge[,2]==k,1]) # parental node
    treeprevnodes <- sapply(prevnode, function(x) which(akt$treenode==x))
    kivesz<-NULL
  }
  
  #column (eg gain ) for each node / node length / parental node copy number
  ratios<-abs(akt[treenodes,colsz])/akt[treenodes,"edge.length"]
  ac<-akt[treeprevnodes,"Copynum"][[1]]
  ac[ac==0]<-1   # If the parental node copynum=0 then we should divide by 1
  ac[kivesz]<-1
  ratios2<-ratios/ac
  if(showprevcop){akt$prevcopynum<-ac}
  akt[,ratename]<-as.numeric(unlist(ratios2))
  return(akt)
}

LadderizeTree <- function(tree, temp_file = "temp", orientation = "left"){
  if(file.exists(paste0("./", temp_file))){
    stop("The chosen temporary file exists! Please choose an other temp_file name")
  }
  if(orientation == "left"){
    right <- FALSE
  }else{
    right <- TRUE
  }
  tree_temp <- ladderize(tree, right = right)
  write.tree(tree_temp, file = paste0("./", temp_file, ".tre"))
  tree_lad <- read.tree(paste0("./", temp_file, ".tre"))
  file.remove(paste0("./", temp_file, ".tre"))
  return(tree_lad)
}

# species part of a proteinID
spresz<-function(s,nametype=c("prefix","sufix","JGI")){
  s2<-str_split(s,"_")
  switch(nametype,
         prefix={
           s3<-sapply(s2, "[", 1)
         },
         sufix={
           s3<-sapply(s2, "[", 2)
         },
         JGI={
           s3<-sapply(1:length(s2), function(x) paste0(s2[[x]][1:(length(s2[[x]])-1)],collapse = "_") )
         })
  return(s3)
}

# character to integer
chrtoint<-function(tabla,oszl=c(1:1)){
  orignames<-colnames(tabla)
  for(h in oszl){
    tabla[,h]<-as.integer(tabla[,h][[1]])
  }
  colnames(tabla)<-orignames
  return(tabla)
}

chrelo<-function(TABLE,x){TABLE[,x]<-as.character(TABLE[,x][[1]])}


dupdont<-function(unn,x,bt,kezdo,st,nametype){
  kell<-bt$edge[which(bt$edge[,1]==unn[x]),2]
  if(kell[1]<=kezdo){ap<-bt$tip.label[kell[1]]}else{ap<-st[[kell[1]-kezdo]]$tip.label} # adott node A irany 1. nodeja fele faj v node van, mindegy kiextrahaljuk  
  if(kell[2]<=kezdo){bp<-bt$tip.label[kell[2]]}else{bp<-st[[kell[2]-kezdo]]$tip.label} # adott node B irany 1. nodeja fele faj v node van, mindegy kiextrahaljuk
  A<-unique(spresz(ap,nametype=nametype))
  B<-unique(spresz(bp,nametype=nametype))
  er<-length(intersect(A,B))>0
  la<-length(A)
  lb<-length(B)
  er2<-c(er,paste0(A,collapse = "|"),paste0(B,collapse = "|"),la,lb,ifelse(la>=lb,kell[1],kell[2]),kell[1],kell[2]) # ha ugyanannyi fajszamra az A es B ag akkor az A-t valasztja!! ennel szofisztikaltabb mod?
  return(er2)
}

# placement of losses
losscost<-function(spt,avx,abrazol){
  nev=unique(paste0(avx$cluster,"_OG",avx$OG))
  exist<- which(spt$tip.label %in% avx$spec)                   # species in the cluster
  lossok<- c(1:length(spt$tip.label))[-exist]                  # missing species from the cluster
  kezd<-length(spt$tip.label)+1                                # first node
  emrca<-getMRCA(spt,exist)                                    # cluster MRCA
  EK<-as_tibble(spt$edge);EK$erd<-1                            # edge table
  elotte<-which(EK<emrca & EK>=kezd)
  if(length(elotte)==0){
    EK$erd<-2
  }else{
    csakaz<-Descendants(spt,node = emrca,type = "all")
    csakaz<-c(emrca,csakaz)
    EK[which(EK$V1 %in% csakaz & EK$V2 %in% csakaz),"erd"]<-2         # nodes before the MRCA of the cluster
  }
  kiv<-which.edge(spt,exist)                                   # edges between the species in the cluster
  EK[kiv,"erd"]<-EK[kiv,"erd"]+1                               # add +1 if a species is member of the clustera => 2->3 => member=3 loss=2
  ecol<-c("darkgrey","blue","red")[EK$erd];table(ecol)         # colouring darkgrey=1 blue=2(loss) red=3(cluster membership)
  gerinc<-unique(EK[kiv,1][[1]])                               # backbone among cluster memberships "gerinc"
  gl<-sapply(gerinc,function(g) sz(EK[EK$V1==g,"erd"][[1]]))   # Along the backbone 33 (continued) or 23 (loss) 
  gcol<-c("red","blue")[gl];table(gcol)                        # colouring of edges red=cluster membership, blue=loss
  err<-as.data.frame(table(gl))
  mmo<-as_tibble(cbind(gerinc,gl))                             # backbone nodes loss(2) non loss(1)
  mmo$cluster<-unique(avx$cluster)
  mmo$OG<-unique(avx$OG)
  losnodes<-EK[EK$erd==2,]  
  losnodes<-losnodes[losnodes$V1 %in% EK[kiv,"V1"][[1]],]
  mmo$lossnode<-losnodes[match(mmo$gerinc,losnodes$V1),"V2"][[1]]
  while(F){
    plot.phylo(spt,cex=0.8,edge.width = 1, label.offset = 0.2, edge.color = ecol)
    tiplabels(pch=1,tip = exist, col="red", cex=1)
    tiplabels(pch=1,tip = lossok, col="blue", cex=1)
    nodelabels(pch=16,node=gerinc,cex=1, col="green3")
    nodelabels(pch=16,node=gerinc,cex=1, col=gcol)
    nodelabels(cex=0.6,frame="none", col="red")
    nodelabels(mmo$lossnode,node=mmo$lossnode, bg="yellow")
  }
  # visualisation =T
  if(abrazol){
    labsize<-log10(length(spt$tip.label))/6
    File2<-paste0("./",fold2,"/lossFig_",nev,".png")
    png (file=File2,width=12,height=18,units="in",res=200)
    plot.phylo(spt,cex=labsize,edge.width = 2, label.offset = 0.05, edge.color = ecol)
    tiplabels(pch=1,tip = exist, col="red", cex=2*labsize)
    tiplabels(pch=1,tip = lossok, col="blue", cex=2*labsize)
    nodelabels(pch=16,node=gerinc,cex=2*labsize, col=gcol)
    nodelabels(mmo$lossnode,node=mmo$lossnode, bg="yellow")
    dev.off()
  }
  return(list(err,mmo))
}

compaRe2<-function(x,spt=spt,col_vector=col_vector,clusters=clusters,fold2=fold2,fold1=fold1,nametype=nametype,abrazol=F){
  OGC<-NULL
  
  orthologcoding<-function(dupe=dupe,vag=vag,osibbek=osibbek,kezdo=kezdo, bt.=bt,st=st){
    kihova<-as_tibble(bt$tip.label);kihova$OG<-0    # default protein tag is 0
    q=0;kovspec<-as.character(vag);btwspecnode=T    # default values
    while(btwspecnode){                             # defining first speciation nodes 
      for(megvizs in kovspec){                      # all potential speciation nodes
        bekerul<-NULL;inspecnode=T;kivgy=NULL       # iteration default values
        while(inspecnode){                          # iteration till we have nodes
          for(i in megvizs){
            # i = as.character(megvizs) =12
            ake<-dupe[dupe$NODE==i,];cat(i,"\t")    # ake variable is one row from the table dupnodeos
            if(ake$dupnode){                        # If the node is a Duplication
              ez<-ake$more                          # checking that direction where we have more species
              bekerul<-c(bekerul,ez[ez<=kezdo])     # if it is tip nodes
              megvizs<-c(megvizs,ez[ez>kezdo])      # internal node => we can analyse further
            }else{                                  # dupnode=F => Speciation=> both direction are important
              ez<-c(ake$Anode,ake$Bnode)            # node A and B seems to be important
              bekerul<-c(bekerul,ez[ez<=kezdo])     # if it is tip nodes
              megvizs<-c(megvizs,ez[ez>kezdo])      # internal node => we can analyse further
            }
            megvizs<-setdiff(megvizs,i)             # we filtered out the analysed nodes
            kivesz<-dupe[which(dupe$NODE %in% megvizs),] # havent analysed yet nodes
            kivesz<-kivesz[which(kivesz$Asp==kivesz$Bsp & kivesz$An==1),"NODE"] # Terminal Paralogs
            megvizs<-setdiff(megvizs,kivesz[[1]])        # filtered out terminal paralogs
            megvizs
            kivgy<-c(kivgy,kivesz[[1]])             # monospecific nodes => terminal paralogs
            osibbek <- data.frame(lapply(osibbek, function(x) {gsub(i, NA, x)})) # filtered out the anaysed nodes
          } 
          if(length(megvizs)<1){inspecnode=F;q=q+1} # if we have no more nodes 
        } # end of inspecnode 
        kihova[which(kihova$value %in% bt$tip.label[bekerul]),"OG"]<-q
        if(length(kivgy)>0){                        # monospecific clade:
          for(jk in kivgy){                         # nearest protein will be the member of the first orthogroup
            facska<-st[names(st)==jk][[1]]          # kivagjuk az elso monospecifikus kladot es vizsgaljuk
            pathok<-sapply(1:facska$Ntip,function(g) length(nodepath(facska,facska$edge[1,1],g)))                # minimum path to tip
            kihova[which(kihova$value %in% facska$tip.label[c(1:facska$Ntip)[which.min(pathok)[1] ]] ),"OG"]<-q  # first from the minimal paths
          }
        }
        kovspec<-na.omit(as.character(osibbek[is.na(osibbek[,2]),1]))
        if(length(kovspec)==0){btwspecnode=F}      # if we have no more speciation node -> exit
      }
    } # end of btwspecnode 
    kihova$spec<-spresz(kihova$value,nametype=nametype) 
    # ordering of orthogroups based on species richness
    if(max(kihova$OG)>1){
      kihovan<-kihova[which(kihova$OG!=0),]                            
      origord<-data.frame(table(unique(kihovan[,c("OG","spec")])$OG))
      origord$Var1<-as.character(origord$Var1)
      ujord<-origord[order(origord$Freq,decreasing = T),]              # oredering of orthogroups by size (decreasing)
      ujord$Var2<-c(1:nrow(ujord))                                     # new names for orthogroups
      ujord<-rbind(ujord,data.frame(Var1=0,Freq=0,Var2=0))
      kihova$OG <- ujord[match(kihova$OG,ujord$Var1),"Var2"]
    } 
    # new names for orthogroups caused by duplications
    if(nrow(kihova[kihova$OG==0,"OG"])>0)(kihova[kihova$OG==0,"OG"]<-max(kihova$OG)+c(1:length(which(kihova$OG==0)))) 
    return(kihova)  
  }
  
  #-----------------
  setwd(fold1)
  aneve<-sapply(str_split(x,"\\."), "[", 1)
  cat(paste0("./",clusters,"/",x,"\n"))
  bt <- read.tree(paste0("./",clusters,"/",x)) #"/Cluster",v,".tree"
  cat(nametype,"\n")
  
  if(sz(setdiff(spresz(bt$tip.label,nametype=nametype),spt$tip.label))!=0){warning("unexpected species in gene tree")}
  
  et<-bt$edge                                # EdgeTable of bt
  kezdo<-length(bt$tip.label)                # number of proteins
  cat("\n",aneve,": ",kezdo)
  allnode<-c((kezdo+1):(kezdo+bt$Nnode))     # number of nodes
  unn<-unique(bt$edge[,1])                  
  
  # INDENTIFICATION OF SPECIATION / DUPLICATION NODES 
  stv<-Sys.time()
  st<-subtrees(bt) # all possible subtree from gene tree
  Sys.time()-stv
  names(st)<-allnode
  # identification of speciation and duplication ------------------------------------------------------------------------------------------------------
  dupe<-lapply(1:length(unn), function(x) dupdont(unn,x,bt,kezdo,st,nametype=nametype)) # identification of speciation and duplication
  dupe<-as_tibble(data.frame(do.call(rbind,dupe)))
  dupe<-as.tibble(sapply(1:ncol(dupe), function(v) chrelo(dupe,v)))
  dupe<-chrtoint(dupe,c(4:8))
  colnames(dupe)<-c("dupnode","Asp","Bsp","An","Bn","more","Anode","Bnode")
  dupe$dupnode<-as.logical(dupe$dupnode)
  dupe$NODE<-unn
  write.table(dupe, file=paste0('./',fold2,'/dupe',aneve,'.csv'), quote=FALSE, sep=',', row.names = F)
  
  dupnode<-unn[dupe$dupnode]
  specnode<-setdiff(allnode,dupnode)
  meret<-length(bt$tip.label)
  
  ### Ortholog coding in different cases -----------------------------------------------     
  if(length(specnode)==0){                                          # without speciation node
    kihova<-tibble(value=bt$tip.label);kihova$OG<-0
    kihova[kihova$OG==0,"OG"]<-max(kihova$OG)+c(1:length(which(kihova$OG==0)))
    kihova$spec<-spresz(kihova$value,nametype = nametype)
    OGC<-kihova;abj=F                                               
  }else{
    snfak<-st[which(names(st) %in% specnode)]                       # speciation subtrees (all)
    snfak4<-snfak            
    ancnodeok<-Ancestors(bt,as.numeric(names(snfak4)))              # Ancestors of speciation nodes
    if(length(ancnodeok)<1){                                        # without ancestor => it is the most ancient
      osibbek<-tibble('.id'=names(snfak4),'1'=NA)
      vag<-osibbek[is.na(osibbek$'1'),1]                            # only those speciation nodes which has no ancestors (start nodes)
      vagfak4<-snfak4[names(snfak4) %in% vag]
    }else{
      names(ancnodeok)<-names(snfak4)         
      osibbek<-plyr::ldply(sapply(ancnodeok, function(g) intersect(g,names(snfak4))),rbind) 
      if(ncol(osibbek)==1){osibbek$'1'<-NA}
      vag<-osibbek[is.na(osibbek$'1'),1]                            # only those speciation nodes which has no ancestors (start nodes)
      vagfak4<-snfak4[names(snfak4) %in% vag]
    }
    # Ortholog coding ---------------------------------------------------------------------------------------------------        
    OGC<-orthologcoding(dupe=dupe,vag=vag,osibbek=osibbek,kezdo=kezdo,st=st) # Ortholog coding
    abj=T
  }
  
  OGC<-dplyr::arrange(OGC,OG)
  OGC$cols<-col_vector[OGC$OG]
  OGC$cluster<-aneve;OGC<-OGC[,c(5,2,3,1,4)]
  write.table(OGC, file=paste0('./',fold2,'/OGC_',aneve,'.csv'), quote=FALSE, sep=',', row.names = F)
  # ------------------------------------------------------------------------------------- 
  
  if(abrazol){ # visualisation
    EK<-as_tibble(et);EK$cols<-"#000000"
    for(co in unique(OGC$cols)){
      EK[which.edge(bt,OGC[which(OGC$cols==co),"value"][[1]]),"cols"]<-co
    }
    bt$edge.length<-NULL
    fasize<-10*log10(meret)
    File2<-paste0("./",fold2,"/Fig_ell_",aneve,".png")
    png (file=File2,width=14,height=fasize,units="in",res=100)
    plot.phylo(bt,cex=1.2,edge.width = 6, label.offset = 6, edge.color = EK$cols)
    tiplabels(OGC[match(bt$tip.label, OGC$value),"OG"][[1]],adj = c(-1,0), bg=OGC[match(bt$tip.label, OGC$value),"cols"][[1]],cex=1.3)
    #tiplabels(tip=unlist(BEK),adj = c(-1,0), bg="green")
    if(length(snfak4)>0){nodelabels(pch=16,node=as.numeric(names(snfak4)), frame="none", col="lightblue",cex=4)}
    if(length(vagfak4)>0){nodelabels(pch=16,node=as.numeric(names(vagfak4)), frame="none", col="yellow",cex=4)}
    if(length(dupnode)>0){nodelabels(dupnode,node=dupnode,bg=NULL, frame="none", col="red",cex=1.3)}
    if(length(specnode)>0){nodelabels(specnode,node=specnode, frame="none", col="green3",cex=1.3)}
    dev.off()
  } # nezo
  
  # calculation of LOSScost for only those orthogroups which have at least 2 species
  lossok<-NULL;st1<-Sys.time()
  for(l in names(table(OGC$OG)[table(OGC$OG)>1])){
    cat(l)
    avx<-OGC[which(OGC$OG==l),]
    LC<-losscost(spt = spt,avx=avx,abrazol=abrazol) # abrazol=T
    lossok<-rbind(lossok,LC[[2]])
  }
  Sys.time()-st1
  write.table(lossok, file=paste0('./',fold2,'/lossok',aneve,'.csv'), quote=FALSE, sep=',', row.names = F)
  
  vg<-lossok[which(lossok$gl==1),]                           # true gains
  vge<-vg[!duplicated(vg$OG),]                               # true gains start of orthogroups
  vgt<-as_tibble(data.frame(table(vge$gerinc),stringsAsFactors = F)) # true gain table 
  vgt$Var1<-as.character(vgt$Var1)
  vl<-lossok[which(lossok$gl==2),]                           # ture losses
  vlt<-as_tibble(data.frame(table(vl$lossnode)))             # ture losses table
  losam<-OGC[which(OGC$OG %in% names(table(OGC$OG)[table(OGC$OG)==1])),] # count terminal duplications
  losam$sptip<-as.integer(match(losam$spec,spt$tip.label))
  sgt<-as.data.frame(table(losam$sptip))                     # singleton gain table
  
  gtermS<-losam[,c("cluster","OG","sptip")]                  
  lnek<-c("cluster","OG","lossnode")
  gbnek<-c("cluster","OG","gerinc")
  # true losses
  if(length(vl)==0){lossS<-as_tibble(data.frame(matrix(nrow=0,ncol=length(lnek))));colnames(lossS)<-lnek;lossS$cluster<-as.character(lossS$cluster)}else{lossS<-vl[,c("cluster","OG","lossnode")]}
  # true gains
  if(length(vge)==0){gbirthS<-as_tibble(data.frame(matrix(nrow=0,ncol=length(gbnek))));colnames(gbirthS)<-gbnek;gbirthS$cluster<-as.character(gbirthS$cluster)}else{gbirthS<-vge[,c("cluster","OG","gerinc")]}
  
  # merging the gain loss tables
  ee<-dplyr::full_join(gbirthS,lossS)
  eee<-dplyr::full_join(ee,gtermS)
  if(nrow(vgt)==0){eee$dupnum<-0}else{eee$dupnum<-vgt[match(eee$gerinc,vgt$Var1),2][[1]]}
  eee[!is.na(eee$sptip),"dupnum"]<-sgt[match(eee[!is.na(eee$sptip),"sptip"][[1]],sgt$Var1),"Freq"] # ahol van sptip oda beirjuk a singleton gaint
  if(nrow(vlt)>0){eee$lossnum<-vlt[match(eee$lossnode,vlt$Var1),2][[1]]}else{eee$lossnum<-0}       # number of losses
  
  # Duplication and Loss Table DLT 
  DLT<-tibble(treenode=1:(length(spt$tip.label)+spt$Nnode))                                        # DuplicationLossTable: diear column conrins the node names of the species tree
  DLT$Terminal_gains<-sgt[match(DLT$treenode,sgt$Var1),"Freq"]                                     # Terminal duplications (in tips)
  if(nrow(vgt)>0){DLT$Orthogroup_gains<-vgt[match(DLT$treenode,vgt$Var1),2][[1]]}else{DLT$Orthogroup_gains<-0}  # true gains (internal node gain) 
  if(nrow(vlt)>0){DLT$Losses<-vlt[match(DLT$treenode,vlt$Var1),2][[1]]}else{DLT$Losses<-0}         # losses
  DLT<- DLT %>% replace(is.na(.), 0)                                                               # NA to 0
  DLT$Gains<-DLT$Terminal_gains+DLT$Orthogroup_gains                                               # Gain = terminal gain + true gain
  DLT$Losses<- c(-1*DLT$Losses)                                                                    # Loss is a negative number
  DLT$Net_gains<-DLT$Gains+DLT$Losses                                                              # Net_Gain = Gain + Loss

  # calculation of Copy numbers
  if(nrow(vgt)>0){
    clelsonode<-min(as.numeric(vgt$Var1))                        # birth of the cluster will be the first true gain
    relnodes<-Descendants(spt,clelsonode,type ="all")            # relevant nodes
    relcp<-as_tibble(cbind(relnodes,sapply(relnodes, function(s) # adding Net_gains throughout the relevant nodes
      sum(DLT[ which(DLT$treenode %in% nodepath(spt,clelsonode,s) ) , "Net_gains"]))))
    DLT$Copynum<-relcp[match(DLT$treenode,relcp$relnodes),2][[1]]
    DLT<- DLT %>% replace(is.na(.), 0);DLT[which(DLT$treenode==clelsonode),"Copynum"]<-DLT[which(DLT$treenode==clelsonode),"Net_gains"]
  }else{DLT$Copynum<-DLT$Net_gains}
  DLT$Cluster<-aneve
  write.table(DLT, file=paste0('./',fold2,'/DLT',aneve,'.csv'), quote=FALSE, sep=',', row.names = F)
  save(list = ls(),file = paste0("./",fold2,"/",aneve,".saved"))
  
  # mapping
  if(F){ # visualisation
    tree<-spt;tree$edge.length <- NULL;size<-3
    for(h in 4:7){
      mit<-DLT[,h][[1]]
      mi<-colnames(DLT)[h]
      
      File2<-paste0("./",fold2,"/",aneve,"_",mi,".png")
      png (file=File2,width=12,height=18,units="in",res=100)
      plot.phylo(tree, main=paste0(aneve," ",mi), label.offset = 0.2, node.depth = 2, show.tip.label = T, cex = 1, edge.width=2)
      nodelabels(pch = 16, col =c("blue","red")[(mit<0*1)+1], node = DLT$treenode, cex = abs(mit) * size)
      nodelabels(text = mit, node = DLT$treenode, adj = c(0.5,0.5), bg = NULL, frame = "none", cex = 1, col="white")
      dev.off()
    }
  } # end of mapping
  return(DLT)
} # end of function

# single species clusters
SP1comp<-function(x,spt=spt,fold2=fold2){
  DLT<-tibble(treenode=c(1:(length(spt$tip.label)+spt$Nnode)),Terminal_gains=0,Orthogroup_gains=0,Losses=0,Gains=0,Net_gains=0,Copynum=0)
  aneve<-gsub("mcl_","Cluster",unique(x$mcl))
  DLT$Cluster<-aneve
  DLT[which(spt$tip.label==unique(x$spec)),"Terminal_gains"]<-nrow(x)
  DLT[which(spt$tip.label==unique(x$spec)),"Gains"]<-nrow(x)
  DLT[which(spt$tip.label==unique(x$spec)),"Net_gains"]<-nrow(x)
  DLT[which(spt$tip.label==unique(x$spec)),"Copynum"]<-nrow(x)
  write.table(DLT, file=paste0('./',fold2,'/DLT',aneve,'.csv'), quote=FALSE, sep=',', row.names = F)
  save(list = c("DLT","aneve","spt","fold2"),file = paste0("./",fold2,"/",aneve,".saved"))
  return(DLT)
}
# >1 species but <4 species clusters
SPTcomp<-function(x,spt=spt,col_vector=col_vector,clusters=clusters,fold2=fold2,fold1=fold1,nametype=nametype,abrazol=F){
  OGC<-NULL
  x$OG<-0
  got3 <- data.table(x)
  got3<-got3[ , Index := 1:.N , by = c("spec")]
  got3<-unite(got3,Index,sep="_",col=indid)
  kihova<-dplyr::arrange(got3,indid)
  kihova[which(kihova$indid==1),"OG"]<-1
  if(length(kihova[kihova$OG==0,"OG"])>0){kihova[kihova$OG==0,"OG"]<-max(kihova$OG)+c(1:length(which(kihova$OG==0)))}
  aneve<-gsub("mcl_","Cluster",unique(x$mcl))
  kihova$mcl<-aneve;colnames(kihova)[c(1,3)]<-c("cluster","value")
  OGC<-kihova[,c(1,4,2,3,5)]
  write.table(OGC, file=paste0('./',fold2,'/OGC_',aneve,'.csv'), quote=FALSE, sep=',', row.names = F)
  # ------------------------------------------------------------------------------------- 
  
  # calculation of LOSScost  
  lossok<-NULL;st1<-Sys.time()
  for(l in names(table(OGC$OG)[table(OGC$OG)>1])){
    cat(l)
    avx<-OGC[which(OGC$OG==l),]
    LC<-losscost(spt = spt,avx=avx,abrazol=abrazol) #abrazol=T
    lossok<-rbind(lossok,LC[[2]])
  }
  Sys.time()-st1
  write.table(lossok, file=paste0('./',fold2,'/lossok',aneve,'.csv'), quote=FALSE, sep=',', row.names = F)
  
  vg<-lossok[which(lossok$gl==1),]                          
  vge<-vg[!duplicated(vg$OG),]                               
  vgt<-as_tibble(data.frame(table(vge$gerinc),stringsAsFactors = F)) 
  vgt$Var1<-as.character(vgt$Var1)
  vl<-lossok[which(lossok$gl==2),]                           
  vlt<-as_tibble(data.frame(table(vl$lossnode)))          
  losam<-OGC[which(OGC$OG %in% names(table(OGC$OG)[table(OGC$OG)==1])),] 
  losam$sptip<-as.integer(match(losam$spec,spt$tip.label))
  sgt<-as.data.frame(table(losam$sptip))                     
  
  gtermS<-losam[,c("cluster","OG","sptip")]                 
  lnek<-c("cluster","OG","lossnode")
  gbnek<-c("cluster","OG","gerinc")

  if(length(vl)==0){lossS<-as_tibble(data.frame(matrix(nrow=0,ncol=length(lnek))));colnames(lossS)<-lnek;lossS$cluster<-as.character(lossS$cluster)}else{lossS<-vl[,c("cluster","OG","lossnode")]}

  if(length(vge)==0){gbirthS<-as_tibble(data.frame(matrix(nrow=0,ncol=length(gbnek))));colnames(gbirthS)<-gbnek;gbirthS$cluster<-as.character(gbirthS$cluster)}else{gbirthS<-vge[,c("cluster","OG","gerinc")]}
  

  ee<-dplyr::full_join(gbirthS,lossS)
  eee<-dplyr::full_join(ee,gtermS)
  if(nrow(vgt)==0){eee$dupnum<-0}else{eee$dupnum<-vgt[match(eee$gerinc,vgt$Var1),2][[1]]}
  eee[!is.na(eee$sptip),"dupnum"]<-sgt[match(eee[!is.na(eee$sptip),"sptip"][[1]],sgt$Var1),"Freq"] 
  if(nrow(vlt)>0){eee$lossnum<-vlt[match(eee$lossnode,vlt$Var1),2][[1]]}else{eee$lossnum<-0}       
  
  # Duplication Loss Tabla DLT elkeszitese
  DLT<-tibble(treenode=1:(length(spt$tip.label)+spt$Nnode))                                        
  DLT$Terminal_gains<-sgt[match(DLT$treenode,sgt$Var1),"Freq"]                                     
  if(nrow(vgt)>0){DLT$Orthogroup_gains<-vgt[match(DLT$treenode,vgt$Var1),2][[1]]}else{DLT$Orthogroup_gains<-0}  
  if(nrow(vlt)>0){DLT$Losses<-vlt[match(DLT$treenode,vlt$Var1),2][[1]]}else{DLT$Losses<-0}       
  DLT<- DLT %>% replace(is.na(.), 0)                                                              
  DLT$Gains<-DLT$Terminal_gains+DLT$Orthogroup_gains                                            
  DLT$Losses<- c(-1*DLT$Losses)                                                              
  DLT$Net_gains<-DLT$Gains+DLT$Losses                                                              # Net_Gain = Gain + Loss
   
  #Calculation of Copy numbers
  if(nrow(vgt)>0){
    clelsonode<-min(as.numeric(vgt$Var1)) 
    relnodes<-Descendants(spt,clelsonode,type ="all") 
    relcp<-as_tibble(cbind(relnodes,sapply(relnodes, function(s)   
      sum(DLT[ which(DLT$treenode %in% nodepath(spt,clelsonode,s) ) , "Net_gains"]))))
    DLT$Copynum<-relcp[match(DLT$treenode,relcp$relnodes),2][[1]]
    DLT<- DLT %>% replace(is.na(.), 0);DLT[which(DLT$treenode==clelsonode),"Copynum"]<-DLT[which(DLT$treenode==clelsonode),"Net_gains"]
  }else{DLT$Copynum<-DLT$Net_gains}
  DLT$Cluster<-aneve
  write.table(DLT, file=paste0('./',fold2,'/DLT',aneve,'.csv'), quote=FALSE, sep=',', row.names = F)
  save(list = ls(),file = paste0("./",fold2,"/",aneve,".saved"))
  
  # mapping
  if(F){ 
    tree<-spt;tree$edge.length <- NULL;size<-3
    for(h in 4:7){
      mit<-DLT[,h][[1]]
      mi<-colnames(DLT)[h]
      
      File2<-paste0("./",fold2,"/",aneve,"_",mi,".png")
      png (file=File2,width=12,height=18,units="in",res=100)
      plot.phylo(tree, main=paste0(aneve," ",mi), label.offset = 0.2, node.depth = 2, show.tip.label = T, cex = 1, edge.width=2)
      nodelabels(pch = 16, col =c("blue","red")[(mit<0*1)+1], node = DLT$treenode, cex = abs(mit) * size)
       nodelabels(text = mit, node = DLT$treenode, adj = c(0.5,0.5), bg = NULL, frame = "none", cex = 1, col="white")
        dev.off()
    }
  } # end of mapping
  return(DLT)
} # end of function

# intro for small clusters-----------------------------------------------------------
# setwd("D:/3_Laska/PL109/compaRe_PL109/")
# cls<-read_csv("cls_PL109m12CQH8020NQ_I2.csv")
# 
# #elovalogato: lehetne az alapjan h mi nem fert be a fak koze v csak a V4SP komplementere
# cls<-dplyr::arrange(cls,mcl,spec)
# clfr<-as_tibble(as.data.frame.matrix(table(cls$mcl,cls$spec)))
# clfrt<-clfr
# clfrt$mcl<-unique(cls$mcl)
# clfrt$protnum<-rowSums(clfr)
# clfrt$spnum<-rowSums((clfr>0)*1)
# smallcl<-clfrt[which(clfrt$protnum < 4 | clfrt$spnum < 2),] # kevesebb mint 4 feherje VAGY kevesebb mint 2 faj clusterek (ami nem V4SP)
# write.table(smallcl, file=paste0('smallclek',SD,'.csv'), quote=FALSE, sep=',', col.names = T, row.names = F)
# 
#   clsm<-cls[which(cls$mcl %in% smallcl$mcl),]
# table(clsm$spec) #ell nem csak 1 fajunk van
# 
# cls_1SP<-cls[which(cls$mcl %in% smallcl[which(smallcl$spnum==1),"mcl"][[1]]),]
# cls_TSP<-cls[which(cls$mcl %in% smallcl[which(smallcl$spnum>1),"mcl"][[1]]),]
# 
# ir <- cls_1SP %>% group_by(mcl)
# cldb<-group_split(ir)

#x<-"D:/3_Laska/PL109/compaRe_PL109/Cluster2264.saved"
DiffMRCA2<-function(x){
  load(x)
  ogc<-OGC
  LDUP<-function(y){
    sor<-unlist(lapply(nevkat, function(h) any(h==y))) # all duplications from the birth of the cluster till the proteins
    asp<-spresz(y,"JGI")
    h<-nevkat[[1]]
    sormas<-unlist(lapply(nevkat, function(h) sz(which(spresz(h,"JGI")==asp))>1 )) # duplication where a species remained as duplicated in the recent species
    sor<-as.numeric(names(which(sor)))
    sormas<-as.numeric(names(which(sormas)))
    sor<-unique( c(sor[[1]], intersect(sor,sormas)))
    sorp<-paste0(sor,collapse = "|")
    lastdup<-max(sor)
    return(c(sorp,lastdup))
  }
  
  if(nrow(ogc)<4){
    ir <- ogc %>% group_by(OG)
    ogdb<-group_split(ir)
    length(ogdb)  
    mrcak<-sapply(ogdb,function(v) getMRCA(spt,unique(v$spec)))
    mrcak<-sapply(mrcak,function(v) ifelse(is.null(v),0,v))
    mrcak<-data.frame(cbind(unique(ogc$OG),unlist(mrcak)))
    ogc$mrca<-mrcak[match(ogc$OG,mrcak$X1),"X2"]
    kerd<-filter(ogc,ogc$mrca==0)
    ogc[which(ogc$mrca==0),"mrca"]<-match(kerd$spec,spt$tip.label)
    colnames(ogc)[4]<-"value"
    ogc$cols<-"#FFFFFF"
    ogc$GTlastdup<-NA
    if(table(ogc$spec)>1){
      ogc[ogc$spec == names(which(table(ogc$spec)>1)),"STlastdup"]<- min(ogc[ogc$spec == names(which(table(ogc$spec)>1)),"mrca"])
      ogc[is.na(ogc$STlastdup),"STlastdup"]<-ogc[is.na(ogc$STlastdup),"mrca"]
    }else{
      ogc$STlastdup<-ogc$mrca
    }
    ogc$minmrca<-ogc$mrca
    ogc$maxmrca<-ogc$mrca
    ogc$ancOG<-1
    
  }else{
    
    ir <- ogc %>% group_by(OG)       # split the clusters by OGs (orthogroups)
    ogdb<-group_split(ir)
    length(ogdb)
    
    mrcak<-sapply(ogdb,function(v) getMRCA(spt,unique(v$spec))) # identification of the MRCAs of the OGs (only >1 species)
    mrcak<-sapply(mrcak,function(v) ifelse(is.null(v),0,v))
    mrcak<-data.frame(cbind(unique(ogc$OG),unlist(mrcak)))
    ogc$mrca<-mrcak[match(ogc$OG,mrcak$X1),"X2"]
    kerd<-filter(ogc,ogc$mrca==0)
    ogc[which(ogc$mrca==0),"mrca"]<-match(kerd$spec,spt$tip.label) # MRCA of the singleton OGs
    
    if(length(dupnode)==0 | sz(ogc$spec)==1){ # if we have not duplication in the cluster
      
      ogc[,"minmrca"]<-ogc[,"mrca"]
      ogc[,"maxmrca"]<-ogc[,"mrca"]
      ogc[,"ancOG"]<-ogc[,"OG"]
      
    }else{                  # If we have duplication node in the cluster
      
      potmrca<-unique(c(dupnode,min(dupnode,specnode)))                       # all dupnode and the smallest specnode
      dupfak<-st[which(names(st) %in% potmrca)]                
      nevkat<-lapply(dupfak,function(x) x$tip.label)                          # portein names
      ogbes<-lapply(nevkat,function(x) ogc[match(x,ogc$value),"OG"][[1]])     # OG membership
      mrcabes<-lapply(nevkat,function(x) ogc[match(x,ogc$value),"mrca"][[1]]) # MRCA based on OG membership
    
      LDE<-lapply(ogc$value,function(y) LDUP(y))
      ogc<-bind_cols(ogc,data.frame(do.call(rbind,LDE)))
      ogc$X2<-as.numeric(as.character(ogc$X2))
      colnames(ogc)[c(7,8)]<-c("dupsor","GTlastdup")

      LDMRCA<-function(g){
        spku<-unique(spresz(g,"JGI"))
        LM<-getMRCA(spt, spku)
        if(is.null(LM)){
          LM<-match(spku,spt$tip.label)
        }
        return(LM)
      }
      
      DupMRCAk<-lapply(nevkat,function(g) LDMRCA(g))
      DupMRCAk<-data.frame(do.call(rbind,DupMRCAk))
      DupMRCAk$GTnode<-rownames(DupMRCAk)
      colnames(DupMRCAk)<-c("STnode","GTnode")
      ogc$STlastdup<-DupMRCAk[match(ogc$GTlastdup,DupMRCAk$GTnode),"STnode"] 

      kihol<-lapply(vagfak4,function(x) x$tip.label)
      for(ii in 1:length(kihol)){
        ogc[which(ogc$value %in% kihol[[ii]]),"GTspnode"]<-as.numeric(names(kihol)[ii]) #First Spec node
      }
      
      ogc$minmrca<-as.numeric(NA) 
      ogc$maxmrca<-as.numeric(NA) 
      
      for(x in 1:length(nevkat)){ 
        dupnnev<-as.numeric(names(nevkat)[x])
        a<-nevkat[[x]]
        b<-ogbes[[x]]
        c<-mrcabes[[x]]
          if(sz(spresz(a,"JGI"))==1 & sz(c)>1){                     # terminal duplications 1 species more duplications
          cm<-ifelse(c>length(spt$tip.label),c,(c+nrow(spt$edge)))
          ogc[match(a,ogc$value),"minmrca"]<-c[which.min(cm)]
          ogc[match(a,ogc$value),"maxmrca"]<-c[which.max(cm)]
          ogc[match(a,ogc$value),"ancOG"]<-min(b)
        }else{                                
          ogc[match(a,ogc$value),"minmrca"]<-NA
          ogc[match(a,ogc$value),"maxmrca"]<-NA
          ogc[which(ogc$value %in% a & ogc$GTspnode<dupnnev),"ancOG"]<-min(b) # dupnode must be higher than the MRCA of the first OG 
        }
      }
      
      ogc[is.na(ogc$minmrca),"minmrca"]<-ogc[is.na(ogc$minmrca),"mrca"]
      ogc[is.na(ogc$maxmrca),"maxmrca"]<-ogc[is.na(ogc$maxmrca),"mrca"]
      ogc[is.na(ogc$ancOG),"ancOG"]<-ogc[is.na(ogc$ancOG),"OG"]
      
    }
  }
  
  write.table(ogc, file=paste0('mrcaogc_',aneve,'.csv'), quote=FALSE, sep=',', col.names = T, row.names = F)
  return(ogc)
}