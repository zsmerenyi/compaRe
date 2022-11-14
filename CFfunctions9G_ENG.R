while(F){
  abrazol=T
  covnorm=T
  spnorm=T
  sajat=T
  covnorm=sc
  spnorm=sn
  sajat=sa
  y=vclsizes2$ALAP[11000]
  y="allMillcstl/mcl_12011M.fas"

 
  bt<-read.tree("species.tree")
  bt<-LadderizeTree(bt)
  
}
  
# Amino Acid Distance calculation of an aligned AAset
AADIS<-function(y){
  alio<-Biostrings::readAAMultipleAlignment(y)
  dat <- phangorn::as.phyDat(alio, type="AA")
  AAdist<-dist.ml(dat,model = "WAG",exclude = "none")
  AAdist<-as.matrix(AAdist)
  AAdist[AAdist==1e-8]<-10       # sequences without overlap 
  diag(AAdist)<-NA
  meanAAD<-mean(AAdist,na.rm=T)
  return(meanAAD)
}
  

filter<-function(y,abrazol=F,covnorm=T,spnorm=T,sajat=T){
    specnorm<-function(x,spk.=spk){  # AAdist specnorm fuggvenye fajtav.=buscotav2
      tab<-fajtav[fajtav[,1]==x,]
      ert<-tab[match(spk,tab[,2][[1]]),"norm"]
      return(ert)
    }
  
    add<-a1<-paste0(fold,"CFC9",substr(sajat,1,1),substr(covnorm,1,1),substr(spnorm,1,1))                                    
    sall<-Sys.time() 
    nevresz<-str_locate(y,"mcl");nevresz2<-str_locate(y,"M.fas")-1;clnev<-str_sub(y,nevresz[1],nevresz2[1]) # name of clusters
    cat("\n\n",clnev,"\n")
    alio<-readAAStringSet(y)                               # original AA alignment
    names(alio)<-sapply(str_split(names(alio)," "), "[", 1) ######################################
    SumTo<-tibble::enframe(names(alio), name=NULL)         # summary table
    SumTo$seqlen<-width(alio)-alphabetFrequency(alio)[,28] # length of seqs
    fejlec<-setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("no","fragment","covAA","AADist"))
    SumTl<-SumTo
    SumTl$spec<-spresz(SumTl$value)
    SumTl$mcl<-clnev
    
    SumTo<-mrcadec(SumTl,bt)
    
    #------------------ < 3 proteins
  # if(length(alio)<=3){    # min 3feh
  #    fejlec[1,"all"]<-length(alio);fejlec[1,"comment"]<-"too small";fejlec[1,"tip"]<-add;fejlec[1,"mcl"]<-clnev
  #    cat("too small");kimenet<-fejlec[,c(ncol(fejlec),1:(ncol(fejlec)-1))]
  #    write.table(SumTo, file=paste0('SumT_',clnev,"_",add,'.csv'), quote=FALSE, sep=',', col.names = T, row.names = F)
  #    kimenet<- kimenet %>% replace(is.na(.), 0);return(kimenet)
  #  }else{
      #----------------- >= 3 proteins in the cluster
      alim<-alio                              
    
    # COVERAGE distancematrix 
      covdist<-covszr(ill=alim)                 # Coverage Distance
      covdisto<-covdist                         # original covdist
      mean_cov<-1-mean(covdisto,trim = 0.2)     # average of pairwise coverage values among all protein pairs of the cluster
      cov_prot<-1-rowMeans(covdist)             # average of pairwise coverage values of each and all other proteins in the cluster
      SumTo$AvgCov<-cov_prot
      SumTo$group<-"core"
      suspic<-names(cov_prot[which(mean_cov/cov_prot >= 2 )]) # suspicious proteins with low coverage
      SumTo[which(SumTo$value %in% suspic),"group"]<-"suspic" # fragmented proteins
      core<-SumTo[which(SumTo$group=="core"),"value"][[1]]
      SumTo$noise<-ifelse(SumTo$group=="core","no","fragment") # fragmented prots which would change to covAA or AAdist
      table(SumTo$noise)
      
    # AA distancematrix
      alio<-Biostrings::readAAMultipleAlignment(y)
      dat <- phangorn::as.phyDat(alio, type="AA")
      AAdist<-dist.ml(dat,model = "WAG",exclude = "none")
      AAdist<-as.matrix(AAdist)
      AAdist<-AAdist[match(SumTo$value,rownames(AAdist)),match(SumTo$value,colnames(AAdist))] #same order than SumTo
      AAdist[AAdist==1e-8]<-10       
      AAdistorig<-AAdist
      SumTo$AvgAAdist<-apply(AAdistorig, 1, function(x) mean(x))
      SumTo$value
      
    # Coverage fragment/covAA with original AAdist
      if(length(suspic)>1 & length(core)>1){
        ineqM<-AAdist[which(rownames(AAdist) %in% suspic),which(colnames(AAdist) %in% core)] # inequal AAdist Matrix between suspicious and core
        minAADistie<-apply(ineqM, 1, function(x) min(x))
        SumTo[SumTo$value %in% names(which(minAADistie/mean(AAdist,na.rm=T) > 2 )),"noise"]<-"covAA"        # noise
        
      }else{
        if(length(suspic)==1){
          ineqM<-AAdist[which(rownames(AAdist) %in% suspic),which(colnames(AAdist) %in% core)] # inequal AAdist Matrix between suspicious and core 
          # a core atlagtol kell 2x rosszabbnak lennie, ez konzervalt feherjeknel szigoru lesz, masoknal lazabb
          coreavg<-mean(AAdist[which(colnames(AAdist) %in% core),which(colnames(AAdist) %in% core)],na.rm = T)
          minAADistie<-ineqM[which.min(ineqM)]
          zaj<-ifelse(minAADistie/coreavg > 2,suspic,0)
          SumTo[SumTo$value %in% zaj,"noise"]<-"covAA"        # noise
        }
      }
      table(SumTo$noise)
      
      
    # AADist NORMALISATION-------------
      #######
      marcsak<-SumTo[SumTo$noise=="no","value"][[1]]
      AAdist<-AAdist[which(rownames(AAdist) %in% marcsak),which(colnames(AAdist) %in% marcsak)]
      covdist<-covdist[which(rownames(covdist) %in% marcsak),which(colnames(covdist) %in% marcsak)]
      
      #######
      if(covnorm){AAdist<-AAdist*(covdist*covdist)}         # coverage normalisated AAdist
      rn<-rownames(AAdist)
       
      if(spnorm){
        spk<-spresz(rownames(AAdist))
        fajtav=buscotav2
        normspec<-sapply(spk,function(x) 1-specnorm(x,spk=spk)) # 1- similarity
        normspec<-data.frame(plyr::ldply(normspec,rbind))
        normspec<-normspec[,-1]
        normspec<-normspec %>% mutate_all(as.numeric)
        AAdistm<-AAdist*normspec
        AAdist<-AAdistm  
        rownames(AAdist)<-colnames(AAdist)<-rn
        AAdist<-as.matrix(AAdist)
      }
      
      diag(AAdist)<-NA 
      if(sajat){    
        # 1
        avgAAdist<-mean(AAdist,na.rm=T,trim=0.1)  #20% trimm => upper 20% of the dataset
        protAAdist<-apply(AAdist, 1, function(x) mean(x,na.rm=T))
        SumTo[SumTo$value %in% names(which(protAAdist/avgAAdist >= 2 & protAAdist > 0.2)),"noise"]<-"AADist"  # THRESHOLD 0.2-------
      }else{
        # 2
        #medAA<-mean(AAdist,na.rm=T,trim=0.5)  
        medAA<-median(AAdist,na.rm=T)  
        
        medprotAA<-apply(AAdist, 1, function(x) mean(x,na.rm=T))
        protIQR<-(medprotAA-medAA)/IQR(AAdist,na.rm=T)
        SumTo[SumTo$value %in% names(which(protIQR >= 2)),"noise"]<-"AADist"
      }
      
      table(SumTo$noise)
      # Filter out those ones which has already erased based on AAdist, from the core
      core<-intersect(core,SumTo[SumTo$noise=="no","value"][[1]])
      add<-a1 
      
      # Final similarities
      AAdist<-AAdistorig
      if( length(which(SumTo$noise=="no"))>1 ){
        coredist<-AAdist[,which(colnames(AAdist) %in% SumTo[which(SumTo$group=="core" & SumTo$noise=="no"),"value"][[1]]) ]
        minAAdistc<-apply(coredist, 1, function(x) min(x,na.rm=T))
        minAAdistcnev<-apply(coredist, 1, function(x) which.min(x))
        SumTo$coreminAAd<-minAAdistc
      }else{
        if(length(which(SumTo$noise=="no"))==1){
          coredist<-AAdist[,which(colnames(AAdist) %in% SumTo[which(SumTo$group=="core" & SumTo$noise=="no"),"value"][[1]]) ]
          SumTo$coreminAAd<-coredist
          }else{SumTo$coreminAAd<-1}
        }
      SumTo$mcl<-clnev
      SumTov<-SumTo[,c(4,1,3,2,5:(ncol(SumTo)))]
      return(SumTov)
    
  } # End of function

mrcadec<-function(x,bt){
  avx<-x
  nev<-unique(avx$mcl)[1]
  exist<- which(bt$tip.label %in% avx$spec)                   # species in the cluster
  emrcaorig<-getMRCA(bt,exist)                                # cluster MRCA
  
  mrcavalt<-function(x){
    avxakt<-avx[-x,]
    exist<- which(bt$tip.label %in% avxakt$spec)              # species in the cluster
    emrca<-getMRCA(bt,exist)
    if(is.null(emrca)){emrca<-exist}
    return(emrca)# cluster MRCA
  }
  
  protmrca<-sapply(1:nrow(avx),function(x) mrcavalt(x))
  pathok<-sapply(1:nrow(avx),function(g) length(nodepath(bt,emrcaorig,protmrca[g])) )
  
  avx$Omrca<-emrcaorig 
  avx$protmrca<-protmrca 
  avx$Dmrca<-pathok-1 
  return(avx)
}  


lossdb<-function(x,bt,abrazol){
  avx<-x
  nev<-unique(avx$mcl)[1]
  exist<- which(bt$tip.label %in% avx$spec)                   # species in the cluster
  lossok<- c(1:length(bt$tip.label))[-exist]                  # missing species
  kezd<-length(bt$tip.label)+1                                # first node
  emrca<-getMRCA(bt,exist)                                    # MRCA of the cluster
  EK<-as_tibble(bt$edge);EK$erd<-1                            # edge table
  elotte<-which(EK<emrca & EK>=kezd)
  if(length(elotte)==0){
    EK$erd<-2
  }else{
    EK[-which(EK<emrca & EK>=kezd),"erd"]<-2                  # edges before the cluster MRCA (non relevant)
  }
  kiv<-which.edge(bt,exist)                                   # edges between the members of the cluster
  EK[kiv,"erd"]<-EK[kiv,"erd"]+1                              # add +1 if a species is member of the clustera => 2->3 => member=3 loss=2
  # if erd = 2 =>loss erd=3 =>presence of a gene
  ecol<-c("darkgrey","blue","red")[EK$erd];table(ecol)        # colouring darkgrey=1 blue=2(loss) red=3(cluster membership)
  gerinc<-unique(EK[kiv,1][[1]])                              # backbone among cluster memberships "gerinc"
  gl<-sapply(gerinc,function(g) sz(EK[EK$V1==g,"erd"][[1]]))  # Along the backbone 33 (continued) or 23 (loss) 
  gcol<-c("red","blue")[gl];table(gcol)                       # colouring of edges red=cluster membership, blue=loss
  err<-as.data.frame(table(gl))
  err$mcl<-nev
  err$spnum<-sz(avx$spec)
  err$protnum<-nrow(avx)
  err$kezdo<-kezd
  mmo<-as_tibble(cbind(gerinc,gl))                            # backbone nodes loss(2) non loss(1)
  
  while(F){
    plot.phylo(bt,cex=0.8,edge.width = 1, label.offset = 0.2, edge.color = ecol)
    tiplabels(pch=1,tip = exist, col="red", cex=1)
    tiplabels(pch=1,tip = lossok, col="blue", cex=1)
    nodelabels(pch=16,node=gerinc,cex=1, col=gcol)
  }
  
  if(abrazol){  
    File2<-paste0("lossFig_",nev,".pdf")
    pdf(file = File2, width=16,height=16, pointsize = 18, family = "sans", bg = "white")
    plot.phylo(bt,cex=1.5,edge.width = 2, label.offset = 0.05, edge.color = ecol)
    tiplabels(pch=1,tip = exist, col="red", cex=2)
    tiplabels(pch=1,tip = lossok, col="blue", cex=2)
    nodelabels(pch=16,node=gerinc,cex=2, col=gcol)
    dev.off()
  }
  return(err)
}


sz<-function(vmi){length(unique(vmi))}
 
  
  covszr<-function(ill){  
    alib<-ill  
    alibs<-as.character(alib)
    alibs0<-str_replace_all(alibs,"-","0")
    alibs0<-str_replace_all(alibs0,"\\*","0")
    alibs1<-str_replace_all(alibs0,"[:alpha:]","1")
    alibst<-str_split(alibs1,"|")
    alib3<-do.call(rbind,alibst)
    adt<-as.data.frame(alib3)
    adt2<-adt[,c(2:(ncol(adt)-1))]
    rownames(adt2)<-names(alib)
    #st<-Sys.time() 
    adt3<-as_tibble(adt2)
    adt4<-as_tibble(sapply(1:ncol(adt3), function(v) chrelo(adt3,v)))
    adt4<- adt4 %>% replace(is.na(.), 0)
    B <- Matrix(as.matrix(adt4), sparse = TRUE)
    allsite<-dist.matrix(B, method = "minkowski", p=0)
    diffsite<-dist.matrix(B, method = "minkowski", p=1) 
    covdist<- (1-((allsite-diffsite)/allsite))
    colnames(covdist)<-names(alib)
    rownames(covdist)<-names(alib)
    #cat("\t",Sys.time()-st,"\n") 
    return(covdist)
  }
 
  chrelo<-function(TABLE,x){TABLE[,x]<-as.numeric(as.character(TABLE[,x][[1]]))}
  
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
  
  MMO<-function(bt,mmo,avx){       
    if(nrow(mmo)>1){
      exist<- which(bt$tip.label %in% avx$spec)                                                       # species in the cluster
      allc<-data.frame(combn(exist,2))                                                                # all species combination (half matrix)
      utak<-lapply(allc,function(hh) which.edge(bt,hh))                                               # paths btw species combinations
      lossutak<-lapply(utak,function(hl) table(mmo[which(mmo$gerinc %in% unique(bt$edge[hl])),"gl"])) # species pairs: nothing(1)/loss(2) crosstable
      lossdbok<-lapply(lossutak,function(hx) hx[names(hx)==2]) 
      lossdbok<-plyr::ldply(lossdbok,rbind)                                                           # number of losses for each species pair 
      mmer<-cbind(t(allc),lossdbok);mmer$A<-bt$tip.label[mmer$'1'];mmer$B<-bt$tip.label[mmer$'2']     # summary table about species, losses
      colnames(mmer)<-c("An","Bn","id","loss","Asp","Bsp");mmer[is.na(mmer$loss),"loss"]<-0           # NA to 0
      mmer2<-mmer[,c(1:4,6,5)];colnames(mmer2)[c(5,6)]<-c("Asp","Bsp");mmer3<-rbind(mmer,mmer2)       # duplication of the table
      splosscont<-plyr::ddply(mmer3,c("Asp"),summarise,mean=mean(loss),sd=sd(loss), count=sum(loss))  # avg sd sum of losses for each species
      splosscont$msd<-splosscont$mean/splosscont$sd                                                   # 1/RSD (relative standard deviation) 
      tqu<-quantile(mmer3$loss)[4]                                                                    # qantile (75%) 
      outl<-table(unlist(mmer3[which(mmer3$loss>tqu),c("Asp","Bsp")]))                                
      splosscont$outl<-outl[match(splosscont$Asp,names(outl))]
    }else{splosscont<-tibble( "Asp"=0,"mean"=0,"sd"=0,"count"=0,"msd"=0)}
    return(splosscont)
  }

  # species names
  spresz<-function(s){
    s2<-str_split(s,"_")
    s3<-sapply(1:length(s2), function(x) paste0(s2[[x]][1:(length(s2[[x]])-1)],collapse = "_") )
    return(s3)
  }
  
  chrtonum<-function(tabla,oszl=c(1:1)){
    orignames<-colnames(tabla)
    for(h in oszl){
      tabla[,h]<-as.numeric(tabla[,h][[1]])
    }
    colnames(tabla)<-orignames
    return(tabla)
  }
  
  