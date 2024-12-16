library(RNetCDF)
library(glmnet)
library(car) #vif
library(igraph)
library(FactoMineR)
library(mgcv)


DRLPC=function(onedata, oneSNPinfo, vifcut=20, CLQcut=0.5, pccut=0.8, Yvar="PHENOTYPE", itermax=1000, Klim=300){
  
  SNPsinLPC=NULL
  tagvarlist=NULL
  group=NULL
  grpnames=NULL
  lowr2ind=NULL
  RPCind=NULL
  
  #1. aliased SNP removal
  allSNPs=setdiff(colnames(onedata),Yvar)   
  genoup=onedata[,allSNPs]
  remainSNPs=allSNPs
  rmvSNPs=NULL
  rmvSNPs2=NULL
  
  
  K=length(allSNPs)

  if(K<Klim){
    
    alremout=alias.remove(onedata, allSNPs, Yvar) 
    remainSNPs=alremout$remainSNPs
    rmvSNPs=alremout$rmvSNPs
    
  }
  
  if(K >= Klim){
    
    #alias removal step 1 : partition and removal of aliased SNPs in each partition
    numset=floor((K-1)/Klim)+1
    alremSNPlist=vector("list",numset)
    alrmvSNPlist=vector("list",numset)
    alremsizes=rep(NA,numset)
    for(k in 1:numset){
      start=Klim*(k-1)+1
      end=K
      if(k<numset) end=Klim*k
      alremout=alias.remove(onedata,allSNPs[start:end],Yvar)
      alremSNPlist[[k]]=alremout$remainSNPs
      if(length(alremout$rmvSNPs)>0) alrmvSNPlist[[k]]=alremout$rmvSNPs
      alremsizes[k]=length(alremout$remainSNPs)
    }
    
    remSNPin=alremSNPlist
    rmvSNPin=alrmvSNPlist
    remsizesin=alremsizes
    
    #repeat removal cycle until no change (alias removal step 6) 
    cycle=1
    k=0
    while(cycle>0){
      k=k+1
      cyout=removalcycle(onedata,remSNPin, rmvSNPin, remsizesin,Klim, Yvar) #alias removal step 2~5
      remSNPin=cyout$remSNPlist
      rmvSNPin=cyout$rmvSNPlist
      remsizesin=cyout$remsizes
      cycle=cyout$change
    }
    
    remainSNPs=unlist(remSNPin)
    rmvSNPs=unlist(rmvSNPin)
    #sorting
    remainSNPs=allSNPs[sort(match(remainSNPs, allSNPs))]
    rmvSNPs=allSNPs[sort(match(rmvSNPs, allSNPs))]
    
    #print(length(remainSNPs))
    
  }
  
  
  if(length(remainSNPs)==1){
    removedSNPs=rmvSNPs
    dataup=data.frame(onedata[,c(Yvar,remainSNPs)]) 
    return(list(vdata=dataup,LPCinfo=SNPsinLPC, taglist=tagvarlist, aliasremoved=rmvSNPs, removed=removedSNPs, group=group, remainedgrps=grpnames, RPCind=lowr2ind))
    
    
  }
  
  #2.Using CLQD using a CLQcut, find clique clusters
  genoin=onedata[,remainSNPs]
  varnames=colnames(genoin)
  
  CLQout=CLQD(genoin, oneSNPinfo, CLQcut=CLQcut, clstgap=40000, hrstType="fast", hrstParam=200,CLQmode="maximal", LD="r2")
 
  
  
  #CLQout: NA is singleton,
  numclst=0
  if(!all(is.na(CLQout)))   numclst=max(CLQout,na.rm=T) 
  
  #3. Replace each clique clusters by PC1 (pccut=0.8, 0.9)
  
  group=vector("list",numclst)
  PCvars=NULL
  if(numclst>0){
    for(i in 1:numclst){
      SNPinCLQ=which(CLQout==i)
      subset=genoin[,SNPinCLQ]
      PCout=PCscores(subset,pccut=pccut, type="PC1") 
      PCvars=cbind(PCvars, PCout)
      group[[i]]=varnames[SNPinCLQ]
    }
    grpnames=paste("grp",1:numclst,sep='')
    colnames(PCvars)=grpnames
    
    
    rmInd=which(!is.na(CLQout))
    rmvSNPs2=varnames[rmInd]
    remainSNPs2=varnames[-rmInd]   
    dataup=data.frame(onedata[,c(Yvar,remainSNPs2)]) #phenotype combined
    if(dim(dataup)[2]==1) colnames(dataup)=Yvar
    dataup=cbind(dataup,PCvars)
  }
  if(numclst==0){
    dataup=data.frame(onedata[,c(Yvar,remainSNPs)]) 
  }
  
  #4~5. Vif removal process
  #4~5 for partitions if the number of variables bigger than the limit!
  varnames=setdiff(colnames(dataup),Yvar)  
  KK=length(varnames)
  cycle=1
  k=0
  while(cycle>0){
    k=k+1
    vrout=vifremovalcycle(dataup, grpnames, Yvar, Klim, itermax)
    dataup=vrout$dataup
    cycle=vrout$change
    remainvars=vrout$remainvars
    grpnames=vrout$grpnames
  }
  
  
  
  #6. Find the tag for removed variables among remaining variables 
  
  if(numclst>0) allgeno=cbind(onedata[,allSNPs],PCvars)  
  if(numclst==0) allgeno=onedata[,allSNPs]
  
  #rmvSNPs - aliased/rmvSNPs-replaced by groups / rmvvars-group removed
  rmvSNPs3=setdiff(setdiff(varnames,remainvars), colnames(PCvars)) #single SNPs among removed SNPs by step 4,5
  removed=c(rmvSNPs,rmvSNPs2, rmvSNPs3)  
  removed=allSNPs[sort(match(removed,allSNPs))]
  
  tagvarlist=NULL
  tagcorsf=NULL
  
  if(length(removed)>0){
    
    cornew=cor(allgeno)
    for(i in 1:length(removed)){
      corvec=cornew[removed[i],]
      corvec=corvec[which(names(corvec) %in% remainvars)]
      maxcor=max(corvec)
      tagvars=names(corvec)[which(corvec==maxcor)]
      tagvarlist=c(tagvarlist,tagvars[1])
      tagcorsf=c(tagcorsf,maxcor)
    }
    names(tagvarlist)=removed
  }
  
  
  #optional 7.Obtain PCs of removed variables and examine R2 after regression on the final variables
  rmvvars=NULL
  rmvvarsf=NULL
  rmvgroups=setdiff(colnames(PCvars),remainvars)
  rmvgrind=match(rmvgroups,colnames(PCvars))
  rmvSNPs2_1=unlist(group[rmvgrind])  #SNPs in removed groups only (rmvSNPs2 is all grouped SNPs)
  removedSNPs=c(rmvSNPs,rmvSNPs2_1,rmvSNPs3)  #
  lowr2ind=NULL
  if(length(removedSNPs)>0){
    
    rmvgeno=onedata[,removedSNPs]
    rmvPCs=as.matrix(rmvgeno)
    if(length(removedSNPs)>1) rmvPCs=PCscores(rmvgeno,pccut=pccut, type="PCmulti") 
    if(dim(rmvPCs)[2]==1) colnames(rmvPCs)="PC1" 
    
    
    regdata=cbind(dataup,rmvPCs)
    rsquared=rep(0,dim(rmvPCs)[2])
    for(i in 1:dim(rmvPCs)[2]){
      model=paste(colnames(rmvPCs)[i], "~", paste(remainvars,collapse=" + "))
      rout=summary(lm(model, data=regdata))
      rsquared[i]=rout$r.squared
      
    }
 
    lowr2ind=which(rsquared<(1-1/vifcut))
    
    
    
    if(length(lowr2ind)>0){
      
      rmvPCs=rmvPCs[,lowr2ind,drop=FALSE] 
      regdata=cbind(dataup,rmvPCs) 
      
      #alias check if there are multiple RPCs expected
      if(length(lowr2ind)>1){
        model=paste(Yvar," ~", paste(setdiff(colnames(regdata),Yvar),collapse=" + "))
        outo=alias(lm(model,data=regdata))
        if(!is.null(outo$Complete)) {
          remvars=colnames(outo$Complete)[-1] #remained variables after alias removal
          rmvvars=rownames(outo$Complete)   #removed variables by alias removal
          rmvvarsf=rmvvars
          for(i in 1:length(rmvvars)){
            if(!( rmvvars[i] %in% colnames(rmvPCs) )){
                  rmvvarsf[i]=intersect(colnames(rmvPCs),remvars)[1]
            }
          }
          rmvvarsf=unique(rmvvarsf)
        }
      }
      remPCs=setdiff(colnames(rmvPCs),rmvvarsf)
      RPCind=match(remPCs,colnames(rmvPCs))
      #RPC1,RPC2... variables for final PCS
      RPCvars=paste("RPC",1:length(remPCs),sep='')
      rmvPCs_selected=rmvPCs[,remPCs,drop=FALSE]
      colnames(rmvPCs_selected)=RPCvars
      dataup=cbind(dataup, rmvPCs_selected)
      
    }
    
  }

  
  #final output
  #LPC variables
  LPCvars=NULL
  SNPsinLPC=NULL
  result=NULL

  
  if(length(grpnames)>0){ 
    LPCvars=paste("LPC",1:length(grpnames), sep='')
    grpnamesf=colnames(dataup)[match(grpnames,colnames(dataup))]
    colnames(dataup)[match(grpnames,colnames(dataup))]=LPCvars
    
    remaingrind=setdiff(1:numclst,rmvgrind)
    for(i in 1:length(LPCvars)){
      SNPsinLPC=c(SNPsinLPC, paste(group[[remaingrind[i]]],collapse="-"))
    }
    names(SNPsinLPC)=LPCvars
    
    #update taglist names (grp -> LPC)
    result <- vector("numeric",length(LPCvars))
    for(i in 1:length(LPCvars)){
      tagvarlist[which(tagvarlist==grpnamesf[i])]=LPCvars[i]
      #result=NULL
      
     geno1<-genodata[group[[remaingrind[i]]]]
      cov1<-cov(geno1)
      e1<-eigen(cov1)
      max.eigen<-max(e1$values)
      sumLPC<-sum(e1$values)
      VarEachLPC<-max.eigen/sumLPC
      result[i]<-VarEachLPC
      

    }
    result
    
    } 
    
  return(list(vdata=dataup,LPCinfo=SNPsinLPC,VarEachLPC=result, taglist=tagvarlist, aliasremoved=rmvSNPs, removed=removedSNPs, group=group, remainedgrps=grpnames, RPCind=RPCind))
}


###VIF removal iteration until no change###################################
vifremovalcycle=function(dataup, grpnames, Yvar, Klim, itermax){
  
  varnames=setdiff(colnames(dataup),Yvar)  
  KK=length(varnames)
  partition=0
  change=0
  remainvars=NULL
  rmvvars=NULL
  rmvgrpnames=NULL
  if(KK< Klim){
    vrout=vif.remove(dataup, grpnames, Yvar, Klim, itermax)

    remainvars=vrout$remainvars
    rmvvars=vrout$rmvvars
    rmvgrpnames=vrout$rmvgrpnames
  }
  if(KK>=Klim){
    numset2=floor((KK-1)/Klim)+1
    eachsize=floor((KK+2)/numset2)
    
    for(k in 1:numset2){
      start=eachsize*(k-1)+1
      end=KK
      if(k<numset2) end=eachsize*k
      dataupin=data.frame(dataup[,c(Yvar,varnames[start:end])])
      vrout=vif.remove(dataupin, grpnames, Yvar, Klim, itermax)
      remainvars=c(remainvars,vrout$remainvars)
      rmvvars=c(rmvvars,vrout$rmvvars)
      rmvgrpnames=c(rmvgrpnames, vrout$rmvgrpnames)
    }
    
    partition=1
    if(partition==1 && length(rmvvars)>0) change=1
  }
  grpnames=setdiff(grpnames, rmvgrpnames)
  dataupout=data.frame(dataup[,c(Yvar,remainvars)]) 
  return(list(dataup=dataupout, remainvars=remainvars, change=change, grpnames=grpnames)) 
}



###VIF removal using partitioning 
vif.remove=function(dataup, grpnames, Yvar, Klim, itermax, highnumlim=60){
  varnames=setdiff(colnames(dataup),Yvar)  
  rmvvars=NULL
  remainvars=setdiff(colnames(dataup),Yvar)
  orggrpnames=grpnames
  
  if(length(varnames)> Klim ){
    print("step 3: too many variables for vif calculation. No vif removal")
    return(list(remainvars=remainvars, rmvvars=rmvvars, dataup=dataup))
    
  }
  
  #preprocessing if any aliased SNPs exists
  alremout=alias.remove(dataup, varnames, Yvar) 
  remainvars=alremout$remainSNPs
  rmvvars=alremout$rmvSNPs
  grpnames=setdiff(grpnames,rmvvars)
  if(length(rmvvars)>0) dataup=data.frame(dataup[,c(Yvar,remainvars)])  
  #######
  
  vifs=1
  if(length(remainvars)>1 && length(remainvars)<250){
    model=paste(Yvar, "~", paste(remainvars,collapse=" + "))
    vifs=vif(lm(model, data=dataup))
  }
  if(length(remainvars)>=250){
    vifs=VIFlarge(dataup, Yvar)
  }

  
  #5. Remove high VIF variables iteratively. 
  
  maxvif=max(vifs)
  maxind=which(vifs==max(vifs))[1]
  var1=remainvars[maxind]
  
  k=0
  while(maxvif>vifcut){
    k=k+1
   
    highvifs=which(vifs>vifcut)
    highnum=length(highvifs)
    if(highnum>highnumlim){ 
      numrmv=floor(highnum*0.1)
      cut=vifs[which(rank(-vifs)==numrmv)]
      tobermv=which(vifs>=cut)
    }
    if(highnum<=highnumlim){
      tobermv=which(names(dataup) == var1)[1]
    }
    
    rmvvars=c(rmvvars, names(dataup)[tobermv])
    dataup=dataup[,-tobermv]  
    remainvars=setdiff(colnames(dataup),Yvar)  

    grpnames=setdiff(grpnames,rmvvars)
    
    if(length(remainvars)==1) break
    
    #calculate VIF
    if(length(remainvars)<250){
      model=paste(Yvar, "~", paste(remainvars,collapse=" + "))
      vifs=vif(lm(model, data=dataup))
    }
    if(length(remainvars)>=250){
      vifs=VIFlarge(dataup, Yvar)
    }
    
    maxvif=max(vifs)
    maxind=which(vifs==max(vifs))[1]
    var1=remainvars[maxind]

    if(k>itermax) {
      print("Iteration reached the max")
      break
    }
  }
  rmvgrpnames=setdiff(orggrpnames, grpnames)
  
  return(list(remainvars=remainvars, rmvvars=rmvvars, dataup=dataup, rmvgrpnames=rmvgrpnames))
  
}



###partition update + additional removal comparing between partitions (alias removal step 5~6)
removalcycle=function(onedata,remSNPin, rmvSNPin, remsizesin,Klim, Yvar){
  change=0
  partition=1
  #alias step 4: repeat alias step 2-3 until no combining the partition performed
  i=0
  while(partition>0){
    i=i+1
    pupout=part.update(onedata,remSNPin, rmvSNPin, remsizesin,Klim, Yvar) #step 2 & 3
    remSNPin=pupout$remSNPlist
    rmvSNPin=pupout$rmvSNPlist
    remsizesin=pupout$remsizes
    partition=pupout$change
  }
  
  
  #step 5 : dependency between different matrix
  for(i in 1:length(remsizesin)){
    for(j in 1:length(remsizesin)){
      if(i!=j){
        ones=rep(1,dim(onedata)[1]) #For intercept 
        dpnlist=mgcv::fixDependence(cbind(ones,onedata[,remSNPin[[j]]]),onedata[,remSNPin[[i]]], tol=1e-8, strict=TRUE)
        
        #update rmv rem SNPs data
        if(length(dpnlist)>0){
          rmvSNPin[[i]]=c(rmvSNPin[[i]],remSNPin[[i]][dpnlist])
          remSNPin[[i]]=remSNPin[[i]][-dpnlist]
          remsizesin[i]=remsizesin[i]-length(dpnlist)
          change=1
        }
      }
    }
  }
  
  return(list(remSNPlist=remSNPin, rmvSNPlist=rmvSNPin, remsizes=remsizesin, change=change)) 
  
}

#partition update (alias removal step 2~3)
part.update=function(onedata,alremSNPlist, alrmvSNPlist, alremsizes,Klim, Yvar ){
  change=0
  newremsizes=NULL
  #alias step 2
  combine=list()
  sumsubset=findsumsubset(alremsizes,Klim)
  while(length(sumsubset)>1){
    combine=c(combine,list(sumsubset))
    sumsubset=findsumsubset(alremsizes,Klim, sort(unlist(combine)))
  }
  
  #alias step 3
  if(length(combine)>0){
    change=1
    alremSNPlist2=vector("list",length(combine))
    alrmvSNPlist2=vector("list",length(combine))
    for(k in 1:length(combine)){
      remSNPs=unlist(alremSNPlist[combine[[k]]])
      prermvSNPs=unlist(alrmvSNPlist[combine[[k]]])
      
      alremout=alias.remove(onedata,remSNPs,Yvar)
      alremSNPlist2[[k]]=alremout$remainSNPs
      alrmvSNPlist2[[k]]=c(prermvSNPs, alremout$rmvSNPs)
      newremsizes=c(newremsizes, length(alremout$remainSNPs))
    }
    #add partitions not combined
    alremSNPlist2=c(alremSNPlist2, alremSNPlist[-unlist(combine)])
    alrmvSNPlist2=c(alrmvSNPlist2, alrmvSNPlist[-unlist(combine)])
    newremsizes=c(newremsizes,alremsizes[-unlist(combine)])
    
  }
  if(length(combine)==0){
    alremSNPlist2=alremSNPlist
    alrmvSNPlist2=alrmvSNPlist
    newremsizes=alremsizes
  }
  
  return(list(remSNPlist=alremSNPlist2, rmvSNPlist=alrmvSNPlist2, remsizes=newremsizes, change=change)) 
}


#find list of group of patition
findsumsubset=function(alremsizes,Klim, exclude=NULL){
  setranks=rank(alremsizes,ties.method="first")
  setorder=order(alremsizes)
  alremsizesex=alremsizes
  alremsizesex[exclude]=0
  cumsizes=cumsum(alremsizesex[setorder])
  numcum=max(which(cumsizes<Klim))
  subset=sort(setdiff(setorder[1:numcum],exclude))
  return(subset)
}

#alias removal 
alias.remove=function(onedata, SNPnames, Yvar){
  model=paste(Yvar," ~", paste(SNPnames,collapse=" + "))
  out=alias(lm(model,data=onedata))
  remSNPs=SNPnames
  rmSNPs=NULL
  if(!is.null(out$Complete)) {
    remSNPs=colnames(out$Complete)[-1] #remained SNPs after alias removal
    #genoup=onedata[,remainSNPs]
    rmSNPs=rownames(out$Complete)   #removed SNPs by alias removal
  }
  return(list(rmvSNPs=rmSNPs, remainSNPs=remSNPs))
}


#VIF calculation for large data
#input: data of Y and X, no missing data
VIFlarge=function(onedata, Yvar=NULL){
  X=onedata[,!(colnames(onedata) %in% Yvar)]
  s=dim(X)[2]
  vfs=rep(0,s)
  for(i in 1:s){
    model=paste( paste(colnames(X)[i], "~") , paste(colnames(X)[-i], collapse=" + ")) 
    out=summary(lm(model, data=X))
    vfs[i]=1/(1-out$r.squared)
  }
  
  return(vfs)
}

#PC80 scores (PCA in FactoMineR package) 
#<input>
#indata: input data of X variables only
#type : PCmulti=choose PCs according to pccut, PC1=choose first PC
#pccut: threhold for selecting PC variables (minimum PCs explaining variances at least >pccut)
#<output>
#New data with PC variables

PCscores=function(indata,pccut=0.8,type=c("PCmulti","PC1")){
  indmax=dim(indata)[2]
  fRmodel=colnames(indata)
  mean_onedata=apply(indata,2,mean)
  max_onedata=apply(indata,2,max)
  sub=c(1:length(mean_onedata))[max_onedata>mean_onedata] #remove no variation data
  pc=FactoMineR::PCA(indata[,sub], ncp=indmax, graph=F)
  
  corone=cor(indata[,sub])
  ei=eigen(corone)
  
  cumei=NULL
  for(m in 1:length(ei$values)){
    cumei[m]=sum(ei$values[1:m])
  }
  cumei=cumei/cumei[m]
  
  ind=1  
  if(type=="PCmulti"){
    indlength=min(sum(cumei<pccut)+1, indmax)
    ind=1:indlength
  }
  
  
  #print(ind)
  scores=as.matrix(pc$ind$coord[,ind])
  colnames(scores)=paste("PC",1:dim(scores)[2],sep="")
  return(scores)
}





#PC80 scores 
#<input>
#indata: input data of X variables only
#type : PCmulti=choose PCs according to pccut, PC1=choose first PC
#pccut: threhold for selecting PC variables (minimum PCs explaining variances at least >pccut)
#<output>
#New data with PC variables

PCscores_prcomp=function(indata,pccut=0.8,type=c("PCmulti","PC1")){
  fRmodel=colnames(indata)
  mean_onedata=apply(indata,2,mean)
  max_onedata=apply(indata,2,max)
  sub=c(1:length(mean_onedata))[max_onedata>mean_onedata] #remove no variation data
  pc=prcomp(indata[,sub], scale.=T, tol=1e-6, rank.=4)
  
  corone=cor(indata[,sub])
  ei=eigen(corone)
  
  cumei=NULL
  for(m in 1:length(ei$values)){
    cumei[m]=sum(ei$values[1:m])
  }
  cumei=cumei/cumei[m]
  
  ind=1
  if(type=="PCmulti"){
    if(cumei[1]<pccut) { 
      ind=c(1:length(sub))[cumei<pccut]  
      ind=c(ind,length(ind)+1)   
    }
    
  }
  
  
  #print(ind)
  return(as.matrix(pc$x[,ind]))
}



CLQD <- function(geno, SNPinfo, CLQcut=0.5, clstgap=40000, hrstType=c("near-nonhrst", "fast", "nonhrst"), hrstParam=200,
                 CLQmode=c("density", "maximal"), LD=c("r2", "Dprime"))
{
  ########################################################################################################
  skipRatio <- 0
  CLQmode <- match.arg(CLQmode)
  LD <- match.arg(LD)
  hrstType <- match.arg(hrstType)
  geno <- as.matrix(geno)
  # filtering all NA SNPs
  allNASNPs = apply(geno, 2, function(x) all(is.na(x)))
  geno<-geno[, !allNASNPs]
  SNPinfo<-SNPinfo[!allNASNPs, ]
  # Main Function
  SNPbps <- as.integer(SNPinfo[, 3])
  if(LD == "r2"){
    OCM <- suppressWarnings(cor(geno, use="pairwise.complete.obs"))
    diag(OCM) <- 0
    OCM[is.na(OCM)]<-0
    OCM[abs(OCM) < CLQcut] <- 0
    OCM[abs(OCM) >= CLQcut] <- 1
  }else if(LD == "Dprime"){
    OCM <- genoDp(geno)
    OCM[(OCM==Inf)|(OCM==-Inf)] <- 0
  }
  binvector <- rep(NA, dim(OCM)[2])
  binnum <- 1
  # re.SNPbps <- SNPbps
  if(all(OCM==0)) return(binvector)
  
  # test graph complexity
  OCM1 <- OCM
  # bothhigh <- (degs>=400 & cores>=400)
  if(hrstType == "nonhrst"){
    heuristic <- FALSE
  }else if(hrstType == "near-nonhrst"){
    heuristic <- TRUE
    heuristicNum <- 0
    while(heuristic == TRUE){
      g <- igraph::graph_from_adjacency_matrix(OCM1, mode="undirected", weighted=TRUE, diag=FALSE, add.colnames=NA)
      cores <- igraph::coreness(g)
      highcore <- sum(cores>=hrstParam)
      local_cores <- table(cores[cores>=(hrstParam)])
      local_cores <- local_cores[local_cores>=(hrstParam)]
      if(length(local_cores>0)){
        if(heuristicNum==0){
          message("Use near-CLQ heuristic procedure!!")
          heuristicNum <- heuristicNum+1
          message(paste("hueristic loop", heuristicNum))
        }else{
          heuristicNum <- heuristicNum+1
          message(paste("hueristic loop", heuristicNum))
        }
        # find dense region
        # local_cores <- table(cores[cores>=(hrstParam)])
        # local_cores <- local_cores[local_cores>=(hrstParam)]
        local_hrstParam <- names(local_cores[which(as.integer(names(local_cores))==max(as.integer(names(local_cores))))])
        local_hrstParam <- as.numeric(local_hrstParam)
        bothhighSNPs <- which(cores == local_hrstParam)
        SNPset1 <- which(is.na(binvector))[bothhighSNPs]
        nowOCM <- OCM[SNPset1, SNPset1]
        heuristicBins <- heuristicCLQ(nowOCM, hrstParam)
        binvector[SNPset1[heuristicBins]]<-binnum
        binnum <- binnum+1
        OCM1 <- OCM[is.na(binvector), is.na(binvector)]
      }else{
        # if(heuristicNum ==0){
        #   # message("We do not need heuristic procedure")
        # }else{
        #   message("End new heuristic procedure")
        # }
        heuristic <- FALSE
      }
    } #end while
  }else if(hrstType == "fast"){
    checkLargest <- TRUE
    r2Mat <- OCM
    re.SNPbps <- as.integer(SNPinfo[,3])
    # firstTerm = T
    while(checkLargest == TRUE){
      g <- igraph::graph_from_adjacency_matrix(r2Mat, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA)
      compo = igraph::components(g)
      componum = which(compo$csize==max(compo$csize))[1]
      compov = which(compo$membership==componum)
      compadjM = OCM1[compov, compov]
      cg = igraph::graph_from_adjacency_matrix(compadjM, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA)
      if((median(igraph::coreness(cg))>80 & max(igraph::coreness(cg))>100)| (quantile(igraph::coreness(cg), 0.75)>100 & max(igraph::coreness(cg))>100)){
        message("use fast heuristic procedure!")
        # if(quantile(degrees, 0.7) == 1) break
        # message(head((quantile(degrees, 0.7))))
        degrees = apply(r2Mat, 1, sum)
        maxdegv = which(degrees >=max(quantile(degrees, 0.7), 80))
        # if(length(maxdegv)>=1){
        maxdegvs = maxdegv
        edgeDens = NULL
        for(maxdegv in maxdegvs){
          Bignbds = which(r2Mat[maxdegv,, drop = FALSE]>0, arr.ind = TRUE)
          Bignbds.c = unique(Bignbds[,2])
          newr2Mat = r2Mat[Bignbds.c,Bignbds.c]
          EdgeDen = sum(newr2Mat)/((dim(newr2Mat)[1])*(dim(newr2Mat)[1]-1))
          edgeDens = c(edgeDens, EdgeDen)
        }
        maxdegvs = maxdegvs[order(edgeDens, decreasing = TRUE)]
        edgeDens = edgeDens[order(edgeDens, decreasing = TRUE)]
        degv = maxdegvs[1]
        edgeD = edgeDens[1]
        Bignbds = which(r2Mat[degv,, drop = FALSE]>0, arr.ind = TRUE)
        Bignbds.c = unique(Bignbds[,2])
        # maxiC = maximal.cliques(g, min = dim(OCM1)[1]*0.9)
        # largestOneRange = range(Bignbds.c)
        # largestSNPn = diff(largestOneRange)
        # largestCsize = length(Bignbds.c)
        nowSNPsbp = re.SNPbps[Bignbds.c]
        nowSNPsbploca = match(nowSNPsbp, SNPbps)
        binvector[nowSNPsbploca] <- binnum
        binnum = binnum + 1
        r2Mat <- r2Mat[-Bignbds.c, -Bignbds.c, drop = FALSE]
        OCM1 <- OCM1[-Bignbds.c, -Bignbds.c, drop = FALSE]
        re.SNPbps <- re.SNPbps[-Bignbds.c]
        # message("case2")
        checkLargest = TRUE
        if(length(re.SNPbps)<500)  checkLargest = FALSE
      }else{
        checkLargest = FALSE
      }
      
    }
  }
  # binvector splitting by clstgap
  if(all(is.na(binvector))==FALSE){
    binveclist <- lapply(seq_len(max(binvector, na.rm = TRUE)),
                         function(x) SNPinfo[,3][which(binvector == x)])
    binvecbplist <- newSplitCliques(binveclist, clstgap)
    binvecbpindex <- lapply(binvecbplist, function(x) match(x, SNPinfo[,3]))
    binvecbpindex <- binvecbpindex[vapply(binvecbpindex, length, c(1))>1]
    binvector <- rep(NA, length(binvector))
    for(i in seq_len(length(binvecbpindex))){
      binvector[binvecbpindex[[i]]]<-i
    }
  }
  # take Toooo Big block First!
  r2Mat <- OCM[is.na(binvector), is.na(binvector)]
  if(sum(is.na(binvector))<=1 | (sum(r2Mat)==0)) return(binvector)
  g <- igraph::graph_from_adjacency_matrix(r2Mat, mode="undirected", weighted=TRUE, diag=FALSE, add.colnames=NA)
  max.cliques <- igraph::max_cliques(g, min=2)
  if(length(max.cliques)==0) stop("max.cliques is empty")
  re.SNPbps <- SNPbps[is.na(binvector)]
  bp.cliques <- lapply(max.cliques, function(x) re.SNPbps[x])
  split.bp.cliques <- newSplitCliques(bp.cliques, clstgap)
  # reduce the number of maximal clique? (candidate) or
  # modify density function. narrow SNP distance, so small cliques are chosen preferencely.
  repeat {
    # message(binnum)
    if (all(is.na(binvector) == FALSE)) {
      break
    }
    if(length(split.bp.cliques)==0) break
    now.split.bp.cliques <- split.bp.cliques
    if(CLQmode =="density"){
      density.v <- vapply(now.split.bp.cliques, function(x) ((length(x)))/(diff(range(x))/1000), 1)
    }else{
      density.v <- vapply(now.split.bp.cliques, length, 1)
    }
    
    max.d <- which(density.v == max(density.v))
    max.cluster <- now.split.bp.cliques[max.d]
    if (length(max.cluster) > 1) {
      # if there are two bins of same density, then we choose the bigger one.
      max.cluster <- max.cluster[order(vapply(max.cluster, length, 1), decreasing=TRUE)]
    }
    max.cluster <- max.cluster[[1]]
    max.cluster.od <- match(max.cluster, re.SNPbps)

    ## excluding all SNPs in max.cluster from re.SNPbps
    split.bp.cliques <- lapply(split.bp.cliques, function(x) setdiff(x, max.cluster))
    split.bp.cliques <- unique(split.bp.cliques)
    split.bp.cliques <- split.bp.cliques[which(vapply(split.bp.cliques, length, 1) > 1)]
    binvector[match(max.cluster, SNPbps)] <- binnum
    binnum=binnum + 1
    
    if (length(re.SNPbps) < 2) {
      break
    }
    if(length(split.bp.cliques)==0) break
  }  ##end repeat
  return(binvector)
}

newSplitCliques <- function(cliques.bp, gapdist)
{
  nowlist <- lapply(cliques.bp, sort)
  fixlist <- NULL
  repeat {
    need.split <- which(lapply(nowlist, function(x) max(diff(x)) > gapdist) == TRUE)
    need.fix <- which(lapply(nowlist, function(x) max(diff(x)) > gapdist) == FALSE)
    addlist <- nowlist[need.fix]
    fixlist <- c(fixlist, addlist)
    if (length(need.split) == 0) {
      break
    }
    nowlist <- nowlist[need.split]
    nowlength <- length(nowlist)
    newlist <- as.list(rep(NA, nowlength))
    for (i in seq_len(nowlength)) {
      gap <- diff(nowlist[[i]])
      frontpart <- nowlist[[i]][seq_len(min(which(gap > gapdist)))]
      restpart <- nowlist[[i]][-(seq_len(min(which(gap > gapdist))))]
      nowlist[[i]] <- frontpart
      newlist[[i]] <- restpart
    }
    addlist <- nowlist[vapply(nowlist, function(x) length(x) > 1, TRUE)]
    fixlist <- c(fixlist, addlist)
    nowlist <- newlist[vapply(newlist, function(x) length(x) > 1, TRUE)]
  }
  return(fixlist)
}
#2
heuristicCLQ <- function(subOCM, hrstParam)
{
  degs <- apply(subOCM, 1, sum)
  top5pctDegVs <- which(degs>=quantile(degs, 0.95))
  top5pctDegVs <- top5pctDegVs[order(degs[top5pctDegVs],decreasing=TRUE)]
  heuristicBinbasket <- NULL
  total_density <- NULL
  v_dens_all <- c()
  find_opt <- FALSE
  for(v in top5pctDegVs){
    v_nbds <- as.integer(which(subOCM[v,]!=0))
    candibin <- c(v, v_nbds)
    v_density <- sum(subOCM[candibin,candibin])/(length(candibin)^2 - length(candibin))
    v_dens <- c(v_density, candibin)
    if(v_density>=0.95) {
      return(candibin)
      find_opt <- TRUE
      break
    }else{
      new_v_density <- 1
      while((v_density+0.0001)<new_v_density){
        newdeg <- apply(subOCM[v_nbds,], 1, sum)
        v_nbds <- v_nbds[-(which(newdeg==min(newdeg)))]
        candibin <- c(v, v_nbds)
        new_v_density <- sum(subOCM[candibin,candibin])/(length(candibin)^2 - length(candibin))
        if(new_v_density>=0.95){
          return(candibin)
        }else{
          v_dens_all <- c(v_dens_all, new_v_density)
          heuristicBinbasket<- c(heuristicBinbasket, list(candibin))
        }
        if(length(v_nbds) < (2*hrstParam/3)) break
      }#end while
    }#end else
  }#end for
  ## return results
  max_den_bin = heuristicBinbasket[which(v_dens_all == max(v_dens_all))]
  if(length(max_den_bin) == 1) {
    return(max_den_bin[[1]])
  }else{
    max_den_bin_leng = vapply(max_den_bin, length, 1)
    candibin = max_den_bin[[which(max_den_bin_leng == max(max_den_bin_leng))[1]]]
    return(candibin)
  }
}

