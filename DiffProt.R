DiffProt = function(Pijk,PCT,muPct,muTau,a,b,vartau,alpha){
  source('BayesDiffExp.R')
  source('BFDR.R')
  if(missing(Pijk)){   
    stop("Pijk missing")     
  }
  
  if(missing(PCT)){   
    stop("PCT missing")     
  }
  
  
  if(missing(muPct)){
    warning('mean of muPct. Using default parameter')
    muPct = 0
  }
  
  if(missing(muTau)){
    warning('mean of muTau. Using default parameter')
    muTau = 0
  }
  
  
  if(missing(a)){
    warning('hyperparameter missing- using default value')
    a = rep(1,dim(PCT)[1])
  }
  
  if(missing(b)){
    warning('hyperparameter missing- using default value')
    b = rep(1,dim(PCT)[1])
  }
  
  if(missing(alpha)){   
    warning("alpha missing-using default value")  
    alpha = 0.05
  }
  
  
  #####    Filtering ################
  
  locz.PCT = apply(PCT,1,function(x) sum(x==0))
  locz.Pijk = apply(Pijk,1,function(x) sum(x==0))
  locz = intersect(which(locz.PCT==dim(PCT)[2]),which(locz.Pijk==dim(Pijk)[2]))
  PCT = PCT[-locz,]
  Pijk = Pijk[-locz,]
  
  xxx = function(x){
    
    x[x!=0]=(x[x!=0] - mean(x[x!=0]))/sd(x[x!=0])
    return(x)
  }
  
  ######## Standardising ###########################
  
  xscale = t(sapply(1:dim(cbind(PCT,Pijk))[1],function(x,i) xxx(x[i,]), x=cbind(PCT,Pijk) ))
  xscale[is.nan(xscale)] =0
  PCT.trans = xscale[,1:dim(PCT)[2]]#PCT/apply(PCT,1,function(x) mean(x!=0))
  Pijk.trans = xscale[,(dim(PCT)[2]+1):dim(xscale)[2]]#Pijk/apply(Pijk,1,function(x) mean(x!=0))
  
  ##### Priors #############
  ZPCT = apply(PCT,1,function(x) sum(x==0))/dim(PCT)[2]
  ZPijK = apply(Pijk,1,function(x) sum(x==0))/dim(Pijk)[2]
  p0.H1 = 0.5 + abs(-ZPijK + ZPCT)/2
  #################### Bayes Factor ###########################################
  tryCatch({
    post.H1 = sapply(1:dim(PCT)[1], function(i,PCT,Pijk,muPct,muTau,pr_gam_ij,a,b,Vbstar) BayesDiffExp(Pijk[i,],PCT[i,],muPct,muTau,pr_gam_ij[i],a[i],b[i],Vbstar) , PCT = PCT.trans, Pijk = Pijk.trans, muPct= muPct, muTau = muTau, pr_gam_ij = p0.H1,a=a,b=b ,Vbstar=vartau)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  FDR = BFDR(post.H1,alpha)
  
  loc = which(post.H1>FDR)
  BF.proteins = row.names(PCT)[loc]
  Results = matrix(NA,sum(post.H1>=FDR,na.rm=T),4)
  colnames(Results) = c('Proteins','Probability','Fold-change','Score')
 
  Results = data.frame(Results)
  
  Results[,1] = BF.proteins
  Results[,2] = post.H1[which(post.H1>=FDR)]
  
  ########################################
  ###         Calculating Fold Change ####
  ########################################
  Results$Fold.change = as.numeric(Results$Fold.change)
  for(i in 1:length(BF.proteins)){
    loc = which(row.names(PCT)%in%BF.proteins[i])
    if(sum(PCT[loc,])==0 || sum(Pijk[loc,])==0) {
      Results[i,3] = Inf
    }else{
      Results[i,3] =  log2(sum(Pijk[loc,])/sum(PCT[loc,]))
    }
  }
  Results$Score = Results$Probability*abs(Results$Fold.change)
  Results = Results[order(Results$Score,decreasing = T),]
  row.names(Results) = 1:dim(Results)[1]
  
  list(Results= Results)
  
  
}