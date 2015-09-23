BFDR = function(p,a){
  if(missing(p)){   
    stop("posterior probabilities missing")     
  }
  
  if(missing(a)){   
    stop("alpha missing")     
  }
  
  thresholds=seq(from = 0.01, by = 0.01, to = 0.99)
  fdr = rep(NA,length(thresholds))
  pth = 1
  for (i in 1:length(thresholds)){
    Ip=(p>=thresholds[i])
    fdr[i]=sum((1-p)*Ip,na.rm = T)/(sum(Ip,na.rm = T)+1)
  }
  
  #fdr[is.na(fdr)]=0
  
  pth =min(thresholds[which(fdr<=a)])
Ip=(p>a)
fdr=sum((1-p)*Ip)/sum(Ip)
  return(pth)
}