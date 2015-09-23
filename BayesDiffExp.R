BayesDiffExp=function(Pijk,PCT,muPct,muTau,pr_gam_ij,a,b,vartau)
{
  #################################################################
  ####################      Checks    #############################
  #################################################################
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
  
  if(missing(pr_gam_ij)){
    stop('priors p0 missing')
  }
  
  if(missing(a)){
    warning('hyperparameter missing- using default value')
    a = 1
  }
  
  if(missing(b)){
    warning('hyperparameter missing- using default value')
    b = 1
  }
  
  
  nc=sum(PCT!=0)
  nt=sum(Pijk!=0)
  n=nt+nc
  
  if(missing(vartau)){
    warning('prior hyperparameter missing- using default value')
    vartau = n
  }
  
  
  
  Pijk_nz=Pijk[Pijk!=0]
  PCT_nz=PCT[PCT!=0]
 ############################################################
 ################## Model 1 #################################
 ############################################################
 
 Py_g1 = -(nc+nt)*log(2*pi) + a*log(b) - lgamma(a) - 0.5 * log(nc+nt+1)
 astar = a + 0.5*(nc+nt)
 bstar = b + 0.5*( muPct^2 + sum(PCT_nz^2) + sum(Pijk_nz^2) - (sum(PCT_nz) + sum(Pijk_nz) + muPct)^2/(1+nc+nt+1))
 Py_g1 = Py_g1 + lgamma(astar) - astar*log(bstar)
  Py_g1 =exp(Py_g1)
############################################################
################## Model 2 #################################
############################################################

sigmat = 1/vartau + nt - nt^2/(nc+nt+1)
mt  = muTau/vartau + sum(Pijk_nz) -nt/(nc+nt+1)*(sum(Pijk_nz)+muPct+sum(PCT_nz))
bstar = b + 0.5*(muPct^2 + muTau^2/vartau + sum(PCT_nz^2) + sum(Pijk_nz^2) - (sum(PCT_nz) + muPct+ sum(Pijk_nz))^2/(nc+nt+1) - mt^2/sigmat )
                   
 Py_g2 =  -(nc+nt)*log(2*pi) + a*log(b) - lgamma(a) - 0.5 * log(nc+nt+1) -0.5*log(sigmat)- 0.5*log(vartau)
Py_g2 = Py_g2 + lgamma(astar) - astar*log(bstar)
  Py_g2=exp(Py_g2)
pr_gam_ij_0=1-pr_gam_ij
p=(Py_g2*pr_gam_ij)/(Py_g2*pr_gam_ij+Py_g1*pr_gam_ij_0)
  return(p)
}