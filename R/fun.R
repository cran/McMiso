
#' @keywords internal
getDoses = function(X,y) {
  DOSES = as.matrix(dplyr::distinct(data.frame(X)))
  K = nrow(DOSES)
  rownames(DOSES) = paste("condition ",1:K,sep="")
  D = ncol(X)
  M = rep(0,K)
  S = rep(0,K)        
  for (k in 1:K) {
    vec = DOSES[k,]
    A = t(t(X)==vec)
    pos = which(rowSums(A)==D)
    M[k] = length(pos)
    S[k] = sum(y[pos])
  }
  list(alldoses=DOSES, M=M, S=S) 
}

#' @keywords internal
getOrd = function(alldoses) {
  K = nrow(alldoses)
  D = ncol(alldoses)
  PartOrd = matrix(rep(0,K*K), nrow=K)
  for (i in 1:(K-1)) {
    curdose = alldoses[i,]
    A = t(t(alldoses)>=curdose)
    posNeg = intersect((i+1):K, which(rowSums(A)==D))
    B = t(t(alldoses)<=curdose)
    posPos = intersect((i+1):K,which(rowSums(B)==D))
    PartOrd[i,posNeg] = -1
    PartOrd[posNeg,i] = 1
    PartOrd[i,posPos] = 1
    PartOrd[posPos,i] = -1
  }
  rownames(PartOrd) = colnames(PartOrd) = apply(round(alldoses,digits=2),1,paste,collapse=",")
  PartOrd
}


#' @keywords internal
getScenesV3 = function(alldoses) {
  K = nrow(alldoses)
  PartOrd = getOrd(alldoses)
  A = t(as.matrix(expand.grid(rep(list(0:1), K))))
  SCENES = NULL
  for (i in 1:ncol(A)) {
    newscene = A[,i]    
    for (k in 1:K) {
      violate = F
      curgamm = newscene[k]
      if (curgamm == 0) { 
        potpos = which(PartOrd[k,]==1)
        if (length(potpos)>0) {
          if (any(newscene[potpos]==1)) { violate = T; break; }
        }
      }
      else {
        potpos = which(PartOrd[k,]== -1)
        if (length(potpos)>0) {
          if (any(newscene[potpos]==0)) { violate = T; break; }
        }
      }
    }
    if (!violate) SCENES = cbind(SCENES, newscene)
  }
  list(PartOrd=PartOrd,SCENES=SCENES)
}

#' Probabilistic Bayesian classifier
#'
#' @param X numeric matrix of doses
#' @param y numeric response vector
#' @param method character, either "DU" or "UD"
#' @param a0 numeric, prior alpha
#' @param b0 numeric, prior beta
#' @param t0 numeric, threshold
#' @return A list with class "pbc"
#' @examples
#' A <- as.matrix(expand.grid(rep(list(0:1), 6)))
#' set.seed(2025)
#' X <- A[sample(nrow(A),size=500, replace = TRUE),]
#' y <- as.numeric(rowSums(X)>=3)
#' PBclassifier(X,y)
#' @references
#' Cheung YK, Diaz KM. Monotone response surface of multi-factor condition: estimation and Bayes classifiers. 
#' *J R Stat Soc Series B Stat Methodol.* 2023 Apr;85(2):497-522. 
#' doi: 10.1093/jrsssb/qkad014. Epub 2023 Mar 22. PMID: 38464683; PMCID: PMC10919322.
#' @export
PBclassifier = function(X,y, method="DU", a0 = 0.25, b0 = 0.25, t0=0.5) {
  if (ncol(X)==1) stop("variable dimension should be greater than 1")
  dosefit = getDoses(X,y)
  M = dosefit$M
  S = dosefit$S
  alldoses = dosefit$alldoses
  pt = 1 - stats::pbeta(t0, a0+S, b0+M-S)
  pt = pmin( pmax(pt, 1e-5), 1-1e-5)
  if (method=="DU") output = QuickDownTog(alldoses,pt,M)
  else if (method=="UD") output = QuickUpTog(alldoses,pt,M)
  else { stop("method not recognized") }
  val = list(alldoses = alldoses, M = M, S = S, yhat = output$combhat, pt=pt, logH=output$logH)
  class(val) = "pbc"
  return(val)
}

#' @keywords internal
QuickDownTog = function(alldoses0, pt0, M0, eps=0.5, K0=10, depth=1) {
  ### Initialization
  K = nrow(alldoses0)
  dep = depth
  totX = apply(alldoses0,1,sum)
  o = order(totX,pt0,decreasing=T)
  alldoses = as.matrix(alldoses0[o,])
  M = M0[o]
  pt = pt0[o]
  bayeshat = rep(0,K)
  bayeshat[pt > (1-eps)] = 1  
  PartOrdAll = getOrd(alldoses)

  ### Algorithm starts  
  sset = 1:min(K0,K)
  dosesSset = as.matrix(alldoses[sset,])
  SCENESsset = getScenesV3(dosesSset)$SCENES
  logprob1 = rep(NA, ncol(SCENESsset))
  for (iter in 1:ncol(SCENESsset)) {
    curscene = SCENESsset[,iter]
    logprob1[iter] = sum(M[sset]*curscene*log(pt[sset]*eps) + (M[sset]*(1-curscene)*log((1-pt[sset])*(1-eps)))) 
  }
  posBest1 = which(logprob1==max(logprob1))[1]
  hat0 = SCENESsset[,posBest1]

  isize = esize = rep(0,K)
  if (K<=K0) { hatcomb = hat0; }
  else {
#    if (dep==1) cat("\nDU Size",K,"iteration: ")
    for (k in (K0+1):K) {
#      if (dep==1) cat(round(k/K,digits=3),"-->")
      if (!bayeshat[k]) { hat1 = c(hat0, 0); } 
      else {
        hat1a = c(hat0,0)
        logH1a = sum(M[1:k]*hat1a*log(pt[1:k]*eps) + (M[1:k]*(1-hat1a)*log((1-pt[1:k])*(1-eps))))
        
        ### Calculate the second candidate
        U1set = which(PartOrdAll[k,1:(k-1)]==-1 & hat0==1)
        U0set = which(PartOrdAll[k,1:(k-1)]==-1 & hat0==0)	
        E1set = which(PartOrdAll[k,1:(k-1)]==0 & hat0==1)
        E0set = which(PartOrdAll[k,1:(k-1)]==0 & hat0==0)	
        hat1b = rep(NA,k)
        hat1b[k] = 1
        hat1b[U1set] = hat1b[E1set] = hat1b[U0set] = 1
        
        if (length(E0set)>0) {
          esize[k] = length(E0set)
          if (length(U0set)==0) { hat1b[E0set] = 0; }
          else if (length(E0set)==1) { hat1b[E0set] = bayeshat[E0set]; isize[k] = 1; } 
          else {  
            I0 = U0set
            if (length(I0)==1) {
              LI0 = intersect(which(PartOrdAll[I0,1:k]==1), E0set)
              UI0 = intersect(which(PartOrdAll[I0,1:k]==-1), E0set)
            }
            else { 
              LI0 = intersect(which(apply(PartOrdAll[I0,1:k],2,max)==1), E0set)
              UI0 = intersect(which(apply(PartOrdAll[I0,1:k],2,min)==-1), E0set)
            }
            I1 = unique(c(I0,LI0,UI0))
            while (T) {   ## Evaluate I_{\infty}
              I0 = I1
              if (length(I0)==1) {
                LI0 = intersect(which(PartOrdAll[I0,1:k]==1), E0set)
                UI0 = intersect(which(PartOrdAll[I0,1:k]==-1), E0set)
              }
              else { 
                LI0 = intersect(which(apply(PartOrdAll[I0,1:k],2,max)==1), E0set)
                UI0 = intersect(which(apply(PartOrdAll[I0,1:k],2,min)==-1), E0set)
              }
              I1 = unique(c(I0,LI0,UI0))				
              if (length(I0)==length(I1)) { if (all(sort(I0)==sort(I1))) break; }
            }    
            E0setB = setdiff(I0,U0set)
            E0setA = setdiff(E0set,E0setB)
            hat1b[E0setA] = hat0[E0setA]
            isize[k] = length(E0setB)     	
            if (isize[k]>0) {
              if (isize[k]==1) { hat1b[E0setB] = bayeshat[E0setB]; }
              else if (isize[k]<=K0) {	
                dosesEset = as.matrix(alldoses[E0setB,])
                SCENESeset = getScenesV3(dosesEset)$SCENES
                logprob1 = rep(NA, ncol(SCENESeset))
                for (iter in 1:ncol(SCENESeset)) {
                  curscene = SCENESeset[,iter]
                  logprob1[iter] = sum(M[E0setB]*curscene*log(pt[E0setB]*eps) + (M[E0setB]*(1-curscene)*log( (1-pt[E0setB])*(1-eps))))  
                }
                posBest1 = which(logprob1==max(logprob1))[1]
                hat1b[E0setB] = SCENESeset[,posBest1]
              }
              else {
                recOBJ = QuickUpTog(alldoses[E0setB,], pt[E0setB], M[E0setB],eps,depth=dep+1)
                hat1b[E0setB] = recOBJ$combhat
              }
            }
          }	
        }  # End calculationg hat1b
        logH1b = sum(M[1:k]*hat1b*log(pt[1:k]) + (M[1:k]*(1-hat1b)*log(1-pt[1:k]))) + sum(M[1:k]*hat1b*log(eps)) + sum(M[1:k]*(1-hat1b)*log(1-eps))
        if (logH1a>=logH1b) hat1 = hat1a
        else hat1 = hat1b
      }    # End calculating hat1
      hat0 = hat1
    }
    hatcomb = hat0
  }
  hatcomb0 =hatcomb
  hatcomb0[o] = hatcomb
  logH = sum(M*hatcomb*log(pt*eps) + (M*(1-hatcomb)*log((1-pt)*(1-eps))))
  obj = list(combhat = hatcomb0, alldoses=alldoses0,pt=pt0,M=M0, esize=esize, isize=isize, curdepth = dep,logH = logH)
}


#' @keywords internal
QuickUpTog = function(alldoses0, pt0, M0, eps=0.5, K0=10, depth=1) {
  ### Initialization
  K = nrow(alldoses0)
  dep = depth
  totX = apply(alldoses0,1,sum)
  o = order(totX,pt0,decreasing=F)
  alldoses = as.matrix(alldoses0[o,])
  
  M = M0[o]
  pt = pt0[o]
  bayeshat = rep(0,K)
  bayeshat[pt > (1-eps)] = 1  
  PartOrdAll = getOrd(alldoses)

  ### Algorithm starts  
  sset = 1:min(K0,K)
  dosesSset = as.matrix(alldoses[sset,])
  SCENESsset = getScenesV3(dosesSset)$SCENES
  logprob1 = rep(NA, ncol(SCENESsset))
  for (iter in 1:ncol(SCENESsset)) {
    curscene = SCENESsset[,iter]
    logprob1[iter] = sum(M[sset]*curscene*log(pt[sset]*eps) + (M[sset]*(1-curscene)*log((1-pt[sset])*(1-eps)))) 
  }
  posBest1 = which(logprob1==max(logprob1))[1]
  hat0 = SCENESsset[,posBest1]

  isize = esize = rep(0,K)
  if (K<=K0) { hatcomb = hat0; }
  else {
#    if (dep==1) cat("\nUD Size",K,"iteration: ")
    for (k in (K0+1):K) {
#      if (dep==1) cat(round(k/K,digits=3),"-->")
      if (bayeshat[k]) { hat1 = c(hat0, 1); } 
      else {
        hat1a = c(hat0,1)
        logH1a = sum(M[1:k]*hat1a*log(pt[1:k]*eps) + (M[1:k]*(1-hat1a)*log((1-pt[1:k])*(1-eps))))
        
        ### Calculate the second candidate
        L1set = which(PartOrdAll[k,1:(k-1)]==1 & hat0==1)
        L0set = which(PartOrdAll[k,1:(k-1)]==1 & hat0==0)
        E1set = which(PartOrdAll[k,1:(k-1)]==0 & hat0==1)
        E0set = which(PartOrdAll[k,1:(k-1)]==0 & hat0==0)	
        hat1b = rep(NA,k)
        hat1b[k] = 0
        hat1b[L1set] = hat1b[E0set] = hat1b[L0set] = 0
        
        if (length(E1set)>0) {
          esize[k] = length(E1set)
          if (length(L1set)==0) { hat1b[E1set] = 1; }
          else if (length(E1set)==1) { hat1b[E1set] = bayeshat[E1set]; isize[k] = 1; } 
          else {  
            I0 = L1set
            if (length(I0)==1) {
              LI0 = intersect(which(PartOrdAll[I0,1:k]==1), E1set)
              UI0 = intersect(which(PartOrdAll[I0,1:k]==-1), E1set)
            }
            else { 
              LI0 = intersect(which(apply(PartOrdAll[I0,1:k],2,max)==1), E1set)
              UI0 = intersect(which(apply(PartOrdAll[I0,1:k],2,min)==-1), E1set)
            }
            I1 = unique(c(I0,LI0,UI0))
            while (T) {   ## Evaluate I_{\infty}
              I0 = I1
              if (length(I0)==1) {
                LI0 = intersect(which(PartOrdAll[I0,1:k]==1), E1set)
                UI0 = intersect(which(PartOrdAll[I0,1:k]==-1), E1set)
              }
              else { 
                LI0 = intersect(which(apply(PartOrdAll[I0,1:k],2,max)==1), E1set)
                UI0 = intersect(which(apply(PartOrdAll[I0,1:k],2,min)==-1), E1set)
              }
              I1 = unique(c(I0,LI0,UI0))				
              if (length(I0)==length(I1)) { if (all(sort(I0)==sort(I1))) break; }
            }    
            E1setB = setdiff(I0,L1set)
            E1setA = setdiff(E1set,E1setB)
            hat1b[E1setA] = hat0[E1setA]
            isize[k] = length(E1setB)     	
            if (isize[k]>0) {
              if (isize[k]==1) { hat1b[E1setB] = bayeshat[E1setB]; }
              else if (isize[k]<=K0) {	
                dosesEset = as.matrix(alldoses[E1setB,])
                SCENESeset = getScenesV3(dosesEset)$SCENES
                logprob1 = rep(NA, ncol(SCENESeset))
                for (iter in 1:ncol(SCENESeset)) {
                  curscene = SCENESeset[,iter]
                  logprob1[iter] = sum(M[E1setB]*curscene*log(pt[E1setB]*eps) + (M[E1setB]*(1-curscene)*log( (1-pt[E1setB])*(1-eps))))  
                }
                posBest1 = which(logprob1==max(logprob1))[1]
                hat1b[E1setB] = SCENESeset[,posBest1]
              }
              else {
                recOBJ = QuickDownTog(alldoses[E1setB,], pt[E1setB], M[E1setB],eps,depth=dep+1)
                hat1b[E1setB] = recOBJ$combhat
              }
            }
          }	
        }  # End calculationg hat1b
        logH1b = sum(M[1:k]*hat1b*log(pt[1:k]) + (M[1:k]*(1-hat1b)*log(1-pt[1:k]))) + sum(M[1:k]*hat1b*log(eps)) + sum(M[1:k]*(1-hat1b)*log(1-eps))
        if (logH1a>=logH1b) hat1 = hat1a
        else hat1 = hat1b
      }    # End calculating hat1
      hat0 = hat1
    }
    hatcomb = hat0
  }
  hatcomb0 =hatcomb
  hatcomb0[o] = hatcomb
  logH = sum(M*hatcomb*log(pt*eps) + (M*(1-hatcomb)*log((1-pt)*(1-eps))))
  obj = list(combhat = hatcomb0, alldoses=alldoses0,pt=pt0,M=M0, esize=esize, isize=isize, curdepth = dep,logH = logH)
}


#' Fit Bayesian misclassification model (binary)
#'
#' @param X numeric matrix
#' @param y numeric response vector
#' @param incr numeric, increment for threshold grid
#' @return A list containing fitted parameters
#' @examples
#' A <- as.matrix(expand.grid(rep(list(0:1), 6)))
#' set.seed(2025)
#' X <- A[sample(nrow(A),size=500, replace = TRUE),]
#' y <- as.numeric(rowSums(X)>=3)
#' miso(X,y)
#' @references
#' Cheung YK, Diaz KM. Monotone response surface of multi-factor condition: estimation and Bayes classifiers. 
#' *J R Stat Soc Series B Stat Methodol.* 2023 Apr;85(2):497-522. 
#' doi: 10.1093/jrsssb/qkad014. Epub 2023 Mar 22. PMID: 38464683; PMCID: PMC10919322.
#' @export
miso = function(X,y, incr = 0.01) {
  dosefit = getDoses(X,y)
  M = dosefit$M
  S = dosefit$S
  alldoses = dosefit$alldoses
  nt = round(1/incr) + 1
  output = SweepCombTogBinom(alldoses,M,S,nt=nt)
  list(alldoses = alldoses, M = M, S = S, thetahat = output$thetahat, nt = nt, logH = output$logH)
}


#' @keywords internal
SweepCombTogBinom = function(alldoses,M,S,eps=0.5,a0=0.25,b0=0.25, nt=11, K0=10, PartOrdAll=NULL, sweepoutLowE=NULL, sweepoutHighE=NULL,showprogress=FALSE) {
  aPost = a0 + S
  bPost = b0 + M - S
  K = nrow(alldoses)
  if (is.null(PartOrdAll)) PartOrdAll = getOrd(alldoses)
  sweepout = matrix(rep(NA, nt*K), nrow=K)
  logH = ctsize = lztsize = ztsize = dtsize = rep(NA,nt)
  thor = seq(0,1,length=nt)
  comb = rep(NA,nt)
  for (i in 1:nt) {
    thres = thor[i]
    pt = 1 - stats::pbeta(thres, aPost,bPost)
    pt = pmin( pmax(pt, 1e-5), 1-1e-5)
    if (all(pt>(1-eps))) {
      sweepout[,i] = rep(1,K)
      ctsize[i] = lztsize[i] = ztsize[i] = dtsize[i] = 0
    }
    else {
#     	cat("Threshold, t:",thres,"\n")
      zset = which(pt<=(1-eps))
      ztsize[i] = nz = length(zset)
      Lset = NULL
      for (j in 1:nz) {
        curRow = zset[j]
        Lset = c(Lset, which(PartOrdAll[curRow,]==1))
      }
      cset = unique(c(zset,Lset))
      lztsize[i] = length(cset)
      sweepout[-cset,i] = rep(1,(K-length(cset)))
      knownPos = which(sweepout[cset,i]==0)
      dtsize[i] = length(knownPos)
      if (length(knownPos)>0) { cset = cset[-knownPos] }
      ctsize[i] = length(cset)
#      if (showprogress) cat(thres, ctsize[i],"\n")
      if (length(cset)>0) {
        if (length(cset)==1) sweepout[cset,i] = 0
        else if (length(cset)<=K0){
          dosesCset = as.matrix(alldoses[cset,])
          SCENEScset = getScenesV3(dosesCset)$SCENES
          logprob1 = rep(NA, ncol(SCENEScset))
          for (iter in 1:ncol(SCENEScset)) {
            curscene = SCENEScset[,iter]
            logprob1[iter] = sum(M[cset]*curscene*log(pt[cset]) + (M[cset]*(1-curscene)*log(1-pt[cset]))) + 
              sum(M[cset]*curscene*log(eps)) + sum(M[cset]*(1-curscene)*log(1-eps))
          }
          posBest1 = which(logprob1==max(logprob1))[1]
          sceneHatCset = SCENEScset[,posBest1]
          sweepout[cset,i] = sceneHatCset		
        }
        else {		
          mpt = mean(pt[cset])
          if (mpt<0.5) { obj = QuickUpTog(as.matrix(alldoses[cset,]),pt[cset],M[cset],eps,K0=K0); comb[i] = "up"; }
          else { obj = QuickDownTog(as.matrix(alldoses[cset,]),pt[cset],M[cset],eps,K0=K0); comb[i] = "down"; }
          sweepout[cset,i] = obj$combhat
        }
        pos = which(sweepout[,i]==0)
        if (i<nt) sweepout[pos,(i+1):nt] = 0 		
      }
    }
    logH[i] = sum(M*sweepout[,i]*log(pt) + (M*(1-sweepout[,i])*log(1-pt))) + sum(M*sweepout[,i]*log(eps)) + sum(M*(1-sweepout[,i])*log(1-eps))
  }
  thetahat = rep(NA,K)
  for (k in 1:K) {
  	pos = which(sweepout[k,]==1)
	  thetahat[k] = thor[max(pos)]
  }
  obj = list(alldoses=alldoses,M=M,S=S,thetahat=thetahat,sweepout=sweepout, eps=eps,a0=a0,b0=b0,tgrid=thor, ztsize=ztsize, ctsize=ctsize, lztsize=lztsize,dtsize=dtsize,logH=logH,comb.directions=comb)
}


#' S3 predict method for class "pbc"
#'
#' @param object object of class "pbc"
#' @param Xnew numeric matrix of inputs
#' @param ... additional arguments (not used)
#' @return List containing predictions
#' @examples
#' A <- as.matrix(expand.grid(rep(list(0:1), 6)))
#' set.seed(2025)
#' X <- A[sample(nrow(A),size=500, replace = TRUE),]
#' y <- as.numeric(rowSums(X)>=3)
#' fit <- PBclassifier(X,y)
#' predict(fit,X)
#' @exportS3Method predict pbc
predict.pbc = function(object, Xnew, ...) {
  X1 = Xnew
  X0 = object$alldoses
  decimals = 6
  idx <- match(
    interaction(round(as.data.frame(X1), decimals), drop = TRUE),
    interaction(round(as.data.frame(X0), decimals), drop = TRUE)
  )
  fit1 = as.numeric(object$yhat[idx])
  msg = rep("Unsure - need more training data",nrow(Xnew))
  msg[!is.na(fit1)] = "exact match"
  
  if (any(is.na(fit1))) {
    pos = which(is.na(fit1))
    for (i in 1:length(pos)) {
      v = as.numeric(X1[pos[i],])
      resL = rowSums(t(t(X0)<=v))==ncol(X0)
      resU = rowSums(t(t(X0)>=v))==ncol(X0)
      lb = max(object$yhat[resL])
      ub = min(object$yhat[resU])
      if (lb==1) {
        fit1[pos[i]] = 1
        msg[pos[i]] = "bound from below"
      }
      else if (ub==0) {
        fit1[pos[i]] = 0
        msg[pos[i]] = "bound from above"
      }
    }
  }
  list(fit=fit1, msg=msg)
}

#' @keywords internal
smallest1 = function(obj) {
 alldoses = obj$alldoses
 partOrd = getOrd(alldoses)
 K = nrow(alldoses)
 val = rep(0, K)
 for (i in 1:K) {
  if (obj$yhat[i]==1) {
    Lpos = which(partOrd[i]==1)
    if (all(obj$yhat[Lpos]==0)) val[i] = 1
  }   
 }
 alldoses[which(val==1),]
}

#' @keywords internal
mcmiso = function(X,y, incr = 0.01) {
  dosefit = getDoses(X,y)
  M = dosefit$M
  S = dosefit$S
  alldoses = dosefit$alldoses
  nt = round(1/incr) + 1
  output = SweepMcCombBinom(alldoses,M,S,nt=nt)
  list(alldoses = alldoses, M = M, S = S, thetahat = output$thetahat, nt = nt, logH = output$logH)
}

#' @keywords internal
SweepMcCombBinom = function(alldoses,M,S,eps=0.5,a0=0.25,b0=0.25, nt=11, K0=10, PartOrdAll=NULL, sweepoutLowE=NULL, sweepoutHighE=NULL,showprogress=FALSE) {
  aPost = a0 + S
  bPost = b0 + M - S
  K = nrow(alldoses)
  if (is.null(PartOrdAll)) PartOrdAll = getOrd(alldoses)
  sweepout = matrix(rep(NA, nt*K), nrow=K)
  logH = ctsize = lztsize = ztsize = dtsize = rep(NA,nt)
  thor = seq(0,1,length=nt)
  comb = rep(NA,nt)
  for (i in 1:nt) {
    thres = thor[i]
    pt = 1 - stats::pbeta(thres, aPost,bPost)
    pt = pmin( pmax(pt, 1e-5), 1-1e-5)
    if (all(pt>(1-eps))) {
      sweepout[,i] = rep(1,K)
      ctsize[i] = lztsize[i] = ztsize[i] = dtsize[i] = 0
    }
    else {
#      cat("Threshold, t:",thres,"\n")
      zset = which(pt<=(1-eps))
      ztsize[i] = nz = length(zset)
      Lset = NULL
      for (j in 1:nz) {
        curRow = zset[j]
        Lset = c(Lset, which(PartOrdAll[curRow,]==1))
      }
      cset = unique(c(zset,Lset))
      lztsize[i] = length(cset)
      sweepout[-cset,i] = rep(1,(K-length(cset)))
      knownPos = which(sweepout[cset,i]==0)
      dtsize[i] = length(knownPos)
      if (length(knownPos)>0) { cset = cset[-knownPos] }
      ctsize[i] = length(cset)
#      if (showprogress) cat(thres, ctsize[i],"\n")
      if (length(cset)>0) {
        if (length(cset)==1) sweepout[cset,i] = 0
        else if (length(cset)<=K0){
          #          dosesCset = alldoses[cset,]
          dosesCset = as.matrix(alldoses[cset,])
          SCENEScset = getScenesV3(dosesCset)$SCENES
          logprob1 = rep(NA, ncol(SCENEScset))
          for (iter in 1:ncol(SCENEScset)) {
            curscene = SCENEScset[,iter]
            logprob1[iter] = sum(M[cset]*curscene*log(pt[cset]) + (M[cset]*(1-curscene)*log(1-pt[cset]))) + 
              sum(M[cset]*curscene*log(eps)) + sum(M[cset]*(1-curscene)*log(1-eps))
          }
          posBest1 = which(logprob1==max(logprob1))[1]
          sceneHatCset = SCENEScset[,posBest1]
          sweepout[cset,i] = sceneHatCset		
        }
        else {
          tmpobj = mcComb(as.matrix(alldoses[cset,]), pt[cset], M[cset])
          obj = tmpobj$val
          comb[i] = tmpobj$winner
          sweepout[cset,i] = obj$combhat
        }
        pos = which(sweepout[,i]==0)
        if (i<nt) sweepout[pos,(i+1):nt] = 0 		
      }
    }
    logH[i] = sum(M*sweepout[,i]*log(pt) + (M*(1-sweepout[,i])*log(1-pt))) + sum(M*sweepout[,i]*log(eps)) + sum(M*(1-sweepout[,i])*log(1-eps))
  }
  thetahat = rep(NA,K)
  for (k in 1:K) {
    pos = which(sweepout[k,]==1)
    thetahat[k] = thor[max(pos)]
  }
  obj = list(alldoses=alldoses,M=M,S=S,thetahat=thetahat,sweepout=sweepout, eps=eps,a0=a0,b0=b0,tgrid=thor, ztsize=ztsize, ctsize=ctsize, lztsize=lztsize,dtsize=dtsize,logH=logH,comb.directions=comb)
}


#' @keywords internal
mcComb = function(alldoses,pt,M) {
  future::plan("multicore")
  #  plan(multisession)
  f1 = future::future( { list(val = QuickDownTog(alldoses,pt,M), winner = "DU") })
  f2 = future::future( { list(val = QuickUpTog(alldoses,pt,M), winner = "UD") })
  
  repeat {
    ready <- c(future::resolved(f1), future::resolved(f2))
    if (any(ready)) {
      winner <- if (ready[1]) f1 else f2
      result <- future::value(winner)
      if (!ready[1]) try(future::cancel(f1), silent=TRUE)
      if (!ready[2]) try(future::cancel(f2), silent=TRUE)
      return(result)
    }
    Sys.sleep(0.01)  # small pause to avoid busy waiting
  }
}






