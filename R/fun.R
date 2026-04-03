#' @keywords internal
getDoses = function(X,y) {
  DOSES = as.matrix(unique(data.frame(X)))
  K = nrow(DOSES)
  rownames(DOSES) = paste("condition ",1:K,sep="")
  D = ncol(DOSES)
  M = rep(0,K)
  S = rep(0,K)        ## binomial sufficient stats
  MY = VY = rep(NA,K)  ## normal sufficient stats
  for (k in 1:K) {
    vec = DOSES[k,]
    A = t(t(X)==vec)
    pos = which(rowSums(A)==D)
    M[k] = length(pos)
    S[k] = sum(y[pos])
    MY[k] = mean(y[pos])
    if (length(pos)>1) VY[k] = stats::var(y[pos])
  }
  list(alldoses=DOSES, M=M, S=S, MY=MY, VY=VY)
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

  lower = lapply(1:K, function(k) which(PartOrd[k,] ==  1))  # doses lower than k
  higher = lapply(1:K, function(k) which(PartOrd[k,] == -1)) # doses higher than k
  
  SCENES = vector("list", 1000)
  count = 0
  assignment = integer(K)
  
  recurse = function(k) {
    if (k > K) {
      count <<- count + 1
      if (count > length(SCENES)) SCENES <<- c(SCENES, vector("list", 1000))
      SCENES[[count]] <<- assignment
      return()
    }
    for (val in 0:1) {
      assignment[k] <<- val
      valid = TRUE
      if (val == 0) {
        # If k=0, no LOWER dose (already assigned) can be 1
        already_assigned_lower = lower[[k]][lower[[k]] < k]
        if (length(already_assigned_lower) > 0 &&
            any(assignment[already_assigned_lower] == 1)) valid = FALSE
      } else {
        # If k=1, no HIGHER dose (already assigned) can be 0
        already_assigned_higher = higher[[k]][higher[[k]] < k]
        if (length(already_assigned_higher) > 0 &&
            any(assignment[already_assigned_higher] == 0)) valid = FALSE
      }
      if (valid) recurse(k + 1)
    }
  }
  recurse(1)
  SCENES_mat = do.call(cbind, SCENES[1:count])
  list(PartOrd = PartOrd, SCENES = SCENES_mat)
}

#' Multivariable Isotonic Classification using Projective Bayes
#'
#' Estimates a monotone binary classification rule for multivariate features
#' using a projective Bayes classifier. The classifier is obtained by projecting
#' an unconstrained nonparametric Bayes estimator onto the partial ordering
#' subspace defined by the assumption that the outcome probability is
#' nondecreasing in each feature. The projection is computed using a recursive
#' sequential update algorithm that yields the exact Bayes solution maximizing
#' the posterior gain. Both a down-up (\code{"DU"}) and an up-down (\code{"UD"})
#' algorithm are available.
#' @param X a numeric matrix of observed feature combinations, one row per
#'   observation, where repeated rows are expected. Each column represents a
#'   feature (e.g., a dose component or experimental factor) and each row
#'   represents the feature combination observed for one unit.
#' @param y a binary numeric vector of length \code{nrow(X)} indicating the
#'   observed outcome for each observation (1 = event, 0 = no event).
#' @param method a character string specifying the search strategy for finding
#'   the optimal monotone classification, either \code{"DU"} (down-up) or
#'   \code{"UD"} (up-down). Defaults to \code{"DU"}.
#' @param a0 a positive numeric value specifying the shape1 hyperparameter of
#'   the Beta prior in the Beta-Binomial conjugate model. Defaults to 0.5
#'   (Jeffreys prior).
#' @param b0 a positive numeric value specifying the shape2 hyperparameter of
#'   the Beta prior in the Beta-Binomial conjugate model. Defaults to 0.5
#'   (Jeffreys prior).
#' @param t0 a numeric value in (0,1) specifying the threshold on the response
#'   rate used to classify each feature combination as event or no-event.
#'   Defaults to 0.5.
#' @return A list of class \code{"pbc"} containing the following components:
#'   \describe{
#'     \item{alldoses}{a numeric matrix of unique feature combinations}
#'     \item{M}{a numeric vector of observation counts at each feature combination}
#'     \item{S}{a numeric vector of event counts at each feature combination}
#'     \item{yhat}{a binary numeric vector of the optimal monotone classification
#'       for each feature combination (1 = event, 0 = no event)}
#'     \item{pt}{a numeric vector of posterior probabilities that the true response
#'       rate exceeds \code{t0} at each feature combination}
#'     \item{logH}{the log posterior probability of the optimal classification}
#'   }
#' @examples
#' A <- as.matrix(expand.grid(rep(list(0:1), 6)))
#' set.seed(2025)
#' X <- A[sample(nrow(A),size=500, replace = TRUE),]
#' y <- as.numeric(rowSums(X)>=3)
#' PBclassifier(X,y)
#' @references
#' Cheung YK, Diaz KM. Monotone response surface of multi-factor condition: estimation and Bayes classifiers.
#' \emph{J R Stat Soc Series B Stat Methodol.} 2023 Apr;85(2):497-522.
#' doi: 10.1093/jrsssb/qkad014. Epub 2023 Mar 22. PMID: 38464683; PMCID: PMC10919322.
#'
#' Cheung YK, Kuhn L. Evaluating multiplex diagnostic test using partially ordered Bayes classifier.
#' \emph{Ann Appl Stat.} In press.
#' @export
PBclassifier = function(X, y, method="DU", a0=0.5, b0=0.5, t0=0.5) {
  
  # Input validation
  if (!is.matrix(X)) stop("'X' must be a matrix")
  if (!is.numeric(X)) stop("'X' must be numeric")
  if (!is.numeric(y)) stop("'y' must be a numeric vector")
  if (length(y) != nrow(X)) stop("length of 'y' must equal nrow(X)")
  if (!all(y %in% c(0, 1))) stop("'y' must be binary (0 or 1)")
  if (anyNA(X)) stop("'X' must not contain missing values")
  if (anyNA(y)) stop("'y' must not contain missing values")
  if (!method %in% c("DU", "UD")) stop("'method' must be either 'DU' or 'UD'")
  if (length(a0) != 1 || !is.numeric(a0) || a0 <= 0) stop("'a0' must be a single positive number")
  if (length(b0) != 1 || !is.numeric(b0) || b0 <= 0) stop("'b0' must be a single positive number")
  if (length(t0) != 1 || !is.numeric(t0) || t0 <= 0 || t0 >= 1) stop("'t0' must be a single number in (0,1)")
  
  dosefit = getDoses(X, y)
  M = dosefit$M
  S = dosefit$S
  alldoses = dosefit$alldoses
  pt = 1 - stats::pbeta(t0, a0+S, b0+M-S)
  pt = pmin(pmax(pt, 1e-5), 1-1e-5)
  if (method=="DU") output = QuickDownTog(alldoses, pt, M)
  else if (method=="UD") output = QuickUpTog(alldoses, pt, M)
  else stop("method not recognized")
  val = list(alldoses=alldoses, M=M, S=S, yhat=output$combhat, pt=pt, logH=output$logH)
  class(val) = "pbc"
  return(val)
}

#' @keywords internal
QuickDownTog = function(alldoses0, pt0, M0, eps=0.5, K0=10, depth=1) {
  ### Initialization
  K = nrow(alldoses0)
  dep = depth
  totX = apply(alldoses0,1,sum)
  o = order(totX,pt0,decreasing=TRUE)
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
    for (k in (K0+1):K) {
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
            while (TRUE) {   ## Evaluate I_{\infty}
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
  return(list(combhat = hatcomb0, alldoses=alldoses0,pt=pt0,M=M0, esize=esize, isize=isize, curdepth = dep,logH = logH))
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
    for (k in (K0+1):K) {
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
            while (TRUE) {   ## Evaluate I_{\infty}
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
  return(list(combhat = hatcomb0, alldoses=alldoses0,pt=pt0,M=M0, esize=esize, isize=isize, curdepth = dep,logH = logH))
}

#' Multivariable Isotonic Regression for Binary Data using Inverse Projective Bayes
#'
#' Estimates the underlying response probability for multivariate features
#' using an inverse projective Bayes approach. The estimator is obtained by
#' inverting the projective Bayes classifier across a grid of thresholds,
#' yielding a monotone nonparametric estimate of the response probability
#' that is nondecreasing in each feature. A Beta-Binomial conjugate model
#' is used to compute posterior probabilities at each threshold.
#'
#' @param X a numeric matrix of observed feature combinations, one row per
#'   observation, where repeated rows are expected. Each column represents a
#'   feature (e.g., a dose component or experimental factor) and each row
#'   represents the feature combination observed for one unit.
#' @param y a binary numeric vector of length \code{nrow(X)} indicating the
#'   observed outcome for each observation (1 = event, 0 = no event).
#' @param incr a numeric value in (0,1) specifying the increment between
#'   threshold grid points used to invert the classifier. Smaller values
#'   yield finer resolution at the cost of increased computation time.
#'   Defaults to 0.01.
#' @return A list containing the following components:
#'   \describe{
#'     \item{alldoses}{a numeric matrix of unique feature combinations observed
#'       in the training data}
#'     \item{M}{a numeric vector of observation counts at each feature combination}
#'     \item{S}{a numeric vector of event counts at each feature combination}
#'     \item{thetahat}{a numeric vector of estimated response probabilities at
#'       each unique feature combination, monotone nondecreasing with respect
#'       to the partial ordering of the features}
#'     \item{nt}{an integer giving the number of threshold grid points used,
#'       equal to \code{round(1/incr) + 1}}
#'     \item{logH}{a numeric vector of length \code{nt} giving the log posterior
#'       gain of the optimal classification at each threshold grid point}
#'   }
#' @examples
#' A <- as.matrix(expand.grid(rep(list(0:1), 6)))
#' set.seed(2025)
#' X <- A[sample(nrow(A), size=500, replace=TRUE),]
#' y <- as.numeric(rowSums(X) >= 3)
#' miso(X, y)
#' @references
#' Cheung YK, Diaz KM. Monotone response surface of multi-factor condition:
#' estimation and Bayes classifiers.
#' \emph{J R Stat Soc Series B Stat Methodol.} 2023 Apr;85(2):497-522.
#' doi: 10.1093/jrsssb/qkad014. Epub 2023 Mar 22. PMID: 38464683; PMCID: PMC10919322.
#' @export
miso = function(X, y, incr=0.01) {
  # Input validation
  if (!is.matrix(X)) stop("'X' must be a matrix")
  if (!is.numeric(X)) stop("'X' must be numeric")
  if (!is.numeric(y)) stop("'y' must be a numeric vector")
  if (length(y) != nrow(X)) stop("length of 'y' must equal nrow(X)")
  if (!all(y %in% c(0,1))) stop("'y' must be binary (0 or 1)")
  if (anyNA(X)) stop("'X' must not contain missing values")
  if (anyNA(y)) stop("'y' must not contain missing values")
  if (length(incr)!=1 || !is.numeric(incr) || incr<=0 || incr>=1) stop("'incr' must be a single number in (0,1)")
  
  dosefit = getDoses(X, y)
  M = dosefit$M
  S = dosefit$S
  alldoses = dosefit$alldoses
  nt = round(1/incr) + 1
  output = SweepCombTogBinom(alldoses, M, S, nt=nt)
  val = list(alldoses=alldoses, M=M, S=S, thetahat=output$thetahat, nt=nt, logH=output$logH)
  class(val) = "miso"
  return(val)
}


#' Multivariable Isotonic Regression for Binary Data using Inverse Projective Bayes
#' with Parallel Computing
#'
#' A parallel computing wrapper for \code{\link{miso}} that uses
#' \code{mcComb} to run the down-up (\code{"DU"}) and up-down (\code{"UD"})
#' algorithms simultaneously at each threshold grid point, returning the
#' result of whichever finishes first. While the two algorithms may produce
#' different classification rules at each threshold, they are guaranteed to
#' achieve the same maximum log posterior gain (Cheung and Kuhn, in press).
#' This approach reduces elapsed time without sacrificing accuracy and is
#' most effective for large datasets where the number of unique feature
#' combinations \code{K} is large (typically \code{K > 500}).
#'
#' Before calling this function, the user must set up a parallel plan using
#' \code{future::plan}. The recommended plan is \code{future::multisession}
#' which works across all platforms including Windows. The \code{future}
#' package must be installed separately (\code{install.packages("future")}).
#' After the analysis is complete, it is good practice to restore the default
#' sequential plan using \code{future::plan(future::sequential)}.
#'
#' Note that parallel overhead may outweigh the benefit for small datasets.
#' When run from RStudio, \code{future::multisession} is used automatically
#' but incurs higher overhead due to session launching costs compared to
#' \code{future::multicore} available in a terminal.
#'
#' @param X a numeric matrix of observed feature combinations, one row per
#'   observation, where repeated rows are expected. Each column represents a
#'   feature (e.g., a dose component or experimental factor) and each row
#'   represents the feature combination observed for one unit.
#' @param y a binary numeric vector of length \code{nrow(X)} indicating the
#'   observed outcome for each observation (1 = event, 0 = no event).
#' @param incr a numeric value in (0,1) specifying the increment between
#'   threshold grid points used to invert the classifier. Smaller values
#'   yield finer resolution at the cost of increased computation time.
#'   Defaults to 0.01.
#' @return A list containing the same components as \code{\link{miso}},
#'   based on whichever of the \code{"DU"} or \code{"UD"} algorithm finishes
#'   first at each threshold grid point:
#'   \describe{
#'     \item{alldoses}{a numeric matrix of unique feature combinations observed
#'       in the training data}
#'     \item{M}{a numeric vector of observation counts at each feature combination}
#'     \item{S}{a numeric vector of event counts at each feature combination}
#'     \item{thetahat}{a numeric vector of estimated response probabilities at
#'       each unique feature combination, monotone nondecreasing with respect
#'       to the partial ordering of the features}
#'     \item{nt}{an integer giving the number of threshold grid points used,
#'       equal to \code{round(1/incr) + 1}}
#'     \item{logH}{a numeric vector of length \code{nt} giving the log posterior
#'       gain of the optimal classification at each threshold grid point}
#'   }
#' @examples
#' \dontrun{
#' # install.packages("future")  # if not already installed
#' future::plan(future::multisession)  # set up parallel plan first
#' A <- as.matrix(expand.grid(rep(list(0:1), 6)))
#' set.seed(2025)
#' X <- A[sample(nrow(A), size=500, replace=TRUE),]
#' y <- as.numeric(rowSums(X) >= 3)
#' fit <- mcmiso(X, y)
#' future::plan(future::sequential)  # restore default plan when done
#' }
#' @references
#' Cheung YK, Diaz KM. Monotone response surface of multi-factor condition:
#' estimation and Bayes classifiers.
#' \emph{J R Stat Soc Series B Stat Methodol.} 2023 Apr;85(2):497-522.
#' doi: 10.1093/jrsssb/qkad014. Epub 2023 Mar 22. PMID: 38464683; PMCID: PMC10919322.
#'
#' Cheung YK, Kuhn L. Evaluating multiplex diagnostic test using partially
#' ordered Bayes classifier. \emph{Ann Appl Stat.} In press.
#' @export
mcmiso = function(X, y, incr=0.01) {
  
  # Input validation
  if (!is.matrix(X)) stop("'X' must be a matrix")
  if (!is.numeric(X)) stop("'X' must be numeric")
  if (!is.numeric(y)) stop("'y' must be a numeric vector")
  if (length(y) != nrow(X)) stop("length of 'y' must equal nrow(X)")
  if (!all(y %in% c(0,1))) stop("'y' must be binary (0 or 1)")
  if (anyNA(X)) stop("'X' must not contain missing values")
  if (anyNA(y)) stop("'y' must not contain missing values")
  if (length(incr)!=1 || !is.numeric(incr) || incr<=0 || incr>=1) stop("'incr' must be a single number in (0,1)")
  if (!requireNamespace("future", quietly=TRUE))
    stop("package 'future' is required. Please install it with install.packages('future')")
  
  dosefit = getDoses(X, y)
  M = dosefit$M
  S = dosefit$S
  alldoses = dosefit$alldoses
  nt = round(1/incr) + 1
  output = SweepMcCombBinom(alldoses, M, S, nt=nt)
  val = list(alldoses=alldoses, M=M, S=S,
             thetahat=output$thetahat, nt=nt, logH=output$logH)
  class(val) = "miso"
  return(val)
}

#' @keywords internal
SweepMcCombBinom = function(alldoses, M, S, eps=0.5, a0=0.5, b0=0.5,
                            nt=101, K0=10, PartOrdAll=NULL) {
  aPost = a0 + S
  bPost = b0 + M - S
  K = nrow(alldoses)
  if (is.null(PartOrdAll)) PartOrdAll = getOrd(alldoses)
  sweepout = matrix(rep(NA, nt*K), nrow=K)
  logH = ctsize = lztsize = ztsize = dtsize = rep(NA,nt)
  thor = seq(0, 1, length=nt)
  comb = rep(NA,nt)
  for (i in 1:nt) {
    thres = thor[i]
    pt = 1 - stats::pbeta(thres, aPost, bPost)
    pt = pmin(pmax(pt, 1e-5), 1-1e-5)
    if (all(pt>(1-eps))) {
      sweepout[,i] = rep(1,K)
      ctsize[i] = lztsize[i] = ztsize[i] = dtsize[i] = 0
    } else {
      zset = which(pt<=(1-eps))
      ztsize[i] = length(zset)
      Lset = which(apply(PartOrdAll[zset, , drop=FALSE], 2, function(x) any(x==1)))
      cset = unique(c(zset, Lset))
      lztsize[i] = length(cset)
      sweepout[-cset,i] = rep(1,(K-length(cset)))
      knownPos = which(sweepout[cset,i]==0)
      dtsize[i] = length(knownPos)
      if (length(knownPos)>0) { cset = cset[-knownPos] }
      ctsize[i] = length(cset)
      if (length(cset)>0) {
        if (length(cset)==1) {
          sweepout[cset,i] = 0
        } else if (length(cset)<=K0) {
          dosesCset = as.matrix(alldoses[cset,])
          SCENEScset = getScenesV3(dosesCset)$SCENES
          logprob1 = rep(NA, ncol(SCENEScset))
          for (iter in 1:ncol(SCENEScset)) {
            curscene = SCENEScset[,iter]
            logprob1[iter] = sum(M[cset]*curscene*log(pt[cset]) +
                                   (M[cset]*(1-curscene)*log(1-pt[cset]))) +
              sum(M[cset]*curscene*log(eps)) +
              sum(M[cset]*(1-curscene)*log(1-eps))
          }
          posBest1 = which(logprob1==max(logprob1))[1]
          sweepout[cset,i] = SCENEScset[,posBest1]
        } else {
          tmpobj = mcComb(as.matrix(alldoses[cset,]), pt[cset], M[cset])
          comb[i] = tmpobj$winner
          sweepout[cset,i] = tmpobj$val$combhat
        }
        pos = which(sweepout[,i]==0)
        if (i<nt) sweepout[pos,(i+1):nt] = 0
      }
    }
    logH[i] = sum(M*sweepout[,i]*log(pt) +
                    (M*(1-sweepout[,i])*log(1-pt))) +
      sum(M*sweepout[,i]*log(eps)) +
      sum(M*(1-sweepout[,i])*log(1-eps))
  }
  thetahat = apply(sweepout, 1, function(row) {
    pos = which(row==1)
    if (length(pos)==0) return(NA)
    thor[max(pos)]
  })
  return(list(alldoses=alldoses, M=M, S=S, thetahat=thetahat,
              sweepout=sweepout, eps=eps, a0=a0, b0=b0, tgrid=thor,
              ztsize=ztsize, ctsize=ctsize, lztsize=lztsize,
              dtsize=dtsize, logH=logH, comb.directions=comb))
}

#' Predict Method for Projective Bayes Classifier
#'
#' Generates predictions from a fitted projective Bayes classifier for new
#' feature combinations. Exact matches to training feature combinations are
#' classified directly from the fitted classifier. For unobserved feature
#' combinations, monotonicity constraints are used to impute the classification
#' when the new combination is bounded from above or below by training
#' combinations with a determinate classification. Combinations that cannot
#' be classified by monotonicity constraints alone are flagged as indeterminate.
#'
#' @param object a fitted object of class \code{"pbc"}, produced by
#'   \code{\link{PBclassifier}}.
#' @param Xnew a numeric matrix of new feature combinations to classify, with
#'   the same number of columns as the training data.
#' @param ... additional arguments (not used).
#' @return A list containing the following components:
#'   \describe{
#'     \item{fit}{a numeric vector of predicted classifications (1 = event,
#'       0 = no event, \code{NA} = indeterminate) for each row of \code{Xnew}}
#'     \item{msg}{a character vector describing how each prediction was obtained:
#'       \code{"exact match"} for training combinations, \code{"bound from below"}
#'       or \code{"bound from above"} for combinations classified by monotonicity
#'       constraints, and \code{"Unsure - need more training data"} for
#'       indeterminate cases}
#'   }
#' @examples
#' A <- as.matrix(expand.grid(rep(list(0:1), 6)))
#' set.seed(2025)
#' X <- A[sample(nrow(A), size=500, replace=TRUE),]
#' y <- as.numeric(rowSums(X) >= 3)
#' fit <- PBclassifier(X, y)
#' predict(fit, X)
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

#' Multivariable Isotonic Classification using Projective Bayes with Parallel Computing
#'
#' A parallel computing wrapper for \code{\link{PBclassifier}} that runs the
#' down-up (\code{"DU"}) and up-down (\code{"UD"}) algorithms simultaneously
#' in parallel and returns the result of whichever finishes first. Since both
#' algorithms are guaranteed to achieve the same maximum log posterior gain,
#' i.e. the same \code{logH} (Cheung and Kuhn in press), this approach reduces
#' elapsed time without sacrificing optimality. 
#'
#' Before calling this function, the user must set up a parallel plan using
#' \code{future::plan}. The recommended plan is \code{future::multisession}
#' which works across all platforms including Windows. The \code{future}
#' package must be installed separately (\code{install.packages("future")}).
#' After the analysis is complete, it is good practice to restore the default
#' sequential plan using \code{future::plan(future::sequential)}.
#'
#' Note that parallel overhead may outweigh the benefit for small datasets.
#' This function is most effective when the number of unique feature
#' combinations \code{K} is large (typically \code{K > 500}) and when
#' run from a terminal rather than RStudio, where \code{future::multicore}
#' (forking) is available. In RStudio, \code{future::multisession} is used
#' automatically but incurs higher overhead due to session launching costs.
#' 
#' @param X a numeric matrix of observed feature combinations, one row per
#'   observation, where repeated rows are expected. Each column represents a
#'   feature (e.g., a dose component or experimental factor) and each row
#'   represents the feature combination observed for one unit.
#' @param y a binary numeric vector of length \code{nrow(X)} indicating the
#'   observed outcome for each observation (1 = event, 0 = no event).
#' @param a0 a positive numeric value specifying the shape1 hyperparameter of
#'   the Beta prior in the Beta-Binomial conjugate model. Defaults to 0.5
#'   (Jeffreys prior).
#' @param b0 a positive numeric value specifying the shape2 hyperparameter of
#'   the Beta prior in the Beta-Binomial conjugate model. Defaults to 0.5
#'   (Jeffreys prior).
#' @param t0 a numeric value in (0,1) specifying the threshold on the response
#'   rate used to classify each feature combination as event or no-event.
#'   Defaults to 0.5.
#' @return A list of class \code{"pbc"} containing the same components as
#'   \code{\link{PBclassifier}}, based on whichever of the \code{"DU"} or
#'   \code{"UD"} algorithm finishes first.
#' @examples
#' \dontrun{
#' # install.packages("future")  # if not already installed
#' future::plan(future::multisession)  # set up parallel plan first
#' A <- as.matrix(expand.grid(rep(list(0:1), 6)))
#' set.seed(2025)
#' X <- A[sample(nrow(A), size=500, replace=TRUE),]
#' y <- as.numeric(rowSums(X) >= 3)
#' fit <- mcPBclassifier(X, y)
#' future::plan(future::sequential)  # restore default plan when done
#' }
#' @references
#' Cheung YK, Diaz KM. Monotone response surface of multi-factor condition:
#' estimation and Bayes classifiers.
#' \emph{J R Stat Soc Series B Stat Methodol.} 2023 Apr;85(2):497-522.
#' doi: 10.1093/jrsssb/qkad014. Epub 2023 Mar 22. PMID: 38464683; PMCID: PMC10919322.
#'
#' Cheung YK, Kuhn L. Evaluating multiplex diagnostic test using partially
#' ordered Bayes classifier. \emph{Ann Appl Stat.} In press.
#' @export
mcPBclassifier = function(X, y, a0=0.5, b0=0.5, t0=0.5) {
  
  # Input validation
  if (!is.matrix(X)) stop("'X' must be a matrix")
  if (!is.numeric(X)) stop("'X' must be numeric")
  if (!is.numeric(y)) stop("'y' must be a numeric vector")
  if (length(y) != nrow(X)) stop("length of 'y' must equal nrow(X)")
  if (!all(y %in% c(0,1))) stop("'y' must be binary (0 or 1)")
  if (anyNA(X)) stop("'X' must not contain missing values")
  if (anyNA(y)) stop("'y' must not contain missing values")
  if (length(a0)!=1 || !is.numeric(a0) || a0<=0) stop("'a0' must be a single positive number")
  if (length(b0)!=1 || !is.numeric(b0) || b0<=0) stop("'b0' must be a single positive number")
  if (length(t0)!=1 || !is.numeric(t0) || t0<=0 || t0>=1) stop("'t0' must be a single number in (0,1)")
  
  if (!requireNamespace("future", quietly=TRUE))
    stop("package 'future' is required. Please install it with install.packages('future')")
  
  dosefit = getDoses(X, y)
  M = dosefit$M
  S = dosefit$S
  alldoses = dosefit$alldoses
  pt = 1 - stats::pbeta(t0, a0+S, b0+M-S)
  pt = pmin(pmax(pt, 1e-5), 1-1e-5)
  
  output = mcComb(alldoses, pt, M)
  
  val = list(alldoses=alldoses, M=M, S=S, yhat=output$val$combhat,
             pt=pt, logH=output$val$logH)
  class(val) = "pbc"
  return(val)
}

#' @keywords internal
mcComb = function(alldoses, pt, M) {
  # NOTE: future::plan must be set by the caller before using this function
  # DU and UD produce identical logH - run in parallel and take first to finish
  f1 = future::future({ list(val=QuickDownTog(alldoses, pt, M), winner="DU") }, seed=TRUE)
  f2 = future::future({ list(val=QuickUpTog(alldoses, pt, M), winner="UD") }, seed=TRUE)
  
  repeat {
    ready = c(future::resolved(f1), future::resolved(f2))
    if (any(ready)) {
      result = future::value(if (ready[1]) f1 else f2)
      # attempt to cancel the slower future
      return(result)
    }
    Sys.sleep(0.01)
  }
}


#' Multivariable Isotonic Regression for Continuous Data using Inverse Projective Bayes
#'
#' Estimates the underlying mean response for multivariate features using an
#' inverse projective Bayes approach. The estimator is obtained by inverting
#' the projective Bayes classifier across a grid of thresholds, yielding a
#' monotone nonparametric estimate of the mean response that is nondecreasing
#' in each feature. A Normal-Inverse-Chi-Squared conjugate model is used to
#' compute posterior probabilities at each threshold. The threshold grid is
#' determined automatically from the data range and \code{nt} controls the
#' resolution of the grid.
#'
#' The prior distribution assumes that the mean response \eqn{\mu} at each
#' feature combination follows a Normal-Inverse-Chi-Squared model. Specifically,
#' the conditional prior on \eqn{\mu} given variance \eqn{\sigma^2} is
#' \eqn{\mu | \sigma^2 \sim N(\mu_0, \sigma^2 / \kappa_0)}, and the marginal
#' prior on \eqn{\mu} is a t-distribution centered at \code{mu0} with
#' \code{nu0} degrees of freedom and scale \code{sig0}. The default values
#' \code{kap0 = nu0 = 0.01} approximate a Jeffreys non-informative prior,
#' minimizing the influence of the prior on the posterior especially for
#' feature combinations with few observations (\code{M = 1}). A more
#' informative prior can be specified by increasing \code{kap0} and \code{nu0}.
#'
#' @param X a numeric matrix of observed feature combinations, one row per
#'   observation, where repeated rows are expected. Each column represents a
#'   feature (e.g., a dose component or experimental factor) and each row
#'   represents the feature combination observed for one unit.
#' @param y a numeric vector of length \code{nrow(X)} containing the
#'   continuous outcome for each observation.
#' @param nt a positive integer specifying the number of threshold grid points
#'   used to invert the classifier. The grid range is determined automatically
#'   from the observed data. Larger values yield finer resolution at the cost
#'   of increased computation time. Defaults to 101, which is approximately
#'   equivalent to the default \code{incr = 0.01} used in \code{\link{miso}}.
#' @param mu0 a numeric value specifying the prior mean of the response.
#'   Defaults to 0.
#' @param sig0 a positive numeric value specifying the prior scale parameter,
#'   interpreted as the prior standard deviation of the response. Defaults to
#'   100, yielding a diffuse prior.
#' @param kap0 a positive numeric value specifying the prior pseudo sample size
#'   for the mean. Smaller values yield a more diffuse prior on the mean.
#'   Defaults to 0.01.
#' @param nu0 a positive numeric value specifying the prior degrees of freedom
#'   for the variance. Smaller values yield a more diffuse prior on the
#'   variance. Defaults to 0.01.
#' @return A list containing the following components:
#'   \describe{
#'     \item{alldoses}{a numeric matrix of unique feature combinations observed
#'       in the training data}
#'     \item{M}{a numeric vector of observation counts at each feature combination}
#'     \item{MY}{a numeric vector of sample means of the outcome at each feature
#'       combination (\code{NA} if \code{M = 0})}
#'     \item{VY}{a numeric vector of sample variances of the outcome at each
#'       feature combination (\code{NA} if \code{M <= 1})}
#'     \item{thetahat}{a numeric vector of estimated mean responses at each
#'       unique feature combination, monotone nondecreasing with respect to
#'       the partial ordering of the features}
#'     \item{nt}{an integer giving the number of threshold grid points used}
#'     \item{logH}{a numeric vector of length \code{nt} giving the log posterior
#'       gain of the optimal classification at each threshold grid point}
#'   }
#' @examples
#' A <- as.matrix(expand.grid(rep(list(0:1), 6)))
#' set.seed(2025)
#' X <- A[sample(nrow(A), size=500, replace=TRUE),]
#' y <- rowSums(X) + rnorm(500)
#' misoN(X, y)
#' @references
#' Cheung YK, Diaz KM. Monotone response surface of multi-factor condition:
#' estimation and Bayes classifiers.
#' \emph{J R Stat Soc Series B Stat Methodol.} 2023 Apr;85(2):497-522.
#' doi: 10.1093/jrsssb/qkad014. Epub 2023 Mar 22. PMID: 38464683; PMCID: PMC10919322.
#' @export
misoN = function(X, y, nt=101, mu0=0, sig0=100, kap0=0.01, nu0=0.01) {
  # Input validation
  if (!is.matrix(X)) stop("'X' must be a matrix")
  if (!is.numeric(X)) stop("'X' must be numeric")
  if (!is.numeric(y)) stop("'y' must be a numeric vector")
  if (length(y) != nrow(X)) stop("length of 'y' must equal nrow(X)")
  if (anyNA(X)) stop("'X' must not contain missing values")
  if (anyNA(y)) stop("'y' must not contain missing values")
  if (length(nt)!=1 || !is.numeric(nt) || nt<2 || nt!=round(nt)) stop("'nt' must be a single integer >= 2")
  if (length(mu0)!=1 || !is.numeric(mu0)) stop("'mu0' must be a single numeric value")
  if (length(sig0)!=1 || !is.numeric(sig0) || sig0<=0) stop("'sig0' must be a single positive number")
  if (length(kap0)!=1 || !is.numeric(kap0) || kap0<=0) stop("'kap0' must be a single positive number")
  if (length(nu0)!=1 || !is.numeric(nu0) || nu0<=0) stop("'nu0' must be a single positive number")
  
  dosefit = getDoses(X, y)
  M  = dosefit$M
  MY = dosefit$MY
  VY = dosefit$VY
  alldoses = dosefit$alldoses
  output = SweepCombTogNorm(alldoses, M, MY, VY, nt=nt,
                            mu0=mu0, sig0=sig0, kap0=kap0, nu0=nu0)
  val = list(alldoses=alldoses, M=M, MY=MY, VY=VY, thetahat=output$thetahat, nt=nt, logH=output$logH)
  class(val) = "misoN"
  return(val)
}

#' @keywords internal
getPTnorm = function(t0, M, MY, VY, mu0=0, sig0 = 100, kap0 = 0.01, nu0 = 0.01) {
  # pt = posterior probability of mu>t0
  # get pt for all conditions (vectorized)
  # allow M=0 or 1 for which MY and VY has a value of NA
  kapn = kap0 + M
  nun = nu0 + M
  mun = (kap0*mu0 + M*MY) / kapn
  mun[M==0] = mu0
  vy = VY
  vy[M==1] = 0
  sig2n = ( nu0*sig0^2 + (M-1)*vy + kap0*M/kapn*(MY-mu0)^2 ) / nun
  sig2n[M==0] = sig0^2
  sca = sqrt(sig2n/kapn)
  tscale = (t0-mun)/sca
  PT = 1 - stats::pt(tscale,df=nun)
  pmin( pmax(PT, 1e-5), 1-1e-5)
}

#' @keywords internal
SweepCombTogNorm = function(alldoses,M,MY,VY,eps=0.5,mu0=0,sig0=100, kap0=0.01,nu0=0.01, nt=101, K0=10, PartOrdAll=NULL) {
  K = nrow(alldoses)
  if (is.null(PartOrdAll)) PartOrdAll = getOrd(alldoses)
  sweepout = matrix(rep(NA, nt*K), nrow=K)
  logH = ctsize = lztsize = ztsize = dtsize = rep(NA,nt)
  tlb = min(MY,na.rm=TRUE)*2 - max(MY,na.rm=TRUE)
  while (TRUE) {
    pt = getPTnorm(tlb,M,MY,VY,mu0=mu0,sig0=sig0,kap0=kap0,nu0=nu0)
    if (all(pt>(1-eps))) break
    tlb = tlb - (max(MY, na.rm = TRUE)-min(MY, na.rm=TRUE))
  }
  tub = 2*max(MY, na.rm = TRUE)-min(MY, na.rm= TRUE)
  while (TRUE) {
    pt = getPTnorm(tub,M,MY,VY,mu0=mu0,sig0=sig0,kap0=kap0,nu0=nu0)
    if (all(pt<(1-eps))) break
    tub = tub + (max(MY, na.rm = TRUE)-min(MY, na.rm=TRUE))
  }
  thor = seq(tlb,tub,length=nt)
  comb = rep(NA,nt)
  for (i in 1:nt) {
    thres = thor[i]
    pt = getPTnorm(thres,M,MY,VY,mu0=mu0,sig0=sig0,kap0=kap0,nu0=nu0)
    if (all(pt>(1-eps))) {
      sweepout[,i] = rep(1,K)
      ctsize[i] = lztsize[i] = ztsize[i] = dtsize[i] = 0
    }
    else {
      zset = which(pt<=(1-eps))
      ztsize[i] = length(zset)
      Lset = which(apply(PartOrdAll[zset, , drop=FALSE], 2, function(x) any(x==1)))
      cset = unique(c(zset,Lset))
      lztsize[i] = length(cset)
      sweepout[-cset,i] = rep(1,(K-length(cset)))
      knownPos = which(sweepout[cset,i]==0)
      dtsize[i] = length(knownPos)
      if (length(knownPos)>0) { cset = cset[-knownPos] }
      ctsize[i] = length(cset)
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
  thetahat = apply(sweepout, 1, function(row) {
    pos = which(row==1)
    if (length(pos)==0) return(NA)
    thor[max(pos)]
  })
  return(list(alldoses=alldoses,thetahat=thetahat,sweepout=sweepout, eps=eps,tgrid=thor, ztsize=ztsize, ctsize=ctsize, lztsize=lztsize,dtsize=dtsize,logH=logH,comb.directions=comb))
}

#' Multivariable Isotonic Regression for Continuous Data using Inverse Projective Bayes
#' with Parallel Computing
#'
#' A parallel computing wrapper for \code{\link{misoN}} that uses
#' \code{mcComb} to run the down-up (\code{"DU"}) and up-down (\code{"UD"})
#' algorithms simultaneously at each threshold grid point, returning the
#' result of whichever finishes first. While the two algorithms may produce
#' different classification rules at each threshold, they are guaranteed to
#' achieve the same maximum log posterior gain (Cheung and Kuhn, in press).
#' This approach reduces elapsed time without sacrificing accuracy and is
#' most effective for large datasets where the number of unique feature
#' combinations \code{K} is large (typically \code{K > 500}).
#'
#' Before calling this function, the user must set up a parallel plan using
#' \code{future::plan}. The recommended plan is \code{future::multisession}
#' which works across all platforms including Windows. The \code{future}
#' package must be installed separately (\code{install.packages("future")}).
#' After the analysis is complete, it is good practice to restore the default
#' sequential plan using \code{future::plan(future::sequential)}.
#'
#' Note that parallel overhead may outweigh the benefit for small datasets.
#' When run from RStudio, \code{future::multisession} is used automatically
#' but incurs higher overhead due to session launching costs compared to
#' \code{future::multicore} available in a terminal.
#'
#' @param X a numeric matrix of observed feature combinations, one row per
#'   observation, where repeated rows are expected. Each column represents a
#'   feature (e.g., a dose component or experimental factor) and each row
#'   represents the feature combination observed for one unit.
#' @param y a numeric vector of length \code{nrow(X)} containing the
#'   continuous outcome for each observation.
#' @param nt a positive integer specifying the number of threshold grid points
#'   used to invert the classifier. The grid range is determined automatically
#'   from the observed data. Larger values yield finer resolution at the cost
#'   of increased computation time. Defaults to 101, which is approximately
#'   equivalent to the default \code{incr = 0.01} used in \code{\link{mcmiso}}.
#' @param mu0 a numeric value specifying the prior mean of the response.
#'   Defaults to 0.
#' @param sig0 a positive numeric value specifying the prior scale parameter,
#'   interpreted as the prior standard deviation of the response. Defaults to
#'   100, yielding a diffuse prior on the variance.
#' @param kap0 a positive numeric value specifying the prior pseudo sample size
#'   for the mean. Smaller values yield a more diffuse prior on the mean and
#'   reduce prior shrinkage of the posterior mean toward \code{mu0}. Defaults
#'   to 0.01, approximating a Jeffreys non-informative prior.
#' @param nu0 a positive numeric value specifying the prior degrees of freedom
#'   for the variance. Smaller values yield a more diffuse prior on the
#'   variance. Defaults to 0.01, approximating a Jeffreys non-informative prior.
#' @return A list containing the same components as \code{\link{misoN}},
#'   based on whichever of the \code{"DU"} or \code{"UD"} algorithm finishes
#'   first at each threshold grid point:
#'   \describe{
#'     \item{alldoses}{a numeric matrix of unique feature combinations observed
#'       in the training data}
#'     \item{M}{a numeric vector of observation counts at each feature combination}
#'     \item{MY}{a numeric vector of sample means of the outcome at each feature
#'       combination (\code{NA} if \code{M = 0})}
#'     \item{VY}{a numeric vector of sample variances of the outcome at each
#'       feature combination (\code{NA} if \code{M <= 1})}
#'     \item{thetahat}{a numeric vector of estimated mean responses at each
#'       unique feature combination, monotone nondecreasing with respect to
#'       the partial ordering of the features}
#'     \item{nt}{an integer giving the number of threshold grid points used}
#'     \item{logH}{a numeric vector of length \code{nt} giving the log posterior
#'       gain of the optimal classification at each threshold grid point}
#'   }
#' @examples
#' \dontrun{
#' # install.packages("future")  # if not already installed
#' future::plan(future::multisession)  # set up parallel plan first
#' A <- as.matrix(expand.grid(rep(list(0:1), 6)))
#' set.seed(2025)
#' X <- A[sample(nrow(A), size=500, replace=TRUE),]
#' y <- rowSums(X) + rnorm(500)
#' fit <- mcmisoN(X, y)
#' future::plan(future::sequential)  # restore default plan when done
#' }
#' @references
#' Cheung YK, Diaz KM. Monotone response surface of multi-factor condition:
#' estimation and Bayes classifiers.
#' \emph{J R Stat Soc Series B Stat Methodol.} 2023 Apr;85(2):497-522.
#' doi: 10.1093/jrsssb/qkad014. Epub 2023 Mar 22. PMID: 38464683; PMCID: PMC10919322.
#'
#' Cheung YK, Kuhn L. Evaluating multiplex diagnostic test using partially
#' ordered Bayes classifier. \emph{Ann Appl Stat.} In press.
#' @export
mcmisoN = function(X, y, nt=101, mu0=0, sig0=100, kap0=0.01, nu0=0.01) {
  
  # Input validation
  if (!is.matrix(X)) stop("'X' must be a matrix")
  if (!is.numeric(X)) stop("'X' must be numeric")
  if (!is.numeric(y)) stop("'y' must be a numeric vector")
  if (length(y) != nrow(X)) stop("length of 'y' must equal nrow(X)")
  if (anyNA(X)) stop("'X' must not contain missing values")
  if (anyNA(y)) stop("'y' must not contain missing values")
  if (length(nt)!=1 || !is.numeric(nt) || nt<2 || nt!=round(nt)) stop("'nt' must be a single integer >= 2")
  if (length(mu0)!=1 || !is.numeric(mu0)) stop("'mu0' must be a single numeric value")
  if (length(sig0)!=1 || !is.numeric(sig0) || sig0<=0) stop("'sig0' must be a single positive number")
  if (length(kap0)!=1 || !is.numeric(kap0) || kap0<=0) stop("'kap0' must be a single positive number")
  if (length(nu0)!=1 || !is.numeric(nu0) || nu0<=0) stop("'nu0' must be a single positive number")
  if (!requireNamespace("future", quietly=TRUE))
    stop("package 'future' is required. Please install it with install.packages('future')")
  
  dosefit = getDoses(X, y)
  M  = dosefit$M
  MY = dosefit$MY
  VY = dosefit$VY
  alldoses = dosefit$alldoses
  output = SweepMcCombNorm(alldoses, M, MY, VY, nt=nt,
                           mu0=mu0, sig0=sig0, kap0=kap0, nu0=nu0)
  val = list(alldoses=alldoses, M=M, MY=MY, VY=VY,
             thetahat=output$thetahat, nt=nt, logH=output$logH)
  class(val) = "misoN"
  return(val)
}

#' @keywords internal
SweepMcCombNorm = function(alldoses, M, MY, VY, eps=0.5, mu0=0, sig0=100,
                           kap0=0.01, nu0=0.01, nt=101, K0=10, PartOrdAll=NULL) {
  K = nrow(alldoses)
  if (is.null(PartOrdAll)) PartOrdAll = getOrd(alldoses)
  sweepout = matrix(rep(NA, nt*K), nrow=K)
  logH = ctsize = lztsize = ztsize = dtsize = rep(NA,nt)
  
  # determine data-driven threshold range
  tlb = min(MY, na.rm=TRUE)*2 - max(MY, na.rm=TRUE)
  while (TRUE) {
    pt = getPTnorm(tlb, M, MY, VY, mu0=mu0, sig0=sig0, kap0=kap0, nu0=nu0)
    if (all(pt>(1-eps))) break
    tlb = tlb - (max(MY, na.rm=TRUE) - min(MY, na.rm=TRUE))
  }
  tub = 2*max(MY, na.rm=TRUE) - min(MY, na.rm=TRUE)
  while (TRUE) {
    pt = getPTnorm(tub, M, MY, VY, mu0=mu0, sig0=sig0, kap0=kap0, nu0=nu0)
    if (all(pt<(1-eps))) break
    tub = tub + (max(MY, na.rm=TRUE) - min(MY, na.rm=TRUE))
  }
  thor = seq(tlb, tub, length=nt)
  comb = rep(NA,nt)
  
  for (i in 1:nt) {
    thres = thor[i]
    pt = getPTnorm(thres, M, MY, VY, mu0=mu0, sig0=sig0, kap0=kap0, nu0=nu0)
    if (all(pt>(1-eps))) {
      sweepout[,i] = rep(1,K)
      ctsize[i] = lztsize[i] = ztsize[i] = dtsize[i] = 0
    } else {
      zset = which(pt<=(1-eps))
      ztsize[i] = length(zset)
      Lset = which(apply(PartOrdAll[zset, , drop=FALSE], 2, function(x) any(x==1)))
      cset = unique(c(zset, Lset))
      lztsize[i] = length(cset)
      sweepout[-cset,i] = rep(1,(K-length(cset)))
      knownPos = which(sweepout[cset,i]==0)
      dtsize[i] = length(knownPos)
      if (length(knownPos)>0) { cset = cset[-knownPos] }
      ctsize[i] = length(cset)
      if (length(cset)>0) {
        if (length(cset)==1) {
          sweepout[cset,i] = 0
        } else if (length(cset)<=K0) {
          dosesCset = as.matrix(alldoses[cset,])
          SCENEScset = getScenesV3(dosesCset)$SCENES
          logprob1 = rep(NA, ncol(SCENEScset))
          for (iter in 1:ncol(SCENEScset)) {
            curscene = SCENEScset[,iter]
            logprob1[iter] = sum(M[cset]*curscene*log(pt[cset]) +
                                   (M[cset]*(1-curscene)*log(1-pt[cset]))) +
              sum(M[cset]*curscene*log(eps)) +
              sum(M[cset]*(1-curscene)*log(1-eps))
          }
          posBest1 = which(logprob1==max(logprob1))[1]
          sweepout[cset,i] = SCENEScset[,posBest1]
        } else {
          tmpobj = mcComb(as.matrix(alldoses[cset,]), pt[cset], M[cset])
          comb[i] = tmpobj$winner
          sweepout[cset,i] = tmpobj$val$combhat
        }
        pos = which(sweepout[,i]==0)
        if (i<nt) sweepout[pos,(i+1):nt] = 0
      }
    }
    logH[i] = sum(M*sweepout[,i]*log(pt) +
                    (M*(1-sweepout[,i])*log(1-pt))) +
      sum(M*sweepout[,i]*log(eps)) +
      sum(M*(1-sweepout[,i])*log(1-eps))
  }
  thetahat = apply(sweepout, 1, function(row) {
    pos = which(row==1)
    if (length(pos)==0) return(NA)
    thor[max(pos)]
  })
  return(list(alldoses=alldoses, M=M, MY=MY, VY=VY,
              thetahat=thetahat, sweepout=sweepout, eps=eps,
              mu0=mu0, sig0=sig0, kap0=kap0, nu0=nu0,
              tgrid=thor, ztsize=ztsize, ctsize=ctsize,
              lztsize=lztsize, dtsize=dtsize, logH=logH,
              comb.directions=comb))
}

#' @keywords internal
SweepCombTogBinom = function(alldoses, M, S, eps=0.5, a0=0.5, b0=0.5,
                             nt=101, K0=10, PartOrdAll=NULL) {
  aPost = a0 + S
  bPost = b0 + M - S
  K = nrow(alldoses)
  if (is.null(PartOrdAll)) PartOrdAll = getOrd(alldoses)
  sweepout = matrix(rep(NA, nt*K), nrow=K)
  logH = ctsize = lztsize = ztsize = dtsize = rep(NA,nt)
  thor = seq(0, 1, length=nt)
  comb = rep(NA,nt)
  for (i in 1:nt) {
    thres = thor[i]
    pt = 1 - stats::pbeta(thres, aPost, bPost)
    pt = pmin(pmax(pt, 1e-5), 1-1e-5)
    if (all(pt>(1-eps))) {
      sweepout[,i] = rep(1,K)
      ctsize[i] = lztsize[i] = ztsize[i] = dtsize[i] = 0
    } else {
      zset = which(pt<=(1-eps))
      ztsize[i] = length(zset)
      Lset = which(apply(PartOrdAll[zset, , drop=FALSE], 2, function(x) any(x==1)))
      cset = unique(c(zset, Lset))
      lztsize[i] = length(cset)
      sweepout[-cset,i] = rep(1,(K-length(cset)))
      knownPos = which(sweepout[cset,i]==0)
      dtsize[i] = length(knownPos)
      if (length(knownPos)>0) { cset = cset[-knownPos] }
      ctsize[i] = length(cset)
      if (length(cset)>0) {
        if (length(cset)==1) {
          sweepout[cset,i] = 0
        } else if (length(cset)<=K0) {
          dosesCset = as.matrix(alldoses[cset,])
          SCENEScset = getScenesV3(dosesCset)$SCENES
          logprob1 = rep(NA, ncol(SCENEScset))
          for (iter in 1:ncol(SCENEScset)) {
            curscene = SCENEScset[,iter]
            logprob1[iter] = sum(M[cset]*curscene*log(pt[cset]) +
                                   (M[cset]*(1-curscene)*log(1-pt[cset]))) +
              sum(M[cset]*curscene*log(eps)) +
              sum(M[cset]*(1-curscene)*log(1-eps))
          }
          posBest1 = which(logprob1==max(logprob1))[1]
          sweepout[cset,i] = SCENEScset[,posBest1]
        } else {
          mpt = mean(pt[cset])
          if (mpt<0.5) {
            obj = QuickUpTog(as.matrix(alldoses[cset,]), pt[cset], M[cset], eps, K0=K0)
            comb[i] = "up"
          } else {
            obj = QuickDownTog(as.matrix(alldoses[cset,]), pt[cset], M[cset], eps, K0=K0)
            comb[i] = "down"
          }
          sweepout[cset,i] = obj$combhat
        }
        pos = which(sweepout[,i]==0)
        if (i<nt) sweepout[pos,(i+1):nt] = 0
      }
    }
    logH[i] = sum(M*sweepout[,i]*log(pt) +
                    (M*(1-sweepout[,i])*log(1-pt))) +
      sum(M*sweepout[,i]*log(eps)) +
      sum(M*(1-sweepout[,i])*log(1-eps))
  }
  thetahat = apply(sweepout, 1, function(row) {
    pos = which(row==1)
    if (length(pos)==0) return(NA)
    thor[max(pos)]
  })
  return(list(alldoses=alldoses, M=M, S=S, thetahat=thetahat,
              sweepout=sweepout, eps=eps, a0=a0, b0=b0, tgrid=thor,
              ztsize=ztsize, ctsize=ctsize, lztsize=lztsize,
              dtsize=dtsize, logH=logH, comb.directions=comb))
}


#' Decision Boundary of a Projective Bayes Classifier
#'
#' Extracts the decision boundary from a fitted projective Bayes classifier.
#' The decision boundary consists of the minimal set of feature combinations
#' that are classified as positive (1), i.e., combinations classified as 1
#' for which no combination lower in the partial ordering is also classified
#' as 1. Any feature combination at or above a boundary combination in the
#' partial ordering is guaranteed to be classified as 1, making the boundary
#' the most compact representation of the classification rule.
#'
#' @param object a fitted object of class \code{"pbc"}, produced by
#'   \code{\link{PBclassifier}} or \code{\link{mcPBclassifier}}.
#' @return A list of class \code{"boundary"} containing the following
#'   components:
#'   \describe{
#'     \item{boundary}{a numeric matrix of feature combinations forming the
#'       decision boundary, where each row is one minimal positive combination.
#'       \code{NULL} if no combinations are classified as 1.}
#'     \item{yhat}{a binary numeric vector giving the classification of each
#'       boundary combination (all values will be 1)}
#'     \item{pt}{a numeric vector of posterior probabilities at each boundary
#'       combination}
#'     \item{n_boundary}{an integer giving the number of boundary combinations}
#'     \item{n_positive}{an integer giving the total number of combinations
#'       classified as 1}
#'     \item{n_total}{an integer giving the total number of unique feature
#'       combinations}
#'   }
#' @examples
#' A <- as.matrix(expand.grid(rep(list(0:1), 6)))
#' set.seed(2025)
#' X <- A[sample(nrow(A), size=500, replace=TRUE),]
#' y <- as.numeric(rowSums(X) >= 3)
#' fit <- PBclassifier(X, y)
#' db <- boundary(fit)
#' print(db)
#' @references
#' Cheung YK, Diaz KM. Monotone response surface of multi-factor condition:
#' estimation and Bayes classifiers.
#' \emph{J R Stat Soc Series B Stat Methodol.} 2023 Apr;85(2):497-522.
#' doi: 10.1093/jrsssb/qkad014. Epub 2023 Mar 22. PMID: 38464683; PMCID: PMC10919322.
#'
#' Cheung YK, Kuhn L. Evaluating multiplex diagnostic test using partially
#' ordered Bayes classifier. \emph{Ann Appl Stat.} In press.
#' @export
boundary = function(object) {
  
  # Input validation
  if (!inherits(object, "pbc")) stop("'object' must be of class 'pbc'")
  
  alldoses = object$alldoses
  K = nrow(alldoses)
  PartOrd = getOrd(alldoses)
  
  # handle edge cases
  if (all(object$yhat == 0)) {
    message("No combinations classified as 1 - decision boundary is empty.")
    return(structure(
      list(boundary   = NULL,
           yhat       = NULL,
           pt         = NULL,
           n_boundary = 0L,
           n_positive = 0L,
           n_total    = K),
      class = "boundary"))
  }
  
  if (all(object$yhat == 1)) {
    message("All combinations classified as 1 - no meaningful boundary exists.")
    return(structure(
      list(boundary   = NULL,
           yhat       = object$yhat,
           pt         = object$pt,
           n_boundary = 0L,
           n_positive = K,
           n_total    = K),
      class = "boundary"))
  }
  
  # find minimal positive combinations
  val = rep(0L, K)
  for (i in 1:K) {
    if (object$yhat[i] == 1) {
      Lpos = which(PartOrd[i,] == 1)
      if (length(Lpos) == 0 || all(object$yhat[Lpos] == 0)) {
        val[i] = 1L
      }
    }
  }
  
  boundary_idx = which(val == 1)
  return(structure(
    list(boundary   = alldoses[boundary_idx, , drop=FALSE],
         yhat       = object$yhat[boundary_idx],
         pt         = object$pt[boundary_idx],
         n_boundary = length(boundary_idx),
         n_positive = sum(object$yhat == 1),
         n_total    = K),
    class = "boundary"))
}

#' @export
print.pbc = function(x, ...) {
  cat("Projective Bayes Classifier\n")
  cat("===========================\n")
  cat(sprintf("Feature combinations : %d\n", nrow(x$alldoses)))
  cat(sprintf("Features (D)         : %d\n", ncol(x$alldoses)))
  cat(sprintf("Total observations   : %d\n", sum(x$M)))
  cat(sprintf("Classified as 1      : %d (%.1f%%)\n",
              sum(x$yhat), 100*mean(x$yhat)))
  cat(sprintf("Log posterior gain   : %.4f\n\n", x$logH))
  db = boundary(x)
  cat(sprintf("Decision boundary    : %d combination(s)\n", db$n_boundary))
  if (!is.null(db$boundary)) {
    out = as.data.frame(db$boundary)
    out$pt = round(db$pt, 3)
    if (nrow(out) > 20) {
      cat(sprintf("(Showing first 20 of %d boundary combinations)\n", nrow(out)))
      print(utils::head(out, 20), row.names=FALSE)
      cat(sprintf("... %d more rows\n", nrow(out)-20))
    } else {
      print(out, row.names=FALSE)
    }
  }
  invisible(x)
}

#' @export
print.boundary = function(x, ...) {
  cat("Decision Boundary of Projective Bayes Classifier\n")
  cat("=================================================\n")
  cat(sprintf("Total feature combinations : %d\n", x$n_total))
  cat(sprintf("Classified as positive (1) : %d (%.1f%%)\n",
              x$n_positive, 100*x$n_positive/x$n_total))
  cat(sprintf("Boundary combinations      : %d\n\n", x$n_boundary))
  if (is.null(x$boundary)) {
    if (x$n_positive == 0) {
      cat("No combinations classified as 1.\n")
    } else {
      cat("All combinations classified as 1 - no boundary.\n")
    }
  } else {
    cat("Boundary combinations (minimal positive set):\n")
    out = as.data.frame(x$boundary)
    out$pt = round(x$pt, 3)
    if (nrow(out) > 20) {
      cat(sprintf("(Showing first 20 of %d boundary combinations)\n", nrow(out)))
      print(utils::head(out, 20), row.names=FALSE)
      cat(sprintf("... %d more rows\n", nrow(out)-20))
    } else {
      print(out, row.names=FALSE)
    }
  }
  invisible(x)
}

#' @export
print.miso = function(x, ...) {
  cat("Multivariable Isotonic Regression (Binary)\n")
  cat("==========================================\n")
  cat(sprintf("Feature combinations : %d\n", nrow(x$alldoses)))
  cat(sprintf("Features (D)         : %d\n", ncol(x$alldoses)))
  cat(sprintf("Total observations   : %d\n", sum(x$M)))
  cat(sprintf("Threshold grid points: %d\n", x$nt))
  cat(sprintf("thetahat range       : [%.4f, %.4f]\n",
              min(x$thetahat, na.rm=TRUE),
              max(x$thetahat, na.rm=TRUE)))
  cat("\nFeature combinations and estimated response probabilities:\n")
  out = as.data.frame(x$alldoses)
  out$M        = x$M
  out$S        = x$S
  out$thetahat = round(x$thetahat, 4)
  if (nrow(out) > 20) {
    cat(sprintf("(Showing first 20 of %d combinations)\n", nrow(out)))
    print(utils::head(out, 20), row.names=FALSE)
    cat(sprintf("... %d more rows\n", nrow(out)-20))
  } else {
    print(out, row.names=FALSE)
  }
  invisible(x)
}

#' @export
print.misoN = function(x, ...) {
  cat("Multivariable Isotonic Regression (Continuous)\n")
  cat("===============================================\n")
  cat(sprintf("Feature combinations : %d\n", nrow(x$alldoses)))
  cat(sprintf("Features (D)         : %d\n", ncol(x$alldoses)))
  cat(sprintf("Total observations   : %d\n", sum(x$M)))
  cat(sprintf("Threshold grid points: %d\n", x$nt))
  cat(sprintf("thetahat range       : [%.4f, %.4f]\n",
              min(x$thetahat, na.rm=TRUE),
              max(x$thetahat, na.rm=TRUE)))
  cat("\nFeature combinations and estimated mean responses:\n")
  out = as.data.frame(x$alldoses)
  out$M        = x$M
  out$MY       = round(x$MY, 4)
  out$thetahat = round(x$thetahat, 4)
  if (nrow(out) > 20) {
    cat(sprintf("(Showing first 20 of %d combinations)\n", nrow(out)))
    print(utils::head(out, 20), row.names=FALSE)
    cat(sprintf("... %d more rows\n", nrow(out)-20))
  } else {
    print(out, row.names=FALSE)
  }
  invisible(x)
}


