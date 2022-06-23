#' Step I: Multiple GDS runs with random interactions
#'
#' @description Runs Gauss Dantzig Selector (GDS) multiple times, each time 
#' with a differently set of randomly selected two-factor interactions.
#' All \code{m} main effects are included in each GDS run. For each set of
#' randomly selected interactions, the best GDS output is chosen among 
#' \code{delta.n} values of \code{delta}. We use kmeans with 2 
#' clusters and BIC to select such best model.
#'
#'
#' @param delta.n a positive integer suggesting the number of delta values
#' to be tried. \code{delta.n} values of \code{delta} will be used 
#' strictly between 0 and \code{max(t(X)y)}. The default value is set to 10.
#' 
#' @param nint a positive integer representing the number of randomly 
#' chosen interactions. The suggested value to use is the ceiling of 20% 
#' of the total number of interactions, that is, for \code{m} factors, we have
#' \code{ceiling(0.2(m choose 2))}.
#' 
#' @param nrep a positive integer representing the number of times GDS should
#' be run. The suggested value is \code{(m choose 2)}.
#' 
#' @param Xmain a \eqn{n \times m}{n x m} matrix of \code{m} main effects.
#' 
#' @param Xint a matrix of \eqn{{m \choose 2}}{\code{m choose 2})} two-factor 
#' interactions.
#' 
#' @param Y a vector of \code{n} responses.
#' 
#' @param opt.heredity a string with either `none`, or `weak`, or `strong`. Denotes
#' whether the effect-heredity (weak or strong) should be embedded in GDS-ARM. 
#' The default value is `none` as suggested in Singh and Stufken (2022).
#' 
#' @return A list containing the (a) matrix of the output of each GDS run with each 
#' row representing the selected effects from the corresponding GDS run, (b) 
#' a vector with the corresponding BIC values of each model.
#'
#' @source Singh, R. and Stufken, J. (2022). Factor selection in screening experiments
#' by aggregation over random models, 1--31. \doi{10.48550/arXiv.2205.13497}
#' 
#' @keywords internal
#'
## needs delta.n, nint, nrep,design as Xmain (main effects), Xint (two factor interactions)
StepI_chooseints <- function(delta.n=10,nint,nrep,Xmain,Xint,Y, opt.heredity = c("none")) {
  n = dim(Xmain)[1]
  m= dim(Xmain)[2]
  
  ncluster = 2
  
  ncol = dim(Xmain)[2]+dim(Xint)[2]
  ListofSelectedEffects=c()
  all.nint=matrix(ncol = nint,nrow = nrep)
  for (i in 1:nrep) {
    all.nint[i,]=sample((dim(Xmain)[2]+1):ncol,nint,replace = F)
  }
  B=matrix(ncol = ncol,nrow = nrep)
  BICReps= numeric(nrep)
  res = matrix(0,ncol = ncol,nrow = nrep)
  for (irep in 1:nrep) {
    Xran=cbind(Xmain,Xint[,all.nint[irep,]-dim(Xmain)[2]])

    delta0 <- max(abs(t(Xran)%*%Y))
    delta = seq(0, delta0, length.out = delta.n+2)
    delta = delta[-1]

    beta = dantzig.delta(Xran, Y, delta, plot=FALSE)
    beta2 = beta;
    beta2[]=0

    #matrix(0,nrow=nrow(beta),ncol=ncol(beta))
    BICafterkmeans = numeric(delta.n-1)
    BICafterkmeans = numeric(delta.n-1)
    for (ct_lam in 1:(delta.n-1)){
      ## doing kmeans
      tempvec = beta[ct_lam,]
      tempvec1 = c(tempvec,rep(0,ncol-length(tempvec)))
      if (tempvec1[order(tempvec1, decreasing = T)][ncluster-1]> 0) {
        kmeansOUT=stats::kmeans(abs(tempvec1), ncluster) #$cluster
        clusternumber = order(kmeansOUT$centers)[ncluster] #pick the cluster
        #with largest cluster mean
        selectcols = tempvec[which(kmeansOUT$cluster==clusternumber)]
        # print(selectcols)
        beta2[ct_lam,names(selectcols)]=1
        beta2[ct_lam,names(selectcols)]=1
		#heredity or not
    namesselectcols = names(selectcols)
   # namesselectcols = (selectcols)
    namesselectints = grep(":", namesselectcols, value = F)
    deleteInt = numeric(length(namesselectints))
    deleteInt[] = 1
    namesselectmain = namesselectcols[setdiff(1:length(namesselectcols),namesselectints)]
   if ((length(namesselectints)) >= 1) {
   for (tt in 1:length(namesselectints)){
    selectint = namesselectcols[namesselectints[tt]]
    selectint1 = unlist(strsplit(selectint,":"))
    if (opt.heredity == "weak" ) {
      if (selectint1[1] %in% namesselectmain | selectint1[2] %in% namesselectmain){
      deleteInt[tt] = 0
      }
    } else if (opt.heredity == "strong" ) {
      if (selectint1[1] %in% namesselectmain & selectint1[2] %in% namesselectmain){
        deleteInt[tt] = 0
      }
    } else {
	 deleteInt[tt] = 0
	}
  }
  namesToUse = c(namesselectmain,namesselectcols[namesselectints[which(deleteInt == 0)]])
  } else {
    namesToUse = namesselectmain  
  }
        datasubset = as.data.frame(cbind(Xran[,namesToUse,drop=FALSE],Y))
        colnames(datasubset) = c(namesToUse,"Y")
		## fit a model to find BIC
        fit = stats::lm(Y~., data=datasubset)
        BICafterkmeans[ct_lam] =stats::AIC(fit, k = log(n))#BIC
      } else {
        BICafterkmeans[ct_lam] = 10^8
      }
    }
    Index_minBIC = which.min(BICafterkmeans)
    SelectedEffectsFromThisRep = colnames(beta2[Index_minBIC,which(beta2[Index_minBIC,]==1),drop=FALSE])
    if (length(SelectedEffectsFromThisRep)> 0) {
      B[irep,1:length(SelectedEffectsFromThisRep)]=SelectedEffectsFromThisRep
    }
    BICReps[irep]=BICafterkmeans[Index_minBIC]
  }

  return(list(output=B,BIC=BICReps))
}

