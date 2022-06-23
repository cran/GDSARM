#' Gauss-Dantzig Selector
#'
#' @description This function runs Gauss-Dantzig Selector on the given columns. 
#' We have two options: either (a) GDS(m) on the \code{m} main 
#' effects, and (b) GDS(m+2fi) on the \code{m} main effects plus two-factor interactions. 
#' For a given \code{delta}, DS minimizes the l1-norm (sum of absolute values)
#' of \code{beta} subject to the constraint that \code{max(t(X)(y-X * beta))} <= \code{delta}.
#' GDS is run for multiple values of \code{delta}. We use kmeans with 2 
#' clusters and BIC to select a best model  for all the \code{delta} values.
#' 
#' @param delta.n a positive integer suggesting the number of deltas
#' to be tried. \code{delta.n} values of \code{delta} will be used 
#' strictly between 0 and \code{max(t(X)y)}. The default value is
#' set to 10.
#' 
#' @param design a \eqn{n \times m}{n x m} matrix of \code{m} two-level main effects. 
#' The columns should have the high-level coded as +1 and the low-level codes as -1.
#' 
#' @param Y a vector of \code{n} responses.
#' 
#' @param which.cols a string with either `main` or `main2fi`. Denotes
#' whether the Gauss-Dantzig Selector should be run on the main effect columns (`main`), or on  
#' all main effects plus all 2 factor interaction columns (`main2fi`). 
#' The default value is `main2fi`.
#' 
#' @return A list returning the effects identified as active as well as the
#' corresponding identified important factors.
#' 
#' @source Cand{\`e}s, E. and Tao, T. (2007). The Dantzig selector: Statistical estimation when p is much
#' larger than n. Annals of Statistics 35 (6), 2313--2351.
#' 
#' @source Dopico-Garc{\' i}a, M.S., Valentao, P., Guerra, L., Andrade, P. B., and Seabra, R. M. (2007).
#' Experimental design for extraction and quantification of phenolic 
#' compounds and organic acids in white "Vinho Verde" grapes 
#' Analytica chimica acta, 583(1): 15--22.
#' 
#' @source Hamada, M. and Wu, C. F. J. (1992). Analysis of designed experiments 
#' with complex aliasing. Journal of Quality Technology 24 (3), 130--137.
#' 
#' @source Hunter, G. B., Hodi, F. S. and Eagar, T. W. (1982). High cycle fatigue of weld repaired
#' cast Ti-6AI-4V. Metallurgical Transactions A 13 (9), 1589--1594.
#' 
#' @source Phoa, F. K., Pan, Y. H. and Xu, H. (2009). Analysis of supersaturated 
#' designs via the Dantzig selector. Journal of Statistical Planning and Inference
#' 139 (7), 2362--2372.
#'
#' @source Singh, R. and Stufken, J. (2022). Factor selection in screening experiments
#' by aggregation over random models, 1--31. \doi{10.48550/arXiv.2205.13497}
#' 
#' @seealso \code{\link{GDSARM}}, \code{\link{dantzig.delta}}
#'
#' @examples
#' data(dataHamadaWu)
#' X = dataHamadaWu[,-8]
#' Y= dataHamadaWu[,8]
#' delta.n=10
#' # GDS on main effects 
#' GDS_givencols(delta.n,design = X, Y=Y, which.cols = "main")
#' 
#' # GDS on main effects and two-factor interactions
#' GDS_givencols(delta.n,design = X, Y=Y)
#' 
#' data(dataCompoundExt)
#' X = dataCompoundExt[,-9]
#' Y= dataCompoundExt[,9]
#' delta.n=10
#' # GDS on main effects
#' GDS_givencols(delta.n,design = X, Y=Y, which.cols = "main")
#' # GDS on main effects and two-factor interactions
#' GDS_givencols(delta.n,design = X, Y=Y, which.cols = "main2fi")
#' @export
#'

GDS_givencols<- function(delta.n=10, design, Y, which.cols = c("main2fi")){
  n = dim(design)[1]
  m = dim(design)[2]
  ncluster = 2
  artificial.response <- stats::rnorm(n,0,1)
  
  Xall <- as.data.frame(stats::model.matrix(stats::lm(artificial.response ~ .^2, design)))[,-1]
  Xall.scale = as.data.frame(base::scale(Xall,center = T,scale = T)) #*(1/sqrt(n-1))*(0.2)^12
  
  if (which.cols == "main") {
    Xran = Xall.scale[,1:m]
  } else {
    Xran = Xall.scale
  }

  ncol = ncol(Xran)
  delta0 <- max(abs(t(Xran)%*%Y))
  delta = seq(0, delta0, length.out = delta.n+2)
  delta = delta[-1]
  
  beta = dantzig.delta(Xran, Y, delta, plot=FALSE)
  beta2 = beta;
  beta2[]=0
  
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
      datasubset = as.data.frame(cbind(Xran[,names(selectcols),drop=FALSE],Y))
      colnames(datasubset) = c(names(selectcols),"Y")
      ## fit a model to find BIC
      fit = stats::lm(Y~., data=datasubset)
      BICafterkmeans[ct_lam] =stats::AIC(fit, k = log(n))#BIC
    } else {
      BICafterkmeans[ct_lam] = 10^8
    }
  }
  Index_minBIC = which.min(BICafterkmeans)
  SelectedEffects = colnames(beta2[Index_minBIC,which(beta2[Index_minBIC,]==1),drop=FALSE])
  
  SelectedFactors1 = strsplit(SelectedEffects,":")
  SelectedFactors = unique(c(unlist(SelectedFactors1)))
  
  return(list(effects = SelectedEffects, factors =SelectedFactors))
}
