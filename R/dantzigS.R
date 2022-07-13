#' Dantzig selector using the lpsolve package
#'
#' @description The Dantzig selector (DS) finds a solution for the model parameters
#' of a linear model, \code{beta} using linear programming. For a given \code{delta},
#' DS minimizes the L_1-norm (sum of absolute values) of \code{beta} subject to the constraint
#' that \code{max(|t(X)(y-X * beta)|)}<= \code{delta}.
#'
#' @param X a design matrix.
#'
#' @param y a vector of  responses.
#'
#' @param delta the specific value of \code{delta} for which the Dantzig Selector
#' optimization needs to be solved
#'
#' @param scale.X a number by which each column of \code{X} should be scaled
#'
#' @return A list containing the (a) \code{opt} (Value of objective function
#' at optimum), (b) \code{status} (Numeric indicator: 0 = success,
#' 2 = no feasible solution), (c) \code{beta} (the estimated values
#' of \code{beta}), (d) \code{delta}
#'
#' @source Cand{\`e}s, E. and Tao, T. (2007). The Dantzig selector: Statistical estimation when p is much
#' larger than n. Annals of Statistics 35 (6), 2313--2351.
#' 
#' @source Phoa, F. K., Pan, Y. H. and Xu, H. (2009). Analysis of supersaturated 
#' designs via the Dantzig selector. Journal of Statistical Planning and Inference
#' 139 (7), 2362--2372.
#' 
#' 
#'
#' @keywords internal

## Dantzig selector solved by linear programming
dantzigS=function(X, y, delta, scale.X=1)
{ # model: y=X beta + epsilon, using lpSolve(), # 2/6/07
  X =as.matrix(X)
  p=ncol(X)
  n=nrow(X)
  Ip=diag(1, p); zerop=diag(0,p)
  X = X / scale.X

  U=t(X) %*% X; V=t(X) %*% y
  Amat=rbind(cbind(U, -U), cbind(-U, U), cbind(2*Ip, -Ip))
  bvec=c(-V-delta, V-delta, rep(0, p))
  cvec = rep(c(1,0), c(p, p))

  lps = lpSolve::lp("min", cvec,  Amat, ">=", bvec)
  if(lps$status != 0)	print(lps)
  beta=lps$sol[(p+1):(2*p)] -lps$sol[1:p]
  names(beta)=dimnames(X)[[2]]
  list(opt=lps$objval, status=lps$status, beta=beta, delta=delta)
}
