#' Dantzig selector with an option to make profile plot
#'
#' @description The Dantzig selector (DS) finds a solution for the model parameters
#' of a linear model, \code{beta} using linear programming. For a given \code{delta},
#' DS minimizes the L_1-norm (sum of absolute values) of \code{beta} subject to the constraint
#' that \code{max(|t(X)(y-X * beta)|) <= delta}.
#'
#' @param X a design matrix.
#'
#' @param y a vector of  responses.
#'
#' @param delta a vector with the values of \code{delta} for which the DS
#' optimization needs to be solved.
#'
#' @param plot a boolean value of either TRUE or FALSE with TRUE
#' indicating that the profile plot should be drawn.
#'
#' @return A matrix of the estimated values of \code{beta} with each 
#' row corresponding to a particular value of \code{delta}.
#'
#' @seealso \code{\link{GDS_givencols}}, \code{\link{GDSARM}}
#'
#' @source Cand{\`e}s, E. and Tao, T. (2007). The Dantzig selector: Statistical estimation when p is much
#' larger than n. Annals of Statistics 35 (6), 2313--2351.
#' 
#' @source Phoa, F. K., Pan, Y. H. and Xu, H. (2009). Analysis of supersaturated 
#' designs via the Dantzig selector. Journal of Statistical Planning and Inference
#' 139 (7), 2362--2372.
#' 
#' @examples
#' data(dataHamadaWu)
#' X = dataHamadaWu[,-8]
#' Y = dataHamadaWu[,8]
#' #scale and center X and y
#' scaleX = base::scale(X, center= TRUE, scale = TRUE)
#' scaleY = base::scale(Y, center= TRUE, scale = FALSE)
#' maxDelta = max(abs(t(scaleX)%*%matrix(scaleY, ncol=1)))
#' # Dantzig Selector on 4 equally spaced delta values between 0 and maxDelta
#' dantzig.delta(scaleX, scaleY, delta = seq(0,maxDelta,length.out=4)) 
#'
#' @export

dantzig.delta = function(X, y, delta, plot=FALSE)
{
  beta = NULL
  for(i in 1:length(delta)) beta = rbind(beta, dantzigS(X, y, delta[i])$beta)
  dimnames(beta)[[1]]=delta
  if(plot) graphics::matplot(delta, beta, type="l")  # profile plot
  beta
}
