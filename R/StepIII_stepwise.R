#' Step III: Stepwise on the consolidated output from different GDS runs
#'
#' @description Runs the stepwise regression on the output received from 
#' top models of the consolidated output of different GDS runs. With
#' \code{n} being the number of runs, the stepwise regression starts
#' with at most \code{(n-3)} selected effects from the previous step. The 
#' remaining effects from the previous step as well as all main effects are
#' given a chance to enter into the model using the forward-backward stepwise
#' regression.
#'
#' @param xstart a vector with effects' names corresponding to the starting model.
#'  
#' @param xremain a vector with effects' names corresponding to the remaining 
#' main effects and other effects that needs to be explored with stepwise regression.
#'  
#' 
#' @param Xmain a \eqn{n \times m}{n x m} matrix of \code{m} main effects.
#' 
#' @param Xint a matrix of \code{m choose 2} two-factor 
#' interactions.
#' 
#' @param Y a vector of \code{n} responses.
#' 
#' @param cri.penter the p-value cutoff for an effect to enter into the 
#' stepwise regression model
#' 
#' @param cri.premove the p-value cutoff for an effect to exit from the 
#' stepwise regression model
#'
#' @param opt.heredity a string with either `none`, or `weak`, or `strong`. Denotes
#' whether the effect-heredity (weak or strong) should be embedded in GDS-ARM. 
#' The default value is `none` as suggested in Singh and Stufken (2022).
#' 
#' @return A list returning the effects identified as active as well as the
#' corresponding identified important factors.
#'
#' @source Singh, R. and Stufken, J. (2022). Factor selection in screening experiments
#' by aggregation over random models, 1--31. \doi{10.48550/arXiv.2205.13497}
#' 
#' @keywords internal
#'
StepIII_stepwise<-function(xstart,xremain,Xmain, Xint, Y,cri.penter=0.01, cri.premove=0.05,opt.heredity="none"){

  nrow = dim(Xmain)[1]
  Xall = cbind(Xmain, Xint)
  
  #dataset needed for stepwise only
  data=as.data.frame(cbind(Xall,Y))
  colnames(data)[dim(data)[2]]="y"
  colnames(data)=gsub(":", "a",colnames(data))
  
  ## 2nd selection -- backward
  modelfin=0
  deadloop=0

  while(length(xstart)+2<nrow & modelfin==0){

    deadloop=deadloop+1
    if(deadloop>30){break}
    movefor=1
    moveback=1

    #forward
    model = paste(xstart,sep = "",collapse ="+")
    if(model==""){model="y~1"}else{
      model = paste("y~1+",model,sep = "")}


    if (length(xremain)>= 1) {
      criteria=vector("numeric",length(xremain))
      for (i in 1:length(xremain)) {
        ##print("hi")
        model2=paste(model,"+",xremain[i],sep = "")
        temp.model2 <- stats::lm(model2, data = data)
        pv=summary(temp.model2)$coefficients[-1,4]
        criteria[i]=pv[length(pv)]
      }
      deleposi=which(criteria==min(criteria)) #r2--max
      if(min(criteria)>cri.penter){movefor=0}else{
        xstart=c(xstart,xremain[deleposi])
        xremain=xremain[-deleposi]
      }
    } else{movefor=0}
    #backward
    model = paste(xstart,sep = "",collapse ="+")
    if(model==""){
      model="y~1"
      max.pvalue=0
    }else{
      model = paste("y~1+",model,sep = "")
      temp.model <- stats::lm(model, data = data)
      pvalue <- base::summary(temp.model)$coefficients[-1,4]
      max.pvalue=max(pvalue)}
    if(max.pvalue<=1 & max.pvalue>cri.premove){
      deleposi = which(xstart==names(pvalue)[which(pvalue==max.pvalue)])
      xremain=c(xremain,xstart[deleposi])
      xstart=xstart[-deleposi]
    }else{moveback=0}

    if(length(xstart)==0){modelfin=1}
    if(length(xremain)==0){modelfin=1}
    if(length(xstart)>=nrow){modelfin=1}
    if(moveback==0 & movefor==0){modelfin=1}
  }

  if(length(xstart)==0){
    result=character(0)
  }else{
    result=matrix(nrow=length(xstart),ncol=2)
    for(i in 1:length(xstart)){
      result[i,]= unlist(strsplit(xstart[i],split="a"))
      if(result[i,1]==result[i,2]){result[i,2]=""}
    }
    result = as.matrix(result)}

  xstarta=gsub("a", ":",xstart)
  xstarta = unique(xstarta)
  #embed heredity or not
  namesselectcols = xstarta
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
  xstarta = unique(namesToUse)
  result1 = unique(unlist(strsplit(xstarta,":")))
  return(list(effects=xstarta, factors=result1) )
}
