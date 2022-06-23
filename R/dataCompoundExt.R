#' Compound Extraction experiment of Dopico-Garc{\' i}a et al. (2007)
#'
#' An analytical experiment conducted by Dopico-Garc{\' i}a et al. (2007) to
#' characterize the chemical composition of white grapes simultaneously
#' determining the most important phenolic compounds and organic acids
#' for the grapes. This example has been further studied in Phoa et al. (2009b)
#' for one phenolic compound, kaempferol-3-Orutinoside + isorhamnetin-3-O glucoside, 
#' which is also what we studied. It is accepted for these data 
#' that fitting a main-effects model suggests that V3 (Factor C), V4 (Factor D), 
#' and inteaction of V1:V4 (A:D) are active effects. 
#' 
#' 
#' @docType data
#'
#' @usage data(dataCompoundExt)
#'
#' @format A data frame with 12 rows and 9 columns:
#' \describe{
#'   \item{V1}{Factor A}
#'   \item{V2}{Factor B}
#'   \item{V3}{Factor C}
#'   \item{V4}{Factor D}
#'   \item{V5}{Factor E}
#'   \item{V6}{Factor F}
#'   \item{V7}{Factor G}   
#'   \item{V8}{Factor H}
#'   \item{V9}{Response}
#' }
#'
#' @keywords datasets
#'
#' @source Dopico-Garc{\' i}a, M.S., Valentao, P., Guerra, L., Andrade, P. B., and Seabra, R. M. (2007).
#' Experimental design for extraction and quantification of phenolic 
#' compounds and organic acids in white "Vinho Verde" grapes 
#' Analytica chimica acta, 583(1): 15--22.
#' 
#' @source Phoa, F. K., Wong, W. K., and Xu, H (2009b). The need of 
#' considering the interactions in the analysis of screening designs. 
#' Journal of Chemometrics: A Journal of the Chemometrics Society, 23(10): 545--553.
#'
#' @examples
#' data(dataCompoundExt)
#' X = dataCompoundExt[,-9]
#' Y= dataCompoundExt[,9]
#'
"dataCompoundExt"
