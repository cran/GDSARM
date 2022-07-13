#' Cast fatigue experiment of Hunter et al. (1982)
#'
#' A cast fatigue experiment with 12 runs and 7 factors was originally studied
#' by Hunter et al. (1982), and was later revisited by Hamada and Wu (1992) 
#' and Phoa et al. (2009), among others. It is widely accepted for these data that 
#' V6 (F) and interaction V6:V7 (F:G) are active effects, with interaction of 
#' V1:V5 (A:E) possibly being active as well.
#' 
#' 
#' @docType data
#'
#' @usage data(dataHamadaWu)
#'
#' @format A data frame with 12 rows and 8 columns:
#' \describe{
#'   \item{V1}{Factor A}
#'   \item{V2}{Factor B}
#'   \item{V3}{Factor C}
#'   \item{V4}{Factor D}
#'   \item{V5}{Factor E}
#'   \item{V6}{Factor F}
#'   \item{V7}{Factor G}
#'   \item{V8}{Response}
#' }
#'
#' @keywords datasets
#'
#' @source Hamada, M. and C. F. J. Wu (1992). Analysis of designed experiments 
#' with complex aliasing. Journal of Quality Technology 24 (3), 130--137.
#' 
#' @source Hunter, G. B., F. S. Hodi, and T. W. Eagar (1982). High cycle fatigue of weld repaired
#' cast Ti-6AI-4V. Metallurgical Transactions A 13 (9), 1589--1594.
#' 
#' @source Phoa, F. K., Y. H. Pan, and H. Xu (2009). Analysis of supersaturated 
#' designs via the Dantzig selector. Journal of Statistical Planning and Inference
#' 139 (7), 2362--2372.
#'
#' @examples
#' data(dataHamadaWu)
#' X = dataHamadaWu[,-8]
#' Y= dataHamadaWu[,8]
#'
"dataHamadaWu"
