###
#' @importFrom Rcpp evalCpp
#' @useDynLib combAUCbc, .registration = TRUE
#'

#' @export
bcAUC <- function(method = c("fi", "msi", "ipw", "spe", "knn"), T, D, V,
                  rhoEst = NULL, piEst = NULL){
  method <- match.arg(method)
  method_name <- toupper(method)
  D.flag <- any(is.na(D))
  if(method == "fi"){
    if(is.null(rhoEst)) stop("The input of argument \"rhoEst\" is needed for ", method_name, " estimator")
    if(length(T) != length(rhoEst$values)) stop(gettextf("arguments imply differing number of observation: %d", length(T), ", %d", length(rhoEst$values)), domain = NA)
    ans <- bc_AUC(T, rhoEst$values)
  }
  else if(method %in% c("msi", "knn")){
    if(is.null(rhoEst)) stop("argument \"rhoEst\" is needed for", method_name, "estimator")
    if(length(T) != length(rhoEst$values)) stop(gettextf("arguments imply differing number of observation: %d", length(T), ", %d", length(rhoEst$values)), domain = NA)
    Dtemp <- D
    if(D.flag) Dtemp[is.na(D)] <- 99
    Dmsi <- Dtemp*V + (1 - V)*rhoEst$values
    ans <- bc_AUC(T, Dmsi)
  }
  else if(method == "ipw"){
    if(is.null(piEst)) stop("argument \"piEst\" is needed for", method_name, "estimator")
    if(length(T) != length(piEst$values)) stop(gettextf("arguments imply differing number of observation: %d", length(T), ", %d", length(piEst$values)), domain = NA)
    Dtemp <- D
    if(D.flag) Dtemp[is.na(D)] <- 99
    Dipw <- Dtemp*V/piEst$values
    ans <- bc_AUC(T, Dipw)
  }
  else{
    if(is.null(rhoEst) | is.null(piEst)) stop("arguments \"rhoEst\" and \"pi.est\" are needed for", method_name, "estimator")
    if(length(T) != length(rhoEst$values)) stop(gettextf("arguments imply differing number of observation: %d", length(T), ", %d", length(rhoEst$values)), domain = NA)
    if(length(T) != length(piEst$values)) stop(gettextf("arguments imply differing number of observation: %d", length(T), ", %d", length(piEst$values)), domain = NA)
    Dtemp <- D
    if(D.flag) Dtemp[is.na(D)] <- 99
    Dspe <- Dtemp*V/piEst$values - (V/piEst$values - 1)*rhoEst$values
    ans <- bc_AUC(T, Dspe)
  }
  return(ans)
}
