################################################################################
# # NCDIF.R
# #
# # Author(s): Victor H. Cervantes
# #
# # Extracting item parameters and covariance matrices from mirt objects
# # Description: Functions to import estimates from mirt fits
# #
# # Inputs: NULL
# #
# # Outputs: functions
# #
# # File history:
# #   20210623: Creation
################################################################################



#' Extract item difficulties and item difficulty variance estimates for Rasch items from a fitted mirt object
#' for one or two groups
#'
#' @param mod       A mirt object containing the fit of unidimensional model.
#' @param focal     Character. Required if mod is MultipleGroupClass, focal should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#' @param reference Character. Required if mod is MultipleGroupClass, reference should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#'
#' @return If mod contains any itemtype == "Rasch", a list with the item parameters and the estimate covariances (if available).
#' If mod is SingleGroupClass, the list contains the item parameters as a matrix and the covariances as a list.
#' If mod is MultipleGroupClass, the list contains the item parameters and covariances for the focal and reference groups only.
#'
#' @examples
#' library(mirt)
#' data <- expand.table(LSAT7)
#' (mod1 <- mirt(data, model = 1, itemtype = "Rasch", SE = TRUE))
#' (DFIT:::ExtractRaschMirt(mod1))
#'
ExtractRaschMirt <- function (mod, focal = NULL, reference = NULL) {

    classMod <- attributes(class(mod))[["package"]]

    if (is.null(classMod) || classMod != "mirt") {
        stop("mod is not a mirt obect")
    }
    if (mod@Model[["model"]] != 1) {
        stop("mod is not a unidimensional model.")
    }

    whichItems <- which(mod@Model[["itemtype"]] == "Rasch")

    if (length(whichItems) == 0) {
        message("None of the items is itemtype Rasch.")
        return(NULL)
    }

    if (class(mod) == "SingleGroupClass") {
        itemPars <- matrix(
            #mirt::coef(mod,
            coef(mod,
                 IRTpars = TRUE,
                 simplify = TRUE)[["items"]][whichItems, "b"],
            ncol = 1)

    } else if (class(mod) == "MultipleGroupClass") {
        if (is.null(focal) || is.null(reference)) {
            stop("focal and reference group names are required.")
        }
        if (!(focal %in% mod@Data$groupNames)) {
            stop("focal does not match any group name.")
        }
        if (!(reference %in% mod@Data$groupNames)) {
            stop("reference does not match any group name.")
        }

        #itemPars <- lapply(mirt::coef(mod,
        itemPars <- lapply(coef(mod,
                                IRTpars = TRUE,
                                simplify = TRUE),
                           function (x) {
                               matrix(x[["items"]][whichItems, "b"],
                                      ncol = 1)
                           }
        )

        itemPars <- itemPars[c(focal, reference)]
        names(itemPars) <- c("focal", "reference")
    }

    if (mod@Options[["SE"]]) {
        if (class(mod) == "SingleGroupClass") {
            #itemCov <- as.list(diag(mirt::vcov(mod))[whichItems])

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)

            whichCov <- ordCoefs[whichItems, "d"]

            whichCov <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichCov
                )

            itemCov <- vcov(mod)[whichCov, whichCov]

        } else if (class(mod) == "MultipleGroupClass") {
            itemCov <- list()

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)
            nVcov     <- nrow(vcovMod)

            nGroups <- length(mod@Data$groupNames)
            nPars   <- nVcov / nGroups
            nParsGroup <- as.numeric(gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)[nPars + 1]) - 1

            whichFocal     <- which(mod@Data$groupNames == focal)
            whichReference <- which(mod@Data$groupNames == reference)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[[1]][["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)

            whichCov <- ordCoefs[whichItems, "d"]

            pickFocal    <- whichCov + (nParsGroup * (whichFocal - 1))
            pickReferece <- whichCov + (nParsGroup * (whichReference - 1))

            pickFocal <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickFocal
                )

            pickReferece <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickReferece
                )

            itemCov[["focal"]]     <- vcovMod[pickFocal, pickFocal]
            itemCov[["reference"]] <- vcovMod[pickReferece, pickReferece]

        }
    } else {
        itemCov <- NULL
    }

    output <- list(itemParameters = itemPars, itemCovariances = itemCov)

    return(output)

}



#' Extract item discrimination and difficulties and estimate covariance estimates for 2PL items from a fitted mirt object
#' for one or two groups
#'
#' @param mod       A mirt object containing the fit of unidimensional model.
#' @param focal     Character. Required if mod is MultipleGroupClass, focal should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#' @param reference Character. Required if mod is MultipleGroupClass, reference should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#'
#' @return If mod contains any itemtype == "2PL", a list with the item parameters and the estimate covariances (if available).
#' If mod is SingleGroupClass, the list contains the item parameters as a matrix and the covariances as a list.
#' If mod is MultipleGroupClass, the list contains the item parameters and covariances for the focal and reference groups only.
#'
#' @examples
#' library(mirt)
#' data <- expand.table(LSAT7)
#' (mod1 <- mirt(data, model = 1, itemtype = "2PL", SE = TRUE))
#' (DFIT:::Extract2PLMirt(mod1))
#'
Extract2PLMirt <- function (mod, focal = NULL, reference = NULL) {

    classMod <- attributes(class(mod))[["package"]]

    if (is.null(classMod) || classMod != "mirt") {
        stop("mod is not a mirt obect")
    }
    if (mod@Model[["model"]] != 1) {
        stop("mod is not a unidimensional model.")
    }

    whichItems <- which(mod@Model[["itemtype"]] == "2PL")

    if (length(whichItems) == 0) {
        message("None of the items is itemtype 2PL.")
        return(NULL)
    }

    if (class(mod) == "SingleGroupClass") {
        itemPars <- as.matrix(
            #mirt::coef(mod,
            coef(mod,
                 IRTpars = TRUE,
                 simplify = TRUE)[["items"]][whichItems, c("a", "b")])

    } else if (class(mod) == "MultipleGroupClass") {
        if (is.null(focal) || is.null(reference)) {
            stop("focal and reference group names are required.")
        }
        if (!(focal %in% mod@Data$groupNames)) {
            stop("focal does not match any group name.")
        }
        if (!(reference %in% mod@Data$groupNames)) {
            stop("reference does not match any group name.")
        }

        #itemPars <- lapply(mirt::coef(mod,
        itemPars <- lapply(coef(mod,
                                IRTpars = TRUE,
                                simplify = TRUE),
                           function (x) {
                               as.matrix(x[["items"]][whichItems, c("a", "b")])
                           }
        )

        itemPars <- itemPars[c(focal, reference)]
        names(itemPars) <- c("focal", "reference")
    }

    if (mod@Options[["SE"]]) {
        if (class(mod) == "SingleGroupClass") {
            #itemCov <- as.list(diag(mirt::vcov(mod))[whichItems])

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)

            whichCov <- ordCoefs[whichItems, c("a1", "d")]
            whichDis <- ordCoefs[whichItems, "a1"]
            whichDif <- ordCoefs[whichItems, "d"]

            whichCov <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichCov
                )

            whichDis <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDis
                )

            whichDif <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDif
                )

            itemCov <- vcov(mod)[whichCov, whichCov]

            transDisc <- paste0("~ x", which(whichCov %in% whichDis))
            transDiff <- paste0("~ -x", which(whichCov %in% whichDif), "/x", whichDif - 1)

            trans2PL <- list()
            trans2PL[which(whichCov %in% whichDis)] <- transDisc
            trans2PL[which(whichCov %in% whichDif)] <- transDiff

#            itemCov <- msm::deltamethod(g = trans2PL,
            itemCov <- deltamethod(g = lapply(trans2PL, as.formula),
                                        mean = as.numeric(t(itemPars)),
                                        cov = itemCov,
                                        ses = FALSE)

        } else if (class(mod) == "MultipleGroupClass") {
            itemCov <- list()

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)
            nVcov     <- nrow(vcovMod)

            nGroups <- length(mod@Data$groupNames)
            nPars   <- nVcov / nGroups
            nParsGroup <- as.numeric(gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)[nPars + 1]) - 1

            whichFocal     <- which(mod@Data$groupNames == focal)
            whichReference <- which(mod@Data$groupNames == reference)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[[1]][["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)

            whichCov <- ordCoefs[whichItems, c("a1", "d")]
            whichDis <- ordCoefs[whichItems, "a1"]
            whichDif <- ordCoefs[whichItems, "d"]

            pickFocal    <- whichCov + (nParsGroup * (whichFocal - 1))
            pickReferece <- whichCov + (nParsGroup * (whichReference - 1))

            whichDis <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDis
                )

            whichDif <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDif
                )

            pickFocal <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickFocal
                )

            pickReferece <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickReferece
                )

            itemCov[["focal"]]     <- vcovMod[pickFocal, pickFocal]
            itemCov[["reference"]] <- vcovMod[pickReferece, pickReferece]

            itemCov <- vcov(mod)[whichCov, whichCov]

            transDisc <- paste0("~ x", which(whichCov %in% whichDis))
            transDiff <- paste0("~ -x", which(whichCov %in% whichDif), "/x", whichDif - 1)

            trans2PL <- list()
            trans2PL[which(whichCov %in% whichDis)] <- transDisc
            trans2PL[which(whichCov %in% whichDif)] <- transDiff

#            itemCov[["focal"]] <- msm::deltamethod(g = trans2PL,
            itemCov[["focal"]] <- deltamethod(g = lapply(trans2PL, as.formula),
                                        mean = as.numeric(t(itemPars)),
                                        cov = itemCov[["focal"]],
                                        ses = FALSE)

            itemCov[["reference"]] <- deltamethod(g = lapply(trans2PL, as.formula),
                                        mean = as.numeric(t(itemPars)),
                                        cov = itemCov[["reference"]],
                                        ses = FALSE)


        }
    } else {
        itemCov <- NULL
    }

    output <- list(itemParameters = itemPars, itemCovariances = itemCov)

    return(output)

}





#' Extract item discrimination, difficulties, and guessing parameters and estimate covariance estimates for 3PL items from a fitted mirt object
#' for one or two groups
#'
#' @param mod       A mirt object containing the fit of unidimensional model.
#' @param focal     Character. Required if mod is MultipleGroupClass, focal should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#' @param reference Character. Required if mod is MultipleGroupClass, reference should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#'
#' @return If mod contains any itemtype == "3PL", a list with the item parameters and the estimate covariances (if available).
#' If mod is SingleGroupClass, the list contains the item parameters as a matrix and the covariances as a list.
#' If mod is MultipleGroupClass, the list contains the item parameters and covariances for the focal and reference groups only.
#'
#' @examples
#' library(mirt)
#' data <- expand.table(LSAT7)
#' (mod1 <- mirt(data, model = 1, itemtype = "3PL", SE = TRUE))
#' (DFIT:::Extract3PLMirt(mod1))
#'
Extract3PLMirt <- function (mod, focal = NULL, reference = NULL) {

    classMod <- attributes(class(mod))[["package"]]

    if (is.null(classMod) || classMod != "mirt") {
        stop("mod is not a mirt obect")
    }
    if (mod@Model[["model"]] != 1) {
        stop("mod is not a unidimensional model.")
    }

    whichItems <- which(mod@Model[["itemtype"]] == "3PL")

    if (length(whichItems) == 0) {
        message("None of the items is itemtype 3PL.")
        return(NULL)
    }

    if (class(mod) == "SingleGroupClass") {
        itemPars <- as.matrix(
            #mirt::coef(mod,
            coef(mod,
                 IRTpars = TRUE,
                 simplify = TRUE)[["items"]][whichItems, c("a", "b", "g")])

    } else if (class(mod) == "MultipleGroupClass") {
        if (is.null(focal) || is.null(reference)) {
            stop("focal and reference group names are required.")
        }
        if (!(focal %in% mod@Data$groupNames)) {
            stop("focal does not match any group name.")
        }
        if (!(reference %in% mod@Data$groupNames)) {
            stop("reference does not match any group name.")
        }

        #itemPars <- lapply(mirt::coef(mod,
        itemPars <- lapply(coef(mod,
                                IRTpars = TRUE,
                                simplify = TRUE),
                           function (x) {
                               as.matrix(x[["items"]][whichItems, c("a", "b", "g")])
                           }
        )

        itemPars <- itemPars[c(focal, reference)]
        names(itemPars) <- c("focal", "reference")
    }

    if (mod@Options[["SE"]]) {
        if (class(mod) == "SingleGroupClass") {
            #itemCov <- as.list(diag(mirt::vcov(mod))[whichItems])

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)

            whichCov <- ordCoefs[whichItems, c("a1", "d", "g")]
            whichDis <- ordCoefs[whichItems, "a1"]
            whichDif <- ordCoefs[whichItems, "d"]
            whichGss <- ordCoefs[whichItems, "g"]

            whichCov <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichCov
                )

            whichDis <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDis
                )

            whichDif <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDif
                )

            whichGss <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichGss
                )

            itemCov <- vcov(mod)[whichCov, whichCov]

            transDisc <- paste0("~ x", which(whichCov %in% whichDis))
            transDiff <- paste0("~ -x", which(whichCov %in% whichDif), "/x", whichDif - 1)
            transGuss <- paste0("~ x", which(whichCov %in% whichGss))

            trans3PL <- list()
            trans3PL[which(whichCov %in% whichDis)] <- transDisc
            trans3PL[which(whichCov %in% whichDif)] <- transDiff
            trans3PL[which(whichCov %in% whichGss)] <- transGuss

#            itemCov <- msm::deltamethod(g = trans3PL,
            itemCov <- deltamethod(g = lapply(trans3PL, as.formula),
                                        mean = as.numeric(t(itemPars)),
                                        cov = itemCov,
                                        ses = FALSE)

        } else if (class(mod) == "MultipleGroupClass") {
            itemCov <- list()

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)
            nVcov     <- nrow(vcovMod)

            nGroups <- length(mod@Data$groupNames)
            nPars   <- nVcov / nGroups
            nParsGroup <- as.numeric(gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)[nPars + 1]) - 1

            whichFocal     <- which(mod@Data$groupNames == focal)
            whichReference <- which(mod@Data$groupNames == reference)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[[1]][["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)

            whichCov <- ordCoefs[whichItems, c("a1", "d", "g")]
            whichDis <- ordCoefs[whichItems, "a1"]
            whichDif <- ordCoefs[whichItems, "d"]
            whichGss <- ordCoefs[whichItems, "g"]

            pickFocal    <- whichCov + (nParsGroup * (whichFocal - 1))
            pickReferece <- whichCov + (nParsGroup * (whichReference - 1))

            whichDis <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDis
                )

            whichDif <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDif
                )

            whichGss <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichGss
                )

            pickFocal <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickFocal
                )

            pickReferece <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickReferece
                )

            itemCov[["focal"]]     <- vcovMod[pickFocal, pickFocal]
            itemCov[["reference"]] <- vcovMod[pickReferece, pickReferece]

            transDisc <- paste0("~ x", which(whichCov %in% whichDis))
            transDiff <- paste0("~ -x", which(whichCov %in% whichDif), "/x", whichDif - 1)
            transGuss <- paste0("~ x", which(whichCov %in% whichGss))

            trans3PL <- list()
            trans3PL[which(whichCov %in% whichDis)] <- transDisc
            trans3PL[which(whichCov %in% whichDif)] <- transDiff
            trans3PL[which(whichCov %in% whichGss)] <- transGuss

#            itemCov[["focal"]] <- msm::deltamethod(g = trans3PL,
            itemCov[["focal"]] <- deltamethod(g = lapply(trans3PL, as.formula),
                                        mean = as.numeric(t(itemPars[["focal"]])),
                                        cov = itemCov[["focal"]],
                                        ses = FALSE)

            itemCov[["reference"]] <- deltamethod(g = lapply(trans3PL, as.formula),
                                        mean = as.numeric(t(itemPars[["reference"]])),
                                        cov = itemCov[["reference"]],
                                        ses = FALSE)


        }
    } else {
        itemCov <- NULL
    }

    output <- list(itemParameters = itemPars, itemCovariances = itemCov)

    return(output)

}




#' Extract item discrimination, difficulties, guessing, and upper asymptote parameters and estimate covariance estimates for 4PL items from a fitted mirt object
#' for one or two groups
#'
#' @param mod       A mirt object containing the fit of unidimensional model.
#' @param focal     Character. Required if mod is MultipleGroupClass, focal should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#' @param reference Character. Required if mod is MultipleGroupClass, reference should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#'
#' @return If mod contains any itemtype == "4PL", a list with the item parameters and the estimate covariances (if available).
#' If mod is SingleGroupClass, the list contains the item parameters as a matrix and the covariances as a list.
#' If mod is MultipleGroupClass, the list contains the item parameters and covariances for the focal and reference groups only.
#'
#' @examples
#' library(mirt)
#' data <- expand.table(LSAT7)
#' (mod1 <- mirt(data, model = 1, itemtype = "4PL", SE = TRUE))
#' (DFIT:::Extract4PLMirt(mod1))
#'
Extract4PLMirt <- function (mod, focal = NULL, reference = NULL) {

    classMod <- attributes(class(mod))[["package"]]

    if (is.null(classMod) || classMod != "mirt") {
        stop("mod is not a mirt obect")
    }
    if (mod@Model[["model"]] != 1) {
        stop("mod is not a unidimensional model.")
    }

    whichItems <- which(mod@Model[["itemtype"]] == "4PL")

    if (length(whichItems) == 0) {
        message("None of the items is itemtype 4PL.")
        return(NULL)
    }

    if (class(mod) == "SingleGroupClass") {
        itemPars <- as.matrix(
            #mirt::coef(mod,
            coef(mod,
                 IRTpars = TRUE,
                 simplify = TRUE)[["items"]][whichItems, c("a", "b", "g", "u")])

    } else if (class(mod) == "MultipleGroupClass") {
        if (is.null(focal) || is.null(reference)) {
            stop("focal and reference group names are required.")
        }
        if (!(focal %in% mod@Data$groupNames)) {
            stop("focal does not match any group name.")
        }
        if (!(reference %in% mod@Data$groupNames)) {
            stop("reference does not match any group name.")
        }

        #itemPars <- lapply(mirt::coef(mod,
        itemPars <- lapply(coef(mod,
                                IRTpars = TRUE,
                                simplify = TRUE),
                           function (x) {
                               as.matrix(x[["items"]][whichItems, c("a", "b", "g", "u")])
                           }
        )

        itemPars <- itemPars[c(focal, reference)]
        names(itemPars) <- c("focal", "reference")
    }

    if (mod@Options[["SE"]]) {
        if (class(mod) == "SingleGroupClass") {
            #itemCov <- as.list(diag(mirt::vcov(mod))[whichItems])

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)

            whichCov <- ordCoefs[whichItems, c("a1", "d", "g", "u")]
            whichDis <- ordCoefs[whichItems, "a1"]
            whichDif <- ordCoefs[whichItems, "d"]
            whichGss <- ordCoefs[whichItems, "g"]
            whichUps <- ordCoefs[whichItems, "u"]

            whichCov <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichCov
                )

            whichDis <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDis
                )

            whichDif <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDif
                )

            whichGss <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichGss
                )

            whichUps <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichUps
                )

            itemCov <- vcov(mod)[whichCov, whichCov]

            transDisc <- paste0("~ x", which(whichCov %in% whichDis))
            transDiff <- paste0("~ -x", which(whichCov %in% whichDif), "/x", whichDif - 1)
            transGuss <- paste0("~ x", which(whichCov %in% whichGss))
            transUpas <- paste0("~ x", which(whichCov %in% whichUps))

            trans4PL <- list()
            trans4PL[which(whichCov %in% whichDis)] <- transDisc
            trans4PL[which(whichCov %in% whichDif)] <- transDiff
            trans4PL[which(whichCov %in% whichGss)] <- transGuss
            trans4PL[which(whichCov %in% whichUps)] <- transUpas

#            itemCov <- msm::deltamethod(g = trans3PL,
            itemCov <- deltamethod(g = lapply(trans4PL, as.formula),
                                        mean = as.numeric(t(itemPars)),
                                        cov = itemCov,
                                        ses = FALSE)

        } else if (class(mod) == "MultipleGroupClass") {
            itemCov <- list()

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)
            nVcov     <- nrow(vcovMod)

            nGroups <- length(mod@Data$groupNames)
            nPars   <- nVcov / nGroups
            nParsGroup <- as.numeric(gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)[nPars + 1]) - 1

            whichFocal     <- which(mod@Data$groupNames == focal)
            whichReference <- which(mod@Data$groupNames == reference)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[[1]][["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)

            whichCov <- ordCoefs[whichItems, c("a1", "d", "g", "u")]
            whichDis <- ordCoefs[whichItems, "a1"]
            whichDif <- ordCoefs[whichItems, "d"]
            whichGss <- ordCoefs[whichItems, "g"]
            whichUps <- ordCoefs[whichItems, "u"]

            pickFocal    <- whichCov + (nParsGroup * (whichFocal - 1))
            pickReferece <- whichCov + (nParsGroup * (whichReference - 1))

            whichDis <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDis
                )

            whichDif <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDif
                )

            whichGss <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichGss
                )

            whichUps <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichUps
                )

            pickFocal <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickFocal
                )

            pickReferece <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickReferece
                )

            itemCov[["focal"]]     <- vcovMod[pickFocal, pickFocal]
            itemCov[["reference"]] <- vcovMod[pickReferece, pickReferece]

            transDisc <- paste0("~ x", which(whichCov %in% whichDis))
            transDiff <- paste0("~ -x", which(whichCov %in% whichDif), "/x", whichDif - 1)
            transGuss <- paste0("~ x", which(whichCov %in% whichGss))
            transUpas <- paste0("~ x", which(whichCov %in% whichUps))

            trans4PL <- list()
            trans4PL[which(whichCov %in% whichDis)] <- transDisc
            trans4PL[which(whichCov %in% whichDif)] <- transDiff
            trans4PL[which(whichCov %in% whichGss)] <- transGuss
            trans4PL[which(whichCov %in% whichUps)] <- transUpas


#            itemCov[["focal"]] <- msm::deltamethod(g = trans4PL,
            itemCov[["focal"]] <- deltamethod(g = lapply(trans4PL, as.formula),
                                        mean = as.numeric(t(itemPars[["focal"]])),
                                        cov = itemCov[["focal"]],
                                        ses = FALSE)

            itemCov[["reference"]] <- deltamethod(g = lapply(trans4PL, as.formula),
                                        mean = as.numeric(t(itemPars[["reference"]])),
                                        cov = itemCov[["reference"]],
                                        ses = FALSE)


        }
    } else {
        itemCov <- NULL
    }

    output <- list(itemParameters = itemPars, itemCovariances = itemCov)

    return(output)

}




#' Extract item discrimination and difficulties and estimate covariance estimates for GRM items from a fitted mirt object
#' for one or two groups
#'
#' @param mod       A mirt object containing the fit of unidimensional model.
#' @param focal     Character. Required if mod is MultipleGroupClass, focal should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#' @param reference Character. Required if mod is MultipleGroupClass, reference should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#'
#' @return If mod contains any itemtype == "graded", a list with the item parameters and the estimate covariances (if available).
#' If mod is SingleGroupClass, the list contains the item parameters as a matrix and the covariances as a list.
#' If mod is MultipleGroupClass, the list contains the item parameters and covariances for the focal and reference groups only.
#'
#' @examples
#' library(mirt)
#' (mod1 <- mirt(Science, model = 1, itemtype = "graded", SE = TRUE))
#' (DFIT:::ExtractGRMMirt(mod1))
#'
ExtractGRMMirt <- function (mod, focal = NULL, reference = NULL) {

    classMod <- attributes(class(mod))[["package"]]

    if (is.null(classMod) || classMod != "mirt") {
        stop("mod is not a mirt obect")
    }
    if (mod@Model[["model"]] != 1) {
        stop("mod is not a unidimensional model.")
    }

    whichItems <- which(mod@Model[["itemtype"]] == "graded")

    if (length(whichItems) == 0) {
        message("None of the items is itemtype GRM (graded).")
        return(NULL)
    }

    if (class(mod) == "SingleGroupClass") {
        itemPars <- as.matrix(
            #mirt::coef(mod,
            coef(mod,
                 IRTpars = TRUE,
                 simplify = TRUE)[["items"]][whichItems, ])
        whichCols <- grep(pattern = "^(a|b)",
                          x = colnames(itemPars))
        itemPars <- itemPars[, whichCols]

    } else if (class(mod) == "MultipleGroupClass") {
        if (is.null(focal) || is.null(reference)) {
            stop("focal and reference group names are required.")
        }
        if (!(focal %in% mod@Data$groupNames)) {
            stop("focal does not match any group name.")
        }
        if (!(reference %in% mod@Data$groupNames)) {
            stop("reference does not match any group name.")
        }

        #itemPars <- lapply(mirt::coef(mod,
        itemPars <- lapply(coef(mod,
                                IRTpars = TRUE,
                                simplify = TRUE),
                           function (x) {
                               as.matrix(x[["items"]][whichItems, ])
                           }
        )

        whichCols <- grep(pattern = "^(a|b)",
                          x = colnames(itemPars[[1]]))
        itemPars <- lapply(itemPars,
                           function (x) x[, whichCols]
        )

        itemPars <- itemPars[c(focal, reference)]
        names(itemPars) <- c("focal", "reference")
    }

    if (mod@Options[["SE"]]) {
        if (class(mod) == "SingleGroupClass") {
            #itemCov <- as.list(diag(mirt::vcov(mod))[whichItems])

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)


            diffNames <- grep("d[1-9]", x = colnames(ordCoefs), value = TRUE)
            whichCov <- ordCoefs[whichItems, c("a1", diffNames)]
            whichDis <- ordCoefs[whichItems, "a1"]
            whichDif <- ordCoefs[whichItems, diffNames]

            whichCov <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichCov
                )

            whichDis <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDis
                )

            whichDif <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDif
                )

            itemCov <- vcov(mod)[whichCov, whichCov]

            transDisc <- paste0("~ x", which(whichCov %in% whichDis))
            transDiff <- paste0("~ -x", which(whichCov %in% whichDif), "/x",
                                rep(which(whichCov %in% whichDis), each = length(whichDif) / length(whichDis)))

            transGRM <- list()
            transGRM[which(whichCov %in% whichDis)] <- transDisc
            transGRM[which(whichCov %in% whichDif)] <- transDiff

#            itemCov <- msm::deltamethod(g = transGRM,
            itemCov <- deltamethod(g = lapply(transGRM, as.formula),
                                        mean = as.numeric(t(itemPars)),
                                        cov = itemCov,
                                        ses = FALSE)

        } else if (class(mod) == "MultipleGroupClass") {
            itemCov <- list()

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)
            nVcov     <- nrow(vcovMod)

            nGroups <- length(mod@Data$groupNames)
            nPars   <- nVcov / nGroups
            nParsGroup <- as.numeric(gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)[nPars + 1]) - 1

            whichFocal     <- which(mod@Data$groupNames == focal)
            whichReference <- which(mod@Data$groupNames == reference)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)


            diffNames <- grep("d[1-9]", x = colnames(ordCoefs), value = TRUE)
            whichCov <- ordCoefs[whichItems, c("a1", diffNames)]
            whichDis <- ordCoefs[whichItems, "a1"]
            whichDif <- ordCoefs[whichItems, diffNames]

            pickFocal    <- whichCov + (nParsGroup * (whichFocal - 1))
            pickReferece <- whichCov + (nParsGroup * (whichReference - 1))

            whichDis <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDis
                )

            whichDif <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDif
                )

            pickFocal <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickFocal
                )

            pickReferece <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickReferece
                )

            itemCov[["focal"]]     <- vcovMod[pickFocal, pickFocal]
            itemCov[["reference"]] <- vcovMod[pickReferece, pickReferece]

            itemCov <- vcov(mod)[whichCov, whichCov]

            transDisc <- paste0("~ x", which(whichCov %in% whichDis))
            transDiff <- paste0("~ -x", which(whichCov %in% whichDif), "/x",
                                rep(which(whichCov %in% whichDis), each = length(whichDif) / length(whichDis)))

            transGRM <- list()
            transGRM[which(whichCov %in% whichDis)] <- transDisc
            transGRM[which(whichCov %in% whichDif)] <- transDiff

#            itemCov[["focal"]] <- msm::deltamethod(g = transGRM,
            itemCov[["focal"]] <- deltamethod(g = lapply(transGRM, as.formula),
                                        mean = as.numeric(t(itemPars)),
                                        cov = itemCov[["focal"]],
                                        ses = FALSE)

            itemCov[["reference"]] <- deltamethod(g = lapply(transGRM, as.formula),
                                        mean = as.numeric(t(itemPars)),
                                        cov = itemCov[["reference"]],
                                        ses = FALSE)

        }
    } else {
        itemCov <- NULL
    }

    output <- list(itemParameters = itemPars, itemCovariances = itemCov)

    return(output)

}



#' Extract item discrimination and difficulties and estimate covariance estimates for GPCM items from a fitted mirt object
#' for one or two groups
#'
#' @param mod       A mirt object containing the fit of unidimensional model.
#' @param focal     Character. Required if mod is MultipleGroupClass, focal should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#' @param reference Character. Required if mod is MultipleGroupClass, reference should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#'
#' @return If mod contains any itemtype == "gpcm", a list with the item parameters and the estimate covariances (if available).
#' If mod is SingleGroupClass, the list contains the item parameters as a matrix and the covariances as a list.
#' If mod is MultipleGroupClass, the list contains the item parameters and covariances for the focal and reference groups only.
#'
#' @examples
#' library(mirt)
#' (mod1 <- mirt(Science, model = 1, itemtype = "gpcm", SE = TRUE))
#' (DFIT:::ExtractGPCMMirt(mod1))
#'
ExtractGPCMMirt <- function (mod, focal = NULL, reference = NULL) {

    classMod <- attributes(class(mod))[["package"]]

    if (is.null(classMod) || classMod != "mirt") {
        stop("mod is not a mirt obect")
    }
    if (mod@Model[["model"]] != 1) {
        stop("mod is not a unidimensional model.")
    }

    whichItems <- which(mod@Model[["itemtype"]] == "gpcm")

    if (length(whichItems) == 0) {
        message("None of the items is itemtype GPCM (gpcm).")
        return(NULL)
    }

    if (class(mod) == "SingleGroupClass") {
        itemPars <- as.matrix(
            #mirt::coef(mod,
            coef(mod,
                 IRTpars = TRUE,
                 simplify = TRUE)[["items"]][whichItems, ])
        whichCols <- grep(pattern = "^(a|b)",
                          x = colnames(itemPars))
        itemPars <- itemPars[, whichCols]

    } else if (class(mod) == "MultipleGroupClass") {
        if (is.null(focal) || is.null(reference)) {
            stop("focal and reference group names are required.")
        }
        if (!(focal %in% mod@Data$groupNames)) {
            stop("focal does not match any group name.")
        }
        if (!(reference %in% mod@Data$groupNames)) {
            stop("reference does not match any group name.")
        }

        #itemPars <- lapply(mirt::coef(mod,
        itemPars <- lapply(coef(mod,
                                IRTpars = TRUE,
                                simplify = TRUE),
                           function (x) {
                               as.matrix(x[["items"]][whichItems, ])
                           }
        )

        whichCols <- grep(pattern = "^(a|b)",
                          x = colnames(itemPars[[1]]))
        itemPars <- lapply(itemPars,
                           function (x) x[, whichCols]
        )

        itemPars <- itemPars[c(focal, reference)]
        names(itemPars) <- c("focal", "reference")
    }

    if (mod@Options[["SE"]]) {
        if (class(mod) == "SingleGroupClass") {
            #itemCov <- as.list(diag(mirt::vcov(mod))[whichItems])

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)


            diffNames <- grep("d[1-9]", x = colnames(ordCoefs), value = TRUE)
            whichCov <- ordCoefs[whichItems, c("a1", diffNames)]
            whichDis <- ordCoefs[whichItems, "a1"]
            whichDif <- ordCoefs[whichItems, diffNames]

            whichCov <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichCov
                )

            whichDis <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDis
                )

            whichDif <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDif
                )

            itemCov <- vcov(mod)[whichCov, whichCov]

            transDisc <- paste0("~ x", which(whichCov %in% whichDis))
            transDiff <- paste0("~ -x", which(whichCov %in% whichDif), "/x",
                                rep(which(whichCov %in% whichDis), each = length(whichDif) / length(whichDis)))

            transGPCM <- list()
            transGPCM[which(whichCov %in% whichDis)] <- transDisc
            transGPCM[which(whichCov %in% whichDif)] <- transDiff

#            itemCov <- msm::deltamethod(g = transGPCM,
            itemCov <- deltamethod(g = lapply(transGPCM, as.formula),
                                        mean = as.numeric(t(itemPars)),
                                        cov = itemCov,
                                        ses = FALSE)

        } else if (class(mod) == "MultipleGroupClass") {
            itemCov <- list()

            vcovMod   <- vcov(mod)
            namesVcov <- rownames(vcovMod)
            nVcov     <- nrow(vcovMod)

            nGroups <- length(mod@Data$groupNames)
            nPars   <- nVcov / nGroups
            nParsGroup <- as.numeric(gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)[nPars + 1]) - 1

            whichFocal     <- which(mod@Data$groupNames == focal)
            whichReference <- which(mod@Data$groupNames == reference)

            strMatrix              <- !is.na(coef(mod, simplify = TRUE)[["items"]])
            strMatrix              <- strMatrix[, sort(colnames(strMatrix))]
            ordCoefs               <- t(strMatrix)
            ordCoefs[t(strMatrix)] <- seq(sum(strMatrix))
            ordCoefs               <- t(ordCoefs)


            diffNames <- grep("d[1-9]", x = colnames(ordCoefs), value = TRUE)
            whichCov <- ordCoefs[whichItems, c("a1", diffNames)]
            whichDis <- ordCoefs[whichItems, "a1"]
            whichDif <- ordCoefs[whichItems, diffNames]

            pickFocal    <- whichCov + (nParsGroup * (whichFocal - 1))
            pickReferece <- whichCov + (nParsGroup * (whichReference - 1))

            whichDis <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDis
                )

            whichDif <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% whichDif
                )

            pickFocal <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickFocal
                )

            pickReferece <- which(as.numeric(
                gsub(pattern = "^.+\\.", replacement = "", x = namesVcov)) %in% pickReferece
                )

            itemCov[["focal"]]     <- vcovMod[pickFocal, pickFocal]
            itemCov[["reference"]] <- vcovMod[pickReferece, pickReferece]

            itemCov <- vcov(mod)[whichCov, whichCov]

            transDisc <- paste0("~ x", which(whichCov %in% whichDis))
            transDiff <- paste0("~ -x", which(whichCov %in% whichDif), "/x",
                                rep(which(whichCov %in% whichDis), each = length(whichDif) / length(whichDis)))

            transGPCM <- list()
            transGPCM[which(whichCov %in% whichDis)] <- transDisc
            transGPCM[which(whichCov %in% whichDif)] <- transDiff

#            itemCov[["focal"]] <- msm::deltamethod(g = transGPCM,
            itemCov[["focal"]] <- deltamethod(g = lapply(transGPCM, as.formula),
                                        mean = as.numeric(t(itemPars)),
                                        cov = itemCov[["focal"]],
                                        ses = FALSE)

            itemCov[["reference"]] <- deltamethod(g = lapply(transGPCM, as.formula),
                                        mean = as.numeric(t(itemPars)),
                                        cov = itemCov[["reference"]],
                                        ses = FALSE)

        }
    } else {
        itemCov <- NULL
    }

    output <- list(itemParameters = itemPars, itemCovariances = itemCov)

    return(output)

}





#' Extracts the item parameters from a unidimensional mirt model
#'
#' @param mod       A mirt object containing the fit of unidimensional model.
#' @param focal     Character. Required if mod is MultipleGroupClass, focal should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#' @param reference Character. Required if mod is MultipleGroupClass, reference should coincide with the label for the focal group. If mod is SingleGroupClass, it is ignored.
#'
#' @return If mod contains any itemtype == "gpcm", a list with the item parameters and the estimate covariances (if available).
#' If mod is SingleGroupClass, the list contains the item parameters as a matrix and the covariances as a list.
#' If mod is MultipleGroupClass, the list contains the item parameters and covariances for the focal and reference groups only.
#'
#' @examples
#' library(mirt)
#' (mod1 <- mirt(Science, model = 1, itemtype = c("graded", "graded", "gpcm", "gpcm"), SE = TRUE))
#' (ExtractMirtPars(mod1))
#'
#' @return
#' @export
#'
#' @examples
ExtractMirtPars <- function (mod, focal = NULL, reference = NULL) {

    modItemTypes <- unique(mod@Model[["itemtype"]])

    mirtParameters <- list()

    for (iiType in modItemTypes) {

        ExtractPars <- switch(iiType,
                            "Rasch"  = ExtractRaschMirt,
                            "2PL"    = Extract2PLMirt,
                            "3PL"    = Extract3PLMirt,
                            "4PL"    = Extract4PLMirt,
                            "graded" = ExtractGRMMirt,
                            "gpcm"   = ExtractGPCMMirt,
                            warning(paste("irtModel ", iiType, " not known or not implemented"))
        )

        mirtParameters[[iiType]] <- ExtractPars(mod = mod, focal = focal, reference = reference)
   }

    return(mirtParameters)
}





