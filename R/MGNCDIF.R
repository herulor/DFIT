################################################################################
# # MGNCDIF.R
# # R Versions: 4.1.2+
# #
# # Author(s): Victor H. Cervantes
# #
# # General DFIT functions for multiple groups
# # Description: Functions for calculating MGNCDIF
# #
# # Inputs: NULL
# #
# # Outputs: functions
# #
# # File history:
# #   20211110: Creation
################################################################################



################################################################################
# # Function MGNcdif: Multi-Group NonCompensatory Differential Item Functioning index
################################################################################

#' Calculates the MG-NCDIF index for an item with given item parameters of two or more groups.
#'
#' @param itemParameters    A list of three or more sets of item parameters.
#'                          Item parameters are assumed to be on the same scale.
#'                          Item parameters for each group should be a matrix with nrow equal to the number of items.
#'                          One set of item parameters must be called "base" and contain the item parameters to which
#'                          all other sets will be compared to; the other group item parameters must have some name but
#'                          it can be arbitrarily given.
#' @param irtModel          A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param focalAbilities    If NULL, NCDIF is calculated by numerical integration of focal distributions.
#'                          If not NULL, it must be a list of numerical vectors containing the abilities for
#'                          the individuals in each group.
#'                          Names of the list must coincide with the names of the itemParameters list that are not the "base" group.
#' @param focalDistribution Either a string or a list of strings stating the distribution name to be used for integrating for each non base group.
#'                          Only used if focalAbilities is NULL.
#'                          If a list, names of the list must coincide with the names of the itemParameters list that are not the "base" group.
#'                          IF a single string, the same distribution will be used for all pairwise computation of NCDIF.
#' @param focalDistrExtra   Extra parameters for the focal group distribution function if needed.
#'                          If focalDistribution is a list, then names(focalDistrExtra) must coincide with names(focalDistribution).
#' @param relSizes          If not NULL, a numeric vector with named components containing the relative sizes of the groups.
#'                          Used to obtain the weighted average of the pairwise NCDIF statistics.
#'                          If NULL and focalAbilities is NULL, defaults to 1/length(focalDistribution)
#' @param subdivisions      A numeric value indicating the number of subdivisions for numerical integration. Only used if focalAbilities is NULL.
#' @param logistic          A logical value stating if the IRT model will use the logistic or the normal metric. Defaults to using the logistic metric by fixing the D constant to 1. If FALSE the constant is set to 1.702 so that the normal metric is used.
#'
#' @return mgncdif            Numeric vector with the MG-NCDIF index value for each item.
#'
#' @examples
#'
#' # Data fom Oshima, Wright and White
#' baseParameters <- matrix(c(0.49, -0.07, 0.19,
#'                            0.92,  0.21, 0.15,
#'                            1.26,  0.54, 0.05,
#'                            0.61, -0.03, 0.18,
#'                            1.74,  0.01, 0.12,
#'                            0.50,  1.96, 0.12,
#'                            0.96,  0.04, 0.13,
#'                            0.59, -0.09, 0.18,
#'                            0.82, -1.16, 0.17,
#'                            1.26,  0.02, 0.11,
#'                            0.82,  0.20, 0.07,
#'                            0.75, -0.43, 0.15,
#'                            1.49, -0.06, 0.09,
#'                            0.97, -0.34, 0.12,
#'                            1.49,  0.05, 0.12,
#'                            0.89, -0.25, 0.15,
#'                            1.45,  0.06, 0.07,
#'                            0.75,  0.31, 0.18,
#'                            1.43,  0.04, 0.08,
#'                            0.60,  0.13, 0.22,
#'                            0.83,  0.52, 0.09,
#'                            0.56, -0.96, 0.19,
#'                            0.67, -0.79, 0.20,
#'                            0.70,  0.37, 0.18,
#'                            1.03, -0.71, 0.14,
#'                            0.89, -0.19, 0.21,
#'                            1.23,  0.74, 0.06,
#'                            0.90, -0.44, 0.18,
#'                            1.23, -0.17, 0.12,
#'                            0.69, 0.53, 0.17),
#'                            byrow = TRUE, ncol = 3)
#'
#'
#' group1Pars <- group2Pars <- group3Pars <- group4Pars <- group5Pars <- baseParameters
#'
#' group1Pars[c(4,  22), 1] <- group1Pars[c(4,  22), 1] + .4
#' group1Pars[c(7,  22), 2] <- group1Pars[c(7,  22), 2] + .7
#'
#' group2Pars[c(13, 28), 1] <- group2Pars[c(13, 28), 1] - .4
#' group2Pars[c(19, 28), 2] <- group2Pars[c(19, 28), 2] - .7
#'
#' group3Pars[c(4,  13), 1] <- group3Pars[c(4,  13), 1] + .4
#' group3Pars[c(4,  22), 2] <- group3Pars[c(4,  22), 2] - .7
#'
#' group4Pars[c(7,  28), 1] <- group4Pars[c(7,  28), 1] - .4
#' group4Pars[c(7,  19), 2] <- group4Pars[c(7,  19), 2] + .7
#'
#' itemParameters <- list(base = baseParameters,
#'                        group1 = group1Pars,
#'                        group2 = group2Pars,
#'                        group3 = group3Pars,
#'                        group4 = group4Pars,
#'                        group5 = group5Pars
#'                        )
#'
#' itemCovariances <- lapply(itemParameters, AseIrt, irtModel = "2pl")
#'
#' relSizes        <- rep(.2, 5)
#' names(relSizes) <- names(itemParameters)[-1]
#'
#' threePlMGNcdif <- MGNcdif(itemParameters = itemParameters, irtModel = '3pl',
#'                       focalAbilities = NULL, focalDistribution = "norm",
#'                       relSizes = relSizes,
#'                       subdivisions = 5000, logistic = TRUE)
#'
#' @references Oshima, T. C., Wright, K., & White, N. (2014). Multiple-Group Noncompensatory Differential Item Functioning in Raju’s Differential Functioning of Items and Tests. International Journal of Testing, 15, 254–273.
#'
#' @author Victor H. Cervantes <vhcervantesb at unal.edu.co>
#'
#' @export
MGNcdif <- function (itemParameters, irtModel = "2pl", focalAbilities = NULL, focalDistribution = "norm", relSizes = NULL,
                     subdivisions = 5000, logistic = TRUE, focalDistrExtra = list(list(mean = 0, sd = 1))) {

  if (!("base" %in% names(itemParameters))) {
    stop("To compute MGNcdif a base group must be defined. The base group is identified as 'base' among itemParameters. See ?MGNcdif")
  }
  if (length(itemParameters) < 3) {
    stop("At least three sets of items parameters are needed to compute MGNcdif: One base group and at least two groups that are compared to the base group.")
  }

  baseParameters <- itemParameters[["base"]]

  groups <- which(names(itemParameters) != "base")

  itemParameters <- itemParameters[groups]

  if (!is.null(focalAbilities)) {

    if (length(focalAbilities) != length(itemParameters)) {
      stop("If focalAbilities is given, it must be a list of the same length as the itemParameters minus one.")
    }
    if (!all(names(focalAbilities) %in% names(itemParameters))) {
      stop("If focalAbilities is given, the names of the list must coincide with the names of the itemParameters list that are not the 'base' group.")
    }
    if (!all(sapply(focalAbilities, is.numeric))) {
      stop("Each element of focalAbilities must be a numeric vector.")
    }

    relSizes <- sapply(focalAbilities, length)
    relSizes <- relSizes / sum(relSizes)
    names(relSizes) <- names(itemParameters)

  } else {
    if (length(focalDistribution) == 1) {
      focalDistribution <- as.list(rep(focalDistribution, length(itemParameters)))
      focalDistrExtra   <- as.list(rep(focalDistrExtra, length(itemParameters)))

      names(focalDistribution) <- names(focalDistrExtra) <- names(itemParameters)
    }
    if (length(focalDistribution) != length(itemParameters)) {
      stop("If focalDistribution is given, it must be a list of the same length as the itemParameters minus one.")
    }
    if (!all(names(focalDistribution) %in% names(itemParameters))) {
      stop("If focalDistribution is given, the names of the list must coincide with the names of the itemParameters list that are not the 'base' group.")
    }
    if (!all(names(focalDistribution) %in% names(focalDistrExtra))) {
      stop("If focalDistribution is given, focalDistrExtra must be a list of lists and names(focalDistrExtra) must coincide names(focalDistribution).")
    }

    if (is.null(relSizes)) {
      relSizes <- rep(1 / length(focalDistribution), length(focalDistribution))
      names(relSizes) <- names(itemParameters)
    }

  }

  mgncdif <- numeric(nrow(baseParameters))

  for (iiGroup in names(itemParameters)) {
    groupParameters <- list(focal     = itemParameters[[iiGroup]],
                            reference = baseParameters)
    mgncdif <- mgncdif + (relSizes[iiGroup] *
                            Ncdif(itemParameters = groupParameters, irtModel = irtModel,
                                  focalAbilities    = focalAbilities[[iiGroup]],
                                  focalDistribution = focalDistribution[[iiGroup]],
                                  focalDistrExtra   = focalDistrExtra[[iiGroup]],
                                  subdivisions = subdivisions, logistic = logistic)
    )
  }

  return(mgncdif)

}






################################################################################
# # Function IprMG: Item parameter replication for multiple groups
################################################################################

#' Item parameter replication
#' @description Generates a sample of item parameters assuming multivariate normality of estimates
#'
#' @param itemParameters    A list of three or more sets of item parameters.
#'                          Item parameters are assumed to be on the same scale.
#'                          Item parameters for each group should be a matrix with nrow equal to the number of items.
#'                          One set of item parameters must be called "base" and contain the item parameters to which
#'                          all other sets will be compared to; the other group item parameters must have some name but
#'                          it can be arbitrarily given.
#' @param itemCovariances A list containing matrices of covariance for item estimates.
#'                        See 'itemParameters' for list structure.
#'                        Each list element may be either a list of covariance matrices for each item or a single matrix of covariance of all parameters.
#' @param nReplicates     A numeric value indicating the number of replications to perform
#'
#' @return itemParameters A list with item parameters for focal and reference groups
#'
#' @importFrom simex diag.block
#'
#' @examples
#'#' # Data fom Oshima, Wright and White
#' baseParameters <- matrix(c(0.49, -0.07, 0.19,
#'                            0.92,  0.21, 0.15,
#'                            1.26,  0.54, 0.05,
#'                            0.61, -0.03, 0.18,
#'                            1.74,  0.01, 0.12,
#'                            0.50,  1.96, 0.12,
#'                            0.96,  0.04, 0.13,
#'                            0.59, -0.09, 0.18,
#'                            0.82, -1.16, 0.17,
#'                            1.26,  0.02, 0.11,
#'                            0.82,  0.20, 0.07,
#'                            0.75, -0.43, 0.15,
#'                            1.49, -0.06, 0.09,
#'                            0.97, -0.34, 0.12,
#'                            1.49,  0.05, 0.12,
#'                            0.89, -0.25, 0.15,
#'                            1.45,  0.06, 0.07,
#'                            0.75,  0.31, 0.18,
#'                            1.43,  0.04, 0.08,
#'                            0.60,  0.13, 0.22,
#'                            0.83,  0.52, 0.09,
#'                            0.56, -0.96, 0.19,
#'                            0.67, -0.79, 0.20,
#'                            0.70,  0.37, 0.18,
#'                            1.03, -0.71, 0.14,
#'                            0.89, -0.19, 0.21,
#'                            1.23,  0.74, 0.06,
#'                            0.90, -0.44, 0.18,
#'                            1.23, -0.17, 0.12,
#'                            0.69, 0.53, 0.17),
#'                            byrow = TRUE, ncol = 3)
#'
#'
#' group1Pars <- group2Pars <- group3Pars <- group4Pars <- group5Pars <- baseParameters
#'
#' group1Pars[c(4,  22), 1] <- group1Pars[c(4,  22), 1] + .4
#' group1Pars[c(7,  22), 2] <- group1Pars[c(7,  22), 2] + .7
#'
#' group2Pars[c(13, 28), 1] <- group2Pars[c(13, 28), 1] - .4
#' group2Pars[c(19, 28), 2] <- group2Pars[c(19, 28), 2] - .7
#'
#' group3Pars[c(4,  13), 1] <- group3Pars[c(4,  13), 1] + .4
#' group3Pars[c(4,  22), 2] <- group3Pars[c(4,  22), 2] - .7
#'
#' group4Pars[c(7,  28), 1] <- group4Pars[c(7,  28), 1] - .4
#' group4Pars[c(7,  19), 2] <- group4Pars[c(7,  19), 2] + .7
#'
#' itemParameters <- list(base = baseParameters,
#'                        group1 = group1Pars,
#'                        group2 = group2Pars,
#'                        group3 = group3Pars,
#'                        group4 = group4Pars,
#'                        group5 = group5Pars
#'                        )
#'
#' itemCovariances <- lapply(itemParameters, AseIrt, irtModel = "3pl", sampleSize = 5000)
#'
#' mgIpr <- IprMG(itemParameters = itemParameters,
#'                itemCovariances = itemCovariances,
#'                nReplicates = 100)
#'
#' @references Oshima, T. C., Wright, K., & White, N. (2014). Multiple-Group Noncompensatory Differential Item Functioning in Raju’s Differential Functioning of Items and Tests. International Journal of Testing, 15, 254–273.
#' @references Oshima, T., Raju, N. & Nanda, A. (2006). A new method for assessing the statistical significance in the Differential Functioning of Items and Tests (DFIT) framework. Journal of educational measurement, 43(1), 1--17. doi:10.1111/j.1745-3984.2006.00001.x
#'
#'
#' @author Victor H. Cervantes <vhcervantesb at unal.edu.co>
#'
#' @export
IprMG <- function (itemParameters, itemCovariances, nReplicates = 5000) {

  # # Data check
  if (all(length(itemParameters[["base"]]) == sapply(itemParameters, length))) {
    nItem <- nrow(itemParameters[["base"]])
  }  else {
    stop("There must be the same number of items for all groups")
  }
  if (all(length(itemCovariances[["base"]]) == sapply(itemCovariances, length))) {
    nCovs <- length(itemCovariances[["base"]])
  }  else {
    stop("There must be the same number of item covariance matrices for all groups")
  }

  if (length(nReplicates) > 1 | !is.numeric(nReplicates)) {
    stop("nReplicates must be a single numeric value")
  }
  if (nItem != nCovs) {
    stop("The number of item parameter vectors must be equal to the number of covariance matrices for each group")
  }

  # # IPR for each group
  itemPars <- as.numeric(t(itemParameters[["base"]]))

  groupNames <- names(itemParameters)

  groupParameters <- list()

  for (iiGroup in groupNames) {
    if (length(itemCovariances[["base"]]) == nrow(itemParameters[[iiGroup]])) {
      itemCovs <- diag.block(itemCovariances[[iiGroup]])
    } else {
      itemCovs <- itemCovariances[[iiGroup]]
    }

    groupParameters[[iiGroup]] <- rmvnorm(n = nReplicates, mean = itemPars, sigma = itemCovs, method = "chol")
    groupParameters[[iiGroup]] <- tapply(groupParameters[[iiGroup]] , row(groupParameters[[iiGroup]]),
                                         function (x) matrix(x, nrow = nItem, byrow = TRUE))
  }

  # # Join the lists for each replication
  itemParameterList <- list()
  for (ii in seq(length(groupParameters[["base"]]))) {
    itemParameterList[[ii]]                <- list()
    for (iiGroup in groupNames) {
      itemParameterList[[ii]][[iiGroup]] <- groupParameters[[iiGroup]][[ii]]
    }
  }

  return(itemParameterList)
}






################################################################################
# # Function IprMGNcdif: NCDIF for Item parameter replication
################################################################################

#' MGNCDIF for Item parameter replication
#' @description Calculates the NCDIF index on a list of item parameters such as those produced by the Ipr function
#'
#' @param itemParameterList A list where each element is a list containing "focal" and "reference" item Parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should be a matrix with nrow equal to the number of items.
#' @param irtModel          A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param focalAbilities    If NULL, NCDIF is calculated by numerical integration of focal distribution. If not NULL, must be a numerical vector containing the abilities for the individuals in the focal group.
#' @param focalDistribution A string stating the distribution name to be used for integrating. Only used if focalAbilities is NULL.
#' @param focalDistrExtra   A list stating the extra parameters needed by the focal distribution function.
#' @param relSizes          If not NULL, a numeric vector with named components containing the relative sizes of the groups.
#'                          Used to obtain the weighted average of the pairwise NCDIF statistics.
#'                          If NULL and focalAbilities is NULL, defaults to 1/length(focalDistribution)
#' @param subdivisions      A numeric value indicating the number of subdivisions for numerical integration. Only used if focalAbilities is NULL.
#' @param logistic          A logical value stating if the IRT model will use the logistic or the normal metric. Defaults to using the logistic metric by fixing the D constant to 1. If FALSE the constant is set to 1.702 so that the normal metric is used.
#'
#' @return mgncdif A numeric matrix with the NCDIF values for all the item parameter in each set of itemParameterList
#'
#' @references Oshima, T., Raju, N. & Nanda, A. (2006). A new method for assessing the statistical significance in the Differential Functioning of Items and Tests (DFIT) framework. Journal of educational measurement, 43(1), 1--17. doi:10.1111/j.1745-3984.2006.00001.x
#'
#' @examples
#'#' # Data fom Oshima, Wright and White
#' baseParameters <- matrix(c(0.49, -0.07, 0.19,
#'                            0.92,  0.21, 0.15,
#'                            1.26,  0.54, 0.05,
#'                            0.61, -0.03, 0.18,
#'                            1.74,  0.01, 0.12,
#'                            0.50,  1.96, 0.12,
#'                            0.96,  0.04, 0.13,
#'                            0.59, -0.09, 0.18,
#'                            0.82, -1.16, 0.17,
#'                            1.26,  0.02, 0.11,
#'                            0.82,  0.20, 0.07,
#'                            0.75, -0.43, 0.15,
#'                            1.49, -0.06, 0.09,
#'                            0.97, -0.34, 0.12,
#'                            1.49,  0.05, 0.12,
#'                            0.89, -0.25, 0.15,
#'                            1.45,  0.06, 0.07,
#'                            0.75,  0.31, 0.18,
#'                            1.43,  0.04, 0.08,
#'                            0.60,  0.13, 0.22,
#'                            0.83,  0.52, 0.09,
#'                            0.56, -0.96, 0.19,
#'                            0.67, -0.79, 0.20,
#'                            0.70,  0.37, 0.18,
#'                            1.03, -0.71, 0.14,
#'                            0.89, -0.19, 0.21,
#'                            1.23,  0.74, 0.06,
#'                            0.90, -0.44, 0.18,
#'                            1.23, -0.17, 0.12,
#'                            0.69, 0.53, 0.17),
#'                            byrow = TRUE, ncol = 3)
#'
#'
#' group1Pars <- group2Pars <- group3Pars <- group4Pars <- group5Pars <- baseParameters
#'
#' group1Pars[c(4,  22), 1] <- group1Pars[c(4,  22), 1] + .4
#' group1Pars[c(7,  22), 2] <- group1Pars[c(7,  22), 2] + .7
#'
#' group2Pars[c(13, 28), 1] <- group2Pars[c(13, 28), 1] - .4
#' group2Pars[c(19, 28), 2] <- group2Pars[c(19, 28), 2] - .7
#'
#' group3Pars[c(4,  13), 1] <- group3Pars[c(4,  13), 1] + .4
#' group3Pars[c(4,  22), 2] <- group3Pars[c(4,  22), 2] - .7
#'
#' group4Pars[c(7,  28), 1] <- group4Pars[c(7,  28), 1] - .4
#' group4Pars[c(7,  19), 2] <- group4Pars[c(7,  19), 2] + .7
#'
#' itemParameters <- list(base = baseParameters,
#'                        group1 = group1Pars,
#'                        group2 = group2Pars,
#'                        group3 = group3Pars,
#'                        group4 = group4Pars,
#'                        group5 = group5Pars
#'                        )
#'
#' itemCovariances <- lapply(itemParameters, AseIrt, irtModel = "3pl", sampleSize = 5000)
#'
#' mgIpr <- IprMG(itemParameters = itemParameters,
#'                itemCovariances = itemCovariances,
#'                nReplicates = 100)
#' mgIpr <- Bound3PlIpr(mgIpr)
#'
#' mgncdifIpr <- IprMGNcdif(itemParameterList = mgIpr, irtModel = "3pl")
#'
#' @author Victor H. Cervantes <vhcervantesb at unal.edu.co>
#'
#' @export
IprMGNcdif <- function (itemParameterList, irtModel = "2pl", focalAbilities = NULL, focalDistribution = "norm",
                        relSizes = NULL,
                        subdivisions = 5000, logistic = TRUE, focalDistrExtra = list(mean = 0, sd = 1)) {

  mgncdif <- sapply(itemParameterList, function (x) MGNcdif(x, irtModel = irtModel, focalAbilities = focalAbilities,
                                                            focalDistribution = focalDistribution,
                                                            focalDistrExtra = focalDistrExtra,
                                                            subdivisions = subdivisions, logistic = logistic))

  if (nrow(itemParameterList[[1]][['base']]) == 1) {
    mgncdif <- matrix(mgncdif, nrow = 1)
  }

  return(mgncdif)
}