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
#' @export
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
MGNcdif <- function (itemParameters, irtModel = "2pl", focalAbilities = NULL, focalDistribution = "norm", relSizes = NULL,
                   subdivisions = 5000, logistic = TRUE, focalDistrExtra = list(list(mean = 0, sd = 1))) {

    if (!("base" %in% names(itemParameters))) {
        stop("To compute MGNcdif a base group must be define. The base group is identified as 'base' among itemParameters. See ?MGNcdif")
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
        names(relSizes) <- names(focalAbilities)

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