% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IPR.R
\name{Ipr}
\alias{Ipr}
\title{Item parameter replication}
\usage{
Ipr(itemParameters, itemCovariances, nReplicates = 5000)
}
\arguments{
\item{itemParameters}{A list containing "focal" and "reference" item parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should me a matrix with nrow equal to the number of items.}

\item{itemCovariances}{A list containing "focal" and "reference" matrices of covariance for item estimates. Each (focal and reference) may be either a list of covariance matrices for each item or a single matrix of covariance of all parameters.}

\item{nReplicates}{A numeric value indicating the number of replications to perform}
}
\value{
itemParameters A list with item parameters for focal and reference groups
}
\description{
Generates a sample of item parameters assuming multivariate normality of estimates
}
\examples{
# # Not run
# #
# # data(dichotomousItemParameters)
# # threePlParameters <- dichotomousItemParameters
# # isNot3Pl          <- ((dichotomousItemParameters[['focal']][, 3] == 0) |
# #                       (dichotomousItemParameters[['reference']][, 3] == 0))
# #
# # threePlParameters[['focal']]          <- threePlParameters[['focal']][!isNot3Pl, ]
# # threePlParameters[['reference']]      <- threePlParameters[['reference']][!isNot3Pl, ]
# # threePlParameters[['focal']][, 3]     <- threePlParameters[['focal']][, 3] + 0.1
# # threePlParameters[['reference']][, 3] <- threePlParameters[['reference']][, 3] + 0.1
# # threePlParameters[['focal']][, 2]     <- threePlParameters[['focal']][, 2] + 1.5
# # threePlParameters[['reference']][, 2] <- threePlParameters[['reference']][, 2] + 1.5
# # threePlParameters[['focal']]          <- threePlParameters[['focal']][-c(12, 16, 28), ]
# # threePlParameters[['reference']]      <- threePlParameters[['reference']][-c(12, 16, 28), ]
# #
# # threePlAse <- list()
# # threePlAse[["focal"]]     <- AseIrt(itemParameters = threePlParameters[["focal"]],
# #                                     logistic = TRUE,
# #                                     sampleSize = 10000,
# #                                     irtModel = "3pl")
# # threePlAse[["reference"]] <- AseIrt(itemParameters = threePlParameters[["reference"]],
# #                                     logistic = TRUE,
# #                                     sampleSize = 15000,
# #                                     irtModel = "3pl")
# #
# # set.seed(41568)
# # threePlIpr <- Ipr(itemParameters = threePlParameters, itemCovariances = threePlAse,
# #                   nReplicates = 100)

}
\references{
Oshima, T., Raju, N. & Nanda, A. (2006). A new method for assessing the statistical significance in the Differential Functioning of Items and Tests (DFIT) framework. Journal of educational measurement, 43(1), 1--17. doi:10.1111/j.1745-3984.2006.00001.x
}
\author{
Victor H. Cervantes <vhcervantesb at unal.edu.co>
}
