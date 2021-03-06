% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NCDIF.R
\name{CalculateGrmExp}
\alias{CalculateGrmExp}
\title{Calculates the expected item score under the GRM model.}
\usage{
CalculateGrmExp(thetaValue, itemParameters, logistic = TRUE)
}
\arguments{
\item{thetaValue}{A numeric value or array for the theta (ability) value(s) where the difference will be calculated}

\item{itemParameters}{A matrix containing the numeric values of item discriminations on the first column and category thresholds on the rest columns where the (column position - 1) indicates the category score or weight.}

\item{logistic}{A logical value stating if the IRT model will use the logistic or the normal metric.}
}
\value{
expectedScore A numeric matrix with the expected score on each thetaValue for each item.
}
\description{
Calculates the expected item score under the GRM model.
}
\references{
de Ayala, R. J., (2009). The theory and practice of item response theory. New York: The Guildford Press

Oshima, T. & Morris, S. (2008). Raju's Differential Functioning of Items and Tests (DFIT). Educational Measurement: Issues and Practice, 27(3), 43--50. doi:10.1111/j.1745-3992.2008.00127.x
}
\author{
Victor H. Cervantes <vhcervantesb at unal.edu.co>
}
