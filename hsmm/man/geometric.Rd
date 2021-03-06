% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dwell_time_distributions.R
\name{geometric}
\alias{geometric}
\alias{shifted_poisson}
\alias{shifted_negbinomial}
\alias{shifted_binomial}
\alias{unstructured_gt}
\alias{unstructured}
\title{Pre-defined dwell time distributions for HSMM models}
\usage{
geometric(theta, n)
}
\arguments{
\item{theta}{parameter vector for dwell time distributions}

\item{n}{number of discrete states for dwell time}
}
\value{
vector of probabilities for dwell times
}
\description{
Various discrete dwell time distributions defined in Zucchini et al second edition.
Note: if you use unstructured without a geometric tail (gt), the number of dwell times must exceed the longest 
observed sequence in the state or -lnl could become undefined (Nan).
}
\references{
Zucchini, W., MacDonald, I.L. and Langrock, R. 2016. Hidden Markov Models for Time Series: 
An introduction using R, 2nd ed. CRC Press.
}
\author{
Jeff Laake
}
\keyword{utility}
