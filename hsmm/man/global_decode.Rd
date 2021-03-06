% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global_decode.R
\name{global_decode}
\alias{global_decode}
\title{Global decoding of HMM}
\usage{
global_decode(object, ddl = NULL, state.names = NULL)
}
\arguments{
\item{object}{fitted crm model (must be an HMM model)}

\item{ddl}{design data list}

\item{state.names}{names for states used to label output; if NULL uses strata.labels + Dead state}
}
\value{
matrix of state predictions
}
\description{
Computes sequence of state predictions for each individual
}
\references{
Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. See page 82.
}
\author{
Jeff Laake
}
\keyword{utility}
