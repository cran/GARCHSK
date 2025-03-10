% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_GJRSK.R
\name{gjrsk_construct}
\alias{gjrsk_construct}
\title{This function constructs GJRSK model of given data and parameters.}
\usage{
gjrsk_construct(params, data)
}
\arguments{
\item{params}{vector of GJRSK model parameters(p1,const2,p2,q2,r2,const3,p3,q3,r3,const4,p4,q4,r4)}

\item{data}{vector time series data}
}
\value{
list of conditional mean(mu), variance(h), skewness(sk) and kurtosis(ku)
}
\description{
This function constructs GJRSK model of given data and parameters.
}
