% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/horseshoe.R
\name{brms.horseshoe}
\alias{brms.horseshoe}
\title{Title Define a model formula with horseshoe prior}
\usage{
brms.horseshoe(shrinkage.term, nonshrinkage.term, ...)
}
\arguments{
\item{shrinkage.term}{term with horseshoe prior.}

\item{nonshrinkage.term}{terms with horseshoe prior.}

\item{...}{arguments in horseshoe() see ?horseshoe}
}
\value{
a list containing formula used for brms and prior.
}
\description{
Title Define a model formula with horseshoe prior
}
