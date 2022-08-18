#' Build the Variable Names or Data Objects for Estimation
#' Modification from builder.r in brmsmargins
#' For horseshoe prior
#'
#' @description
#' These are a set of internal utility functions.
#' They are not intended for general use.
#'
#' @details
#' \itemize{
#'   \item{\code{.namesSD_hs}}{Create the names of random effect standard deviation estimates.}
#'   \item{\code{.buildSD_hs}}{Return matrix of random effect standard deviation estimates. Rows are posterior draws.}
#'   \item{\code{.namesZ_hs}}{Create the names of random effects data for predictions.}
#'   \item{\code{.buildZ_hs}}{Return matrix of data for random effect predictions.}
#' }
#'
#' @param data A data object. For example the result of [make_standata()]
#'   for [.buildZ()], which is a list,
#'   or a dataset of the posterior draws such as from [as_draws_df()]
#'   for [.buildL()] and [.buildSD()].
#' @param ranef A data set with information about the model object random effects.
#'   Only used for \code{.namesSD} and \code{.buildSD}.
#' @param block Which random effect block to use. An integer.
#' @param number The number of elements in that random effect block. An integer.
#' @param dpar Which dpar to use. Does not apply to the L matrix.
#' @return A character vector for all \code{.names} functions or a matrix
#'   for all \code{.build} functions.
#' @keywords internal
#' @name builders
NULL



## make Rcmd check happy
utils::globalVariables(c("group", "coef", "id"))

#' @rdname builders
.namesSD_hs <- function(ranef, block) {
  stopifnot(data.table::is.data.table(ranef))
  n <- ranef[id == block]

  n[, sprintf("sd_%s__%s_%s", group, nlpar, coef)]
}

#' @rdname builders
.buildSD_hs <- function(data, ranef, block) {
  stopifnot(data.table::is.data.table(data))
  n <- .namesSD_hs(ranef, block)
  as.matrix(data[, ..n])
}

## make Rcmd check happy
utils::globalVariables(c("Number"))

#' @rdname builders
.namesZ_hs <- function(block, number, nlpar) {
  n <- expand.grid(Block = block,
                   Number = seq_len(number))
  n <- as.data.table(n)
  n[, sprintf("Z_%d_%s_%d", Block, nlpar, Number)]
}

#' @rdname builders
.buildZ_hs <- function(data, block, number, nlpar) {
  n <- .namesZ_hs(block, number, nlpar)
  as.matrix(do.call(cbind, data[n]))
}

