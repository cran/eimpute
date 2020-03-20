#' @title Search rank magnitude of the best approximating matrix
#'
#' @description Estimate a preferable matrix rank magnitude for fitting a low-rank matrix approximation to a matrix with missing values.
#' The algorithm use BIC/GIC rule to search the rank in a given range, and then fill the missing values with the estimated rank.
#'
#' @aliases r.search
#'
#' @param r.min the start rank for searching. Default \code{r.min = 1}.
#' @param r.max the max rank for searching.
#' @param rule.type a character string indicating the information criterion rule.
#' If \code{rule.type = "gic"}, generalized information criterion rule is used,
#' else if \code{rule.type = "bic"}, bayesian information criterion rule is used.
#' Any unambiguous substring can be given. Default \code{rule.type = "gic"}.
#' @param maxit.rank maximal number of iterations in searching rank. Default \code{maxit.rank = 1}.
#' @inheritParams eimpute
#'
#' @return A list containing the following components
#' \item{\code{x.imp}}{the matrix after completion with the estimated rank.}
#' \item{\code{r.est}}{the rank estimation.}
#' \item{\code{rmse}}{the relative mean square error of matrix completion, i.e., training error.}
#' \item{\code{iter.count}}{the number of iterations.}
#'
#' @rdname r.search
#' @importFrom Rcpp sourceCpp
#' @export
#' @examples
#' ################# Quick Start #################
#' m <- 100
#' n <- 100
#' r <- 10
#' x_na <- incomplete.generator(m, n, r)
#' head(x_na[, 1:6])
#' x_impute <- r.search(x_na, 1, 15, "rsvd", "gic")
#' x_impute[["r.est"]]
#' head(x_impute[["x.imp"]][, 1:6])
#' x_impute[["rmse"]]
r.search <- function(x, r.min = 1, r.max, svd.method = c("tsvd", "rsvd"), rule.type = c("gic", "bic"),
                        thresh = 1e-05, maxit.rank = 1, maxit = 100,override = FALSE, control = list(...), ...){

  m <- nrow(x)
  n <- ncol(x)

  x_sd <- biscale(x, control = control)

  x_train <- x_sd[[1]]
  ind_ob <- which(!is.na(x_train), arr.ind = TRUE)
  x_ob <- x_train[ind_ob]
  ind <- ind_ob - 1

  svdm <- match.arg(svd.method)
  type <- ifelse(svdm == "tsvd", 1, 2)

  rule <- match.arg(rule.type)
  rulet <- ifelse(rule == "bic", "bic", "gic")

  res <- vector()
  rank.seq <- r.min : r.max
  if(rule == "gic"){
    for(i in 1 : length(rank.seq)){
      r = rank.seq[i]
      Z_temp <- kkt_fix(ind, x_ob, m, n, r, maxit.rank, thresh, type)
      Z.fit <- Z_temp[[1]] * (x_sd[[4]] %*% t(x_sd[[5]])) + matrix(rep(x_sd[[2]], n), nrow = m) + t(matrix(rep(x_sd[[3]], m), nrow = n))

      train_error <- sum((x_ob - Z.fit[ind_ob])^2)
      res[i] <- m*log(train_error) + r*log(log(m)) * log(n)
    }
  }else{
    for(i in 1 : length(rank.seq)){
      r = rank.seq[i]
      Z_temp <- kkt_fix(ind, x_ob, m, n, r, maxit.rank, thresh, type)
      Z.fit <- Z_temp[[1]] * (x_sd[[4]] %*% t(x_sd[[5]])) + matrix(rep(x_sd[[2]], n), nrow = m) + t(matrix(rep(x_sd[[3]], m), nrow = n))

      train_error <- sum((x_ob - Z.fit[ind_ob])^2)
      res[i] <- 2*log(train_error) + r*log(m)
    }
  }

  r.est <- rank.seq[which.min(res)]
  if (!override) {
    Z.fit[ind_ob] <- x_ob
  }
  Z_temp <- kkt_fix(ind, x_ob, m, n, r.est, maxit, thresh, type)
  Z.fit <- Z_temp[[1]] * (x_sd[[4]] %*% t(x_sd[[5]])) + matrix(rep(x_sd[[2]], n), nrow = m) + t(matrix(rep(x_sd[[3]], m), nrow = n))

  train_error <- sum((x_ob - Z.fit[ind_ob])^2) / sum(x_ob^2)
  list(x.imp = Z.fit, r.est = r.est, rmse = train_error, iter.count = Z_temp[[2]])
}
