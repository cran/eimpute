#' @title Search rank magnitude of the best approximating matrix
#'
#' @description Estimate a preferable matrix rank magnitude for fitting a low-rank matrix approximation to a matrix with missing values.
#' The algorithm use GIC/CV to search the rank in a given range, and then fill the missing values with the estimated rank.
#'
#' @aliases r.search
#'
#' @param r.min the start rank for searching. Default \code{r.min = 1}.
#' @param r.max the max rank for searching.
#' @param rule.type a character string indicating the information criterion rule.
#' If \code{rule.type = "gic"}, generalized information criterion rule is used,
#' else if \code{rule.type = "cv"}, cross validation is used.
#' Any unambiguous substring can be given. Default \code{rule.type = "gic"}.
#' @param noise.var the variance of noise.
#' @param init if init = FALSE(the default), the missing entries will initialize with mean.
#' @param init.mat the initialization matrix.
#' @param maxit.rank maximal number of iterations in searching rank. Default \code{maxit.rank = 1}.
#' @param nfolds number of folds in cross validation. Default \code{nfolds = 5}.
#' @param thresh convergence threshold, measured as the relative change in the Frobenius norm between two successive estimates.
#' @param maxit	 maximal number of iterations.
#' @param override logical value indicating whether the observed elements in \code{x} should be overwritten by its low-rank approximation.
#' @param control a list of parameters that control details of standard procedure, See \link{biscale.control}.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#'
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

r.search <- function(x, r.min = 1, r.max = "auto",
                     svd.method = c("tsvd", "rsvd"), rule.type = c("gic", "cv"), noise.var = 0,
                     init = FALSE, init.mat = 0,
                     maxit.rank = 1, nfolds = 5, thresh = 1e-05, maxit = 100, override = FALSE, control = list(...), ...){

  m <- nrow(x)
  n <- ncol(x)

  x_sd <- biscale(x, control = control)

  x_train <- x_sd[[1]]
  ind_ob <- which(!is.na(x_train), arr.ind = TRUE)
  x_ob <- x_train[ind_ob]
  ind <- ind_ob - 1

  ind_miss <- which(is.na(x_train), arr.ind = TRUE)

  if(init == TRUE){
    ind_miss <- cbind(ind_miss, as.vector(init.mat[ind_miss]))
  }else{
    ind_miss <- cbind(ind_miss, rnorm(nrow(ind_miss), 0, noise.var))
  }
  ind_miss[,1:2] <- ind_miss[,1:2] - 1
  if(r.max == "auto"){
    r.max = floor((m + n + sqrt((n + m)^2 - 4 * nrow(ind_ob)))/2)
  }

  svdm <- match.arg(svd.method)
  type <- ifelse(svdm == "tsvd", 1, 2)

  rule <- match.arg(rule.type)
  rulet <- ifelse(rule == "cv", "cv", "gic")

  res <- vector()
  res2 <- vector()
  rank.seq <- r.min : r.max

  if(rulet == "cv"){
    cvind <- sample(1 : length(x_ob), length(x_ob))
    res <- cv_rank(ind[cvind,], ind_miss, x_ob[cvind], m, n, r.min, r.max, nfolds, maxit.rank, thresh, type, init)
  }else{
    Z_temp <- ic_rank(ind, ind_miss, x_ob, m, n, r.min, r.max, maxit.rank, thresh, type, init)
    m1 <- max(m, n)
    m2 <- min(m, n)
    res <- m2*log(Z_temp) + rank.seq*log(log(m2)) * log(m1)
    res2 <- log(Z_temp)
  }

  r.est <- rank.seq[which.min(res)]

  # Z_temp <- kkt_fix(ind, ind_miss, x_ob, m, n, r.est, maxit, thresh, type, init)
  # Z.fit <- Z_temp[[1]] * (x_sd[[4]] %*% t(x_sd[[5]])) + matrix(rep(x_sd[[2]], n), nrow = m) + t(matrix(rep(x_sd[[3]], m), nrow = n))
  # if (!override) {
  #   Z.fit[ind_ob] <- x_ob
  # }

  list(r.est = r.est, gicseq = res, errseq = res2)
}




