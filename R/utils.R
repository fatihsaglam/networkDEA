#' @importFrom lpSolve lp
#' @importFrom stats sd
NULL

# Solve a linear program with lpSolve.
# lpSolve restricts all variables to be non-negative, so variables listed in
# `free` are split into a positive and a negative part and merged afterwards.
.solve_lp <- function(direction, obj, A, dir, rhs, free = integer(0)) {
  nv <- length(obj)
  A <- as.matrix(A)
  if (length(free)) {
    A2 <- cbind(A, -A[, free, drop = FALSE])
    obj2 <- c(obj, -obj[free])
  } else {
    A2 <- A
    obj2 <- obj
  }
  res <- lpSolve::lp(direction, obj2, A2, dir, rhs)
  if (res$status != 0) {
    return(list(status = res$status, objval = NA_real_,
                solution = rep(NA_real_, nv)))
  }
  sol <- res$solution
  x <- sol[seq_len(nv)]
  if (length(free)) x[free] <- x[free] - sol[nv + seq_along(free)]
  list(status = 0L, objval = res$objval, solution = x)
}

.prep_matrix <- function(x, name) {
  if (is.null(x)) stop("`", name, "` must not be NULL.", call. = FALSE)
  m <- as.matrix(x)
  if (!is.numeric(m)) stop("`", name, "` must be numeric.", call. = FALSE)
  storage.mode(m) <- "double"
  if (anyNA(m) || any(!is.finite(m))) {
    stop("`", name, "` contains missing or non-finite values.", call. = FALSE)
  }
  if (any(m < 0)) {
    stop("`", name, "` contains negative values, which are not supported.",
         call. = FALSE)
  }
  if (is.null(colnames(m))) colnames(m) <- paste0(name, seq_len(ncol(m)))
  m
}

# Column scale factors. Efficiency scores are invariant to rescaling of the
# variables, so every variable is divided by its mean before solving to keep
# the linear programs numerically well conditioned.
.scale_factors <- function(m) {
  s <- colMeans(m)
  s[!is.finite(s) | s <= 0] <- 1
  s
}

.rts <- function(rts) match.arg(rts, c("crs", "vrs"))
.orient <- function(orientation) match.arg(orientation, c("input", "output"))

.check_status <- function(r, o, what) {
  if (r$status != 0) {
    warning(sprintf("%s: the linear program for DMU %s did not solve (lpSolve status %d).",
                    what, o, r$status), call. = FALSE)
  }
  invisible(NULL)
}

#' Simulate Data for a Two-Stage Process
#'
#' Generates a small artificial data set with inputs, intermediate products and
#' final outputs, useful for examples and tests. Stage 1 turns inputs into
#' intermediates and stage 2 turns intermediates into final outputs, each with
#' random inefficiency.
#'
#' @param n Number of decision making units.
#' @param m,k,s Number of inputs, intermediates and final outputs.
#' @param seed Optional random seed.
#'
#' @return A list with matrices `x`, `z` and `y`.
#' @examples
#' d <- ndea_simulate(n = 20, seed = 1)
#' str(d)
#' @export
ndea_simulate <- function(n = 30, m = 2, k = 2, s = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- matrix(stats::runif(n * m, 10, 100), n, m)
  w1 <- matrix(stats::runif(m * k, 0.2, 1), m, k)
  z <- (x %*% w1) * exp(-abs(stats::rnorm(n * k, 0, 0.25)))
  w2 <- matrix(stats::runif(k * s, 0.2, 1), k, s)
  y <- (z %*% w2) * exp(-abs(stats::rnorm(n * s, 0, 0.25)))
  colnames(x) <- paste0("x", seq_len(m))
  colnames(z) <- paste0("z", seq_len(k))
  colnames(y) <- paste0("y", seq_len(s))
  rownames(x) <- rownames(z) <- rownames(y) <- paste0("DMU", seq_len(n))
  list(x = x, z = z, y = y)
}
