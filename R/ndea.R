#' Generalized Network DEA
#'
#' Efficiency of DMUs whose internal structure is an arbitrary network of
#' nodes (serial, parallel or mixed) linked by intermediate products, following
#' the generalized network DEA model of the MaxDEA Team (2026).
#'
#' The computation has two steps for every DMU. First each node is evaluated on
#' its own, treating its intermediate inputs as inputs and its intermediate
#' outputs as outputs, giving the stand-alone node efficiencies
#' \eqn{E^{(t)*}}. Second, a network model couples the nodes through a single
#' conservation constraint on the total value of the intermediates and
#' measures the *synergy loss* `D*` (zero or negative):
#'
#' * `D* = 0` (**lossless synergy**): the stand-alone node efficiencies can all
#'   be attained simultaneously. The node efficiencies are unique and the system
#'   efficiency is \eqn{E^* = \prod_t (E^{(t)*})^{\omega_t}}.
#' * `D* < 0` (**lossy synergy**): the node efficiencies of the network model are
#'   smaller than the stand-alone ones and are generally **not unique**; the
#'   values reported come from one optimal solution of the linear program. A
#'   typical cause is a zero weight on an intermediate product.
#'
#' The weights \eqn{\omega_t} express the role of each node: all 1 (default)
#' gives the product of node efficiencies, suited to serial structures;
#' \eqn{\omega_t = 1/N} gives the geometric mean, suited to parallel
#' structures. They must satisfy \eqn{1 \le \sum_t \omega_t \le N}.
#'
#' With two nodes (a stage producing intermediates and a stage consuming them)
#' this model gives the same results as [ndea_two_stage()] with
#' `synergy = "flexible"`.
#'
#' Not (yet) provided: reference sets (lambdas) and projections for the
#' network model, and models based on full-dimensional efficient facets
#' (FDEF). Data are internally rescaled by column means; scores do not depend on
#' measurement units.
#'
#' @param data Data frame with one row per DMU and one column per variable.
#' @param network A network created with [network()].
#' @param rts Returns to scale: `"crs"` or `"vrs"`.
#' @param orientation `"input"` or `"output"`.
#' @param omega Numeric vector of node weights \eqn{\omega_t} (one per node,
#'   in the order of the nodes). Default: all 1.
#' @param dmu Optional DMU names. Defaults to the row names of `data`.
#' @param tol Tolerance for declaring `D* = 0`.
#'
#' @return An object of class `"ndea"`, a list with
#'   * `efficiency`: data frame with `system`, one column per node, `D_star`
#'     and `synergy` (`"lossless"` / `"lossy"`);
#'   * `standalone`: matrix of stand-alone node efficiencies;
#'   * `weights`: list (one entry per DMU) of lists (one per node) with the
#'     multipliers (shadow prices) in original units;
#'   * `status`: lpSolve status codes.
#'
#' @references
#' MaxDEA Team (2026). Innovative Two-Stage and Generalized Network DEA.
#' MaxDEA Software Ltd. (Not peer reviewed.)
#'
#' @seealso [node()], [network()], [ndea_two_stage()]
#' @examples
#' d <- ndea_simulate(n = 15, m = 2, k = 1, s = 1, seed = 3)
#' dat <- data.frame(d$x, d$z, d$y)
#' net <- network(
#'   node("stage1", inputs = c("x1", "x2"), z_out = "z1"),
#'   node("stage2", z_in = "z1", outputs = "y1")
#' )
#' res <- ndea(dat, net, rts = "vrs", orientation = "input")
#' res
#' @export
ndea <- function(data, network, rts = c("crs", "vrs"),
                 orientation = c("input", "output"), omega = NULL,
                 dmu = NULL, tol = 1e-6) {
  rts <- .rts(rts)
  orientation <- .orient(orientation)
  if (!inherits(network, "ndea_network")) {
    stop("`network` must be created with `network()`.", call. = FALSE)
  }
  nodes <- network$nodes
  N <- length(nodes)
  data <- as.data.frame(data)
  vars <- unique(unlist(lapply(nodes, function(nd) c(nd$inputs, nd$outputs, nd$z_in, nd$z_out))))
  miss <- setdiff(vars, names(data))
  if (length(miss)) {
    stop("Variable(s) not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  n <- nrow(data)
  if (n < 2) stop("At least two DMUs are required.", call. = FALSE)
  if (is.null(dmu)) {
    dmu <- rownames(data)
    if (is.null(dmu)) dmu <- as.character(seq_len(n))
  }
  if (length(dmu) != n) stop("`dmu` must have one name per row.", call. = FALSE)
  if (is.null(omega)) omega <- rep(1, N)
  if (length(omega) != N || any(!is.finite(omega)) || any(omega <= 0)) {
    stop("`omega` must be a positive vector with one entry per node.", call. = FALSE)
  }
  if (sum(omega) < 1 - 1e-8 || sum(omega) > N + 1e-8) {
    stop("The node weights must satisfy 1 <= sum(omega) <= number of nodes.", call. = FALSE)
  }

  # per-variable scale factors (shared by all nodes that use a variable)
  sc <- vapply(vars, function(v) .scale_factors(.prep_matrix(data[v], v)), numeric(1))
  getm <- function(cols) {
    if (!length(cols)) return(matrix(0, n, 0))
    m <- .prep_matrix(data[cols], "v")
    colnames(m) <- cols
    sweep(m, 2, sc[cols], "/")
  }
  nd_data <- lapply(nodes, function(nd) {
    list(name = nd$name, X = getm(nd$inputs), Zin = getm(nd$z_in),
         Y = getm(nd$outputs), Zout = getm(nd$z_out))
  })

  node_names <- names(nodes)
  E_sa <- matrix(NA_real_, n, N, dimnames = list(dmu, node_names))
  E_net <- E_sa
  D_star <- rep(NA_real_, n)
  synergy <- rep(NA_character_, n)
  status <- integer(n)
  wts <- vector("list", n)
  names(wts) <- dmu

  for (o in seq_len(n)) {
    sa <- vapply(nd_data, function(d) .node_standalone(d, o, rts, orientation), numeric(1))
    E_sa[o, ] <- sa
    if (anyNA(sa)) { status[o] <- 99L; next }
    net <- .net_solve(nd_data, o, rts, orientation, sa)
    status[o] <- net$status
    if (net$status != 0) {
      warning(sprintf("ndea: the network linear program for DMU %s did not solve (lpSolve status %d).",
                      dmu[o], net$status), call. = FALSE)
      next
    }
    D_star[o] <- net$D
    lossless <- abs(net$D) <= tol
    synergy[o] <- if (lossless) "lossless" else "lossy"
    E_net[o, ] <- if (lossless) sa else net$E_node
    # weights back to original units
    w <- net$weights
    for (t in seq_len(N)) {
      nd <- nodes[[t]]
      w[[t]]$v <- w[[t]]$v / sc[nd$inputs]
      w[[t]]$w_in <- w[[t]]$w_in / sc[nd$z_in]
      w[[t]]$u <- w[[t]]$u / sc[nd$outputs]
      w[[t]]$w_out <- w[[t]]$w_out / sc[nd$z_out]
    }
    wts[[o]] <- w
  }

  system <- apply(E_net, 1, function(e) prod(e^omega))
  eff <- data.frame(dmu = dmu, system = unname(system), E_net,
                    D_star = D_star, synergy = synergy,
                    stringsAsFactors = FALSE, check.names = FALSE)
  rownames(eff) <- NULL
  out <- list(efficiency = eff, standalone = E_sa, weights = wts,
              status = status, call = match.call(),
              settings = list(rts = rts, orientation = orientation,
                              omega = stats::setNames(omega, node_names),
                              nodes = node_names, tol = tol),
              network = network)
  class(out) <- "ndea"
  out
}

# Coefficient matrices of a node ---------------------------------------------
# Weight vector layout per node: v (inputs), w_in (z_in), u (outputs), w_out
# (z_out). den_j = v X_j + w_in Zin_j ; num_j = u Y_j + w_out Zout_j.
.node_mats <- function(d, o) {
  n <- nrow(d$X); m <- ncol(d$X); kin <- ncol(d$Zin)
  s <- ncol(d$Y); kout <- ncol(d$Zout)
  den <- cbind(d$X, d$Zin, matrix(0, n, s + kout))
  num <- cbind(matrix(0, n, m + kin), d$Y, d$Zout)
  list(den = den, num = num, den0 = den[o, ], num0 = num[o, ],
       sizes = c(m = m, kin = kin, s = s, kout = kout))
}

# Stand-alone efficiency of one node (MaxDEA manual 9.43). Input orientation
# E = max (num0 - u0)/den0 ; output orientation E = max num0/(den0 + u0), subject
# to den_j - num_j + u0 >= 0; u0 = 0 under constant returns to scale.
.node_standalone <- function(d, o, rts, orientation) {
  vrs <- rts == "vrs"
  M <- .node_mats(d, o)
  n <- nrow(M$den); nb <- ncol(M$den); nv <- nb + vrs
  A <- cbind(M$den - M$num, if (vrs) 1)
  dir <- rep(">=", n); rhs <- rep(0, n)
  if (orientation == "input") {
    obj <- c(M$num0, if (vrs) -1)
    nrm <- c(M$den0, if (vrs) 0)
    sense <- "max"
  } else {
    obj <- c(M$den0, if (vrs) 1)
    nrm <- c(M$num0, if (vrs) 0)
    sense <- "min"
  }
  A <- rbind(A, nrm); dir <- c(dir, "="); rhs <- c(rhs, 1)
  r <- .solve_lp(sense, obj, A, dir, rhs, free = if (vrs) nv else integer(0))
  if (r$status != 0) return(NA_real_)
  if (orientation == "input") r$objval else 1 / r$objval
}

# Network model (MaxDEA manual 9.44 / 9.45 with theta = stand-alone efficiency).
# Returns D* (<= 0) and, for lossy synergy, the node efficiencies implied by the
# optimal weights.
.net_solve <- function(nd_data, o, rts, orientation, E_sa) {
  vrs <- rts == "vrs"
  N <- length(nd_data)
  mats <- lapply(nd_data, .node_mats, o = o)
  sizes <- vapply(mats, function(M) ncol(M$den), numeric(1))
  nvt <- sizes + vrs
  offs <- c(0, cumsum(nvt))
  nv <- offs[N + 1]
  n <- nrow(mats[[1]]$den)
  obj <- numeric(nv)
  cons <- numeric(nv)
  Ablocks <- list(); dirs <- character(); rhs <- numeric()
  u0cols <- integer(0)
  for (t in seq_len(N)) {
    M <- mats[[t]]; d <- nd_data[[t]]
    cols <- offs[t] + seq_len(sizes[t])
    Ab <- matrix(0, n, nv)
    Ab[, cols] <- M$den - M$num
    nrm <- numeric(nv)
    if (orientation == "input") {
      obj[cols] <- M$num0 - E_sa[t] * M$den0
      nrm[cols] <- M$den0
    } else {
      obj[cols] <- M$num0 / E_sa[t] - M$den0
      nrm[cols] <- M$num0
    }
    m <- M$sizes[["m"]]; kin <- M$sizes[["kin"]]; s <- M$sizes[["s"]]; kout <- M$sizes[["kout"]]
    cons[offs[t] + m + seq_len(kin)] <- d$Zin[o, ]
    cons[offs[t] + m + kin + s + seq_len(kout)] <- -d$Zout[o, ]
    if (vrs) {
      u0 <- offs[t] + sizes[t] + 1
      u0cols <- c(u0cols, u0)
      Ab[, u0] <- 1
      obj[u0] <- -1
      # the intercept enters the intermediate conservation constraint of nodes
      # whose intercept sits on the intermediate side (no y under input
      # orientation; no x under output orientation)
      eta <- if (orientation == "input") s == 0 else m == 0
      if (eta) cons[u0] <- 1
    }
    Ablocks <- c(Ablocks, list(Ab, nrm))
    dirs <- c(dirs, rep(">=", n), ">=")
    rhs <- c(rhs, rep(0, n), 1)
  }
  A <- do.call(rbind, c(Ablocks, list(cons)))
  dirs <- c(dirs, "="); rhs <- c(rhs, 0)
  r <- .solve_lp("max", obj, A, dirs, rhs, free = u0cols)
  if (r$status != 0) return(list(status = r$status))
  sol <- r$solution
  E_node <- numeric(N)
  weights <- vector("list", N)
  for (t in seq_len(N)) {
    M <- mats[[t]]
    cols <- offs[t] + seq_len(sizes[t])
    wv <- sol[cols]
    u0v <- if (vrs) sol[offs[t] + sizes[t] + 1] else 0
    den0 <- sum(M$den0 * wv); num0 <- sum(M$num0 * wv)
    E_node[t] <- if (orientation == "input") (num0 - u0v) / den0 else num0 / (den0 + u0v)
    m <- M$sizes[["m"]]; kin <- M$sizes[["kin"]]; s <- M$sizes[["s"]]; kout <- M$sizes[["kout"]]
    weights[[t]] <- list(
      v = stats::setNames(wv[seq_len(m)], colnames(nd_data[[t]]$X)),
      w_in = stats::setNames(wv[m + seq_len(kin)], colnames(nd_data[[t]]$Zin)),
      u = stats::setNames(wv[m + kin + seq_len(s)], colnames(nd_data[[t]]$Y)),
      w_out = stats::setNames(wv[m + kin + s + seq_len(kout)], colnames(nd_data[[t]]$Zout)),
      u0 = u0v)
  }
  names(weights) <- vapply(nd_data, function(d) d$name, character(1))
  list(status = 0L, D = r$objval, E_node = E_node, weights = weights)
}

#' @export
print.ndea <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("Generalized network DEA (", toupper(s$rts), ", ", s$orientation,
      "-oriented), ", length(s$nodes), " nodes: ", paste(s$nodes, collapse = ", "),
      "\n", sep = "")
  cat("Node weights (omega):", paste(format(s$omega, digits = 3), collapse = " "), "\n")
  e <- x$efficiency
  cat("Synergy: ", sum(e$synergy == "lossless", na.rm = TRUE), " lossless, ",
      sum(e$synergy == "lossy", na.rm = TRUE), " lossy DMUs\n\n", sep = "")
  print(format(e, digits = digits), row.names = FALSE)
  invisible(x)
}

#' @export
summary.ndea <- function(object, ...) {
  e <- object$efficiency
  cols <- c("system", object$settings$nodes)
  tab <- t(vapply(cols, function(cn) {
    v <- e[[cn]]
    c(Mean = mean(v, na.rm = TRUE), SD = stats::sd(v, na.rm = TRUE),
      Min = min(v, na.rm = TRUE), Max = max(v, na.rm = TRUE),
      `# efficient` = sum(v > 1 - 1e-6, na.rm = TRUE))
  }, numeric(5)))
  structure(list(table = tab, settings = object$settings, n = nrow(e),
                 synergy = table(e$synergy)), class = "summary.ndea")
}

#' @export
print.summary.ndea <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("Generalized network DEA (", toupper(s$rts), ", ", s$orientation,
      "-oriented), ", x$n, " DMUs\n\n", sep = "")
  print(round(x$table, digits))
  cat("\nSynergy:\n")
  print(x$synergy)
  invisible(x)
}

#' @export
as.data.frame.ndea <- function(x, ...) x$efficiency
