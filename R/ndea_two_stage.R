#' Two-Stage Network DEA
#'
#' Efficiency of decision making units (DMUs) whose production process has two
#' stages: stage 1 turns inputs `x` into intermediate products `z`, and stage 2
#' turns the intermediates into final outputs `y`.
#'
#' Two synergy mechanisms for the intermediate products are available.
#'
#' * `synergy = "flexible"` (default): the model of the MaxDEA Team (2026),
#'   which only requires the *total value* of the intermediates to be the same
#'   when produced by stage 1 and consumed by stage 2
#'   (\eqn{w^{(1)} z_0 = w^{(2)} z_0}). The system efficiency equals the
#'   product of the two stage efficiencies, and the stage efficiencies are the
#'   efficiencies of each stage evaluated on its own. Available for constant
#'   (`"crs"`) and variable (`"vrs"`) returns to scale and both orientations.
#' * `synergy = "rigid"`: the classical centralized model (Kao and Hwang, 2008;
#'   Liang, Cook and Zhu, 2008), which forces the intermediate prices of the two
#'   stages to be equal. Only constant returns to scale is implemented. The
#'   system efficiency is at most the product of the stand-alone stage
#'   efficiencies and the decomposition into stage efficiencies is generally
#'   not unique, so only the stand-alone stage efficiencies are reported.
#'
#' Data are internally rescaled by column means for numerical stability; scores
#' do not depend on measurement units. Projections are radial and no
#' second-phase slack maximisation is performed.
#'
#' @param x,y,z Numeric matrices or data frames (rows = DMUs, columns =
#'   variables) of inputs, final outputs and intermediate products. Values must
#'   be non-negative.
#' @param rts Returns to scale: `"crs"` (constant) or `"vrs"` (variable).
#' @param orientation `"input"` (minimise inputs) or `"output"` (maximise final
#'   outputs).
#' @param synergy `"flexible"` or `"rigid"`; see Details.
#' @param multipliers Logical; if `TRUE`, also solve the multiplier model and
#'   return the weights (shadow prices). Weights are not necessarily unique.
#' @param dmu Optional DMU names. Defaults to the row names of `x`.
#'
#' @return An object of class `"ndea_two_stage"`, a list with
#'   * `efficiency`: data frame with the system efficiency (`system`) and the
#'     stage efficiencies (`stage1`, `stage2`), each in (0, 1];
#'   * `lambda1`, `lambda2`: reference-set intensities of stage 1 and stage 2
#'     (rows = evaluated DMUs, columns = reference DMUs);
#'   * `targets`: radial projections of inputs, intermediates and outputs, in
#'     original units;
#'   * `weights` (only if `multipliers = TRUE`);
#'   * `status`: lpSolve status codes (0 = solved).
#'
#' @references
#' Kao, C. and Hwang, S.-N. (2008). Efficiency decomposition in two-stage data
#' envelopment analysis: An application to non-life insurance companies in
#' Taiwan. *European Journal of Operational Research*, 185(1), 418-429.
#'
#' Liang, L., Cook, W. D. and Zhu, J. (2008). DEA models for two-stage processes:
#' Game approach and efficiency decomposition. *Naval Research Logistics*, 55(7),
#' 643-653.
#'
#' MaxDEA Team (2026). Innovative Two-Stage and Generalized Network DEA.
#' MaxDEA Software Ltd. (Not peer reviewed.)
#'
#' @examples
#' d <- ndea_simulate(n = 20, seed = 1)
#' res <- ndea_two_stage(d$x, d$z, d$y, rts = "vrs", orientation = "input")
#' res
#' summary(res)
#' @export
ndea_two_stage <- function(x, z, y, rts = c("crs", "vrs"),
                           orientation = c("input", "output"),
                           synergy = c("flexible", "rigid"),
                           multipliers = FALSE, dmu = NULL) {
  rts <- .rts(rts)
  orientation <- .orient(orientation)
  synergy <- match.arg(synergy)
  if (synergy == "rigid" && rts == "vrs") {
    stop("synergy = \"rigid\" is only implemented for rts = \"crs\".", call. = FALSE)
  }
  X0 <- .prep_matrix(x, "x"); Z0 <- .prep_matrix(z, "z"); Y0 <- .prep_matrix(y, "y")
  n <- nrow(X0)
  if (nrow(Z0) != n || nrow(Y0) != n) {
    stop("`x`, `z` and `y` must have the same number of rows.", call. = FALSE)
  }
  if (n < 2) stop("At least two DMUs are required.", call. = FALSE)
  if (is.null(dmu)) {
    dmu <- rownames(X0)
    if (is.null(dmu)) dmu <- rownames(Z0)
    if (is.null(dmu)) dmu <- rownames(Y0)
    if (is.null(dmu)) dmu <- as.character(seq_len(n))
  }
  if (length(dmu) != n) stop("`dmu` must have one name per row.", call. = FALSE)

  sx <- .scale_factors(X0); sz <- .scale_factors(Z0); sy <- .scale_factors(Y0)
  X <- sweep(X0, 2, sx, "/"); Z <- sweep(Z0, 2, sz, "/"); Y <- sweep(Y0, 2, sy, "/")

  eff <- data.frame(dmu = dmu, system = NA_real_, stage1 = NA_real_,
                    stage2 = NA_real_, stringsAsFactors = FALSE)
  l1 <- l2 <- matrix(NA_real_, n, n, dimnames = list(dmu, dmu))
  tg_x <- matrix(NA_real_, n, ncol(X), dimnames = list(dmu, colnames(X0)))
  tg_z <- matrix(NA_real_, n, ncol(Z), dimnames = list(dmu, colnames(Z0)))
  tg_y <- matrix(NA_real_, n, ncol(Y), dimnames = list(dmu, colnames(Y0)))
  status <- integer(n)
  wts <- if (multipliers) vector("list", n) else NULL

  for (o in seq_len(n)) {
    s1 <- .radial_env(X, Z, o, rts, orientation)
    s2 <- .radial_env(Z, Y, o, rts, orientation)
    eff$stage1[o] <- s1$eff
    eff$stage2[o] <- s2$eff
    r <- if (synergy == "flexible") .ts_flex_env(X, Z, Y, o, rts, orientation)
         else .ts_rigid_env(X, Z, Y, o, orientation)
    .check_status(r, dmu[o], "ndea_two_stage")
    status[o] <- r$status
    eff$system[o] <- r$eff
    l1[o, ] <- r$lambda1
    l2[o, ] <- r$lambda2
    if (!is.na(r$eff)) {
      tau <- if (synergy == "flexible") r$tau else NA_real_
      if (orientation == "input") {
        tg_x[o, ] <- r$factor * X[o, ] * sx
        tg_y[o, ] <- Y[o, ] * sy
      } else {
        tg_x[o, ] <- X[o, ] * sx
        tg_y[o, ] <- r$factor * Y[o, ] * sy
      }
      tg_z[o, ] <- if (is.na(tau)) NA_real_ else tau * Z[o, ] * sz
    }
    if (multipliers) {
      mm <- .ts_mult(X, Z, Y, o, rts, orientation, synergy)
      w <- mm$weights
      w$v <- w$v / sx; w$w1 <- w$w1 / sz; w$w2 <- w$w2 / sz; w$u <- w$u / sy
      names(w$v) <- colnames(X0); names(w$u) <- colnames(Y0)
      names(w$w1) <- names(w$w2) <- colnames(Z0)
      wts[[o]] <- w
    }
  }
  if (!is.null(wts)) names(wts) <- dmu
  out <- list(efficiency = eff, lambda1 = l1, lambda2 = l2,
              targets = list(x = tg_x, z = tg_z, y = tg_y),
              weights = wts, status = status,
              call = match.call(),
              settings = list(rts = rts, orientation = orientation,
                              synergy = synergy))
  class(out) <- c("ndea_two_stage", "ndea")
  out
}

#' @export
print.ndea_two_stage <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("Two-stage network DEA (", s$synergy, " synergy, ", toupper(s$rts), ", ",
      s$orientation, "-oriented)\n", sep = "")
  cat("DMUs:", nrow(x$efficiency), "\n\n")
  print(format(x$efficiency, digits = digits), row.names = FALSE)
  invisible(x)
}

#' @export
summary.ndea_two_stage <- function(object, ...) {
  e <- object$efficiency
  cols <- c("system", "stage1", "stage2")
  tab <- t(vapply(cols, function(cn) {
    v <- e[[cn]]
    c(Mean = mean(v, na.rm = TRUE), SD = stats::sd(v, na.rm = TRUE),
      Min = min(v, na.rm = TRUE), Max = max(v, na.rm = TRUE),
      `# efficient` = sum(v > 1 - 1e-6, na.rm = TRUE))
  }, numeric(5)))
  structure(list(table = tab, settings = object$settings,
                 n = nrow(e)), class = "summary.ndea_two_stage")
}

#' @export
print.summary.ndea_two_stage <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("Two-stage network DEA (", s$synergy, " synergy, ", toupper(s$rts), ", ",
      s$orientation, "-oriented), ", x$n, " DMUs\n\n", sep = "")
  print(round(x$table, digits))
  invisible(x)
}

#' @export
as.data.frame.ndea_two_stage <- function(x, ...) x$efficiency
