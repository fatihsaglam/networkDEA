# Linear programs for two-stage network DEA ---------------------------------
#
# Notation: X (n x m) inputs, Z (n x k) intermediates, Y (n x s) final outputs;
# stage 1 turns X into Z, stage 2 turns Z into Y. `o` is the evaluated DMU.
# All data are expected to be pre-scaled (see .scale_factors()).

# Stand-alone radial (CCR/BCC) envelopment model for one stage. Returns the
# efficiency score in (0, 1] (input orientation: theta, output orientation:
# 1/phi).
.radial_env <- function(X, Y, o, rts, orientation) {
  n <- nrow(X); m <- ncol(X); s <- ncol(Y)
  obj <- c(1, rep(0, n))
  if (orientation == "input") {
    A <- rbind(cbind(-X[o, ], t(X)), cbind(0, t(Y)))
    dir <- c(rep("<=", m), rep(">=", s))
    rhs <- c(rep(0, m), Y[o, ])
    sense <- "min"
  } else {
    A <- rbind(cbind(0, t(X)), cbind(-Y[o, ], t(Y)))
    dir <- c(rep("<=", m), rep(">=", s))
    rhs <- c(X[o, ], rep(0, s))
    sense <- "max"
  }
  if (rts == "vrs") {
    A <- rbind(A, c(0, rep(1, n)))
    dir <- c(dir, "=")
    rhs <- c(rhs, 1)
  }
  r <- .solve_lp(sense, obj, A, dir, rhs)
  eff <- if (is.na(r$objval)) NA_real_ else if (orientation == "input") r$objval else 1 / r$objval
  list(status = r$status, eff = eff, lambda = r$solution[-1])
}

# Two-stage network DEA, flexible synergy (MaxDEA Team 2026, models 9.29,
# 9.31, 9.36, 9.41), envelopment form.
.ts_flex_env <- function(X, Z, Y, o, rts, orientation) {
  n <- nrow(X); m <- ncol(X); k <- ncol(Z); s <- ncol(Y)
  nv <- 2 + 2 * n
  i_l1 <- 2 + seq_len(n)
  i_l2 <- 2 + n + seq_len(n)
  blk <- function(r) matrix(0, r, nv)
  Ax <- blk(m);  Ax[, i_l1] <- t(X)
  Az1 <- blk(k); Az1[, i_l1] <- t(Z); Az1[, 2] <- -Z[o, ]
  Az2 <- blk(k); Az2[, i_l2] <- t(Z); Az2[, 2] <- -Z[o, ]
  Ay <- blk(s);  Ay[, i_l2] <- t(Y)
  if (orientation == "input") {
    Ax[, 1] <- -X[o, ]
    rhs <- c(rep(0, m + 2 * k), Y[o, ])
    sense <- "min"
  } else {
    Ay[, 1] <- -Y[o, ]
    rhs <- c(X[o, ], rep(0, 2 * k + s))
    sense <- "max"
  }
  A <- rbind(Ax, Az1, Az2, Ay)
  dir <- c(rep("<=", m), rep(">=", k), rep("<=", k), rep(">=", s))
  if (rts == "vrs") {
    r1 <- numeric(nv); r2 <- numeric(nv)
    if (orientation == "input") {
      r1[i_l1] <- 1; r1[2] <- -1; b1 <- 0
      r2[i_l2] <- 1;              b2 <- 1
    } else {
      r1[i_l1] <- 1;              b1 <- 1
      r2[i_l2] <- 1; r2[2] <- -1; b2 <- 0
    }
    A <- rbind(A, r1, r2)
    dir <- c(dir, "=", "=")
    rhs <- c(rhs, b1, b2)
  }
  obj <- c(1, rep(0, nv - 1))
  r <- .solve_lp(sense, obj, A, dir, rhs)
  eff <- if (is.na(r$objval)) NA_real_ else if (orientation == "input") r$objval else 1 / r$objval
  list(status = r$status, eff = eff, factor = r$objval, tau = r$solution[2],
       lambda1 = r$solution[i_l1], lambda2 = r$solution[i_l2])
}

# Two-stage network DEA, classical centralized model with rigid synergy
# (Kao and Hwang 2008; Liang, Cook and Zhu 2008; MaxDEA manual 9.26),
# constant returns to scale. Input and output orientation give the same score.
.ts_rigid_env <- function(X, Z, Y, o, orientation) {
  n <- nrow(X); m <- ncol(X); k <- ncol(Z); s <- ncol(Y)
  nv <- 1 + 2 * n
  i_l1 <- 1 + seq_len(n)
  i_l2 <- 1 + n + seq_len(n)
  Ax <- matrix(0, m, nv);  Ax[, i_l1] <- t(X)
  Az <- matrix(0, k, nv);  Az[, i_l1] <- t(Z); Az[, i_l2] <- -t(Z)
  Ay <- matrix(0, s, nv);  Ay[, i_l2] <- t(Y)
  if (orientation == "input") {
    Ax[, 1] <- -X[o, ]
    rhs <- c(rep(0, m + k), Y[o, ])
    sense <- "min"
  } else {
    Ay[, 1] <- -Y[o, ]
    rhs <- c(X[o, ], rep(0, k + s))
    sense <- "max"
  }
  A <- rbind(Ax, Az, Ay)
  dir <- c(rep("<=", m), rep(">=", k), rep(">=", s))
  obj <- c(1, rep(0, nv - 1))
  r <- .solve_lp(sense, obj, A, dir, rhs)
  eff <- if (is.na(r$objval)) NA_real_ else if (orientation == "input") r$objval else 1 / r$objval
  list(status = r$status, eff = eff, factor = r$objval, tau = NA_real_,
       lambda1 = r$solution[i_l1], lambda2 = r$solution[i_l2])
}

# Multiplier form of the two-stage models. Used to report weights (shadow
# prices) and to verify LP duality against the envelopment forms.
.ts_mult <- function(X, Z, Y, o, rts, orientation, synergy) {
  n <- nrow(X); m <- ncol(X); k <- ncol(Z); s <- ncol(Y)
  vrs <- rts == "vrs"
  if (synergy == "rigid") {
    # variables: v (m), w (k), u (s)
    iv <- seq_len(m); iw <- m + seq_len(k); iu <- m + k + seq_len(s)
    nv <- m + k + s
    A1 <- matrix(0, n, nv); A1[, iv] <- X; A1[, iw] <- -Z
    A2 <- matrix(0, n, nv); A2[, iw] <- Z; A2[, iu] <- -Y
    A <- rbind(A1, A2)
    dir <- rep(">=", 2 * n); rhs <- rep(0, 2 * n)
    obj <- numeric(nv); nrm <- numeric(nv)
    if (orientation == "input") {
      obj[iu] <- Y[o, ]; nrm[iv] <- X[o, ]; sense <- "max"
    } else {
      obj[iv] <- X[o, ]; nrm[iu] <- Y[o, ]; sense <- "min"
    }
    A <- rbind(A, nrm); dir <- c(dir, "="); rhs <- c(rhs, 1)
    r <- .solve_lp(sense, obj, A, dir, rhs)
    w <- list(v = r$solution[iv], w1 = r$solution[iw], w2 = r$solution[iw],
              u = r$solution[iu], u01 = NA_real_, u02 = NA_real_)
  } else {
    # variables: v (m), w1 (k), w2 (k), u (s), [u01, u02]
    iv <- seq_len(m); iw1 <- m + seq_len(k); iw2 <- m + k + seq_len(k)
    iu <- m + 2 * k + seq_len(s)
    nv <- m + 2 * k + s + if (vrs) 2 else 0
    i01 <- if (vrs) m + 2 * k + s + 1 else NA
    i02 <- if (vrs) m + 2 * k + s + 2 else NA
    A1 <- matrix(0, n, nv); A1[, iv] <- X; A1[, iw1] <- -Z
    A2 <- matrix(0, n, nv); A2[, iw2] <- Z; A2[, iu] <- -Y
    if (vrs) { A1[, i01] <- 1; A2[, i02] <- 1 }
    A <- rbind(A1, A2)
    dir <- rep(">=", 2 * n); rhs <- rep(0, 2 * n)
    link <- numeric(nv); link[iw1] <- Z[o, ]; link[iw2] <- -Z[o, ]
    obj <- numeric(nv); nrm <- numeric(nv)
    if (orientation == "input") {
      obj[iu] <- Y[o, ]; if (vrs) obj[i02] <- -1
      nrm[iv] <- X[o, ]
      if (vrs) link[i01] <- -1
      sense <- "max"
    } else {
      obj[iv] <- X[o, ]; if (vrs) obj[i01] <- 1
      nrm[iu] <- Y[o, ]
      if (vrs) link[i02] <- -1
      sense <- "min"
    }
    A <- rbind(A, nrm, link); dir <- c(dir, "=", "="); rhs <- c(rhs, 1, 0)
    r <- .solve_lp(sense, obj, A, dir, rhs, free = if (vrs) c(i01, i02) else integer(0))
    w <- list(v = r$solution[iv], w1 = r$solution[iw1], w2 = r$solution[iw2],
              u = r$solution[iu],
              u01 = if (vrs) r$solution[i01] else NA_real_,
              u02 = if (vrs) r$solution[i02] else NA_real_)
  }
  eff <- if (is.na(r$objval)) NA_real_ else if (orientation == "input") r$objval else 1 / r$objval
  list(status = r$status, eff = eff, weights = w)
}
