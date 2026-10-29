configs <- expand.grid(rts = c("crs", "vrs"), orientation = c("input", "output"),
                       stringsAsFactors = FALSE)

test_that("hand-computed example: product form of efficiency (MaxDEA manual 9.1.5 setup)", {
  # A: 100 -> 100 -> 100 ; B: 100 -> 10 -> 9
  x <- matrix(c(100, 100), ncol = 1)
  z <- matrix(c(100, 10), ncol = 1)
  y <- matrix(c(100, 9), ncol = 1)
  for (ori in c("input", "output")) {
    r <- ndea_two_stage(x, z, y, rts = "crs", orientation = ori)
    expect_equal(r$efficiency$system, c(1, 0.09), tolerance = 1e-8)
    expect_equal(r$efficiency$stage1, c(1, 0.1), tolerance = 1e-8)
    expect_equal(r$efficiency$stage2, c(1, 0.9), tolerance = 1e-8)
    rr <- ndea_two_stage(x, z, y, rts = "crs", orientation = ori, synergy = "rigid")
    expect_equal(rr$efficiency$system, c(1, 0.09), tolerance = 1e-8)
  }
})

test_that("single input/output reduces to the textbook CCR ratio", {
  x <- matrix(c(2, 3, 4, 5, 6), ncol = 1)
  y <- matrix(c(2, 5, 6, 6, 9), ncol = 1)
  for (ori in c("input", "output")) {
    e <- vapply(1:5, function(o) .radial_env(x, y, o, "crs", ori)$eff, numeric(1))
    expect_equal(e, (y / x)[, 1] / max(y / x), tolerance = 1e-9)
  }
})

test_that("flexible synergy: system efficiency equals product of stage efficiencies", {
  d <- ndea_simulate(n = 25, seed = 42)
  for (i in seq_len(nrow(configs))) {
    r <- ndea_two_stage(d$x, d$z, d$y, rts = configs$rts[i],
                        orientation = configs$orientation[i])
    e <- r$efficiency
    expect_true(all(r$status == 0))
    expect_equal(e$system, e$stage1 * e$stage2, tolerance = 1e-7)
    expect_true(all(e$system > 0 & e$system <= 1 + 1e-9))
    expect_true(any(abs(e$system - 1) < 1e-7))
  }
})

test_that("envelopment and multiplier forms agree (LP duality)", {
  set.seed(7)
  for (rep in 1:3) {
    d <- ndea_simulate(n = 15, m = sample(1:3, 1), k = sample(1:3, 1), s = sample(1:3, 1))
    sc <- function(m) sweep(m, 2, colMeans(m), "/")
    X <- sc(d$x); Z <- sc(d$z); Y <- sc(d$y)
    for (syn in c("flexible", "rigid")) for (i in seq_len(nrow(configs))) {
      rts <- configs$rts[i]; ori <- configs$orientation[i]
      if (syn == "rigid" && rts == "vrs") next
      for (o in c(1, 8, 15)) {
        env <- if (syn == "flexible") .ts_flex_env(X, Z, Y, o, rts, ori)$eff
               else .ts_rigid_env(X, Z, Y, o, ori)$eff
        mul <- .ts_mult(X, Z, Y, o, rts, ori, syn)$eff
        expect_equal(env, mul, tolerance = 1e-7)
      }
    }
  }
})

test_that("rigid synergy is bounded by flexible synergy and by the stage product", {
  d <- ndea_simulate(n = 25, seed = 3)
  f <- ndea_two_stage(d$x, d$z, d$y, rts = "crs", synergy = "flexible")
  r <- ndea_two_stage(d$x, d$z, d$y, rts = "crs", synergy = "rigid")
  expect_true(all(r$efficiency$system <= f$efficiency$system + 1e-8))
  expect_true(all(r$efficiency$system <=
                    r$efficiency$stage1 * r$efficiency$stage2 + 1e-8))
  # CRS: input and output orientation coincide
  ro <- ndea_two_stage(d$x, d$z, d$y, rts = "crs", orientation = "output", synergy = "rigid")
  expect_equal(r$efficiency$system, ro$efficiency$system, tolerance = 1e-7)
})

test_that("scores are invariant to measurement units", {
  d <- ndea_simulate(n = 20, seed = 9)
  for (i in seq_len(nrow(configs))) {
    a <- ndea_two_stage(d$x, d$z, d$y, rts = configs$rts[i], orientation = configs$orientation[i])
    b <- ndea_two_stage(sweep(d$x, 2, c(1000, 0.01), "*"), sweep(d$z, 2, c(5, 1e4), "*"),
                        sweep(d$y, 2, c(1e-3, 7), "*"), rts = configs$rts[i],
                        orientation = configs$orientation[i])
    expect_equal(a$efficiency$system, b$efficiency$system, tolerance = 1e-6)
  }
})

test_that("reference sets and targets are feasible and consistent", {
  d <- ndea_simulate(n = 20, seed = 5)
  r <- ndea_two_stage(d$x, d$z, d$y, rts = "vrs", orientation = "input")
  e <- r$efficiency$system
  for (o in 1:20) {
    expect_equal(sum(r$lambda2[o, ]), 1, tolerance = 1e-7)
    xr <- as.vector(t(d$x) %*% r$lambda1[o, ])
    expect_true(all(xr <= r$targets$x[o, ] + 1e-6))
    yr <- as.vector(t(d$y) %*% r$lambda2[o, ])
    expect_true(all(yr >= d$y[o, ] - 1e-6))
    expect_equal(r$targets$x[o, ], e[o] * d$x[o, ], tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("multipliers are returned in original units and satisfy normalisation", {
  d <- ndea_simulate(n = 12, seed = 8)
  r <- ndea_two_stage(d$x, d$z, d$y, rts = "crs", orientation = "input", multipliers = TRUE)
  w <- r$weights[[3]]
  expect_equal(sum(w$v * d$x[3, ]), 1, tolerance = 1e-7)
  expect_equal(sum(w$u * d$y[3, ]), r$efficiency$system[3], tolerance = 1e-7)
  expect_equal(sum(w$w1 * d$z[3, ]), sum(w$w2 * d$z[3, ]), tolerance = 1e-7)
})

test_that("input validation", {
  d <- ndea_simulate(n = 10, seed = 1)
  expect_error(ndea_two_stage(d$x, d$z, d$y, rts = "vrs", synergy = "rigid"), "crs")
  expect_error(ndea_two_stage(d$x[-1, ], d$z, d$y), "same number of rows")
  neg <- d$x; neg[1, 1] <- -1
  expect_error(ndea_two_stage(neg, d$z, d$y), "negative")
  na <- d$y; na[2, 1] <- NA
  expect_error(ndea_two_stage(d$x, d$z, na), "missing")
  expect_error(ndea_two_stage(d$x[1, , drop = FALSE], d$z[1, , drop = FALSE],
                              d$y[1, , drop = FALSE]), "two DMUs")
})

test_that("methods work", {
  d <- ndea_simulate(n = 10, seed = 1)
  r <- ndea_two_stage(d$x, d$z, d$y, rts = "vrs", orientation = "output")
  expect_output(print(r), "Two-stage network DEA")
  expect_s3_class(as.data.frame(r), "data.frame")
  expect_output(print(summary(r)), "system")
  expect_equal(r$efficiency$dmu, rownames(d$x))
})
