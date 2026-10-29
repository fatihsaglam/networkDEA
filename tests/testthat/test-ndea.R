make_twostage <- function(seed = 11) {
  d <- ndea_simulate(n = 25, m = 2, k = 2, s = 2, seed = seed)
  list(d = d, dat = data.frame(d$x, d$z, d$y),
       net = network(node("s1", inputs = c("x1", "x2"), z_out = c("z1", "z2")),
                     node("s2", z_in = c("z1", "z2"), outputs = c("y1", "y2"))))
}

test_that("two-node generalized model reproduces two-stage flexible synergy", {
  m <- make_twostage()
  for (rts in c("crs", "vrs")) for (ori in c("input", "output")) {
    g <- ndea(m$dat, m$net, rts = rts, orientation = ori)
    t2 <- ndea_two_stage(m$d$x, m$d$z, m$d$y, rts = rts, orientation = ori)
    expect_true(all(g$status == 0))
    expect_equal(g$efficiency$system, t2$efficiency$system, tolerance = 1e-7)
    expect_equal(g$efficiency$s1, t2$efficiency$stage1, tolerance = 1e-7)
    expect_equal(g$efficiency$s2, t2$efficiency$stage2, tolerance = 1e-7)
    expect_true(all(g$efficiency$synergy == "lossless"))
    expect_equal(g$efficiency$D_star, rep(0, 25), tolerance = 1e-7)
  }
})

test_that("stand-alone node efficiency matches the radial envelopment model", {
  m <- make_twostage()
  sc <- function(v) v / mean(v)
  X <- sapply(m$d$x[, 1:2, drop = FALSE] |> as.data.frame(), sc)
  Y <- sapply(m$d$y |> as.data.frame(), sc)
  nd <- list(X = X, Zin = matrix(0, 25, 0), Y = Y, Zout = matrix(0, 25, 0), name = "n")
  for (rts in c("crs", "vrs")) for (ori in c("input", "output")) {
    a <- vapply(1:25, function(o) .node_standalone(nd, o, rts, ori), numeric(1))
    b <- vapply(1:25, function(o) .radial_env(X, Y, o, rts, ori)$eff, numeric(1))
    expect_equal(a, b, tolerance = 1e-7)
  }
})

test_that("serial three-node network: product form and node weights", {
  set.seed(5); n <- 30
  x1 <- runif(n, 10, 100); x2 <- runif(n, 10, 100)
  ineff <- function() exp(-abs(rnorm(n, 0, .25)))
  za <- (0.6 * x1 + 0.2 * x2) * ineff(); zb <- 0.5 * za * ineff()
  y1 <- 0.8 * zb * ineff(); y2 <- (0.3 * zb + 0.1 * za) * ineff()
  dat <- data.frame(x1, x2, za, zb, y1, y2)
  net <- network(node("A", inputs = c("x1", "x2"), z_out = "za"),
                 node("B", z_in = "za", z_out = "zb"),
                 node("C", z_in = "zb", outputs = c("y1", "y2")))
  g <- ndea(dat, net, rts = "vrs", orientation = "input")
  e <- g$efficiency
  expect_equal(e$system, e$A * e$B * e$C, tolerance = 1e-9)
  g2 <- ndea(dat, net, rts = "vrs", orientation = "input", omega = rep(1 / 3, 3))
  expect_equal(g2$efficiency$system, (e$A * e$B * e$C)^(1 / 3), tolerance = 1e-9)
  # node order does not matter
  net2 <- network(node("C", z_in = "zb", outputs = c("y1", "y2")),
                  node("A", inputs = c("x1", "x2"), z_out = "za"),
                  node("B", z_in = "za", z_out = "zb"))
  g3 <- ndea(dat, net2, rts = "vrs", orientation = "input")
  expect_equal(g3$efficiency$system, e$system, tolerance = 1e-8)
})

test_that("lossy synergy is detected, bounded by stand-alone efficiencies", {
  set.seed(2); n <- 25
  x1 <- runif(n, 10, 100); x2 <- runif(n, 10, 100)
  yA <- 0.7 * x1 * exp(-abs(rnorm(n, 0, .2)))
  za <- runif(n, 1, 5)
  y1 <- (0.5 * za + 0.01 * x2) * exp(-abs(rnorm(n, 0, .2)))
  dat <- data.frame(x1, x2, yA, za, y1)
  net <- network(node("A", inputs = "x1", outputs = "yA", z_out = "za"),
                 node("B", inputs = "x2", z_in = "za", outputs = "y1"))
  for (rts in c("crs", "vrs")) for (ori in c("input", "output")) {
    g <- ndea(dat, net, rts = rts, orientation = ori)
    e <- g$efficiency
    expect_true(all(g$status == 0))
    expect_true(any(e$synergy == "lossy"))
    expect_true(all(e$D_star <= 1e-7))
    expect_true(all(as.matrix(e[c("A", "B")]) <= g$standalone + 1e-7))
    expect_true(all(e$system > 0 & e$system <= 1 + 1e-8))
    # lossless DMUs keep their stand-alone node efficiencies
    ok <- e$synergy == "lossless"
    expect_equal(as.matrix(e[ok, c("A", "B")]), unname(g$standalone[ok, ]),
                 ignore_attr = TRUE, tolerance = 1e-7)
  }
})

test_that("scores are invariant to measurement units", {
  m <- make_twostage(3)
  dat2 <- m$dat
  dat2$x1 <- dat2$x1 * 1e3; dat2$z2 <- dat2$z2 * 1e-3; dat2$y2 <- dat2$y2 * 55
  for (rts in c("crs", "vrs")) {
    a <- ndea(m$dat, m$net, rts = rts); b <- ndea(dat2, m$net, rts = rts)
    expect_equal(a$efficiency$system, b$efficiency$system, tolerance = 1e-6)
  }
})

test_that("weights in original units satisfy the node normalisation", {
  m <- make_twostage(4)
  g <- ndea(m$dat, m$net, rts = "crs", orientation = "input")
  w <- g$weights[[2]]
  expect_gte(sum(w$s1$v * unlist(m$dat[2, c("x1", "x2")])), 1 - 1e-7)
  expect_gte(sum(w$s2$w_in * unlist(m$dat[2, c("z1", "z2")])), 1 - 1e-7)
})

test_that("network and node validation", {
  expect_error(node("A"), "at least one input")
  expect_error(node("A", inputs = "x"), "at least one output")
  expect_error(network(node("A", inputs = "x", outputs = "y")), "at least two nodes")
  a <- node("A", inputs = "x", z_out = "z"); b <- node("B", z_in = "z", outputs = "y")
  expect_error(network(a, node("A", z_in = "z", outputs = "y")), "unique")
  expect_error(network(a, node("B", inputs = "x", z_in = "z", outputs = "y")), "more than one node")
  expect_warning(network(a, node("B", inputs = "w", outputs = "y")), "not both produced")
  net <- network(a, b)
  expect_output(print(net), "Network with 2 nodes")
  dat <- data.frame(x = 1:5 + 1, z = 1:5 + 2, y = 1:5 + 3)
  expect_error(ndea(dat[, c("x", "z")], net), "not found")
  expect_error(ndea(dat, net, omega = c(1, 5)), "omega")
  expect_error(ndea(dat, net, omega = c(0.2, 0.2)), "sum")
  expect_error(ndea(dat, "net"), "network")
})

test_that("methods work", {
  m <- make_twostage()
  g <- ndea(m$dat, m$net)
  expect_output(print(g), "Generalized network DEA")
  expect_output(print(summary(g)), "Synergy")
  expect_s3_class(as.data.frame(g), "data.frame")
})
