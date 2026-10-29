# networkDEA

Network data envelopment analysis (DEA) in R: efficiency of decision making
units whose internal structure consists of connected sub-processes (nodes)
linked by intermediate products.

## Installation

```r
# install.packages("remotes")
remotes::install_github("fatihsaglam/networkDEA")
```

## What is implemented

| Function | Model |
|---|---|
| `ndea_two_stage()` | Two-stage network DEA. Flexible synergy (CRS/VRS, input/output oriented) and classical rigid synergy (CRS). |
| `ndea()` + `network()` + `node()` | Generalized network DEA for any serial/parallel/mixed node structure; reports stand-alone node efficiencies, synergy loss `D*` and system efficiency `prod(E_t^omega_t)`. |

Linear programs are solved with `lpSolve`.

## Example

```r
library(networkDEA)

d <- ndea_simulate(n = 20, seed = 1)

# two-stage, variable returns to scale, input oriented
res <- ndea_two_stage(d$x, d$z, d$y, rts = "vrs", orientation = "input")
summary(res)

# generalized: define the network yourself
dat <- data.frame(d$x, d$z, d$y)
net <- network(
  node("stage1", inputs = c("x1", "x2"), z_out = c("z1", "z2")),
  node("stage2", z_in = c("z1", "z2"), outputs = c("y1", "y2"))
)
ndea(dat, net, rts = "vrs", orientation = "input")
```

## Sources and attribution

This is an independent implementation of published methods; it is not affiliated
with MaxDEA Software Ltd. and contains no MaxDEA code, text or data.

* Kao, C. and Hwang, S.-N. (2008), *European Journal of Operational Research* 185(1), 418-429.
* Liang, L., Cook, W. D. and Zhu, J. (2008), *Naval Research Logistics* 55(7), 643-653.
* MaxDEA Team (2026). *Innovative Two-Stage and Generalized Network DEA*. MaxDEA Software Ltd.
  (flexible-synergy two-stage model and generalized network model; not peer reviewed).

## Status

Version 0.1.0. The results of the flexible-synergy and generalized models have been
verified against the properties stated in the sources (product decomposition, LP
duality, equivalence of the two-node generalized model with the two-stage model)
but not yet against an independent implementation. Please treat results with care.

## License

MIT
